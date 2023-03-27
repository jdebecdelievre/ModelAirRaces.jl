function fixedpropulsion()
    propulsion = InertiaBuildUp(:fuselage)
    x_cg_propulsion   = envt.length_propulsion/2 # propeller at the nose
    propeller       = InertialElement(0.04, x_cg_propulsion) # kg based on Walrus prop + spinner
    @addnode propulsion propeller
    motor           = InertialElement(0.07, x_cg_propulsion)
    @addnode propulsion motor
    sumall!(propulsion)
    return propulsion
end

function fixedfuselage()
    """
    Fuselage and propulsion system constants assume the following lay out:
        ------------------------------------------------------> x
        | wing |           tailboom           | htail
    propulsion | battery | elec | 
    <-------------------------->
            plywood box surrounded by foam to make cylinder
    Here we evaluate the weight, length, wetted area, mass and inertia of the plywood box:
    We add 15 % length in the front and back to taper. 
    The propulsion system fits into the front taper.
    """
    fuselageInertia = InertiaBuildUp(:fuselage)
    spacing = 0.01 # spacing between components
    # Battery
    # https://hobbyking.com/en_us/turnigy-nano-tech-2200mah-3s-45-90c-lipo-pack.html?queryID=54916626f7f0b6271392ef525792077b&objectID=30641&indexName=hbk_live_products_analytics
    # values below are measured from ordered item.
    # Matek: http://www.mateksys.com/?portfolio=f722-wing
    length_battery = 0.112 * 1.3 # battery length x 1.3
    x_cg_battery   = spacing + length_battery/2
    M_battery      = 0.225 * 1.1 # kg
    battery        = InertialElement(M_battery, x_cg_battery)
    @addnode fuselageInertia battery
    
    length_elec       = 0.055 * 3   # 3x Matek length
    x_cg_elec         = spacing + length_battery + spacing + length_elec / 2
    matek             = InertialElement(0.04926, x_cg_elec) # kg
    GPS               = InertialElement(0.00929, x_cg_elec) # kg Mateksys m8Q
    RF_receiver       = InertialElement(0.00580, x_cg_elec) # d4r2 frsky receiver
    telemetry_antenna = InertialElement(0.02000, x_cg_elec) # kg SiK Holybro antenna
    misc_elec         = InertialElement(0.02000, x_cg_elec) # kg various wires mass, and tape to hold elec in place
    ESC               = InertialElement(0.05250, x_cg_elec) # based on Walrus 30A rated ESC
    electronics       = @addbranch fuselageInertia electronics(matek, GPS, RF_receiver, telemetry_antenna, misc_elec)


    ρ_plywood     = 640.7*1.3 # kg/m^3, 40 lb/ft^3 for good quality birch ply wood, x1.3 safety factor
    ply_thickness = 1/16 * 0.0254 # 1/16th of an inch
    ply_width     = 0.05 # box of 5cmx5cm to fit Matek
    ply_area      = ply_width^2 - (ply_width-2*ply_thickness)^2
    ply_length    = spacing + length_battery + spacing + length_elec + spacing
    M_ply         = (ply_length * ply_area + 2 * ply_width^2*(ply_thickness)) * ρ_plywood # weight of fuselage plywood box per unit length
    x_cg_ply      = ply_length / 2
    I_ply         = 1/12 * M_ply * ((ply_length^2 + ply_width^2) - ((ply_length-ply_thickness)^2 + (ply_width-ply_thickness)^2)) ## yy and zz are the same, xx is 0
    plybox        = InertialElement(M_ply, x_cg_ply, 0., 0., 0., I_ply, I_ply)
    motor_mount   = InertialElement(0.02, envt.length_propulsion) # based on Walrus plywood mount
    @addbranch fuselageInertia frame(plybox, motor_mount)
    
    # Combine
    (; Iyy, Izz) = fuselageInertia.value
    @assert Iyy == Izz

    # Mass of foam to take square plywood box into cylinder
    cylinder_radius = ply_width / sqrt(2)
    M_foam = MaterialProperties.ρ_foam * ply_length * (π*cylinder_radius^2-ply_width^2)
    fuselageInertia[:foam] = InertialElement(M_foam, x_cg_ply)

    # Move cg backward by the length of the nose
    tot_width = 2 * cylinder_radius
    fineness_ratio_nose = 1.5
    fineness_ratio_rear = 2.0
    l_nose = fineness_ratio_nose*tot_width
    l_rear = fineness_ratio_rear*tot_width
    fuselageInertia = OptimUtils.scalar_oper(f->OptimUtils.translate(f, Δx=l_nose), fuselageInertia)

    # Add nose weight and CG
    M_nose = V_ellipsoid(cylinder_radius,cylinder_radius,l_nose) / 2 * MaterialProperties.ρ_foam
    xcg_nose = l_nose - centroid_half_ellipsoid(l_nose) 
    fuselageInertia[:nose] = InertialElement(M_nose, xcg_nose)
    
    # Add rear weight and CG
    M_rear = V_parabola(cylinder_radius,l_rear) * MaterialProperties.ρ_foam
    xcg_rear = l_nose + ply_length + centroid_parabola(cylinder_radius, l_rear)
    fuselageInertia[:rear] = InertialElement(M_rear, xcg_rear)
    
    # Tally mass and inertia
    sumall!(fuselageInertia)

    # Wetted area
    tot_length = ply_length+l_nose+l_rear
    Sw_nose = Sw_ellipsoid(cylinder_radius,cylinder_radius,l_nose) / 2
    Sw_rear = Sw_parabola(cylinder_radius, l_rear)
    Sw_fuselage = (2π * cylinder_radius * ply_length) + Sw_nose + Sw_rear
    
    # Aerodynamic center at end of tapering area
    aerodynamic_center = l_nose
    
    fixedPropulsion = fixedpropulsion()
    return (; fuselageInertia, fixedPropulsion, Sw_fuselage, tot_length, tot_width, cylinder_radius, aerodynamic_center, l_nose, l_rear)
end

Sw_ellipsoid(r_x, r_y, r_z, p=1.6075) = 4π * (((r_x*r_y)^p+(r_x*r_z)^p+(r_z*r_y)^p)/3)^(1/p) # paraboloid wetted area
V_ellipsoid(r_x, r_y, r_z) = 4π/3 * r_x * r_y * r_z # paraboloid volume
centroid_half_ellipsoid(h) = 3*h/8 # taken from base

"""
Returns the arclength between the parabola vertex and a point located at (x,h(x)).
Assumes h(x) = a*x^2+c. Source: https://en.wikipedia.org/wiki/Parabola#Arc_length
"""
function arclength_parabola(a, x)
    f = 1/(4*a) # arclength
    h = x / 2
    q = sqrt(f^2+h^2)
    return h*q/f + f*log((h+q)/f)
end

"""
Returns the wetted area of a parabola of revolution that starts non-zero radius with a zero slope
and reaches a radius of 0 in the provided length. The required parabola is h(x) = radius*(length^2-x^2)/length^2
``Sw = 2π \\int_0^length radius*(length^2-x^2)/length^2 dx``
"""
Sw_parabola(radius, length) = 4π/3*radius*length

"""
Returns the volume of a parabola of revolution that starts non-zero radius with a zero slope
and reaches a radius of 0 in the provided length. The required parabola is h(x) = radius*(length^2-x^2)/length^2
``V = π \\int_0^length (radius*(length^2-x^2)/length^2)^2 dx``
"""
V_parabola(radius, length) = π*radius^2*length*8/15

"""
Returns the centroid of a parabola of revolution that starts non-zero radius with a zero slope
and reaches a radius of 0 in the provided length. The required parabola is h(x) = radius*(length^2-x^2)/length^2
``xcg = π \\int_0^length x(radius*(length^2-x^2)/length^2)^2 dx / Volume``
"""
centroid_parabola(radius, length) = length*5/16

function fixedfuselage_shapes()
    fuselage = fixedfuselage()
    (; tot_length, cylinder_radius, l_nose, l_rear) = fuselage

    # Nose cone
    
    # Propeller

    # Fuselage
    s = LinRange(0.,1.,10)
    y_parab = reverse(map(s->sqrt(1-s^2),s) * cylinder_radius)
    y_ellip = map(s->(1-s^2),s) * cylinder_radius
    y = [y_parab; y_ellip]
    x = [l_nose*s; reverse(tot_length .- l_rear*s)]
    x = [x; reverse(x)]
    y = [y; -reverse(y)]
    return ((x,y),)
end
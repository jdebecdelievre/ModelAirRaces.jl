const MaterialProperties = (;
    ## Constants for weight estimation
    ρ_balsa = 2e2, # kg/m^3
    tail_tc = 0.031, # 1/4th of an inch of balsa for horizontal and vertical tail
    E_balsa = 4e9,  # N/m**2 
    G_balsa = 0.3e9,  # N/m**2 
    MaxStressBalsa = 1e7,  # N/m**2 
    # Source: http://www.utc.fr/~hagegebe/UV/MQ12/CORRECTIONS_TD/%5BASHBY99%5D%20-%20Materials%20Selection%20In%20Mechanical%20Design%202Ed.pdf
    # page 347
    # Material selection in mechanical design, Michael Ashby

    # Constants used in calculation
    ρ_foam   = 28  , # kg/m^3
    ρ_carbon = 1600, # kg/m^3
    
    ρ_epoxy         = 1300   , # kg/m^3 to estimate the mass used to glue wing spar
    epoxy_thickness = 1e-3/4 , # m      1/4 mm of epoxy for wing and tail glueå
    
    M_servo        = 0.008, # assuming 4 Hitec HS-55
    M_wing_housing = 0.030, # kg
    M_tail_housing = 0.01,  # kg

    # Assuming Std fabric, 0/90/0 http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
    E_carbon        = 70e9 , # N/m**2 
    G_carbon        = 5e9 , # N/m**2
    MaxStressCarbon = 600e6, # N/m**2
    MaxShearCarbon  = 90e6, # N/m**2
)
abstract type BeamCrossSection{T} end
crosssectionarea!(::Vector, ::CS) where {CS <:BeamCrossSection} = raise("Type $CS did not implement crosssection area")
momentofarea!(::Vector, ::CS) where {CS <:BeamCrossSection} = raise("Type $CS did not implement crosssection area")
polarmoment!(::Vector, ::CS) where {CS <:BeamCrossSection} = raise("Type $CS did not implement crosssection area")
zmax(::CS, ::Int64) where {CS <:BeamCrossSection} = raise("Type $CS did not implement zmax")
struct CircularCrossSection{T}<:BeamCrossSection{T}
    r::Vector{T} # radius
    t::Vector{T} # thickness
    function CircularCrossSection(T, npan)
        return new{T}(zeros(T,npan),zeros(T,npan))
    end
end
crosssectionarea!(A::Vector, cs::CircularCrossSection) = map!((r,t)->π*(r^2-(r-t)^2), A, cs.r, cs.t)
momentofarea!(I::Vector, cs::CircularCrossSection) = map!((r,t)->π/4*(r^4-(r-t)^4), I, cs.r, cs.t)
polarmoment!(J::Vector, cs::CircularCrossSection) = map!((r,t)->π/4*(r^4-(r-t)^4), J, cs.r, cs.t)
zmax(cs::CircularCrossSection, i::Int64) = cs.r[i]
struct RectangularCrossSection{T}<:BeamCrossSection{T}
    w_i::Vector{T} # inner width (0 for filled-in beam)
    h_i::Vector{T} # inner height (0 for filled-in beam)
    w_o::Vector{T} # outer width
    h_o::Vector{T} # outer height
    function RectangularCrossSection(T, npan)
        return new{T}(zeros(T,npan),zeros(T,npan),zeros(T,npan),zeros(T,npan))
    end
end
crosssectionarea!(A::Vector, cs::RectangularCrossSection) = map!((w_o, h_o, w_i, h_i)->w_o*h_o-w_i*h_i, A, cs.w_o, cs.h_o, cs.w_i, cs.h_i)
momentofarea!(I::Vector, cs::RectangularCrossSection) = map!((w_o, h_o, w_i, h_i)->(w_o*h_o^3-w_i*h_i^3)/12, I, cs.w_o, cs.h_o, cs.w_i, cs.h_i)
polarmoment!(J::Vector, cs::RectangularCrossSection) = map!((w_o, h_o, w_i, h_i)->(w_o*h_o*(w_o^2+h_o^2)-w_i*h_i*(w_i^2+h_i^2))/12, J, cs.w_o, cs.h_o, cs.w_i, cs.h_i)
zmax(cs::RectangularCrossSection, i::Int64) = cs.h_o[i] / 2
struct StructMesh{T,CS<:BeamCrossSection}
    """
    ^z
    |_______________
    |_______________
    |-----------> y
    """
    # Geometry
    npan::Int64
    cs::CS        # cross section
    y::Vector{T}  # spanwise coordinate
    A::Vector{T}  # cross section area

    # Bending
    I::Vector{T}   # moment of area
    q::Vector{T}   # load
    S::Vector{T}   # shear force
    M::Vector{T}   # moment
    ϕ::Vector{T}   # deflection angle
    w::Vector{T}   # deflection distance
    σyy::Vector{T} # axial stress
    
    # Torsion
    J::Vector{T}   # polar moment
    tl::Vector{T}  # torsional load
    T::Vector{T}   # torque
    θ::Vector{T}   # torsion angle
    τxz::Vector{T} # shear stress
    function StructMesh(T, npan, cs::CS=CircularCrossSection(T,npan)) where {CS<:BeamCrossSection}
        return new{T,CS}(npan, cs, (zeros(T, npan) for _=1:14)...)
    end
end
StructMesh() = StructMesh(Float64, 10)
npanels(s::StructMesh{T}) where {T} = s.npan
crosssectionarea!(mesh::StructMesh) = crosssectionarea!(mesh.A, mesh.cs)
momentofarea!(mesh::StructMesh) = momentofarea!(mesh.I, mesh.cs)
polarmoment!(mesh::StructMesh) = polarmoment!(mesh.J, mesh.cs)

function updatecrosssection!(mesh::StructMesh)
    crosssectionarea!(mesh)
    momentofarea!(mesh)
    polarmoment!(mesh)
end


function cantilevered_beam_bending!(mesh, E)
# Ref: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10.pdf

    # load_4g = 4 * W / halfspan * (a1*np.sin(mesh['theta_c']) + a3*np.sin(3*mesh['theta_c']))
    N = length(mesh.y)
    S, M, w, ϕ, I, q, y = mesh.S, mesh.M, mesh.w, mesh.ϕ, mesh.I, mesh.q, mesh.y

    # Shear and moment
    S[N] = M[N] = 0.
    for i=N-1:-1:1
        S[i] = S[i+1] - (q[i+1] + q[i]) * (y[i+1] - y[i]) / 2
        M[i] = M[i+1] - (S[i+1] + S[i]) * (y[i+1] - y[i]) / 2
    end

    # Deflections
    w[1] = ϕ[1] = 0.
    for i=2:N
        ϕ[i] = ϕ[i-1] + (M[i] / I[i] + M[i-1] / I[i-1]) * (y[i] - y[i-1]) / 2 / E
        w[i] = w[i-1] - (ϕ[i] + ϕ[i-1]) * (y[i] - y[i-1]) / 2
    end

    # Axial stress
    map!(i-> mesh.M[i] / mesh.I[i] * zmax(mesh.cs,i), mesh.σyy, 1:N)
end

function cantilevered_beam_torsion!(mesh, G)
    # Ref: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10.pdf
    
        N = length(mesh.y)
        T, tl, θ, J, y = mesh.T, mesh.tl, mesh.θ, mesh.J, mesh.y
    
        # Torque
        T[N] = 0.
        for i=N-1:-1:1
            T[i] = T[i+1] - (tl[i+1] + tl[i]) * (y[i+1] - y[i]) / 2
        end
    
        # Torsion angle
        θ[1] = 0.
        for i=2:N
            θ[i] = θ[i-1] + (T[i] / J[i] + T[i-1] / J[i-1]) * (y[i] - y[i-1]) / 2 / G
        end

        # Shear stress
        map!(i-> mesh.T[i] / mesh.J[i] * zmax(mesh.cs,i), mesh.τxz, 1:N)
    end

function mass_estimation(mesh::StructMesh, ρ::Real)
    M = 0.
    for i=1:npanels(mesh)-1
        dy = mesh.y[i+1] - mesh.y[i]
        M += (mesh.A[i+1]+mesh.A[i])/2 * ρ * dy
    end
    return M
end

function inertia_estimation(mesh::StructMesh{T}, ρ::Real) where T
    inertia = InertialElement(T(0.))
    for i=1:npanels(mesh)-1
        dy = mesh.y[i+1] - mesh.y[i]
        yi = (mesh.y[i+1] + mesh.y[i]) / 2
        Ai = (mesh.A[i+1] + mesh.A[i]) / 2
        M = Ai * ρ * dy
        inertia += InertialElement(M, 0., yi)
    end
    return inertia
end

function structuralanalysis1(sharedconstants)
    ## Constants shared with other disciplines
    (; airfoil, fuselageInertia, max_rfuse) = sharedconstants
    airfoil_tc = geom(airfoil[1]).t_c

    ## Local and global variables
    struct_local = (;
        twing = Var(ini=.005, lb=.001, ub=.01, N=2),
        tfuse = Var(ini=.002, lb=.001, ub=.005, N=2),
        rwing = Var(ini=.008, lb=.004, ub=.015, N=2),
        rfuse = Var(ini=.005, lb=.001, ub= max_rfuse, N=2),
        xelec = Var(ini=0.1,  lb=  0., ub= 0.3, N=1),
    )

    struct_global = (; # bounds will be overriden by system level bounds
        cwing  = Var(ini=.10,lb= .05, ub=0.4, N=2),
        xwing  = Var(ini=.3,lb=  0., ub= 1.),
        bwing  = Var(ini= 1.4,lb=  .6, ub= 1.6),
        chtail = Var(ini=0.1,lb=0.05, ub=0.15),
        xhtail = Var(ini= 1.,lb= 0.6, ub= 1.3),
        xcg    = Var(ini= 0.25,lb= 0., ub= 5.),
        a1     = Var(ini= 0.,lb= -4., ub=4.),
        a3     = Var(ini= 0.,lb= -5., ub=5.),
        mass   = Var(ini=1.2,lb=.6, ub= 2.),
    )
    variables = merge(struct_global, struct_local)
    idx = indexbyname(variables)
    x = ini_scaled(variables)
    Nx = len(variables)

    # Cache structure
    chunksize = Nx
    TF        = Float64
    DTF       = ForwardDiff.Dual{Nothing,TF, chunksize}
    wmesh     = StructMesh(TF, 10)
    fmesh     = StructMesh(TF, 5)
    wmesh_d   = StructMesh(DTF, 10)
    fmesh_d   = StructMesh(DTF, 5)

    # Output
    eq = (; ((k=>Var(:eq)) for k=[:mass, :xcg])...)
    ineq = (
        sparradius     = Var(:ineq, N=2),
        wsparthickness = Var(:ineq, N=2),
        fsparthickness = Var(:ineq, N=2),
        wingstress     = Var(:ineq, N=npanels(wmesh)),
        wingdisp       = Var(:ineq, N=npanels(wmesh)),
        fusestress     = Var(:ineq, N=npanels(fmesh)),
        fuseangle      = Var(:ineq, N=npanels(fmesh)),
    )
 
    output     = merge(eq,ineq)
    idy        = indexbyname(output)
    
    c          = ini(output)
    cache      = (wmesh, fmesh)

    function structuralanalysis(x::Vector{T}, inertiaLog::BuildUp{TF}=NoBuildUp, c::Vector{T}=c, 
                                cache::Tuple{StructMesh{T},StructMesh{T}}=cache) where {T,TF}
        """
        Simple structures model. Models wing and tail boom as beams and performs a bending analysis.
        Computes mass and CG location.
        Ahead of the wing, the fuselage is a ply wood box that contains the electronics and motor casing.
        Assumes 0 sweep.
        ------------------------------------------------------> x
                                | wing | tailboom | htail
        propulsion | battery | elec | 
        <-------------------------->
                plywood box
        Tail assumed to go from wing leading to tail trailing edge.
        """
        # Extract & Unscale
        c.=0.
        v = unscale_unpack(x, idx, variables)
        wmesh, fmesh = cache

        # Fill-in mesh with input
        wmesh.y .= LinRange(0., v.bwing/2, npanels(wmesh))
        wmesh.r .= LinRange(v.rwing[1], v.rwing[2], npanels(wmesh))
        wmesh.t .= LinRange(v.twing[1], v.twing[2], npanels(wmesh))
        fmesh.y .= LinRange(0., v.xhtail+v.chtail - v.xwing, npanels(fmesh))
        fmesh.r .= LinRange(v.rfuse[1], v.rfuse[2], npanels(fmesh))
        fmesh.t .= LinRange(v.tfuse[1], v.tfuse[2], npanels(fmesh))
        wmesh.q .= 0.

        # Constrain wing spar to fit comfortably into wing foam core
        @. c[idy.sparradius] = (v.rwing - 0.4 * (v.cwing * airfoil_tc)) / (v.cwing * airfoil_tc);
        
        # Constrain spar thickness to fit into spar radius
        @. c[idy.wsparthickness] = (v.twing - v.rwing) / v.twing;
        @. c[idy.fsparthickness] = (v.tfuse - v.rfuse) / v.tfuse;
        
        # Mass estimation
        inertia = mass_estimation(v,sharedconstants,wmesh,fmesh,fuselageInertia,inertiaLog)
        c[idy.mass] = (inertia.mass - v.mass) / inertia.mass;
        c[idy.xcg] = (inertia.xcg - v.xcg) / inertia.xcg;

        # Evaluate loads on wing
        W = v.mass * envt.g
        a1 = 4 * W * v.a1 # elliptic loading that creates 4g of lift
        a3 = 4 * W * v.a3
        N = length(wmesh.y)
        for i=2:N
            θ = acos(2 * wmesh.y[i] / v.bwing) # Could be computed once and for all
            dy = (wmesh.y[i]-wmesh.y[i-1])
            wmesh.q[i] += (a1 * sin(θ) + a3 * sin(3*θ))*dy # unroll Fourrier transform
        end

        # Beam bending on wing
        @. wmesh.I = pi / 4 * (wmesh.r^4 - (wmesh.r - wmesh.t)^4) # circular beam
        cantilevered_beam_bending!(wmesh, E_carbon)
        @. c[idy.wingstress] = (abs(wmesh.M * wmesh.t/2 / wmesh.I) - MaxStressCarbon*.5)/ MaxStressCarbon; # stress
        @. c[idy.wingdisp]   = (abs(wmesh.w) - MaxTipDeflection*v.bwing) / (MaxTipDeflection*v.bwing); # deflection

        # # Loads on fuselage
        # lift_wing = sum(wmesh.q) + wing.mass * envt.g # lift on wing = loads + inertial relief
        # lift_tail = 4*W - lift_wing
        # for i=eachindex(fmesh.y)
        #     fmesh.q[i] = fmesh.y[i] * lift_tail
        # end

        # # Beam bending on fuselage
        # @. fmesh.I = pi / 4 * (fmesh.r^4 - (fmesh.r - fmesh.t)^4) # circular beam
        # cantilevered_beam_bending!(fmesh, E_carbon)
        # @. c[idy.fusestress] = (abs(fmesh.M * fmesh.t/2 / fmesh.I) - MaxStressCarbon*.8)/MaxStressCarbon; # stress
        # @. c[idy.fuseangle]  = (abs(fmesh.w) - MaxTailAngle)/MaxTailAngle; # angle

        return c, (wmesh, fmesh) # return cache
    end

    # Analysis and gradient together
    FDresults = DiffResults.JacobianResult(c, x);
    fval = FDresults.value
    dfval = FDresults.derivs[1];
    cfg = ForwardDiff.JacobianConfig(structuralanalysis, c, ForwardDiff.Chunk{chunksize}(), nothing);
    c_d = DTF.(c)
    f(x) = structuralanalysis(x, NoBuildUp, c_d, (wmesh_d, fmesh_d))[1]
    function structuralanalysis_wgrad(x::Vector)
        ForwardDiff.jacobian!(FDresults, f, x, cfg,  Val{false}());
        return fval, dfval
    end

    # struct_local = NamedTuple((:prout, Var(1., 1., 1,N=10)) for k=1:10)
    return structuralanalysis, structuralanalysis_wgrad, struct_global, struct_local, output
end
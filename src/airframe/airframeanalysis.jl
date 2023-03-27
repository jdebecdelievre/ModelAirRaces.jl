const sharpedge_gustspeed = 10 * 0.3048 # 30ft/s
 sharpedge_gustloadfactor(mass_si, Sref_si, V_si, CLα) = 1 + 1/2*envt.rho_atm*V_si^2*Sref_si*(CLα * sharpedge_gustspeed / V_si) / mass_si / envt.g

function smooth_gustloadfactor(mass_si, Sref_si, cref_si, V_si, CLα)
    weight_lb = mass_si * envt.g * 0.2248090795
    Sref_ft = Sref_si / (0.3048)^2
    cref_ft = cref_si / 0.3048
    WS  = weight_lb / Sref_ft
    g_imp = 32.174 # ft/s^2
    ρ_imp = 0.002377 # slugs/ft^3
    μg  = 2 * WS / (ρ_imp * cref_ft * CLα * g_imp)
    Kg  = 0.88*μg / (5.3+μg)
    Ude = 25 # ft/s
    V_kt = V_si / 0.514444
    return 1 + (Kg*Ude*V_kt*CLα)/(498*WS)
end

function mass_estimation(v, options, wmesh::StructMesh, fmesh::StructMesh, tmesh::StructMesh,
    fuselageInertia::InertialElement{T} where T, inertiaLog::BuildUp=NoBuildUp)
    airfoil = options.airfoil
    (; tot_length, tot_width) = options.fuselage
    (; ρ_balsa, ρ_foam, ρ_carbon, M_servo, M_wing_housing, M_tail_housing) = MaterialProperties

    # Compute wing mass, CG, moments of inertia, section moments of area
    x_spar       = v.xwing + v.cwing[1]/4
    spar         = InertialElement(0., x_spar)
    foam         = InertialElement(0., x_spar)
    wing_housing = InertialElement(M_wing_housing, x_spar)
    servo        = (InertialElement(M_servo, x_spar, v.bwing*3/8)
                  + InertialElement(M_servo, x_spar,-v.bwing*3/8))
    wmesh.q .= 0.
    (t_c, s_c) = geom(airfoil[1])
    
    for i=1:npanels(wmesh)-1
        dy = wmesh.y[i+1] - wmesh.y[i]
        yi = (wmesh.y[i+1] + wmesh.y[i]) / 2
        Ai = (wmesh.A[i+1] + wmesh.A[i]) / 2

        # Compute local chord and thickness
        chord = v.cwing[1] + yi * 2 / v.bwing * (v.cwing[2] - v.cwing[1])
        thickness = t_c * chord
        area = thickness * chord * s_c

        # Compute foam weight
        M_foam = area * ρ_foam * dy

        # wmesh.q[i+1] -= M_foam * envt.g / dy # wing inertial relief
        if yi>tot_width/2
            foam += InertialElement(M_foam, x_spar, yi)
            foam += InertialElement(M_foam, x_spar,-yi)
        end

        # Compute carbon weight
        M_spar = Ai * ρ_carbon * dy
        # wmesh.q[i+1] -= M_spar * envt.g / dy # wing inertial relief
        spar += InertialElement(M_spar, x_spar, yi)
        spar += InertialElement(M_spar, x_spar,-yi)
    end
    rightspar = inertia_estimation(wmesh, ρ_carbon)
    leftspar  = OptimUtils.rotate(rightspar, SA[1., 0., 0.], 180.)
    spar = OptimUtils.translate(rightspar+leftspar, Δx=x_spar)
    wing = spar+foam+wing_housing+servo
    @addbranch inertiaLog wing(spar,foam,wing_housing,servo)
    # Compute tailboom mass, inertia, CG
    tailboom = OptimUtils.rotate(inertia_estimation(fmesh, ρ_carbon), SA[0., 0., 1.], -90.)
    @addnode inertiaLog tailboom

    # Estimate mass and CG of tail
    (t_c, s_c, x_c) = geom(airfoil[3])
    # foam         = InertialElement(t_c * s_c * v.chtail^2 * ρ_foam * v.bhtail, v.xhtail + v.chtail )
    # carbon       = InertialElement(2 * 0.008 * 0.0008 * ρ_carbon * v.bhtail, v.xhtail + v.chtail * x_c ) # 8mm wide, .8mm thick, top and bottom 
    htail        = OptimUtils.translate(inertia_estimation(tmesh, ρ_balsa), Δx=v.xhtail)
    htail       += OptimUtils.rotate(htail, SA[1., 0., 0.], 180.)
    vtail        = InertialElement(htail.mass * 0.6, htail.xcg) # vertical tail assumed to be 60% of horizontal tail
    servo        = InertialElement(M_servo*2, tot_length)
    tail_housing = InertialElement(M_tail_housing, htail.xcg)
    tail         = htail + vtail + servo + tail_housing
    @addbranch inertiaLog tail(htail, vtail, servo, tail_housing)

    # Finish mass and CG estimation
    icalc  = wing + tail + fuselageInertia + tailboom

    sumall!(inertiaLog)
    return icalc
end

function build_aero_inp(v, δtail, fuselage)
    θwing = haskey(v,:θwing) ? v.θwing : 0.0

    bhtail = v.bwing * v.bhtail
    c(y) = v.cwing[1]  + (v.cwing[2] - v.cwing[1]) * 2 * y / v.bwing
    tw(y) = θwing + (v.wtwist- θwing) * 2 * y / v.bwing
    y = SA[0., fuselage.tot_width/2, bhtail/2, v.bwing/2]
    chord = c.(y)
    twist = tw.(y)
    
    idx_rt = SA[1,2,3]
    idx_tp = SA[2,3,4]
    inputwing = (;
        xrle = (v.xwing + c(0.)/4) .- chord[idx_rt]/4,
        xrte = (v.xwing + c(0.)/4) .+ 3*chord[idx_rt]/4,
        xtle = (v.xwing + c(0.)/4) .- chord[idx_tp]/4,
        xtte = (v.xwing + c(0.)/4) .+ 3*chord[idx_tp]/4,
        yrle = y[idx_rt],
        ytle = y[idx_tp],
        rooti= twist[idx_rt],
        tipi = twist[idx_tp],
        signsymcontrol = SA[1,1,-1],
        controls = SA[0., 0., 0.]
        )
    
    inputtail = (;
        xrle = SA[v.xhtail,v.xhtail],
        xrte = SA[v.xhtail+v.chtail,v.xhtail+v.chtail],
        xtle = SA[v.xhtail,v.xhtail],
        xtte = SA[v.xhtail+v.chtail,v.xhtail+v.chtail],
        yrle = SA[y[1],y[2]],
        ytle = SA[y[2],y[3]],
        rooti = SA[0., 0.],
        tipi = SA[0., 0.],
        signsymcontrol=SA[1,1],
        controls = SA[δtail, δtail]
    )

    nl = fuselage.l_nose
    inputfuse = (;
        xrle = SA[0.,v.xwing+v.cwing[1]],
        xrte = SA[v.xwing,fuselage.tot_length],
        xtle = SA[0.,v.xwing+v.cwing[1]],
        xtte = SA[v.xwing,fuselage.tot_length],
        yrle = SA[y[1],y[1]],
        ytle = SA[y[2],y[2]],
        rooti = SA[0.,0.],
        tipi = SA[0.,0.],
        signsymcontrol=SA[1,1],
        controls = SA[0.,0.]
    )
    inp = map((a,b,c)->SA[a...,b...,c...], inputwing, inputtail, inputfuse)
    return inp
end

function aeroanalysis(v, V::Real, xcg::Real, 
            α::Real, δtail::Real, 
            cache::NamedTuple, 
            options::NamedTuple, 
            name::Symbol, 
            memory::NamedTuple;
            q::Real=0.0
            )
    # get cache
    (; wres, wmesh, fmesh, Δy) = cache
    T = promote_type(typeof(V), typeof(xcg), typeof(α), typeof(δtail),)
    wresT = get_results(wres, T)
    (; Sw_fuselage, tot_width) = options.fuselage
    tot_length = options.fuselage.tot_length
    (; airfoil, panels_tot_width) = options
    nwing = sum(wres.Nelempan[i] for i=1:3)


    # Geometry
    inp = build_aero_inp(v, δtail, options.fuselage)

    bref = v.bwing
    cref = sum(v.cwing) / length(v.cwing)
    Sref = bref * cref
    
    # controls
    omega_stability = SA[0., 0., 0.]
    
    # Apply fuselage downwash on horizontal tail
    if options.fuselage_dw isa Val{true}
        VX = (x,y,sa,ca) -> ((y < options.fuselage.tot_width/2) && (x > 0.7*options.fuselage.tot_length)) ? SA[1-ca,sa,zero(ca),zero(ca),zero(ca)] : SA[zero(ca),zero(ca),zero(ca),zero(ca),zero(ca)] # 0 angle of attack behind 
        VZ = (x,y,sa,ca) -> ((y < options.fuselage.tot_width/2) && (x > 0.7*options.fuselage.tot_length)) ? SA[-sa,-ca,zero(ca),zero(ca),zero(ca)] : SA[zero(ca),zero(ca),zero(ca),zero(ca),zero(ca)] # 0 angle of attack behind 
    else
        VX = VZ = (x,y,sa,ca) -> SA[zero(sa),zero(sa),zero(sa),zero(sa),zero(sa)]
    end
    Drag_trn = addnode(memory.CD, name)
    addnode(memory.speed, name, V)
    Lift_trn = addnode(memory.CL, name)

    # First Weissinger call solve mostly to fill in wres geometry and factorize AIC (Re dependent)
    refvars = SA[xcg, Sref, cref, bref, V]
    weissinger!(wres, inp.xrle, inp.xrte, inp.xtle, inp.xtte, inp.yrle, inp.ytle, inp.rooti, inp.tipi, inp.controls,
    α, omega_stability, airfoil, inp.signsymcontrol, refvars, Δy=Δy, VX=VX, VZ=VZ);
    chord = get_array(wres, :chord, T)
    paneli = get_array(wres, :paneli, T)

    # Apply beam torsion on wing and update wing incidence
    if options.torsion isa Val{true}
        wmesh.tl .= 0.
        striptheory(momentcoefficient, wmesh.tl, V, airfoil, inp.controls, wres, rightwing=Val{true}(), elems=1:3)
        # Translate forward to tip to account for unloaded panels inside fuselage
        for i=nwing:-1:1
            wmesh.tl[i+panels_tot_width] = 1/2 * wmesh.tl[i] * (1/2 * envt.rho_atm * V^2 * chord[i]^2)
        end
        @. wmesh.tl[1:panels_tot_width] = 0.
        # Handle small difference in panels
        for i=nwing+1:-1:2
            wmesh.tl[i+panels_tot_width] += wmesh.tl[i+panels_tot_width-1]
        end
        cantilevered_beam_torsion!(wmesh, MaterialProperties.G_carbon)
        for i=1:nwing
            paneli[i] -= 1/2*(wmesh.θ[panels_tot_width+i]+wmesh.θ[panels_tot_width+i+1])*180/pi
        end
        
        # Re-solve with twisted wing
        weissinger!(wres, inp.xrle, inp.xrte, inp.xtle, inp.xtte, inp.yrle, inp.ytle, inp.rooti, inp.tipi, inp.controls,
            α, omega_stability, airfoil, inp.signsymcontrol, refvars, reuseGeom=Val{true}(), Δy=Δy, VX=VX, VZ=VZ);
    end

    # Add vertical-tail parasitic drag
    CL  = get_array(wres, :CL, T)
    CDp         = get_array(wres, :CDp, T)
    wing, htail = CDp[1][1]+CDp[2][1]+CDp[3][1], CDp[4][1]+CDp[5][1]
    vtail       = htail * 2 / 3 # assume vertical tail has an area of 2/3 of the horizontal tail and the same friction coeft (it's an upper bound)
    d_wing  = dCD_dCL(CDp[1],CL)+dCD_dCL(CDp[2],CL)+dCD_dCL(CDp[3],CL)
    d_htail = dCD_dCL(CDp[4],CL)+dCD_dCL(CDp[5],CL)
    d_vtail = 2/3*d_htail

    # Add fuselage parasitic drag
    Re_fuse      = tot_length * V / envt.nu
    CDp_fuse     = 0.074 * Re_fuse^(-0.2)*1.4 # form factor of 1.4
    fuse         = CDp_fuse * Sw_fuselage / Sref # add fuselage friction drag
    wet_tailboom_length = (fmesh.y[end] - tot_length)
    Re_tailboom  = wet_tailboom_length * (V/5) / envt.nu # assuming that incoming flow if slowed down by fuselage
    CDp_tailboom = 0.074 * Re_tailboom^(-0.2) # form factor of 1 because long and slender
    tailboom     = (wet_tailboom_length * 2π * fmesh.cs.r[1])/Sref * CDp_tailboom 

    # Full parasitic drag
    interference = (wing+htail+vtail+fuse+tailboom) * 0.04 # add 4% interference drag (see Fluid Dynamics Drag, page 8-12, Hoerner)
    CDp          = (wing+htail+vtail+fuse+tailboom+interference)
    @addbranch Drag_trn CDp(wing,htail,vtail,fuse,tailboom,interference) 
    d_interference = (d_wing+d_htail+d_vtail) * 0.04
    d_CDp          = (d_wing+d_htail+d_vtail+d_interference)
    
    # Compute full drag
    # CDi = sum(c[1] for c=get_array(wres, :CDi, T)) * 1.2 # add 20% margin for fuselage effects hurting lift distribution
    CDi = get_array(wres, :CDff, T)[1] * 1.2
    d_CDi = 0.
    # for c=get_array(wres, :CDi, T)
    #     d_CDi += dCD_dCL(c,CL)
    # end
    CD = CDi + CDp
    d_CD = d_CDi + d_CDp
    @addnode Drag_trn CDi

    # # Add fuselage lift and pitching moment
    # # tot_width = 0.07071067811865475
    # (; tot_length, tot_width, aerodynamic_center) = options.fuselage
    # AR_fuse = tot_width / tot_length
    # Sref_fuse = tot_width * tot_length
    # CLa_fuse = π*AR_fuse/2 * Sref_fuse / Sref
    # wresT.CL[1] += CLa_fuse * (α*π/180)
    # wresT.CL[2] += CLa_fuse
    # Cma_fuse = (xcg-aerodynamic_center) / cref * CLa_fuse
    # wresT.Cm[1] += Cma_fuse * (α*π/180)
    # wresT.Cm[2] += Cma_fuse
    
    # Log lift coefficients
    if typeof(memory.CL) != BuildUp{Nothing}
        CLtail, Stail = elementwise_cl(wresT, (4,5))
        CLwing, Swing = elementwise_cl(wresT, (1,2,3))
        CLfuse, Sfuse = elementwise_cl(wresT, (6,7))
        Lift_trn[:tail] = CLtail*Stail/Sref
        Lift_trn[:wing] = CLwing*Swing/Sref
        Lift_trn[:fuse] = CLfuse*Sfuse/Sref#CLa_fuse * (α*π/180)
        memory.speed[name] = V
        memory.cl[name] = wresT.cl[sum(wresT.Nelempan)+1:end,1]
        memory.cla[name] = wresT.cl[sum(wresT.Nelempan)+1:end,2]

        # for k=keys(inp)
        #     memory.weissinger_input[k] = Vector(inp[k])
        # end
    end
    return CD, d_CD
end

function trimmed_aeroanalysis(v, cltarget, V::Real, xcg::Real, 
    α::Real, δ::Real, # initial guesses for alpha and delta
    cache::NamedTuple, 
    options::NamedTuple, 
    name::Symbol, 
    memory::NamedTuple, tol=1e-8, maxsolve = 5)
    # This function assumes that the tail corresponds to the last two wing elements

    CD, d_CD = aeroanalysis(v, V, xcg, α, δ,  cache, options, name, memory)
    CL  = cache.wres.CL[1]
    Cm  = cache.wres.Cm[1]
    i = 5
    while (abs(CL-cltarget)+abs(Cm) > tol) && i>0
        i-=1
        CLα = cache.wres.CL[2]
        CLδ = cache.wres.CL[5+5] + cache.wres.CL[5+4]
        Cmα = cache.wres.Cm[2]
        Cmδ = cache.wres.Cm[5+5] + cache.wres.Cm[5+4]
        δα = (Cmδ*(cltarget - CL)+Cm*CLδ) / (CLα*Cmδ-CLδ*Cmα)
        δ -= (Cm+Cmα*δα)/Cmδ*180/pi
        α += δα*180/pi
        CD, d_CD = aeroanalysis(v, V, xcg, α, δ,  cache, options, name, memory)
        CL  = cache.wres.CL[1]
        Cm  = cache.wres.Cm[1]
    end
    i==0 && println("$(CL-cltarget), $Cm")
    
    return CD, d_CD, α, δ
end

function bendinganalysis(aeroload::Function, wmesh::StructMesh, tmesh::StructMesh, options::NamedTuple,
                        idg::NamedTuple, g::AbstractArray, idx::Int, name::Symbol, memory::NamedTuple)
    (; E_carbon, E_balsa, MaxStressCarbon, MaxStressBalsa, MaxShearCarbon) = MaterialProperties
    MaxTipDeflection = 5π / 180 # per unit half - span. 5 degrees of dihedral
    panels_tot_width = options.panels_tot_width
    nwing = npanels(wmesh)-1
    ntail = npanels(tmesh)-1

    # Wing bending
    wmesh.q .= 0.0
    for i=1:nwing-panels_tot_width
        ld = aeroload(i)
        wmesh.q[panels_tot_width+i] += ld/2
        wmesh.q[panels_tot_width+i+1] += ld/2
    end
    wmesh.q[1] *= 2
    cantilevered_beam_bending!(wmesh, E_carbon)
    g[idg.wingstress[idx]] = ((1.5 * wmesh.σyy[1] / MaxStressCarbon)^2 + (1.5 * wmesh.τxz[1] / MaxShearCarbon)^2 - 1.0) # stress
    g[idg.wingdisp[idx]] = wmesh.w[end] / (MaxTipDeflection * wmesh.y[end]) # deflection

    # Tail bending
    tmesh.q .= 0.0
    for i=1:ntail-panels_tot_width
        ld = aeroload(nwing-panels_tot_width+i)
        tmesh.q[panels_tot_width+i] += ld/2
        tmesh.q[panels_tot_width+i+1] += ld/2
    end
    tmesh.q[1] *= 2
    cantilevered_beam_bending!(tmesh, E_balsa)
    g[idg.tailstress[idx]] = tmesh.σyy[1] / MaxStressBalsa * 1.5
    g[idg.taildisp[idx]] = tmesh.w[end] / (MaxTipDeflection * tmesh.y[end]) # deflection

    # Log
    if !isa(memory.wingΔz, Nothing)
        memory.wingΔz[name] = copy(wmesh.w)
        memory.tailΔz[name] = copy(tmesh.w)
        memory.wingθ[name] = copy(wmesh.θ)
        memory.tailθ[name] = copy(tmesh.θ)
    end
end

function dCD_dCL(CD,CL)
    return sum(CD[i]/CL[i] for i=eachindex(CD)) - CD[1]/CL[1]
end

function airframeanalysis(v::NamedTuple{var}, g::Vector{TF}, idg::NamedTuple{out}, cache, options, memory=NoMemory_Airframe) where {TF,var,out}
    (; Nelempan, fuselage, airfoil, panels_tot_width) = options
    Npanhalfwing = sum(Nelempan)
    (; wres, wmesh, fmesh, tmesh, lift_data, drag_data, speed_data, Δy) = cache
    wresT = get_results(wres, TF)

    ### BASIC GEOMETRIC CONSTRAINTS ###
    Sref = sum(v.cwing) / length(v.cwing) * v.bwing

    # Choose panel spacing on wing for both aero and structures
    bhtail = v.bhtail * v.bwing
    s2, s3, s5 = sum(Δy[2]), sum(Δy[3]), sum(Δy[5])
    @. Δy[2] *= (bhtail-fuselage.tot_width)/2 / s2
    @. Δy[3] *= (v.bwing-bhtail)/2 / s3
    @. Δy[5] *= (bhtail-fuselage.tot_width)/2 / s5
    
    # Fill-in wing mesh with input
    idx = panels_tot_width+1
    for j = 1:3
        for i=eachindex(Δy[j])
            wmesh.y[idx+i] = wmesh.y[idx+i-1]+Δy[j][i]
        end
        idx += length(Δy[j])
    end
    map!(y -> v.rwing[1] + (v.rwing[2] - v.rwing[1]) * 1 / wmesh.y[end] * y, wmesh.cs.r, wmesh.y)
    map!(y -> v.twing[1] + (v.twing[2] - v.twing[1]) * 1 / wmesh.y[end] * y, wmesh.cs.t, wmesh.y)
    updatecrosssection!(wmesh)
    wmesh.q .= 0.0
    wmesh.tl .= 0.0
    
    # Fill-in tailboom mesh with input
    if haskey(v, :masspropulsion)
        fixedInertia = fuselage.fuselageInertia.value + InertialElement(v.masspropulsion, 0.0)
        memory.inertia[:fuselage] = fuselage.fuselageInertia
    else
        fixedInertia = fuselage.fuselageInertia.value + fuselage.fixedPropulsion.value
        memory.inertia[:fuselage] = fixedInertia
    end
    fmesh.y[1] = 0.0
    fmesh.y[end] = v.xhtail + v.chtail
    fmesh.cs.r .= 0.005
    fmesh.cs.t .= 0.001
    updatecrosssection!(fmesh)
    
    # Fill-in tail with w/ input
    foreach(i->tmesh.y[i] = wmesh.y[i], eachindex(tmesh.y))
    @. tmesh.cs.w_i = 0.0
    @. tmesh.cs.h_i = 0.0
    @. tmesh.cs.w_o = v.chtail
    @. tmesh.cs.h_o = tmesh.cs.w_o * 0.031
    updatecrosssection!(tmesh)
    tmesh.q .= 0.0
    tmesh.tl .= 0.0

    # Constrain wing spar to fit comfortably into wing foam core
    wtc = geom(airfoil[1]).t_c
    @. g[idg.sparradius] = (v.rwing - 0.4 * wtc * v.cwing) / 0.005

    # Constrain spar thickness to fit into spar radius
    @. g[idg.wsparthickness] = (v.twing - v.rwing) / 0.005

    # Constrain tail to be behind wing
    @. g[idg.xwtail] = v.xwing + v.cwing - v.xhtail

    # Constrain tail incidence change to remain small (would be for mechanical reasons)
    g[idg.θtailvar] = (v.θtail[2]-v.θtail[4])/10 - 1. # cruise minus take off less than 10 deg

    ### Mass Build-up ####
    inertia = mass_estimation(v, options, wmesh, fmesh, tmesh, fixedInertia, memory.inertia)
    weight = inertia.mass * envt.g
    xcg = inertia.xcg

    #### Aerodynamic calculations at every point ####
    ## Turns
    Vt, loadfactor_t = v.speed_maxload, v.maxload 
    CD, dCD = aeroanalysis(v, Vt, xcg, v.alpha[1], v.θtail[1], cache, options, :turn, memory, q=envt.g/Vt*(loadfactor_t-1/loadfactor_t))
    (; CL, Cm, chord, cl) = wresT
    aeroload(V,i,cl) = (envt.rho_atm/2*chord[i]*V^2*cl[Npanhalfwing+i])

    # Lift, trim and stability constraints
    g[idg.lift[1]] = 1.2*loadfactor_t * weight / (envt.rho_atm / 2 * Vt^2 * Sref) - CL[1]
    g[idg.Cm[1]] = Cm[1]
    c̄ = sum(v.cwing) / length(v.cwing)
    g[idg.staticmargin] = Cm[2]*c̄ + CL[2] * (0.25 * c̄ + 0.03) # 25% static margin + 3cm CG uncertainty
    lift_data[1] = envt.rho_atm / 2 * Sref * CL[1]
    drag_data[1] = envt.rho_atm / 2 * Sref * CD
    speed_data[1] = 1.0 #Vt

    # Stall constraints turning flight
    istall = idg.stall.start - 1
    ielem = 0
    for e = eachindex(Nelempan)
        for ir = 1:Nelempan[e]
            g[istall+ir+ielem] = cl[Npanhalfwing+ir+ielem] / clmaxrange(airfoil[e], 0.0, Vt * chord[ir+ielem] / envt.nu)
        end
        ielem += Nelempan[e]
    end

    # Bending analysis turing turn
    bendinganalysis(i->aeroload(Vt,i,cl), wmesh, tmesh, options, idg, g, 1, :turn, memory)
    
    # Tail torsion at CLmax
    @. tmesh.tl = tmesh.q * v.chtail / 4
    cantilevered_beam_torsion!(tmesh, MaterialProperties.G_balsa)

    ## Solve for straight line
    Vs = v.maxspeed
    CD, dCD = aeroanalysis(v, Vs, xcg, v.alpha[2], v.θtail[2], cache, options, :straight, memory)

    # Lift and trim constraints
    g[idg.lift[2]] = weight / (envt.rho_atm / 2 * Vs^2 * Sref) - CL[1]
    g[idg.Cm[2]] = Cm[1]
    lift_data[2] = envt.rho_atm / 2 * Sref * CL[1]
    drag_data[2] = envt.rho_atm / 2 * Sref * CD
    speed_data[2] = 1.0 #Vs
    
    # Bending analysis at max speed
    bendinganalysis(i->aeroload(Vs,i,cl), wmesh, tmesh, options, idg, g, 3, :straight, memory)

    ## Gust analysis in straight line
    Vg = Vs
    α_g = sharpedge_gustspeed / Vg
    for i=1:2*Npanhalfwing
        cl[i] = cl[i] + cl[i, 2] * α_g
    end

    # Bending analysis during gust
    bendinganalysis(i->aeroload(Vg,i,cl), wmesh, tmesh, options, idg, g, 2, :fastgust, memory)

    # Log coefficients
    if eltype(memory.CL)!=Nothing
        CLtail, Stail = elementwise_cl(wresT, (4, 5))
        CLwing, Swing = elementwise_cl(wresT, (1, 2, 3))
        tail = CLtail * Stail / Sref
        wing = CLwing * Swing / Sref
        @addbranch memory.CL fastgust(tail, wing)
        memory.CD[:fastgust] = NaN
        memory.speed[:fastgust] = Vg
        memory.cl[:fastgust] = cl[Npanhalfwing+1:2*Npanhalfwing]
        memory.α_gust[:fastgust] = α_g*180/pi
        memory.wingΔz[:fastgust] = copy(wmesh.w)
        memory.tailΔz[:fastgust] = copy(tmesh.w)
    end

    # Tip stall constraints in gust
    istall = idg.fastgust_stall.start - 1
    istart = sum(Nelempan[i] for i=1:2)
    for ir = 1:Nelempan[3] # outboard wing section
        iright = Npanhalfwing + istart + ir
        g[istall+ir] = cl[iright] / clmaxrange(airfoil[3], 0.0, Vg * chord[istart+ir] / envt.nu)
    end

    # Tail stall constraints in gust
    istall += Nelempan[3]
    istart = sum(Nelempan[i] for i=1:4)
    for ir = 1:Nelempan[5] # horizontal tail
        iright = Npanhalfwing + istart  + ir
        g[istall+ir] = cl[iright] / clmaxrange(airfoil[5], 0.0, Vg * chord[istart+ir] / envt.nu)
    end

    ## Additional analysis
    loadfactor_c = (1 + v.maxload)/2
    CD, dCD = aeroanalysis(v, Vt, xcg, v.alpha[3], v.θtail[3], cache, options, :extra_eval, memory)
    
    # Lift and trim constraints
    g[idg.lift[3]] = loadfactor_c * weight / (envt.rho_atm / 2 * Vt^2 * Sref) - CL[1]
    g[idg.Cm[3]]  = Cm[1]
    lift_data[3]  = envt.rho_atm / 2 * Sref * CL[1]
    drag_data[3]  = envt.rho_atm / 2 * Sref * CD
    speed_data[3] = 1.0 #Vt

    # Take Off constraints
    Vto = 15.
    CD, dCD = aeroanalysis(v, Vto, xcg, v.alpha[4], v.θtail[4], cache, options, :takeoff, memory)

    # Lift and trim constraints
    g[idg.lift[4]] = weight / (envt.rho_atm / 2 * Vto^2 * Sref) - CL[1]
    g[idg.Cm[4]]  = Cm[1]

    # Stall constraints turning takeoff flight
    istall = idg.takeoff_stall.start - 1
    ielem = 0
    for e = eachindex(Nelempan)
        for ir = 1:Nelempan[e]
            g[istall+ir+ielem] = cl[Npanhalfwing+ir+ielem] / clmaxrange(airfoil[e], 0.0, Vto * chord[ir+ielem] / envt.nu)
        end
        ielem += Nelempan[e]
    end

    # dragmodelfit!(v.dragmodel, lift_data, drag_data, speed_data, g, idg.dragmodelfit)
    dragmodel = estimate_drag_model(lift_data, drag_data, speed_data)
    mass = inertia.mass
    return (; mass, dragmodel)
end
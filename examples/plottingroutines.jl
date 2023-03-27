using ModelAirRaces
using Plots
using LaTeXStrings
using StaticArrays
using Printf
using Statistics

function showplane(wres::WeissingerResults, wmesh::StructMesh, fmesh::StructMesh, xcg::Union{Float64,Nothing}=nothing; xlims=(-0.1, 1.), ylims=(-0.7, 0.7))
    p = plot(aspectratio=:equal, xlabel="x (m)", ylabel="y (m)")

    # # Propulsion
    # width_propulsion = envt.length_propulsion/2
    # xfuse = [0., envt.length_propulsion, envt.length_propulsion, 0., 0.]
    # yfuse = [tot_width/2, tot_width/2, -tot_width/2, -tot_width/2, tot_width/2]
    # p = plot(Shape(xfuse,yfuse), aspectratio=:equal, alpha=.5,label="", xlabel="x (m)", ylabel="y (m)")

    # fuselage
    fuse = ModelAirRaces.fixedfuselage_shapes()
    fuselage = ModelAirRaces.fixedfuselage()
    for xyfuse= fuse
        p = plot!(p, Shape(xyfuse...), alpha=.5, label="")
    end

    # Wing
    nwing = sum(wres.Nelempan[1:3])
    ywing = [wres.yc[1:nwing]; reverse(wres.yc[1:nwing])]
    xwing = [wres.xc[1:nwing] - wres.chord[1:nwing]/4; reverse(wres.xc[1:nwing] + 3*wres.chord[1:nwing]/4)]
    plot!(p, Shape(xwing,ywing), alpha=.5,label="")
    plot!(p, Shape(xwing,-ywing), alpha=.5,label="",color=2)

    # horizontal tail
    ihtail = nwing+1:sum(wres.Nelempan[1:5])
    yhtail = [wres.yc[ihtail]; reverse(wres.yc[ihtail])]
    xhtail = [wres.xc[ihtail] - wres.chord[ihtail]/4; reverse(wres.xc[ihtail] + 3*wres.chord[ihtail]/4)]
    xhtail = [reverse(xhtail); xhtail]
    yhtail = [-reverse(yhtail); yhtail]
    plot!(p, Shape(xhtail,yhtail), aspectratio=:equal, alpha=.5,label="",color=3)

    # tail boom
    xboom = [fmesh.y; reverse(fmesh.y)]
    yboom = [fmesh.cs.r; -reverse(fmesh.cs.r)]
    plot!(p, Shape(xboom,yboom), aspectratio=:equal, alpha=.5, color=:black,label="")
    
    # wing spar
    yspar = [wmesh.y; reverse(wmesh.y)]
    xspar = wres.xr[1] .+ [wmesh.cs.r; -reverse(wmesh.cs.r)]
    yspar = [-reverse(yspar); yspar]
    xspar = [reverse(xspar); xspar]
    plot!(p, Shape(xspar,yspar), aspectratio=:equal, alpha=.5, color=:black,label="", ylims=ylims, xlims=xlims)

    # CG
    if !isa(xcg, Nothing)
        vline!(p, [xcg], color=:gray, label="CG location")
    end
    p
end

function showpolar(a::Airplane,S::Real)
    CL = LinRange(0., 0.9, 100)
    Q = envt.rho_atm / 2 * S
    cd2 = a.dragmodel[3] *Q
    cl_cdmin = - a.dragmodel[2] / 2 / cd2 
    cdmin = a.dragmodel[1]/Q - cd2*cl_cdmin^2
    CD = map(cl->(cdmin+cd2*(cl-cl_cdmin)^2), CL)
    # plot(CL,CD, ylims=(0.,0.1), label="Drag Model", xlabel=L"C_L", ylabel=L"C_D")
    plot(CL,CD, label="Drag Model", xlabel=L"C_L", ylabel=L"C_D")
end

function showplanform!(p, 
    xrle,
    xrte,
    xtle,
    xtte,
    yrle,
    ytle;
    horizontal::Bool=true,
    fuselage=false
    )
    for i=eachindex(xrle)
            x = SA[xrle[i], xtle[i], xtte[i], xrte[i]]
            y = [yrle[i], ytle[i], ytle[i], yrle[i]]
        if horizontal
            plot!(p, x, y, label="", aspectratio=:equal, linewidth=2)
        else
            plot!(p, y, -x, label="", aspectratio=:equal, linewidth=2)
        end
    end

    if fuselage==true
        fuse = ModelAirRaces.fixedfuselage_shapes()
        for xyfuse=fuse
            p = plot!(p, Shape(xyfuse...), alpha=.2, label="")
        end
    end

    return p
end

function showpanels!(p, wres; horizontal::Bool=true)
    i = 1
    xmax = maximum(wres.xctl + wres.chord)
    xmax *= 1.15
    for ielem=eachindex(wres.Nelempan)
    for iy = 1:wres.Nelempan[ielem]
    y = SA[wres.yr[i], wres.yr[i], wres.yt[i], wres.yt[i]]
    x = SA[xmax, wres.xc[i], wres.xc[i], xmax]
    if horizontal
        plot!(p,x,y,label="", alpha = 0.4, linewidth=.5, color=ielem, aspectratio=1)
    else
        plot!(p,y,-x,label="", alpha = 0.4, linewidth=.5, color=ielem, aspectratio=1)
    end
    i += 1
    end
    end
    if horizontal
        xlims!(p,0.,xmax)
        ylims!(p,0.,maximum(wres.yt)*1.02)
    else
        ylims!(p,-xmax,0.)
        xlims!(p,0.,maximum(wres.yt)*1.02)
    end
    return p
end

function optreport(v, memory, dirname, P; showtitles=true)
    if !(P isa NamedTuple)
        cache = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.usercache,P.usercache,P.usercache))
        options = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.options,P.options,P.options))
        idg = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.idg,P.idg,P.idg))
        g = NamedTuple{(:airframe, :simpletraj, :propulsion)}((zeros(len(P.output)), zeros(len(P.output)), zeros(len(P.output))))
    else
        cache = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.airframe.usercache, P.simpletraj.usercache,P.propulsion.usercache))
        options = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.airframe.options, P.simpletraj.options,P.propulsion.options))
        idg = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.airframe.idg, P.simpletraj.idg,P.propulsion.idg))
        g = NamedTuple{(:airframe, :simpletraj, :propulsion)}((P.airframe.g, P.simpletraj.g,P.propulsion.g))
    end
    (; CD, CL, speed, inertia, traj) = memory
    (wres, wmesh, fmesh, tmesh, lift_data, drag_data, speed_data, Δy) = cache.airframe;

    out = AD.airframeanalysis(v, g.airframe, idg.airframe, cache.airframe, options.airframe, memory)
    airplane = Airplane(inertia.value.mass, SVector{3}(out.dragmodel), SA[v.speed_maxload, v.maxload], SA[v.maxspeed, 1.0])
    out = AD.propulsionanalysis(v, g.propulsion, idg.propulsion, cache.propulsion, options.propulsion, memory)
    propulsion = AD.MotorPropeller(out.maxthrustmodel, out.pwrmodel, out.maxpower)
    motor = AD.TMotor(Kv=v.Kv, R=v.R, Pmax=out.maxpower)
    propeller = AD.APC(v.propdiam, v.proppitch*v.propdiam)

    sumall!(memory.CD)
    sumall!(memory.CL)
    sumall!(memory.speed)
    sumall!(memory.inertia)
    sumall!(memory.traj)

    bhtail = v.bhtail * v.bwing
    Sratio = v.chtail * bhtail / (mean(v.cwing) * v.bwing)
    S = v.bwing * mean(v.cwing)
    mass = innertree(inertia, :mass)
    xcg = innertree(inertia, :xcg)
    lift = (envt.rho_atm * S / 2) * speed^2 * CL
    drag = (envt.rho_atm * S / 2) * speed^2 * CD
    loadfactor = lift / (inertia.value.mass * envt.g)
    dragmodel = AD.estimate_drag_model(lift_data, drag_data, speed_data);

    trajint, transition_idx = AD.traj_integral(airplane, propulsion, v.speed_maxload, v.dvdt, v.dh*options.simpletraj.leg_length/100, traj);
    endtime = map(i->trajint.time[i], transition_idx)

    propint = (; torque = copy(trajint.time), volt = copy(trajint.time), RPM = copy(trajint.time), 
    power = copy(trajint.time), eff=copy(trajint.time), effprp=copy(trajint.time), 
    effmtr=copy(trajint.time), I = copy(trajint.time), J = copy(trajint.time))
    for i=eachindex(trajint.time)
        pwr, vlt, I, Q, ω, J = AD.power(motor, propeller, trajint.thrust[i],trajint.speed[i])
        propint.power[i] = pwr
        propint.torque[i] = Q
        propint.I[i] = I
        propint.volt[i] = vlt
        propint.RPM[i] = ω * 60 / 2π
        propint.J[i] = J
        propint.effmtr[i] = Q*ω / pwr
        propint.effprp[i] = trajint.thrust[i]*trajint.speed[i] / (Q*ω)
        propint.eff[i] = propint.effmtr[i]*propint.effprp[i]
    end

    ### Trajectory ###
    trajdirname = joinpath(dirname,"traj")
    mkpath(trajdirname)
    pl = plot(trajint.time, trajint.speed, label="", ylims=(30, 40), xlabel="time (s)", 
                ylabel="Speed (m/s)",
                )
    showtitles && title!("Airspeed During Straight Segments and a Turn")
    spd = map(i->trajint.speed[i], transition_idx)
    scatter!(endtime, spd, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 0.5+spd[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 0.5+spd[i],text("$i", :green, :bottom, 10))
        end
    end
    savefig("$trajdirname/speed_v_time.pdf")
    
    ## Altitude
    pl = plot(trajint.time, trajint.z .+ 1., label="", ylims=(0., 35), xlabel="time (s)", ylabel="Altitude (m)",
    )
    showtitles && title!("Altitude During Straight Segments and a Turn")
    z = map(i->trajint.z[i], transition_idx)
    scatter!(endtime, z, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 1.05*z[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 1.05*z[i],text("$i", :green, :bottom, 10))
        end
    end
    savefig("$trajdirname/alt_v_time.pdf")
    
    ## Load factor
    pl = plot(trajint.time, trajint.n_z, label="", ylims=(0.8, 5.0), xlabel="time (s)", ylabel=L"n_z",
    )
    showtitles && title!("Load Factor During Straight Segments and a Turn")
    nz = map(i->trajint.n_z[i], transition_idx)
    scatter!(endtime, nz, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 1.05*nz[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 1.05*nz[i],text("$i", :green, :bottom, 10))
        end
    end
    savefig("$trajdirname/nz_v_time.pdf")
    
    ## Top view
    x,y,dist, xpylons, ypylons, xlabels, ylabels = AD.showsimpletraj(v.δ/100*envt.marathonlength/120, traj)
    plot(y,x,aspectratio=1,label="Optimized Path")
    scatter!(ypylons, xpylons,label="Pylons")
    scatter!(ylabels, xlabels,label="Segments")
    for i=eachindex(xlabels)
        if i<length(xlabels)
            annotate!(ylabels[i], -2+xlabels[i],text("$i", :green, :top, :left, 10))
        else
            annotate!(ylabels[i], -2+xlabels[i],text("$i", :green, :top, :right, 10))
        end
    end
    xlabel!("x (m)")
    ylabel!("y (m)")
    savefig("$trajdirname/trajtopview.pdf")
    
    ## Thrust
    pl = plot(trajint.time, trajint.thrust, label="", xlabel="time (s)", ylabel=L"T (N)")
    th = map(i->trajint.thrust[i], transition_idx)
    scatter!(endtime, th, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 1.05*th[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 1.05*th[i],text("$i", :green, :bottom, 10))
        end
    end
    showtitles && title!("Thrust During Straight Segments and a Turn")
    savefig("$trajdirname/thrust_v_time.pdf")

    ### Propulsion
    propdirname = joinpath(dirname,"prop")
    mkpath(propdirname)

    ## Propeller efficiency
    J = LinRange(0., 1., 100)
    CT = map(j->AD.ct(propeller,j),J)
    CP = map(j->AD.cp(propeller,j),J)
    plot(J, CT, label="", xlabel="J", ylabel=L"C_T")
    showtitles && title!("Propeller Thrust Coefficient vs. Advance Ratio")
    savefig("$propdirname/ct.pdf")
    plot(J, CP, label="", xlabel="J", ylabel=L"C_P")
    showtitles && title!("Propeller Power Coefficient vs. Advance Ratio")
    savefig("$propdirname/cp.pdf")
    plot(J, CT.*J./CP, label="", xlabel="J", ylabel=L"\eta", ylims=(0.,0.8))
    showtitles && title!("Propeller Efficiency vs. Advance Ratio")
    savefig("$propdirname/etaprp.pdf")

    ## Efficiency contour with max thrust
    spd = LinRange(20., 40., 50)
    eff = (V,T) -> T*V / AD.power(motor, propeller, T, V)[1]
    thr = LinRange(1e-3, out.maxthrustmodel[1] + out.maxthrustmodel[2]*spd[1], 50)
    contourf(spd, thr, eff, colorbar_title="η", 
        ylabel="Thrust (N)", xlabel="Airspeed (m/s)")
    plot!(spd, s->out.maxthrustmodel[1]+s*out.maxthrustmodel[2], color=:red, label="Max. Thrust")
    savefig("$propdirname/etacontour.pdf")
    plot!(trajint.speed, trajint.thrust, color=:black, label="Trajectory")
    savefig("$propdirname/etacontour_wtraj.pdf")

    ## Efficiency contour with max thrust
    spd = LinRange(minimum(cache.propulsion.spd), maximum(cache.propulsion.spd), 50)
    thr = [AD.thrustcalc(motor, propeller, 3*3.0, s)[1] for s=spd] # using smallest voltage encountered
    plot(spd, s->out.maxthrustmodel[1]+s*out.maxthrustmodel[2], label="Model Used by Trajectory Opt.")
    plot!(spd, thr,  ylabel="Max. Thrust (N)", xlabel="Airspeed (m/s)", label="Max. Thrust Computed by Propulsion Opt.")
    # plot!(spd, s->out.maxthrustmodel[1]+s*out.maxthrustmodel[2], color=:red, label="Max. Thrust")
    savefig("$propdirname/tmax.pdf")
    
    
    ## Efficiency contour with max thrust
    eff2 = (V,T) -> T*V / AD.power(propulsion, T, V)[1]
    vmid = (v.speed_maxload + v.maxspeed)/2
    v1 = vmid - options.propulsion.speedbounds
    v2 = vmid + options.propulsion.speedbounds
    spd = LinRange(v1,v2, 50)
    thr = LinRange(1e-3, out.maxthrustmodel[1] + out.maxthrustmodel[2]*spd[1], 50)
    contour(spd, thr, eff, colorbar_title="η", 
            ylabel="Thrust (N)", xlabel="Airspeed (m/s)", label="Propulsive Efficiency")
    contour!(spd, thr, eff2, colorbar_title="η", 
            ylabel="Thrust (N)", xlabel="Airspeed (m/s)", label="Efficiency Model Used by Traj.Opt", linestyle=:dash)
    savefig("$propdirname/etacomparison.pdf")


    ## Evaluate along trajectory
    
    ## Power 
    # trajpower = map((t,v)->AD.power(motor, propeller, t,v)[1], trajint.thrust, trajint.speed)
    # proppower = map((t,v)->AD.power(propulsion, t,v)[1], trajint.thrust, trajint.speed)
    pl = plot(trajint.time, trajint.dEdt, label="Model Used by Trajectory Opt.", xlabel="time (s)", ylabel="Power (W)",)
    plot!(pl, trajint.time, propint.power, label="Power Computed by Propulsion Opt.", xlabel="time (s)",
    color=3, legend=:bottomright)
    showtitles && title!("Power During Straight Segments and a Turn")
    dEdt = map(i->trajint.dEdt[i], transition_idx)
    scatter!(endtime, dEdt, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 1.02*dEdt[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 1.02*dEdt[i],text("$i", :green, :bottom, 10))
        end
    end
    savefig("$propdirname/power.pdf")
    
    ## Efficiency
    pl = plot(trajint.time, propint.eff, label="Total", xlabel="time (s)",ylabel="Efficiency")
    showtitles && title!("Efficiency During Straight Segments and a Turn")
    eff = map(i->propint.eff[i], transition_idx)
    scatter!(endtime, eff, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 1.05*eff[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 1.05*eff[i],text("$i", :green, :bottom, 10))
        end
    end
    plot!(trajint.time, propint.effmtr, label="Motor")
    plot!(trajint.time, propint.effprp, label="Propeller")
    savefig("$propdirname/efficiency.pdf")
    

    ### Airframe
    airdirname = joinpath(dirname,"air")
    mkpath(airdirname)
    
    ## Drag
    drg = similar(trajint.time)
    α = v.alpha[2]
    δ = v.θtail[2]
    for idx=eachindex(drg)
        cltarget = trajint.n_z[idx] * airplane.mass * envt.g / (envt.rho_atm/2*trajint.speed[idx]^2*S)
        drg[idx], _, α, δ = AD.trimmed_aeroanalysis(v, cltarget, trajint.speed[idx], inertia.value.xcg, 
                            α, δ,  cache.airframe, options.airframe, :none, ModelAirRaces.NoMemory_Airframe)
        drg[idx] = drg[idx] * (envt.rho_atm/2*trajint.speed[idx]^2*S)
    end
    pl = plot(trajint.time, trajint.drg, label="Model Used by Trajectory Opt.", xlabel="time (s)",ylabel="Drag (N)",legend=:topleft)
    plot!(pl, trajint.time, drg, label="Drag Computed by Airframe Opt.", xlabel="time (s)")
    showtitles && title!("Drag During Straight Segments and a Turn")
    drg = map(i->trajint.drg[i], transition_idx)
    scatter!(endtime, drg, label="", color=3)
    for i=eachindex(endtime)
        if i<length(endtime)
            annotate!(endtime[i], 1.05*drg[i],text("$i", :green, :bottom, 10))
        else
            annotate!(endtime[i], 1.05*drg[i],text("$i", :green, :bottom, 10))
        end
    end
    savefig("$airdirname/drag.pdf")

    ## Top view
    yl = (0., v.bwing / 2 * 1.05)
    xl = (-0.05, v.xhtail*1.1)
    pl = showplane(wres, wmesh, fmesh, xlims=xl, ylims=yl)
    vline!(pl, [inertia.value.xcg], label="center of gravity")
    savefig(pl,"$airdirname/circleplane.pdf")

    ##
    ilifting = 1:sum(wres.Nelempan[1:5])
    cl = memory.cl
    pl = scatter(wres.yc[ilifting], cl[:straight][ilifting],legend=:bottomright, label=(@sprintf "Max. Speed (1g, %.1fm/s)" v.maxspeed))
    scatter!(pl, wres.yc[ilifting], cl[:turn][ilifting],label=(@sprintf "Max. Load (1.2x%.1fg, %.1fm/s)" (v.maxload) v.speed_maxload))
    scatter!(pl, wres.yc[ilifting], cl[:fastgust][ilifting],label=(@sprintf "Sharp-Edge Gust (%.1fg, %.1fm/s)" (loadfactor[:fastgust].value) speed[:fastgust].value))
    scatter!(pl, wres.yc[ilifting], cl[:takeoff][ilifting],label="Takeoff (1g, 15m/s)")
    # scatter!(pl, wres.yc, cl[:extra_eval],label="")
    showtitles && title!(pl,"Lift Coefficient Distribution Along Semi-Span",)
    ylabel!(pl,L"c_l",)
    xlabel!(pl,L"{2y}/{b_{\sf wing}} \; (m)", legend=:bottomright)
    savefig(pl,"$airdirname/cl.pdf")
    ##
    
    cl = memory.cl
    pl = scatter(wres.yc[ilifting], envt.rho_atm/2*S*speed[:straight].value^2*cl[:straight][ilifting],label=(@sprintf "Max. Speed (1g, %.1fm/s)" v.maxspeed))
    scatter!(pl,wres.yc[ilifting], envt.rho_atm/2*S*speed[:turn].value^2*cl[:turn][ilifting],label=(@sprintf "Max. Load (1.2x%.1fg, %.1fm/s)" (v.maxload) v.speed_maxload))
    # scatter!(pl,wres.yc[ilifting], envt.rho_atm/2*S*speed[:extra_eval].value^2*cl[:extra_eval][ilifting],label="extra polar point")
    scatter!(pl,wres.yc[ilifting], envt.rho_atm/2*S*speed[:fastgust].value^2*cl[:fastgust][ilifting],label=(@sprintf "Sharp-Edge Gust (%.1fg, %.1fm/s)" (loadfactor[:fastgust].value) speed[:fastgust].value))
    scatter!(pl, wres.yc[ilifting], envt.rho_atm/2*S*speed[:takeoff].value^2*cl[:takeoff][ilifting],label="Takeoff (1g, 15m/s)")
    showtitles && title!(pl,"Lift Distribution Along Semi-Span",)
    ylabel!(pl,L"l \; (N/m)",)
    xlabel!(pl,L"{2y}/{b_{\sf wing}} \; (m)", legend=:bottomright)
    savefig("$airdirname/lift.pdf")
    ##
    wingΔz = memory.wingΔz
    pl = scatter(wmesh.y,-wingΔz[:straight],label=(@sprintf "Max. Speed (1g, %.1fm/s)" v.maxspeed),aspectratio=1,alpha=0.8)
    scatter!(pl, wmesh.y,-wingΔz[:turn],label=(@sprintf "Max. Load (1.2x%.1fg, %.1fm/s)" (v.maxload) v.speed_maxload),alpha=0.8)
    scatter!(pl, wmesh.y,-wingΔz[:fastgust],label=(@sprintf "Sharp-Edge Gust (%.1fg, %.1fm/s)" (loadfactor[:fastgust].value) speed[:fastgust].value))
    showtitles && title!(pl,"Wing Deflection",)
    ylabel!(pl,L"Δz",)
    xlabel!(pl,L"{2y}/{b_{\sf wing}} \; (m)")
    savefig(pl,"$airdirname/wingdz.pdf")
    ##
    tailΔz = memory.tailΔz
    pl = scatter(tmesh.y,-tailΔz[:straight],label=(@sprintf "Max. Speed (1g, %.1fm/s)" v.maxspeed),aspectratio=1)
    scatter!(pl, tmesh.y,-tailΔz[:turn],label=(@sprintf "Max. Load (1.2x%.1fg, %.1fm/s)" (v.maxload) v.speed_maxload))
    scatter!(pl, tmesh.y,-tailΔz[:fastgust],label=(@sprintf "Sharp-Edge Gust (%.1fg, %.1fm/s)" (loadfactor[:fastgust].value) speed[:fastgust].value))
    showtitles && title!(pl,"Tail Deflection",)
    ylabel!(pl,L"Δz",)
    xlabel!(pl,L"{2y}/{b_{\sf tail}} \; (m)")
    savefig(pl,"$airdirname/taildz.pdf")
    ##
    θlims = (-10., 2.)
    geo = map(y->v.wtwist*y*2/v.bwing, wmesh.y)
    p = scatter(wmesh.y*2/v.bwing, geo-memory.wingθ[:straight]*180/pi, label=(@sprintf "Max. Speed (1g, %.1fm/s)" v.maxspeed))
    scatter!(wmesh.y*2/v.bwing, geo-memory.wingθ[:turn]*180/pi, label=(@sprintf "Max. Load (%.1fg, %.1fm/s)" (1.2*v.maxload) v.speed_maxload),
    ylabel="Twist Angle (deg)",
    xlabel=L"{2y}/{b_{\sf wing}} \; (m)", ylims=θlims)
    showtitles && title!("Wing Twist")
    scatter!(wmesh.y*2/v.bwing, geo-memory.wingθ[:fastgust]*180/pi, label=(@sprintf "Sharp-Edge Gust (%.1fg, %.1fm/s)" (loadfactor[:fastgust].value) speed[:fastgust].value))
    savefig("$airdirname/torsionwing.pdf")
    p
    ##
    dragmodel = AD.estimate_drag_model(lift_data, drag_data, speed_data)
    airplane = Airplane(inertia.value.mass, dragmodel, SA[v.speed_maxload, v.maxload], SA[v.maxspeed, 1.0])
    pl = showpolar(airplane, S)
    scatter!(pl, [CL[:straight].value], [CD[:straight].value], label=(@sprintf "Max. Speed (1g, %.1fm/s)" v.maxspeed))
    scatter!(pl, [CL[:turn].value], [CD[:turn].value], label=(@sprintf "Max. Load (1.2x%.1fg, %.1fm/s)" (v.maxload) v.speed_maxload))
    scatter!(pl, [CL[:extra_eval].value], [CD[:extra_eval].value], label=(@sprintf "Mid-Domain (%.1fg, %.1fm/s)" (1+v.maxload)/2 speed[:extra_eval].value),
            legend=:bottomright)
    savefig(pl,"$airdirname/dragpolar.pdf")

    ##
    p = plot(wmesh.y*2/v.bwing, -memory.wingθ[:straight]*180/pi, label="torsion",
    ylabel="Twist Angle (deg)",
    xlabel=L"{2y}/{b_{\sf wing}} \; (m)",
    linewidth=2., ylims=θlims)
    geo = map(y->v.wtwist*y*2/v.bwing, wmesh.y)
    showtitles && title!("Wing Twist From Various Sources \n at Max. Load")
    plot!(p,wmesh.y*2/v.bwing, geo, label="jig", linewidth=2)
    plot!(p,wmesh.y*2/v.bwing, geo-memory.wingθ[:straight]*180/pi, label="combined", linewidth=2)
    savefig("$airdirname/torsionstraight.pdf")
    ##
    p = plot(wmesh.y*2/v.bwing, -memory.wingθ[:turn]*180/pi, label="torsion",
    ylabel="Twist Angle (deg)",
    xlabel=L"{2y}/{b_{\sf wing}} \; (m)",
    linewidth=2., ylims=θlims)
    geo = map(y->v.wtwist*y*2/v.bwing, wmesh.y)
    showtitles && title!("Wing Twist From Various Sources \n in Turning Flight")
    plot!(p,wmesh.y*2/v.bwing, geo, label="jig", linewidth=2)
    plot!(p,wmesh.y*2/v.bwing, geo-memory.wingθ[:turn]*180/pi, label="combined", linewidth=2)
    savefig("$airdirname/torsionturn.pdf")
    ##
    p = plot(tmesh.y*2/v.bwing, -memory.tailθ[:turn]*180/pi,
    ylabel="Twist Angle (deg)",
    xlabel=L"{2y}/{b_{\sf wing}} \; (m)",
    linewidth=2., ylims=θlims)
    showtitles && title!("Tail Twist in Turning Flight")
    savefig("$airdirname/torsiontailturn.pdf")
    ##
    inp = AD.build_aero_inp(v, 0., options.airframe.fuselage)
    p = plot()
    showplanform!(p, inp.xrle, inp.xrte, inp.xtle, inp.xtte, inp.yrle, inp.ytle,horizontal=true)
    showpanels!(p, wres,horizontal=true)
    xlabel!("y (m)")
    ylabel!("x (m)")
    # title!("Vortex System at Convergence")
    savefig("$airdirname/vortexsystem.pdf")

    return trajint, propint
end
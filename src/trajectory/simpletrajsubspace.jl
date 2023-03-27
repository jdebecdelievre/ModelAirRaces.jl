struct SimpleTraj <: OptimUtils.AbstractProblem end

makememory(T, ::SimpleTraj) = (; traj=BuildUp(:all, TrajElement(t=zero(T))))
makememory(::Nothing, ::SimpleTraj) = NoBuildUp

# Variables and constraints
function OptimUtils.ProblemCache(prb::SimpleTraj, sharedvariables; variables=(;), output=(;), options=(;))
    
    # Options
    default_options = (
        nstraight = 4, # number of straight segments in  trajopt,
        imaxspeed = 2, # segment whose end is at max speed. Typically nstraight / 2
        speedbounds=3.0,
        leg_length = envt.marathonlength/120,
        href = 10.0,
    )
    options = merge(default_options, options)
    
    # Variables
    nstraight = options.nstraight
    default_variables = (
        dvdt=Var(ini=0.0, lb=-3.0, ub=3.0, N=nstraight), # acceleration during each segment
        δ=Var(ini=25.0, lb=5.0, ub=80.0, N=nstraight+1),
        dh=Var(ini=0., lb=-5.0, ub=5.0, N=nstraight), # percentage of leg_length
    
        # Shared
        mass=sharedvariables.mass,
        time=sharedvariables.time,
        dragmodel=sharedvariables.dragmodel,
        pwrmodel=sharedvariables.pwrmodel,
        maxthrust=sharedvariables.maxthrust,
        maxload=sharedvariables.maxload,
        maxspeed=sharedvariables.maxspeed,
        speed_maxload=sharedvariables.speed_maxload,
        )
    variables = merge(default_variables, variables)
    
    # Output
    default_output = (
        sumδ = Var(lb=-Inf, ub=0.),
        MaxLoad = Var(lb=-Inf, ub=0., N=2),
        turnthrust = Var(lb=-Inf, ub=0., N=2),
        minthrust = Var(lb=-Inf, ub=0., N=2*nstraight),
        maxthrust = Var(lb=-Inf, ub=0., N=2*nstraight),
        minΔE = Var(lb=-Inf, ub=0., N=nstraight),
        Ebat = Var(lb=-Inf, ub=0.),
        MaxSpeed=Var(lb=-Inf, ub=0., N=nstraight+1),
        dvdt=Var(lb=0., ub=0., N=1),
        dh=Var(lb=0., ub=0., N=1),
        glidepath=Var(lb=-Inf, ub=0.0, N=2*nstraight),
        slowdownpath=Var(lb=-Inf, ub=0.0, N=nstraight),
        safety_dv=Var(lb=-Inf, ub=0.0, N=nstraight),

        dT=Var(lb=-Inf, ub=0.),
        # dTdV=Var(lb=-Inf, ub=0.),

        # Shared
        time=Var(lb=-Inf, ub=0.),
    )
    output = merge(default_output, output)

    OptimUtils.ProblemCache(prb, variables, output, options)
end

function OptimUtils.analysis(p::OptimUtils.ProblemCache{SimpleTraj}, g::Vector{TF}, x::Vector{TF}, memory=(;traj=NoBuildUp)) where {TF}
    for i=eachindex(x)
        x[i] = x[i] * (p.x_ub[i] - p.x_lb[i]) + p.x_lb[i]
    end
    v = unpack(x, p.idx)
    g .= 0.0

    idg = p.idg

    # Create propulsion
    T1, T2 = v.maxthrust
    vmid = (v.maxspeed+v.speed_maxload)/2
    dV = p.options.speedbounds
    v1 = vmid - dV
    v2 = vmid + dV
    dTmdx = (T2-T1) / (v2-v1)
    Tm0 = T1 - dTmdx*v1
    maxthrustmodel = SA[Tm0, dTmdx]
    propulsion = MotorPropeller(maxthrustmodel, SVector{6,TF}(v.pwrmodel), TF(0.0))
    g[idg.dT] = T2-T1
    
    # Create Airplane
    airplane = Airplane(v.mass, SVector{3,TF}(v.dragmodel), SA[v.speed_maxload, v.maxload], SA[v.maxspeed, 1.0])
    
    ## Run trajectory analysis
    out = simpletraj(airplane, propulsion, p.options, v, g, idg, memory.traj)
    
    ## Scale output
    g[idg.time] = out.time-v.time

    sumall!(memory.traj)
    for i=eachindex(x)
        x[i] = (x[i] - p.x_lb[i]) / (p.x_ub[i] - p.x_lb[i])
    end
    return g
end

function subspace_objective(p::OptimUtils.ProblemCache{SimpleTraj}, df::AbstractVector{T}, x::Vector{T}, z::Vector{T}, idz::NamedTuple) where T
    (; idx, idg) = p
    g = p.fdresults.value
    dg = p.fdresults.derivs[1]
    df .= 0.

    # Input
    df[idx.maxload] = x[idx.maxload] - z[idz.maxload]
    df[idx.speed_maxload] = x[idx.speed_maxload] - z[idz.speed_maxload]
    df[idx.maxspeed] = x[idx.maxspeed] - z[idz.maxspeed]
    df[idx.mass] = x[idx.mass] - z[idz.mass]
    df[idx.time] = x[idx.time] - z[idz.time]
    for (ix,iz) = zip(idx.dragmodel, idz.dragmodel)
        df[ix] = x[ix] - z[iz]
    end
    for (ix,iz) = zip(idx.pwrmodel, idz.pwrmodel)
        df[ix] = x[ix] - z[iz]
    end
    for (ix,iz) = zip(idx.maxthrust, idz.maxthrust)
        df[ix] = x[ix] - z[iz]
    end
    f = sum(f_^2/2 for f_=df)
    return f
end
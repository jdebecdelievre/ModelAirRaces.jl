struct APC_TMOTOR <: OptimUtils.AbstractProblem end

function OptimUtils.makecache(T::DataType, ::APC_TMOTOR, options)
    thr, spd, pwr, rot, vlt = (zeros(T, 6) for _=1:5)
    # X = ones(T,6,6)
    return (; thr, spd, pwr, rot, vlt)
end

makememory(T::DataType, ::APC_TMOTOR) = (; inertia=BuildUp(:all, InertialElement(zero(T))))
makememory(::Nothing, ::APC_TMOTOR) = (; inertia=NoBuildUp)

# Variables and constraints
function OptimUtils.ProblemCache(prb::APC_TMOTOR, sharedvariables; variables=(;), output=(;), options=(;))
    
    # Options
    default_options = (
        leg_length = envt.marathonlength/120,
        sharedvariables=sharedvariables,
        pwrmodelknots = scaleLHC(LHCoptim(30, 2, 1)[1],[(0.,1.), (0.,1.)]),
        speedbounds=3.0,
    )
    options = merge(default_options, options)
    
    # Variables
    default_variables = (;
        Kv=Var(lb=800., ub=2200., ini=1000.),
        R=Var(lb=0.01, ub=1.0, ini=0.2),
        propdiam=Var(lb=0.15, ub=0.5, ini=10.0*0.0254),
        proppitch=Var(lb=0.5, ub=1., ini=0.75),

        masspropulsion=sharedvariables.masspropulsion,
        maxspeed=sharedvariables.maxspeed,
        speed_maxload=sharedvariables.speed_maxload,
        maxthrust=sharedvariables.maxthrust,
        pwrmodel=sharedvariables.pwrmodel,
    )
    variables = merge(default_variables, variables)    

    # Output
    default_output = (;
        masspropulsion=Var(lb=-Inf, ub=0.),
        maxrot=Var(lb=-Inf, ub=0., N=6),
        maxthrustmodel = Var(lb=0., ub=0., N=2),
        pwrmodelfit = Var(lb=-Inf, ub=0., N=6),#N=size(options.pwrmodelknots,1)),
        dT=Var(lb=-Inf, ub=0.),
        dV=Var(lb=-Inf, ub=0.)
    )
    output = merge(default_output, output)

    OptimUtils.ProblemCache(prb, variables, output, options)
end

function OptimUtils.analysis(p::OptimUtils.ProblemCache{APC_TMOTOR}, g::Vector{TF}, x::Vector{TF}, memory=makememory(nothing, APC_TMOTOR())) where {TF}
    cache = (TF<:ForwardDiff.Dual) ? p.usercache_dual : p.usercache
    for i=eachindex(x)
        x[i] = x[i] * (p.x_ub[i] - p.x_lb[i]) + p.x_lb[i]
    end
    v = unpack(x, p.idx)
    g .= 0.0
    idg = p.idg

    ## Propulsion Analysis
    # out = propulsionanalysis(v, g, idg, cache, p.options, memory)
    # propulsion = MotorPropeller(out.maxthrustmodel, out.pwrmodel, out.maxpower)
    # propulsion_v = MotorPropeller(out.maxthrustmodel, SVector{6,TF}(v.pwrmodel), out.maxpower)
    # Tmax = maximum(cache.thr)
    # Tmin = minimum(cache.thr)
    # Vmax = maximum(cache.spd)
    # Vmin = minimum(cache.spd)
    # for (i,l)=enumerate(eachrow(p.options.pwrmodelknots)) # shared pwrmodel above computed pwrmodel at multiple locations
    #     V = Vmin + (Vmax-Vmin) * l[1]
    #     T = Tmin + (Tmax-Tmin) * l[2]
    #     g[idg.pwrmodelfit[i]] = (power(propulsion, T, V) - power(propulsion_v, T, V)) / Pref
    # end
    out = propulsionanalysis(v, g, idg, cache, p.options, memory)
    for i=eachindex(idg.pwrmodelfit)
        g[idg.pwrmodelfit[i]] = (out.pwrmodel[i]-v.pwrmodel[i])
    end

    for i=eachindex(x)
        x[i] = (x[i] - p.x_lb[i]) / (p.x_ub[i] - p.x_lb[i])
    end
    return g
end

function subspace_objective(p::OptimUtils.ProblemCache{APC_TMOTOR}, df::AbstractVector{T}, x::Vector{T}, z::Vector{T}, idz::NamedTuple) where T
    (; idx, idg) = p
    g = p.fdresults.value
    dg = p.fdresults.derivs[1]
    df .= 0.0

    # Inputs
    df[idx.speed_maxload] = x[idx.speed_maxload] - z[idz.speed_maxload]
    df[idx.maxspeed] = x[idx.maxspeed] - z[idz.maxspeed]
    df[idx.masspropulsion] = x[idx.masspropulsion] - z[idz.masspropulsion]
    for (ix,iz) = zip(idx.maxthrust,idz.maxthrust)
        df[ix] = x[ix] - z[iz]
    end
    for (ix,iz) = zip(idx.pwrmodel,idz.pwrmodel)
        df[ix] = x[ix] - z[iz]
    end
    f = sum(f_^2/2 for f_=df)

    return f    
end

struct Airframe <: OptimUtils.AbstractProblem end

function OptimUtils.makecache(T::DataType, ::Airframe, options)
    (; nVfit, Nelempan, panels_tot_width) = options
    tot_width = options.fuselage.tot_width
    if T <: ForwardDiff.Dual
        wres = WeissingerResultsAD(Float64, T, Nelempan)
    else
        wres = WeissingerResults(T, Nelempan)
    end
    wmesh = StructMesh(T, panels_tot_width + sum(Nelempan[1:3]) + 1)
    if panels_tot_width > 0
        wmesh.y[1:1+panels_tot_width] = LinRange(0., tot_width/2, panels_tot_width+1)
    end
    tmesh = StructMesh(T, panels_tot_width + Nelempan[1] + Nelempan[2] + 1, RectangularCrossSection(T, panels_tot_width + Nelempan[1] + Nelempan[2] + 1))
    fmesh = StructMesh(T, 2)
    lift_data = zeros(T, nVfit)
    drag_data = zeros(T, nVfit)
    speed_data = zeros(T, nVfit)

    Npanhalfwing = sum(Nelempan)
    y = zeros(T, Npanhalfwing)
    ModelAirRaces.cos_wing_tail_yspacing!(y, Nelempan[1:3]..., tot_width/2, 1.0, 1.0)
    Δy = (diff(y[1:Nelempan[1]+1]), diff(y[Nelempan[1]+1:Nelempan[1]+Nelempan[2]+1]), diff(y[Nelempan[1]+Nelempan[2]+1:Nelempan[1]+Nelempan[2]+Nelempan[3]+1]))
    Δy = (Δy...,copy(Δy[1]),copy(Δy[2]),copy(Δy[1]),copy(Δy[1]))
    return (; wres, wmesh, fmesh, tmesh, lift_data, drag_data, speed_data, Δy)
end

function makememory(T, ::Airframe) 
    return (;
    CD      = BuildUp(:all, zero(T)),
    CL      = BuildUp(:all, zero(T)),
    speed   = BuildUp(:all, zero(T)),
    cl      = Dict{Symbol,Vector{T}}(),
    cla     = Dict{Symbol,Vector{T}}(),
    α_gust  = Dict{Symbol,T}(),
    wingΔz  = Dict{Symbol,Vector{T}}(),
    tailΔz  = Dict{Symbol,Vector{T}}(),
    wingθ   = Dict{Symbol,Vector{T}}(),
    tailθ   = Dict{Symbol,Vector{T}}(),
    inertia = BuildUp(:all, InertialElement(zero(T))),
)
end

const NoMemory_Airframe = (;
    CD      = NoBuildUp,
    CL      = NoBuildUp,
    speed   = NoBuildUp,
    cl      = nothing,
    cla     = nothing,
    α_gust  = nothing,
    wingΔz  = nothing,
    tailΔz  = nothing,
    wingθ   = nothing,
    tailθ   = nothing,
    inertia = NoBuildUp,
    )
makememory(::Nothing, ::Airframe) = NoMemory_Airframe

function OptimUtils.ProblemCache(prb::Airframe, sharedvariables; variables=(;), output=(;), options=(;))
    
    # Options
    default_options = (
        airfoil = (SimplerRG15(), SimplerRG15(), SimplerRG15(), # wing
                    SimplerFlatPlate(η=0.2),SimplerFlatPlate(η=0.8), # tail
                    SimplerFlatPlate(),SimplerFlatPlate()), # fuselage
        fuselage = ModelAirRaces.fixedfuselage(),
        Nelempan = [5,12,48,5,12,5,5],
        nVfit = 3, # number of points in fit
        panels_tot_width = 0, # number of panels on fuselage
        takeoff_speed = 15.0,
        fuselage_dw = Val{false}(),
        torsion = Val{true}(),
        dragmodelknots=LinRange(0.,1.,10),
        ) 
    options = merge(default_options, options)

    # Variables
    (; Nelempan, nVfit) = options
    Npanhalfwing = sum(Nelempan)
    default_variables = (
        # Structures
        twing=Var(ini=0.003, lb=0.001, ub=0.005, N=2),
        rwing=Var(ini=0.006, lb=0.001, ub=0.01, N=2),

        # aerodynamics
        cwing=Var(ini=0.15, lb=0.01, ub=0.4, N=2),
        xwing=Var(ini=0.3, lb=0.15, ub=0.35),
        bwing=Var(ini=1.0, lb=0.3, ub=2.5),
        chtail=Var(ini=0.055, lb=0.03, ub=0.12),
        bhtail=Var(ini=0.25, lb=0.2, ub=0.4), # as fraction of tail span
        xhtail=Var(ini=0.8, lb=options.fuselage.tot_length, ub=1.6),
        wtwist=Var(ini=0.0, lb=-15.0, ub=15.0, N=1),

        θtail=Var(ini=0.0, lb=-15.0, ub=5.0, N=options.nVfit+1),
        alpha=Var(ini=3*rand(options.nVfit+1) .+ 2.0, lb=-5.0, ub=15.0, N=options.nVfit+1),

        # Shared variables
        maxload = sharedvariables.maxload,
        speed_maxload = sharedvariables.speed_maxload,
        maxspeed = sharedvariables.maxspeed,
        masspropulsion = sharedvariables.masspropulsion,
        dragmodel = sharedvariables.dragmodel,
        mass = sharedvariables.mass,
    )
    variables = mergevar(default_variables, variables)

    # Output
    default_output = (
        # Aerodynamics
        lift=Var(lb=0.0, ub=0.0, N=nVfit+1),
        stall=Var(lb=-1, ub=1.0, N=Npanhalfwing),
        fastgust_stall=Var(lb=-1, ub=1.0, N=(Nelempan[5] + Nelempan[3]) * (1)),
        takeoff_stall=Var(lb=-1, ub=1.0, N=Npanhalfwing * (1)),
        staticmargin=Var(lb=-Inf, ub=0.0, N=1),
        Cm=Var(lb=0.0, ub=0.0, N=nVfit+1),
        θtailvar=Var(lb=-Inf, ub=0., N=1),

        # Structure
        xwtail=Var(lb=-Inf, ub=0.0, N=2),
        sparradius=Var(lb=-Inf, ub=0.0, N=2),
        wsparthickness=Var(lb=-Inf, ub=0.0, N=2),
        wingstress=Var(lb=-Inf, ub=0., N=nVfit),
        wingdisp=Var(lb=-1., ub=1., N=nVfit),
        tailstress=Var(lb=-1.0, ub=1.0, N=nVfit),
        taildisp=Var(lb=-1, ub=1.0, N=nVfit),

        # Shared variables
        mass = Var(lb=-Inf, ub=0., N=1),
        dragmodelfit = Var(lb=-Inf, ub=0., N=3),
    )
    output = mergevar(default_output, output)
    
    chunksize = (haskey(options,:chunksize)) ? options.chunksize : len(variables)
    return OptimUtils.ProblemCache(prb, variables, output, options, chunksize)
end

function OptimUtils.analysis(p::OptimUtils.ProblemCache{Airframe}, g::Vector{TF}, x::Vector{TF}, memory=NoMemory_Airframe) where {TF}
    cache = (TF<:ForwardDiff.Dual) ? p.usercache_dual : p.usercache
    for i=eachindex(x)
        x[i] = x[i] * (p.x_ub[i] - p.x_lb[i]) + p.x_lb[i]
    end
    v = unpack(x, p.idx)
    g .= 0.0
    idg = p.idg

    # Airframe analysis
    out = airframeanalysis(v, g, idg, cache, p.options, memory)
    g[idg.mass] = (out.mass-v.mass)*2.
    for i=eachindex(idg.dragmodelfit)
        g[idg.dragmodelfit[i]] = (out.dragmodel[i]-v.dragmodel[i])
    end
    # out = airframeanalysis(v, g, idg, cache, p.options, memory)
    # g[idg.mass] = (out.mass-v.mass)*2.
    # dcd0 = (out.dragmodel[1]-v.dragmodel[1])
    # dcd1 = (out.dragmodel[2]-v.dragmodel[2])
    # dcd2 = (out.dragmodel[3]-v.dragmodel[3])
    # lmax = cache.lift_data[1]
    # lmin = cache.lift_data[2]
    # for (w,ig)=zip(p.options.dragmodelknots, idg.dragmodelfit) # shared dragmodel above computed dragmodel at multiple locations
    #     l = lmin + (lmax - lmin) * w
    #     g[ig] = (dcd2 * l^2 + dcd1*l + dcd0)*1e4
    # end

    sumall!(memory.CD)
    sumall!(memory.CL)
    sumall!(memory.speed)
    sumall!(memory.inertia)

    for i=eachindex(x)
        x[i] = (x[i] - p.x_lb[i]) / (p.x_ub[i] - p.x_lb[i])
    end
    return g

end


function subspace_objective(p::OptimUtils.ProblemCache{Airframe}, df::AbstractVector{T}, x::Vector{T}, z::Vector{T}, idz::NamedTuple) where T
    (; idx, idg) = p
    g = p.fdresults.value
    dg = p.fdresults.derivs[1]
    df .= 0.0

    # Inputs
    df[idx.maxload] = x[idx.maxload] - z[idz.maxload]
    df[idx.speed_maxload] = x[idx.speed_maxload] - z[idz.speed_maxload]
    df[idx.maxspeed] = x[idx.maxspeed] - z[idz.maxspeed]
    df[idx.masspropulsion] = x[idx.masspropulsion] - z[idz.masspropulsion]
    df[idx.mass] = x[idx.mass] - z[idz.mass]
    for (ix,iz) = zip(idx.dragmodel, idz.dragmodel)
        df[ix] = x[ix] - z[iz]
    end
    f = sum(f_^2/2 for f_=df)

    return f
end
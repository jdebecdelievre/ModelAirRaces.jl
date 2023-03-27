using Pkg
Pkg.activate(@__DIR__)
using Revise

using Sobol
using JLD2
using LinearAlgebra
using ForwardDiff
using StaticArrays
using ModelAirRaces
const AD = ModelAirRaces
using OptimUtils
using Snopt
using SNOW
using Statistics
using LaTeXStrings
using Printf, Dates
using LatinHypercubeSampling
const T = Float64

include("$(@__DIR__)/plottingroutines.jl")

struct AeroTrajMarathon <: OptimUtils.AbstractProblem end

function OptimUtils.makecache(T::DataType, ::AeroTrajMarathon, options)
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

    ## Propulsion
    thr, spd, pwr, rot, vlt = (zeros(T, 6) for _=1:5)
    propulsion = (; thr, spd, pwr, rot, vlt)
    return (; wres, wmesh, fmesh, tmesh, lift_data, drag_data, speed_data, Δy, propulsion...)
end

const Memory = (;
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
    traj = BuildUp(:all, zero(TrajElement{T})),
)

const NoMemory = (;
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
    traj = NoBuildUp,
    )
##
function ProblemCache(prb::AeroTrajMarathon)
    options = (;
        # airfoil = (ThinAirfoil(), ThinAirfoil(), ThinAirfoil())
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
        nstraight = 4, # number of straight segments in  trajopt,
        imaxspeed = 2, # segment whose end is at max speed. Typically nstraight / 2
        chunksize = 19, # number of chunks for fwd diff
        speedbounds=3.0,
        pwrmodelknots = scaleLHC(LHCoptim(30, 2, 1)[1],[(0.,1.), (0.,1.)]),
        leg_length = envt.marathonlength/120,
        href = 10.0,
    )
    
    # Variables
    (; Nelempan, nVfit, nstraight) = options
    Npanhalfwing = sum(Nelempan)
    variables = (
        # Structures
        twing=Var(ini=0.003, lb=0.001, ub=0.01, N=2),
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

        # shared
        time = Var(lb=0.5, ub=1.1, ini=1., group=(:simpletraj)),
        mass = Var(lb=0.6, ub=1.0, ini=1., group=(:airframe, :simpletraj)),
        dragmodel = Var(lb=[0.0,-0.05,0.5], ub=[0.005, 0.0, 2.], ini=[0.001, -0.005, 1.2], N=3, group=(:airframe, :simpletraj)),
        maxspeed = Var(ini=31.0, lb=25.0, ub=45.0, group=(:airframe, :simpletraj, :propulsion)),
        maxload = Var(ini=3.0, lb=2.0, ub=5.0, group=(:airframe, :simpletraj)),
        speed_maxload = Var(ini=30.0, lb=25.0, ub=40.0, group=(:airframe, :simpletraj, :propulsion)),
        masspropulsion=Var(lb=0.01, ub=0.1, ini=0.05, group=(:airframe, :propulsion)),
        maxthrust=Var(lb=[1.,0.5], ub=[7.,5.], ini=[4.0,4.0], N=2, group=(:propulsion, :simpletraj)),
        pwrmodel = Var(lb=[0.08, -0.4, 0.0, 0.2, 0.1, 0.1], 
                        ub=[0.2,-0.1, 0.1, 0.5, 0.4, 0.2], 
                        ini=[0.15,-0.3,0.05,0.3,0.3,0.15,],
                        N=6, group=(:propulsion, :simpletraj)),

        # Trajectory
        dvdt=Var(ini=0.0, lb=-3.0, ub=3.0, N=nstraight), # acceleration during each segment
        δ=Var(ini=25.0, lb=5.0, ub=80.0, N=nstraight+1),
        dh=Var(ini=0., lb=-5.0, ub=5.0, N=nstraight), # percentage of leg_length
    
        # Propulsion
        Kv=Var(lb=800., ub=2200., ini=1000.),
        R=Var(lb=0.01, ub=1.0, ini=0.2),
        propdiam=Var(lb=0.15, ub=0.5, ini=10.0*0.0254),
        proppitch=Var(lb=0.1, ub=1., ini=0.75),
    )

    # Output
    output = (
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

        mass = Var(lb=-Inf, ub=0., N=1),
        dragmodelfit = Var(lb=-Inf, ub=0., N=3),
        time=Var(ub=0., lb=-Inf),

        # Trajectory
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
        
        # Propulsion
        masspropulsion=Var(lb=-Inf, ub=0.),
        maxrot=Var(lb=-Inf, ub=0., N=6),
        maxthrustmodel = Var(lb=0., ub=0., N=2),
        pwrmodelfit = Var(lb=-Inf, ub=0., N=6),#N=size(options.pwrmodelknots,1)),
        dT=Var(lb=-Inf, ub=0.),
        dV=Var(lb=-Inf, ub=0.)
    )

    return OptimUtils.ProblemCache(prb, variables, output, options, options.chunksize)
end

function OptimUtils.analysis(p::OptimUtils.ProblemCache{AeroTrajMarathon}, g::AbstractVector{TF}, x::AbstractVector{TF}, memory=NoMemory) where TF
    cache = (TF<:ForwardDiff.Dual) ? p.usercache_dual : p.usercache
    
    for i=eachindex(x)
        x[i] = x[i] * (p.x_ub[i] - p.x_lb[i]) + p.x_lb[i]
    end
    v = unpack(x, p.idx)
    g .= 0.0

    idg = p.idg

    ## Airframe
    out = AD.airframeanalysis(v, g, idg, cache, p.options, memory)
    g[idg.mass] = (out.mass-v.mass)*2.
    for i=eachindex(idg.dragmodelfit)
        g[idg.dragmodelfit[i]] = (out.dragmodel[i]-v.dragmodel[i])
    end

    ## Trajectory
    T1, T2 = v.maxthrust
    vmid = (v.maxspeed+v.speed_maxload)/2
    dV = p.options.speedbounds
    v1 = vmid - dV
    v2 = vmid + dV
    dTmdx = (T2-T1) / (v2-v1)
    Tm0 = T1 - dTmdx*v1
    maxthrustmodel = SA[Tm0, dTmdx]
    propulsion = AD.MotorPropeller(maxthrustmodel, SVector{6,TF}(v.pwrmodel), TF(0.0))
    g[idg.dT] = T2-T1

    # Create Airplane
    airplane = AD.Airplane(v.mass, SVector{3,TF}(v.dragmodel), SA[v.speed_maxload, v.maxload], SA[v.maxspeed, 1.0])
    
    # Run trajectory analysis
    out = AD.simpletraj(airplane, propulsion, p.options, v, g, idg, memory.traj)
    
    # Scale output
    g[idg.time] = scale(out.time, P.variables.time)-1.0 #+ (out.time-v.time)^2 # minimize time

    ## Propulsion
    out = AD.propulsionanalysis(v, g, idg, cache, p.options, memory)
    for i=eachindex(idg.pwrmodelfit)
        g[idg.pwrmodelfit[i]] = (out.pwrmodel[i]-v.pwrmodel[i])
    end

    sumall!(memory.CD)
    sumall!(memory.CL)
    sumall!(memory.speed)
    sumall!(memory.inertia)
    sumall!(memory.traj)

    for i=eachindex(x)
        x[i] = (x[i] - p.x_lb[i]) / (p.x_ub[i] - p.x_lb[i])
    end
    return g
end

P = ProblemCache(AeroTrajMarathon())
x0 = ini_scaled(P.variables)
g = ini(P.output)
OptimUtils.analysis(P, P.g, x0, Memory);
v = unscale_unpack(x0, P.idx, P.variables)

##
sbl = SobolSeq(17)
# for i=1:1
i=1
# dirname = joinpath(@__DIR__,Dates.format(now(),"yyyy-mm-dd_HMS"))
dirname = joinpath(@__DIR__,"aao$i")
mkpath(dirname)
snoptions = Dict{Any,Any}("Major iterations limit" => 150,
    "Iterations limit" => 6000,
    "Summary file" => "$dirname/snoptsummary.txt",
    "Print file" => "$dirname/snoptprint.txt",
    "Scale option" => 0,
    "Major optimality tol" => 1e-8,
    # "Verify level"=>1
)
x0 = ini_scaled(P.variables);
x0 .+= rand(size(x0,1)) / 1000.0
# x_ = next!(sbl)
# j=1
# for k=[:time, :mass, :dragmodel, :maxspeed, :maxload, :speed_maxload, :masspropulsion, :maxthrust, :pwrmodel]
#     for idx=P.idx[k]
#         x0[idx] = x_[j]; 
#         j+=1
#     end
# end
ipoptions = Dict{Any,Any}(
    "max_iter" => 250,
    "output_file" => "$dirname/ipoptsummary.txt",
    "print_timing_statistics" => "no",
    "tol"=>1e-6
)
# optimoptions = Options(derivatives=SNOW.UserDeriv(), solver=SNOW.SNOPT(options=snoptions))
optimoptions = Options(derivatives=SNOW.UserDeriv(), solver=IPOPT(ipoptions))

nx = len(P.variables)
lg = lower(P.output)
ug = upper(P.output)
ng = len(P.output)
xopt = copy(x0)
fun = (g_,df_,dg_,x_)->OptimUtils.optfun(P, g_,df_,dg_,x_, :time)
##
# xopt .= x_
xopt, fopt, info, out = SNOW.minimize(fun, xopt, ng, fill(0.0, nx), fill(1.0, nx), lg, ug, optimoptions)
x_ = copy(xopt)
##
g = OptimUtils.analysis(P, P.g, xopt, Memory);
(wres, wmesh, fmesh, tmesh, lift_data, drag_data, speed_data, Δy) = P.usercache;
v = unscale_unpack(xopt, P.idx, P.variables)
P.g .= OptimUtils.analysis(P, P.g, xopt, Memory);
# save("$(dirname)/results.jld2", "xopt", xopt, "Memory", Memory, "P", P, "v", v, "variables", P.variables, "output", P.output, "options", P.options)

(; CD, CL, speed, inertia, traj) = Memory
x_unscaled = OptimUtils.unscale(xopt, P.variables)
c = unpack(P.g, P.idg);

OptimUtils.getactive(P.g, P.output)
OptimUtils.getactive(x_unscaled, P.variables)

##
AR = v.bwing / mean(v.cwing)
taper = v.cwing[2] / v.cwing[1]
bhtail = v.bhtail * v.bwing
ARht = (v.bwing * v.bhtail) / v.chtail
Sratio = v.chtail * bhtail / (mean(v.cwing) * v.bwing)
S = v.bwing * mean(v.cwing)
mass = innertree(inertia, :mass)
xcg = innertree(inertia, :xcg)
lift = (envt.rho_atm * S / 2) * speed^2 * CL
drag = (envt.rho_atm * S / 2) * speed^2 * CD
loadfactor = lift / (inertia.value.mass * envt.g)
dragmodel = AD.estimate_drag_model(lift_data, drag_data, speed_data)
airplane = Airplane(inertia.value.mass, dragmodel, SA[v.speed_maxload, v.maxload], SA[v.maxspeed, 1.0])
dist = innertree(traj, :dist)
ΔE = innertree(traj, :ΔE)
sp = innertree(traj, :speed)
thrust = innertree(traj, :thrust)

CDdata = drag_data ./ (envt.rho_atm/2*S*speed_data.^2 )
CLdata = lift_data ./ (envt.rho_atm/2*S*speed_data.^2 )

Npanhalfwing = sum(wres.Nelempan)
nwing = sum(wres.Nelempan[1:3])
t = OptimUtils.unscale(c.time+1, P.variables.time)*120
##
g = zeros(len(P.output))
df = zeros(len(P.variables))
dg = g .+ df';
OptimUtils.optfun(P, g,df,dg,xopt, :time)
dg_rows = sum(t->t.^2, eachcol(dg))
dg_cols = sum(t->t.^2, eachrow(dg))
##
include(joinpath(@__DIR__,"..","plottingroutines.jl"))
mkpath(joinpath(dirname, "plots"))
trajint, propint = optreport(v, Memory, joinpath(dirname, "plots"), P, showtitles=false);

##
units=(;
    # structures
    twing  = "m",
    rwing  = "m",
    
    # aerodynamics
    cwing  = "m",
    xwing  = "m",
    bwing  = "m",
    chtail = "m",
    bhtail = "\\%",
    xhtail = "m",
    wtwist = "deg",
    
    # trajectory
    θtail  = "deg",
    alpha  = "deg",
    dvdt  = "m/s^2",
    δ = "\\%",
    dh = "rad",
    
    # Propulsion
    Kv="RPM/V",
    R="\\Omega",
    maxpower="W",
    propdiam="m",
    proppitch="m",
    
    # Shared
    masspropulsion="kg",
    time = "s",
    mass = "kg",
    dragmodel = "",
    maxspeed  = "m/s",
    speed_maxload  = "m/s",
    maxload  = "m/s",
    maxthrust="N",
    pwrmodel = ""
)

description=(;
    # structures
    twing  = "Wing spar thickness",
    rwing  = "Wing spar radius",
    
    # aerodynamics
    cwing  = "Wing chord",
    xwing  = "Wing location",
    bwing  = "Wing span",
    wtwist = "Wing tip twist",
    chtail = "Horiz. tail chord",
    bhtail = "Horiz. tail span",
    xhtail = "Horiz. tail location",
    
    # trajectory
    θtail  = "Horiz. tail incidence",
    speed  = "Airspeed",
    alpha  = "Angle of attack",
    
    dvdt  = "Accelerations",
    δ = "Length of segments",
    dh = "Glide angles",

    # Propulsion
    Kv="Motor constant",
    R="Internal motor resistance",
    maxpower="Maximum power consumption",
    propdiam="Propeller diameter",
    proppitch="Propeller pitch",
    
    # Shared
    masspropulsion="Propulsion sys. mass",
    time = "Flight time",
    mass = "Total mass",
    dragmodel = "Drag model",
    maxspeed  = "Max. speed",
    maxthrust = "Max. thrust",
    speed_maxload  = "Speed at max. load factor",
    maxload  = "Max. load factor",
    pwrmodel = "Power model"
)
##
dirname ="/home/adgboost/jdev/ModelAirRaces/examples/simplertraj_aerostruct/aao_mar16_23"
open(joinpath(dirname, "results.tex"),"w") do f
    write(f,OptimUtils.results_summary(P.variables, v, units=units, description=description)) 
end
open(joinpath(dirname, "bounds.tex"),"w") do f
    write(f,OptimUtils.variables_summary(P.variables, units=units, description=description))
end
#
###
xcg = Memory.inertia.value.xcg
xt = (v.xhtail+v.chtail/4)-xcg
xw = (v.xwing+v.cwing[1]/4)-xcg
CL = Memory.CL[:takeoff]
xt*CL[:tail].value + xw*CL[:wing].value

##
include(joinpath(@__DIR__,"..","plottingroutines.jl"))
inp = AD.build_aero_inp(v, 0., P.options.fuselage)
p = plot()
showplanform!(p, inp.xrle, inp.xrte, inp.xtle, inp.xtte, inp.yrle, inp.ytle,horizontal=true, fuselage=true)
showpanels!(p, P.usercache.wres, horizontal=true)
xlabel!("y (m)")
ylabel!("x (m)")
# title!("Vortex System at Convergence")
# savefig("examples/simplertraj_aerostruct/vortexsystem.pdf")
## Custom logging tool for BuildUp package
struct TrajElement{TF}
    t::TF
    dist::TF
    cl::TF
    cd::TF
    λ::TF
    radius::TF
    speed::TF
    lift::TF
    thrust::TF
    ΔE::TF
    Δh::TF
    N::Int64 # useful for averaring operations
end
TrajElement(; t=0., dist=zero(t), cl=zero(t), cd=zero(t), λ=zero(t), radius=zero(t), 
            speed=zero(t), lift=zero(t), thrust=zero(t), ΔE=zero(t), Δh=zero(t), N=1) = TrajElement(promote(t, dist, cl, cd, λ, radius, speed, lift, thrust, ΔE, Δh)...,N)

Base.copy(IE::TrajElement) = TrajElement(IE.t, IE.dist, IE.cl, IE.cd, IE.λ, IE.radius, IE.speed, IE.lift, IE.thrust, IE.ΔE, IE.Δh, IE.N)
Base.zero(::Union{Type{TrajElement{T}},TrajElement{T}}) where T = TrajElement(t=zero(T))
Base.show(io::IO, IE::TrajElement) = Base.show(io, "t, dist, cl, cd, λ, radius, speed, lift, thrust, ΔE, Δh, N")

function Base.:(+)(a::TrajElement{T1} where T1, b::TrajElement{T2} where T2)
    t    = a.t + b.t
    dist = a.dist + b.dist
    λ    = a.λ + b.λ
    ΔE   = a.ΔE + b.ΔE
    Δh   = a.Δh + b.Δh
    N    = a.N + b.N
    
    cl     = (a.cl*a.t + b.cl*b.t) / t
    cd     = (a.cd*a.t + b.cd*b.t) / t
    radius = (a.radius*a.t + b.radius*b.t) / t
    speed  = (a.speed*a.t + b.speed*b.t) / t
    lift  = (a.lift*a.t + b.lift*b.t) / t
    thrust  = (a.thrust*a.t + b.thrust*b.t) / t

    return TrajElement(t, dist, cl, cd, λ, radius, speed, lift, thrust, ΔE, Δh, N)
end

function Base.max(a::TrajElement{T1} where T1, b::TrajElement{T2} where T2)
    t      = max(a.t, b.t)
    dist   = max(a.dist, b.dist)
    λ      = max(a.λ, b.λ)
    ΔE     = max(a.ΔE, b.ΔE)
    Δh     = max(a.Δh, b.Δh)
    N      = max(a.N, b.N)

    cl     = max(a.cl,b.cl)
    cd     = max(a.cd, b.cd)
    speed  = max(a.speed, b.speed)
    lift  = max(a.lift, b.lift)
    thrust  = max(a.thrust, b.thrust)
    radius = max(a.radius, b.radius)

    return TrajElement(t, dist, cl, cd, λ, radius, speed, lift, thrust, ΔE, Δh, N)
end

function Base.min(a::TrajElement{T1} where T1, b::TrajElement{T2} where T2)
    t      = min(a.t, b.t)
    dist   = min(a.dist, b.dist)
    λ      = min(a.λ, b.λ)
    ΔE     = min(a.ΔE, b.ΔE)
    Δh     = min(a.Δh, b.Δh)
    N      = min(a.N, b.N)

    cl     = min(a.cl,b.cl)
    cd     = min(a.cd, b.cd)
    speed  = min(a.speed, b.speed)
    lift  = min(a.lift, b.lift)
    thrust  = min(a.thrust, b.thrust)
    radius = min(a.radius, b.radius)

    return TrajElement(t, dist, cl, cd, λ, radius, speed, lift, thrust, ΔE, Δh, N)
end


## Tools for (1-cos(2π*s)) curvature schedule
"""
Rotation matrix solution of tangential vector equation:
    T(s/sf) = Rs(s/sf,λ) * T(0)
where sf is the total path length, λ the heading change over the turn, T(0) the initial tangential vector.
"""
function Rs(ss,λ)
    return [[cos(λ*ss-sin(2π*ss)*λ/2π), sin(λ*ss-sin(2π*ss)*λ/2π)] [-sin(λ*ss-sin(2π*ss)*λ/2π), cos(λ*ss-sin(2π*ss)*λ/2π)]]
end

"""
Legendre quadrature on [0.,1.]
"""
function nodes_weights_01(n)
    nodes, weights = gausslegendre(n);
    return (nodes .+ 1) /2, weights / 2
end
const nodes, weights = nodes_weights_01(100)
integrate01(f::Function) = sum(f(s)*w for (s,w)=zip(nodes,weights))

"""
Distance between first and last point of turn if the path is of length sf = 1.
If true cross distance cdist is known, get pathlength sf as:
sf = cdist / cross_distance(λ)
"""
function cross_distance(λ)
    nodes, weights = gausslegendre(1000);
    nodes = (nodes .+ 1) /2
    weights = weights / 2
    function helper(s)
        sin_s = sin(2π*s)
        return sincos(λ*(s-sin_s/2π))
    end
    sinInt = 0.
    cosInt = 0.
    for i=eachindex(nodes)
        sc = helper(nodes[i])
        sinInt += sc[1] * weights[i]
        cosInt += sc[2] * weights[i]
    end
    return sqrt(sinInt^2+cosInt^2)
end
function create_cross_distance_fit()
    λ = LinRange(π/10,9π/10,100)
    dis = cross_distance.(λ)
    δ = dis / 2 ./ cos.(λ/2)
    # acd = acos.(dis)
    x = [collect(λ).^2 collect(λ).^4]; ab = (x'x) \ x'*(dis .- 1.)
    return ab
end
const crossdistfit = create_cross_distance_fit()
cross_dist_fit(λ) = 1+ab[1]*λ^2 + ab[2]*λ^4

"""
Distance between the first point and the midpoint of the turn along x direction (assuming T(0)=[0.,1.])
"""
halfturn_x(λ) = integrate01(t->cos(λ*t/2-sin(2π*(t/2))*λ/2π)) / 2
function halfturn_x_fit()
    λ = LinRange(π/10,9π/10,100)
    ht = halfturn_x.(λ)
    dis = ht .* λ ./ sin.(λ / 2) .- 1.
    x = [collect(λ).^2 collect(λ).^4]; ab = (x'x) \ x'*dis
    return ab
end
const halfturn_x_coefs = halfturn_x_fit()
halfturn_x_fit(λ) = (1+halfturn_x_coefs[1]*λ^2 + halfturn_x_coefs[2]*λ^4) * sin(λ/2) / λ

"""
Distance between first point in turn and Bezier control point of turn
if the path is of length sf = 1. If true δ is known, get pathlength sf as:
    sf = δ / δ_dis(λ)
"""
function δ_dis(λ)
    dis = cross_distance(λ)
    δ = dis / 2 / cos(λ/2)
    #δc =  tan(λ/2) ./ λ # circle
    return δ
end

"""
Inverse radius at initial and final point if the path is of length sf = 1. 
If true radius is known, get pathlength sf as:
    sf = radius * r_inv(λ)
Function should be close to linear.
"""
function r_inv(λ) # radius of equivalent circle at initial and final point
    r_inv = 1 ./(cross_distance.(λ)./(2*sin.(λ/2)));
    # r_inv_c = λ (circle case)
    return r_inv
end

function create_r_inv_fit()
    λ = LinRange(π/10,9π/10,100)
    ri = r_inv.(λ)
    x = [ones(100) collect(λ) collect(λ).^2 collect(λ).^3 collect(λ).^4]; 
    ab = (x'x) \ x'*ri
    return ab
end
const r_inv_fit = create_r_inv_fit()
path_length_from_radius(λ, r) = r * (r_inv_fit[1] + r_inv_fit[2]*λ + r_inv_fit[3]*λ^2 + r_inv_fit[4]*λ^3 + r_inv_fit[5]*λ^4)

"""
Safe log(1+x)/x
"""
function log1p_invx(x)
    if x<1e-12
        return 1-x/2+x^2/3 # include information to 2nd order for autodiff
    else
        return log1p(x)/x
    end
end

### Airplane structure with relevant information
struct Airplane{T,n}
    mass::T
    dragmodel::SVector{n,T}
    MaxLoad::SVector{2,T} # (V,n)
    MaxSpeed::SVector{2,T} # (V,n)
end
weight(a::Airplane) = a.mass * envt.g
drag(a::Airplane{T,3}, lift::T, speed::T) where T = a.dragmodel[1] * speed^2 + a.dragmodel[2] * lift + a.dragmodel[3] * lift^2 / speed^2
function clmaxS(a::Airplane) 
    V, n = a.MaxLoad
    W = a.mass * envt.g
    return n*W /(envt.rho_atm/2*V^2)
end
vmax(a::Airplane) = a.MaxSpeed[1]

function cdmodel(a::Airplane, S::Real)
    Q = envt.rho_atm / 2 * S
    cd2 = a.dragmodel[3] *Q
    cl_cdmin = - a.dragmodel[2] / 2 / cd2 
    cdmin = a.dragmodel[1]/Q - cd2*cl_cdmin^2
    return cdmin, cl_cdmin, cd2
end

function dragcoefficient(a::Airplane, S::Real, CL::Real)
    cdmin, cl_cdmin, cd2 = cdmodel(a,S)
    return cdmin + cd2*(CL-cl_cdmin)^2
end

function estimate_drag_model(lift_data::AbstractVector, drag_data::AbstractVector, speed_data::AbstractVector)
    x11 = sum(s^4 for s=speed_data)
    x22 = sum(l^2 for l=lift_data)
    x33 = sum((l^2/s^2)^2 for (s,l)=zip(speed_data,lift_data))
    x12 = sum((l*s^2) for (s,l)=zip(speed_data,lift_data))
    x13 = sum(l^2 for l=lift_data)
    x23 = sum(l^3/s^2 for (s,l)=zip(speed_data,lift_data))
    xtx = SA[x11 x12 x13
             x12 x22 x23
             x13 x23 x33]
    y1 = sum(s^2*d for (s,d)=zip(speed_data,drag_data))
    y2 = sum(l*d for (l,d)=zip(lift_data,drag_data))
    y3 = sum((l^2/s^2)*d for (s,l,d)=zip(speed_data,lift_data,drag_data))
    xty = SA[y1, y2, y3]
    return xtx \ xty
end

function dragmodelfit!(dragmodel::AbstractVector, lift_data::AbstractVector, drag_data::AbstractVector, speed_data::AbstractVector, con, idx=1:3)
    ndata = length(lift_data)
    spd = sum(t->t^2, speed_data) / ndata
    lft = sum(lift_data) / ndata
    lft2_spd2 = sum((l^2/s^2) for (s,l)=zip(speed_data,lift_data)) / ndata
    drg = sum(drag_data) / ndata

    # KKT conditions satified
    drag(lift, speed) = dragmodel[1] * speed^2 + dragmodel[2] * lift + dragmodel[3] * lift^2 / speed^2
    for i=eachindex(lift_data)
        ΔD = (drag_data[i] - drag(lift_data[i], speed_data[i])) / drg / ndata
        con[idx[1]] += ΔD * speed_data[i]^2 / spd
        con[idx[2]] += ΔD * lift_data[i] / lft
        con[idx[3]] += ΔD * (lift_data[i]/speed_data[i])^2 / lft2_spd2
    end
    return con
end

### Propulsion structure with relevant information

### Trajectory Analysis
"""
    straight_leg(A::Airplane, δ,v1,v2,Δh)

Assumes constant glide path angle, airspeed varies linearly with path travelled.
δ: horizontal path length
γ: Glide path angle (positive for altitude gain)
v1: initial speed
v2: final speed
"""
function straight_leg(airplane::Airplane,propulsion::AbstractPropulsionSystem, δ,vi,dvdt,Δh,g,idg,ig)
    W = weight(airplane)

    # Time calculation
    dist = sqrt(δ^2+Δh^2)
    sinγ = Δh / dist
    cosγ = δ / dist
    # g[idg.safety_dv[ig]] = -(vi^2+2*dvdt*dist)
    vf = sqrt(vi^2+2*dvdt*dist)
    t = dist/((vi+vf)/2)

    # Thrust 
    # assumes that min and max are reached at extreme points: not always true because Drag is quadratic in V.
    T1 = W*sinγ + airplane.mass * dvdt + drag(airplane, W*cosγ, vi)
    T2 = W*sinγ + airplane.mass * dvdt + drag(airplane, W*cosγ, vf)
    g[idg.maxthrust[2*ig-1]] = (T1 - Tmax(propulsion, vi)) / 5
    g[idg.maxthrust[2*ig  ]] = (T2 - Tmax(propulsion, vf)) / 5
    g[idg.minthrust[2*ig-1]] = - T1 / 5
    g[idg.minthrust[2*ig  ]] = - T2 / 5 

    # Integrate T*V
    thrust(τ) = W*sinγ + airplane.mass * dvdt + drag(airplane, W*cosγ, vi+dvdt*τ)
    ΔE = sum(power(propulsion, thrust(t*s),vi+dvdt*t*s)*w*t for (s,w)=zip(nodes,weights))

    # Energy constraint
    g[idg.minΔE[ig]] = -ΔE / (envt.batteryenergy*δ/envt.marathonlength)
    

    return TrajElement(speed=dist/t, ΔE=ΔE, t=t, dist=dist, Δh=δ*sinγ/cosγ, thrust=max(T1,T2))
end

function turn_cosinecurvature(airplane::Airplane, propulsion::AbstractPropulsionSystem, δ, λ, speed, γ, g, idg)    
    # Geometry
    horizontal_path = δ / halfturn_x_fit(λ) # path traveled during turn (δ is the path traveled along incoming direction)
    sinγ, cosγ = sincos(γ)
    dist = horizontal_path / cosγ

    # Time calculation
    t = dist / speed
    Q = envt.rho_atm/2*speed^2 # dynamic pressure

    # Stall constraint
    W = weight(airplane)
    max_curvature = 2 * λ / horizontal_path
    L(τ) = sqrt((W*cosγ)^2 + (airplane.mass * speed^2 * λ / horizontal_path * (1-cos(2π*τ/t)))^2)
    max_lift = sqrt((W*cosγ)^2 + (airplane.mass * speed^2 * max_curvature)^2)
    g[idg.MaxLoad[1]] = (max_lift / Q - clmaxS(airplane)) / (0.5*.1) # divide by approx CL*S
    g[idg.MaxLoad[2]] = max_lift / W - airplane.MaxLoad[2]

    # Thrust constraints
    thrust(τ) = W*sinγ + drag(airplane, L(τ), speed)
    thrust_min = thrust(0)
    thrust_max = thrust(t/2)
    g[idg.turnthrust[1]] = (thrust_max - Tmax(propulsion,speed)) / 5.0
    g[idg.turnthrust[2]] = -thrust_min

    # Energy integral
    ΔE = sum(power(propulsion, thrust(s*t), speed)*w*t for (s,w)=zip(nodes,weights))
    l = TrajElement(speed=speed, ΔE=ΔE, t=t, dist=dist, Δh=dist*sinγ,thrust=thrust_max,lift=max_lift)
    return l
end

function turn_loadfactor(λ, δturn, speed)
    sf = δturn / halfturn_x_fit(λ)
    max_curvature = 2 * λ / sf
    return sqrt(1+(speed^2 * max_curvature / envt.g)^2)
end

function simpletraj(airplane::Airplane, propulsion::AbstractPropulsionSystem,
                options, v, g, idg, traj::OptionalBuildUp{TrajElement{T}} where T=NoBuildUp)
    (; leg_length, href, nstraight, imaxspeed) = options
    vturn = v.speed_maxload
    g[idg.sumδ] = 1.0-sum(v.δ)/100

    ## Turn
    δt = v.δ[end]/2*leg_length/100
    L = turn_cosinecurvature(airplane, propulsion, δt, 2π/3, vturn, 0., g, idg)
    addnode(traj, :turn, L)

    ## Straight legs
    trajsl = addnode(traj, :straight)
    vi = vturn
    vmax = vturn
    for i=1:nstraight
        dist = v.δ[i]*leg_length/100
        Δh = v.dh[i]*leg_length/100
        g[idg.slowdownpath[i]] = (-(vi^2+2*v.dvdt[i]*dist) + 1.) / Vref^2
        
        l = straight_leg(airplane, propulsion, dist, vi, v.dvdt[i], Δh, g, idg, i)
        
        vi = vi + l.t*v.dvdt[i]
        vmax = max(vi,vmax)
        g[idg.MaxSpeed[i]] = vi
        g[idg.glidepath[2*i-1]] = (v.dh[i] - v.δ[i] * 0.3) / href # 18 deg max glide up and down
        g[idg.glidepath[2*i]] =  (-v.dh[i] - v.δ[i] * 0.3) / href # 18 deg max glide up and down
        
        L += l
        addnode(trajsl, :leg, l, index=i)
    end
    vmid = g[idg.MaxSpeed[imaxspeed]]
    for i=1:nstraight
        g[idg.MaxSpeed[i]] = (g[idg.MaxSpeed[i]] - vmid)/Vref
    end
    g[idg.MaxSpeed[imaxspeed]] = (v.maxspeed-vmid)/Vref
    g[idg.MaxSpeed[nstraight+1]] = (v.speed_maxload-v.maxspeed)/Vref
    g[idg.dvdt] = (vi - v.speed_maxload)/Vref
    g[idg.dh] = sum(v.dh) / href

    g[idg.Ebat] = envt.marathonlength / leg_length * L.ΔE / (0.8*envt.batteryenergy) - 1.0

    time = envt.marathonlength / leg_length * L.t * (Vref/ envt.marathonlength)
    return (; time)
end

##
function showsimpletraj(dist, traj)
    X = [[0., 0.]]
    ls = sum(dist[1:end-1])
    N = 1000
    ndist = length(dist)

    # Straigth away
    T̄ = [0., 1.]
    for d=dist[1:end-1]
        push!(X, X[end] + d * T̄)
    end

    # Turn
    s = LinRange(0.,traj[:turn].value.dist,N)
    ds = s[2]
    x0 = normalize(X[end])
    for i=eachindex(s)
        push!(X, X[end] + ds * Rs(s[i]/s[end], 2π/3) * x0)
    end
    p1 = X[end-N÷2]
    
    # Straigth away
    T̄ = Rs(1., 2π/3) * normalize(x0)
    for d=dist[1:end-1]
        push!(X, X[end] + d * T̄)
    end
    
    # Turn
    x0 = normalize(X[end]-X[end-1])
    for i=eachindex(s)
        push!(X, X[end] + ds * Rs(s[i]/s[end], 2π/3) * x0)
    end
    p2 = X[end-N÷2]
    
    # Straigth away
    T̄ = Rs(1., 2π/3) * normalize(x0)
    for d=dist[1:end-1]
        push!(X, X[end] + d * T̄)
    end
    
    # Turn
    x0 = normalize(X[end] - X[end-1])
    for i=eachindex(s)
        push!(X, X[end] + ds * Rs(s[i]/s[end], 2π/3) * x0)
    end
    p3 = X[end-N÷2]

    X = -X .+ [[0., ls/2]]
    D = sum(norm, diff(X))
    xpylons = -[p1[1],p2[1],p3[1]]
    ypylons = -[p1[2],p2[2],p3[2]] .+ ls/2

    Nseg = N+ndist
    Xlabels = [X[Nseg:Nseg+ndist-1]..., X[2*Nseg-1]]
    return first.(X), last.(X), D, xpylons, ypylons, first.(Xlabels), last.(Xlabels)
end

function traj_integral(airplane::Airplane, propulsion::AbstractPropulsionSystem, vturn, dvdt, dh, traj)
    N         = 100
    nstraight = length(dvdt)
    speed     = zeros((nstraight+1)*100-nstraight)
    time      = zeros((nstraight+1)*100-nstraight)
    thrust    = zeros((nstraight+1)*100-nstraight)
    drg    = zeros((nstraight+1)*100-nstraight)
    dEdt      = zeros((nstraight+1)*100-nstraight)
    TV      = zeros((nstraight+1)*100-nstraight)
    n_z       = zeros((nstraight+1)*100-nstraight)
    z       = zeros((nstraight+1)*100-nstraight)
    transition_idx = ones(Int64, nstraight+2)
    trajtime = innertree(traj, :t)
    trajdist = innertree(traj, :dist)

    # Straight lines
    W = weight(airplane)
    speed[1] = vturn
    iglb = 1
    for i=1:nstraight
        # time
        tf = trajtime[:straight][Symbol("leg$i")].value
        df = trajdist[:straight][Symbol("leg$i")].value
        transition_idx[i+1] = transition_idx[i] + N-1
        t  = LinRange(0., tf, N)
        
        # thrust
        sinγ = dh[i] / df
        cosγ = sqrt(1-sinγ^2)
        kv = dvdt[i]
        T(v) = W*sinγ + airplane.mass * kv + drag(airplane, W*cosγ, v)

        for j= iglb : iglb+N-1
            τ = t[j-iglb+1]
            time[j] = τ + time[iglb]

            # Speed
            speed[j] = speed[iglb] + kv * τ
            
            # Altitude
            z[j] = z[iglb] + dh[i] * τ/tf
            
            thrust[j] = T(speed[j])
            drg[j] = drag(airplane, W*cosγ, speed[j])

            # dE
            TV[j] = thrust[j] * speed[j]
            dEdt[j] = power(propulsion, thrust[j], speed[j])

            # Load factor
            n_z[j] = cosγ
        end
        # update index
        iglb += N-1
    end

    ## Turn
    tf = trajtime[:turn].value
    t  = LinRange(0., tf, N)
    @. time[iglb : iglb+N-1] = t + time[iglb]
    hpath = trajdist[:turn].value
    λ = 2π/3
    transition_idx[end] = transition_idx[end-1]+N-1

    # Lift and thrust
    cosγ = 1.; sinγ=0.
    L(τ) = sqrt((W*cosγ)^2 + (airplane.mass * vturn^2 * λ / hpath * (1-cos(2π*τ/tf)))^2)
    T_(τ) = W*sinγ + drag(airplane, L(τ), vturn)
    for i=iglb : iglb+N-1
        τ = t[i-iglb+1]
        # Speed
        speed[i] = vturn
        
        # Thrust
        thrust[i] = T_(τ)
        drg[i] = drag(airplane, L(τ), vturn)

        # dE
        TV[i] = thrust[i] * speed[i]
        dEdt[i] = power(propulsion, thrust[i], speed[i])

        # Load factor
        n_z[i] = L(τ)/W
    end

    return (; time, speed, thrust, TV, dEdt, n_z, z, drg), transition_idx
end
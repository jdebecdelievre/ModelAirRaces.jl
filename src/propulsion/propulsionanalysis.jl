## Propeller models
abstract type Propeller end
struct APC{T} <: Propeller
    diam::T
    pitch::T
    ctmodel::NTuple{3,T}
    cpmodel::NTuple{3,T}
end

"""
Maximum RPM based on APC recommendations: https://www.apcprop.com/technical-information/rpm-limits/
"""
maxrpm(prop::APC) = 150000 / (prop.diam/0.0254)

const APC_mass = load(joinpath(@__DIR__,"../../background/propeller/propmassmodel.jld2"), "k")
const APC_ct = load(joinpath(@__DIR__,"../../background/propeller/propaeromodel.jld2"), "ct")
const APC_cp = load(joinpath(@__DIR__,"../../background/propeller/propaeromodel.jld2"), "cp")

mass(prop::APC) = APC_mass[1] + APC_mass[2]*prop.diam
ct(prop::APC, J::Real) = dot(APC_ct, SA[1., J, J^2, prop.pitch/prop.diam, (prop.pitch/prop.diam)^2, prop.pitch/prop.diam*J])
cp(prop::APC, J::Real) = dot(APC_cp, SA[1., J, J^2, prop.pitch/prop.diam, (prop.pitch/prop.diam)^2, prop.pitch/prop.diam*J])

function APC(diam::T, pitch::T) where T
    ctmodel = (APC_ct[1] + APC_ct[4]*pitch/diam + APC_ct[5]*(pitch/diam)^2,
                APC_ct[2] + APC_ct[6]*pitch/diam,
                T(APC_ct[3]))
    cpmodel = (APC_cp[1] + APC_cp[4]*pitch/diam + APC_cp[5]*(pitch/diam)^2,
                APC_cp[2] + APC_cp[6]*pitch/diam,
                T(APC_cp[3]))
    return APC(diam, pitch, ctmodel, cpmodel)
end


"""
Motor model where motor torque varies linearly with RPM, and propeller torque and thrust vary quadratically with RPM.
"""
abstract type Motor end
struct TMotor{T} <: Motor
    Kv::T # RPM/V
    R::T # Ohm
    I0::T # A
    Pmax::T # W
end

const TMOT_I0 = load(joinpath(@__DIR__,"../../background/motor/motor.jld2"), "i0")
const TMOT_M = load(joinpath(@__DIR__,"../../background/motor/motor.jld2"), "m")

function TMotor(; Kv::T=1000.0, R::T=0.100, Pmax::T=100.0) where T
    I0 = 1/(TMOT_I0[1] + TMOT_I0[2] * R)
    return TMotor{T}(Kv, R, I0, Pmax)
end
mass(mot::TMotor) = TMOT_M[1] + TMOT_M[2]* mot.Pmax / (mot.Kv* 2π/60)

function power(m::TMotor, p::APC, thrust::Real, V::Real)
    # RPM from thrust and speed
    a = p.ctmodel[1]*p.diam^4*envt.rho_atm
    b = p.ctmodel[2]*p.diam^3*V*envt.rho_atm
    c = p.ctmodel[3]*p.diam^2*V^2*envt.rho_atm - thrust
    n = (-b + sqrt(b^2-4*a*c)) / (2*a)
    ω = n * 2 * pi

    # Torque from RPM and speed
    J = (V/n/p.diam)
    Q = ((p.cpmodel[1] + p.cpmodel[2]*J + p.cpmodel[3]*J^2) / 2 / pi) * envt.rho_atm * p.diam^5 * n^2

    # Current from Torque
    Kvrad = m.Kv * 2*pi / 60
    i =  Kvrad *Q+m.I0

    # Voltage from current and RPM
    v = n*60/m.Kv+ i*m.R

    # Add losses from battery internal energy
    Rbat = 3 * 0.005 # 5 mΩ per cell for LiPo
    v = v - i * Rbat

    # Add speed controller losses (Modeled as 10% voltage drop)
    v *= 0.9

    return v*i, v,i, Q, ω, J
end

function thrustcalc(m::TMotor, p::APC, volt::Real, V::Real)
    Kvrad = m.Kv * 2*pi / 60

    # RPM from voltage and speed
    a = p.cpmodel[1]/2π*p.diam^5*envt.rho_atm
    b = p.cpmodel[2]/2π*p.diam^4*V*envt.rho_atm + 1/(m.R*Kvrad*m.Kv/60)
    c = p.cpmodel[3]/2π*p.diam^3*V^2*envt.rho_atm - (volt-m.R*m.I0)/(m.R*Kvrad)
    
    n = (-b + sqrt(b^2-4*a*c)) / (2*a)
    ω = n * 2 * pi

    # Torque and thrust from RPM and speed
    J = (V/n/p.diam)
    Q = (p.cpmodel[1] + p.cpmodel[2]*J + p.cpmodel[3]*J^2)/2π * envt.rho_atm * p.diam^5 * n^2
    T = (p.ctmodel[1] + p.ctmodel[2]*J + p.ctmodel[3]*J^2) * envt.rho_atm * p.diam^4 * n^2

    # Current from Torque
    i =  Kvrad *Q+m.I0

    # Add losses from battery internal energy
    Rbat = 3 * 0.005 # 5 mΩ per cell for LiPo
    volt = volt - i * Rbat

    # Add speed controller losses (Modeled as 10% voltage drop)
    volt *= 0.9

    return T, volt*i, i, Q, ω, J
end

function propulsionanalysis(v, g::Vector{TF}, idg, cache, options, memory=NoBuildUp) where TF
    mtr = TMotor(Kv=v.Kv, R=v.R, Pmax=TF(0.))
    prp = APC(v.propdiam, v.proppitch*v.propdiam)

    # Blockage coefficients: induced velocity lower than  expected du to fuselage (est. from Torenbeek)
    # kb = 1-0.329*

    # Make sure max thrust can be reached at min and max speed with nearly discharged battery (3V)
    T1, T2 = v.maxthrust
    vmid = (v.maxspeed+v.speed_maxload)/2
    dV = options.speedbounds
    v1 = vmid - dV
    v2 = vmid + dV
    g[idg.maxthrustmodel[1]] = thrustcalc(mtr, prp, 3*3.0, v1)[1]-T1
    g[idg.maxthrustmodel[2]] = T2 - thrustcalc(mtr, prp, 3*3.0, v2)[1]

    # Build max. thrust model
    dTmdx = (T2-T1) / (v2-v1)
    Tm0 = T1 - dTmdx*v1
    maxthrustmodel = SA[Tm0, dTmdx]
    g[idg.dT] = T2-T1
    g[idg.dV] = (v.speed_maxload-v.maxspeed)

    # Evaluate at 6 points
    (; thr, spd, pwr, rot, vlt) = cache
    spd[1] = vmid
    vlt[1] = 3*3.7
    thr[1], pwr[1], _, _, rot[1] = thrustcalc(mtr, prp, vlt[1], spd[1])
    dT = thr[1] /2 * 0.98

    spd[2] = spd[1] + dV
    thr[2] = thr[1] - dT
    pwr[2], vlt[2], _, _, rot[2] = power(mtr, prp, thr[2], spd[2])
    spd[3] = spd[1] - dV
    thr[3] = thr[1] - dT
    pwr[3], vlt[3], _, _, rot[3] = power(mtr, prp, thr[3], spd[3])
    spd[4] = spd[1] + dV
    thr[4] = thr[1] - 2*dT
    pwr[4], vlt[4], _, _, rot[4] = power(mtr, prp, thr[4], spd[4])
    spd[5] = spd[1] - dV
    thr[5] = thr[1] - 2*dT
    pwr[5], vlt[5], _, _, rot[5] = power(mtr, prp, thr[5], spd[5])
    spd[6] = spd[1]
    thr[6] = thr[1] - 2*dT
    pwr[6], vlt[6], _, _, rot[6] = power(mtr, prp, thr[6], spd[6])

    # Constrain maximum RPM
    Ωmax = maxrpm(prp)
    @. g[idg.maxrot] = (rot*60/2π - Ωmax) / 10000

    # Build reduced order model
    pwrmodel = fitpowermodel(spd, thr, pwr)

    # Evaluate mass
    maxpower = thrustcalc(mtr, prp, 3*4.2, v.maxspeed)[2] # no burning the motor at full throttle, full charge
    mtr = TMotor(Kv=v.Kv, R=v.R, Pmax = maxpower * 1.5)
    propeller = InertialElement(mass(prp), 0.0)
    motor = InertialElement(mass(mtr),0.0)
    propulsionInertia = motor+propeller
    g[idg.masspropulsion] = (propulsionInertia.mass - v.masspropulsion)*10.
    @addbranch memory.inertia propulsion(propeller, motor)

    return (; pwrmodel, maxpower, maxthrustmodel)
end

function fitpowermodel(spd::AbstractVector{TF}, thr::AbstractVector{TF}, pwr::AbstractVector{TF},) where TF
    SPD = SVector{6,TF}(spd) / Vref
    THR = SVector{6,TF}(thr) / Tref
    PWR = SVector{6,TF}(pwr) / Pref
    X = hcat(SVector{6,TF}((1. for _=1:6)...), SPD, THR, SPD .* THR, SPD.^2, THR.^2)
    pwrmodel = X'X \ (X'PWR)
    return pwrmodel
end

const Vref = 25.0
const Pref = 200.0
const Tref = 3.0

function pwrmodelfit!(spd, thr, pwr, v, g, idg)
    n = length(spd)

    power(V,T) = Pref*(v.pwrmodel[1] + v.pwrmodel[2]*V/Vref + v.pwrmodel[3]*T/Tref + v.pwrmodel[4]*T*V/Vref/Tref + v.pwrmodel[5]*V^2/Vref^2 + v.pwrmodel[6]*T^2/Tref^2)

    g[idg.pwrmodelfit[1]] = sum((power(spd[i],thr[i]) - pwr[i]) for i=1:6) / Pref / n
    g[idg.pwrmodelfit[2]] = sum(spd[i]*(power(spd[i],thr[i]) - pwr[i]) for i=1:6) / Pref / n / Vref
    g[idg.pwrmodelfit[3]] = sum(thr[i]*(power(spd[i],thr[i]) - pwr[i]) for i=1:6) / Pref / n / Tref
    g[idg.pwrmodelfit[4]] = sum(thr[i]*spd[i]*(power(spd[i],thr[i]) - pwr[i]) for i=1:6) / Pref / n / Vref / Tref
    g[idg.pwrmodelfit[5]] = sum(spd[i]^2*(power(spd[i],thr[i]) - pwr[i]) for i=1:6) / Pref / n / Vref^2
    g[idg.pwrmodelfit[6]] = sum(thr[i]^2*(power(spd[i],thr[i]) - pwr[i]) for i=1:6) / Pref / n / Tref^2
end

abstract type AbstractPropulsionSystem end
struct MotorPropeller{T} <: AbstractPropulsionSystem
    maxthrust::SVector{2,T}
    powermodel::SVector{6,T}
    maxpower::T
end
Tmax(m::MotorPropeller, V::Real) = m.maxthrust[1] + m.maxthrust[2]*V
power(m::MotorPropeller, T::Real, V::Real) = Pref*(m.powermodel[1] + m.powermodel[2]*V/Vref + m.powermodel[3]*T/Tref + m.powermodel[4]*T/Tref*V/Vref + m.powermodel[5]*V^2/Vref^2 + m.powermodel[6]*T^2/Tref^2)

struct FixedEfficiencyProp{T} <: AbstractPropulsionSystem
    η::T
    maxpower::T
end
efficiency(m::FixedEfficiencyProp, T::Real, V::Real) = m.η 
power(m::FixedEfficiencyProp, T::Real, V::Real) = T*V/efficiency(m,T,V)
Tmax(m::FixedEfficiencyProp, V::Real) = 10.0 #m.max_thrust

## Fixed Efficiency example
const AA246FixedEfficiency = FixedEfficiencyProp(0.375, 150.0)
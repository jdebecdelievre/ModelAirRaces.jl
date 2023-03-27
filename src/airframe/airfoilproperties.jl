#=
General Airfoil type with required methods
=#

abstract type Airfoil end
geom(a::Airfoil) = (; t_c = a.t_c, s_c = a.s_c, x_c = a.x_c)

# Input α and δ in degrees, but derivatives are in per rad.
liftcoefficient(::Airfoil, α::Real, δ::Real, Re::Real) = 0.
dcl_α(::Airfoil, α::Real, δ::Real, Re::Real) = 0.
dcl_δ(::Airfoil, α::Real, δ::Real, Re::Real) = 0.
clmax(::Airfoil, δ, Re) = 0.
clmaxrange(a::Airfoil, δ, Re) = clmax(a, δ, Re) / 1.5

dragcoefficient(::Airfoil, cl::Real, δ::Real, Re::Real) = 0.
dcd_cl(::Airfoil, cl::Real, δ::Real, Re::Real) = 0.
dcd_δ(::Airfoil, cl::Real, δ::Real, Re::Real) = 0.

momentcoefficient(::Airfoil, cl::Real, δ::Real, Re::Real) = 0.
dcm_cl(::Airfoil, cl::Real, δ::Real, Re::Real) = 0.
dcm_δ(::Airfoil, cl::Real, δ::Real, Re::Real) = 0.

#=
Airfoil based on xfoil fitting without flap data
=#
struct XfoilFits{TCL,TCD,TCM} <: Airfoil
    # Lift coefficient fits
    # (clα, αstall, clstar, α0, clmax, cl0)
    clk_cvx::NTuple{6,Int64}    # -1 for concave in loglog plot, 1 for convex
    clk_sgn::NTuple{6,Int64}    # +1 for positive, -1 for negative
    clk_nterms::NTuple{6,Int64} # number of terms in posynomial regression
    clk_fits::TCL
    
    # Drag coefficient fits
    # (cd0, cd1, cd2, cde)
    cdk_cvx::NTuple{4,Int64}    # -1 for concave in loglog plot, 1 for convex
    cdk_sgn::NTuple{4,Int64}    # +1 for positive, -1 for negative
    cdk_nterms::NTuple{4,Int64} # number of terms in posynomial regression
    cdk_fits::TCD
    
    # Moment coefficient fits
    # (cm0, cm1, cm2, cm3)
    cmk_cvx::NTuple{4,Int64}    # -1 for concave in loglog plot, 1 for convex
    cmk_sgn::NTuple{4,Int64}    # +1 for positive, -1 for negative
    cmk_nterms::NTuple{4,Int64} # number of terms in posynomial regression
    cmk_fits::TCM

    # Geom
    t_c::Float64
    s_c::Float64
    x_c::Float64
    symmetric::Bool
end

function XfoilFits(cl_path::String, cd_path::String, cm_path::String, symmetric::Bool;
                t_c::Float64=.1, s_c::Float64=.4, x_c::Float64=.3)
    @load "$cl_path" clk_cvx clk_sgn clk_nterms clk_fits
    @load "$cd_path" cdk_cvx cdk_sgn cdk_nterms cdk_fits
    @load "$cm_path" cmk_cvx cmk_sgn cmk_nterms cmk_fits
    XfoilFits(clk_cvx,clk_sgn,clk_nterms,clk_fits,
            cdk_cvx,cdk_sgn,cdk_nterms,cdk_fits,
            cmk_cvx,cmk_sgn,cmk_nterms,cmk_fits, t_c, s_c, x_c,symmetric)
end

# Possible improvement: use Artifact instead
rel2abs(path) = joinpath(joinpath(@__DIR__, "../../background/airfoils/"), path)

RG15() = XfoilFits(rel2abs("rg15/rg15_clfit.jld2"),
                     rel2abs("rg15/rg15_cdfit.jld2"),
                     rel2abs("rg15/rg15_cmfit.jld2"), false,
                     t_c=0.089, s_c=0.66, x_c=0.3,)

NACA0009() = XfoilFits(rel2abs("naca0009/naca0009_clfit.jld2"),
                     rel2abs("naca0009/naca0009_cdfit.jld2"),
                     rel2abs("naca0009/naca0009_cmfit.jld2"), true,
                     t_c=0.09, s_c=0.6, x_c=0.4,)

function FlatPlate(t_c=0.031)
    s_c = 1.
    x_c = 0.5

    # Lift coefficient fits
    # (clα, αstall, clstar, α0, clmax, cl0)
    clα = 6.25
    clstar = 0.65
    clmax = 0.55
    cl0 = 0.
    αstall = clstar/clα
    clk_cvx = (1,1,1,1,1,1)
    clk_sgn = (1,1,1,1,1,1)
    clk_nterms = (1,1,1,1,1,1)
    clk_fits = ([0.,log(clα)],[0.,log(αstall*180/π)],[0.,log(clstar)], [0.,-100], [0.,log(clmax)],[0.,-100])
    
    # Drag coefficient fits
    # (cd0, cd1, cd2, cde)
    cdk_cvx = (1,1,1,1)
    cdk_sgn = (1,1,1,1)
    cdk_nterms = (1,1,1,1)
    cdk_fits = ([-0.2, log((1+2*t_c) * 2.04 * 0.074)], [0.,-100], [0.,-100], [0.,-100])
    
    # Moment coefficient fits
    # (cm0, cm1, cm2, cm3)
    cmk_cvx = (1,1,1,1)
    cmk_sgn = (1,1,1,1)
    cmk_nterms = (1,1,1,1)
    cmk_fits = ([0.,-100],[0.,-100],[0.,-100],[0.,-100])

    symmetric = true
    return XfoilFits(clk_cvx, clk_sgn, clk_nterms, clk_fits, cdk_cvx, cdk_sgn, cdk_nterms, cdk_fits, cmk_cvx, cmk_sgn, cmk_nterms, cmk_fits,
        t_c, s_c, x_c, symmetric)
end

params(re, cvx, sgn, nterms, fit; logRe=log(re)) = sgn * sum(exp(fit[j]*logRe + fit[nterms+j]) for j=1:nterms)^cvx
params(Re, cvx::Tuple, sgn::Tuple, nterms::Tuple, fit::Tuple; logRe=log(Re)) = ntuple(i->params(Re, cvx[i], sgn[i], nterms[i], fit[i],logRe=logRe),length(cvx))
ftr(c) = 1/(1-c) - 1. - c
dftr(c) = 1/(1-c)^2 - 1.

# Lift coefficient as a function of angle of attack (in deg), and Re number
function liftcoefficient(a::XfoilFits, α::Real, δ::Real, Re::Real)
    # Consider full flying airfoil
    α = α + δ
    
    # Re number variation of parameters
    logRe = log(Re)
    clα    = params(Re, a.clk_cvx[1], a.clk_sgn[1], a.clk_nterms[1], a.clk_fits[1], logRe = logRe)
    clstar = params(Re, a.clk_cvx[3], a.clk_sgn[3], a.clk_nterms[3], a.clk_fits[3], logRe = logRe)
    clmax  = params(Re, a.clk_cvx[5], a.clk_sgn[5], a.clk_nterms[5], a.clk_fits[5], logRe = logRe)
    cl0    = params(Re, a.clk_cvx[6], a.clk_sgn[6], a.clk_nterms[6], a.clk_fits[6], logRe = logRe)
    
    # CL calculation
    f(al) = clstar - sqrt((clα*al*π/180 + cl0 - clstar)^2 + (clmax-clstar)^2)
    if a.symmetric
        α0 = 1/clα * (clstar - sqrt((2*clstar-clmax)*clmax)-cl0)*180/π # very nearly 0 but including for smoothness
        return (α > α0) ? f(α) : f(0.)-f(α0-α)
    else
        return f(α)
    end
    return cl0
end

function dcl_α(a::XfoilFits, α::Real, δ::Real, Re::Real)
    # Consider full flying airfoil
    α = α + δ
    
    # Re number variation of parameters
    logRe = log(Re)
    clα    = params(Re, a.clk_cvx[1], a.clk_sgn[1], a.clk_nterms[1], a.clk_fits[1], logRe=logRe)
    clstar = params(Re, a.clk_cvx[3], a.clk_sgn[3], a.clk_nterms[3], a.clk_fits[3], logRe=logRe)
    clmax  = params(Re, a.clk_cvx[5], a.clk_sgn[5], a.clk_nterms[5], a.clk_fits[5], logRe=logRe)
    cl0    = params(Re, a.clk_cvx[6], a.clk_sgn[6], a.clk_nterms[6], a.clk_fits[6], logRe=logRe)
    # (clα, _, clstar, _, clmax, cl0) = params(Re, a.cdk_cvx, a.cdk_sgn, a.cdk_nterms, a.cdk_fits)
    
    # CL calculation
    cla(α) = clα * (clstar - (clα*α*π/180 + cl0)) / sqrt((clα*α*π/180 + cl0 - clstar)^2 + (clmax-clstar)^2)
    if a.symmetric
        α0 = 1/clα * (clstar - sqrt((2*clstar-clmax)*clmax)-cl0)*180/π # very nearly 0 but including for smoothness
        return (α > α0) ? cla(α) : cla(α0-α)
    else
        return cla(α)
    end
end

dcl_δ(a::XfoilFits, α::Real, δ::Real, Re::Real) = dcl_α(a, α, δ, Re)
clmax(a::XfoilFits, δ::Real, Re::Real) = params(Re, a.clk_cvx[5], a.clk_sgn[5], a.clk_nterms[5], a.clk_fits[5])

# Drag coefficient as a function of C_L, and Re number
function dragcoefficient(a::XfoilFits, cl::Real, δ::Real, Re::Real)
    # Re number variation of parameters
    (cd0, cd1, cd2, cde) = params(Re, a.cdk_cvx, a.cdk_sgn, a.cdk_nterms, a.cdk_fits)

    # CD calculation
    return cd0 + cd1 * cl + cd2 * cl^2 + cde * ftr(cl / clmax(a, δ, Re))
end

function dcd_cl(a::XfoilFits, cl::Real, δ::Real, Re::Real)
    # Re number variation of parameters
    (_, cd1, cd2, cde) = params(Re, a.cdk_cvx, a.cdk_sgn, a.cdk_nterms, a.cdk_fits)
    
    # CD calculation
    return cd1 + 2 * cd2 * cl + cde * dftr(cl / clmax(a, δ, Re)) / clmax(a, δ, Re)
end
dcd_δ(a::XfoilFits, cl::Real, δ::Real, Re::Real) = 0.

# Moment coefficient wrt quarter chord as a function of C_L, and Re number
function momentcoefficient(a::XfoilFits, cl::Real, δ::Real, Re::Real)
    # Re number variation of parameters
    (cm0, cm1, cm2, cm3) = params(Re, a.cmk_cvx, a.cmk_sgn, a.cmk_nterms, a.cmk_fits,)

    # CM calculation
    return cm0 + cm1 * cl + cm2 * cl^2 + cm3 * cl^3
end

function dcm_α(a::XfoilFits, cl::Real, δ::Real, Re::Real)
    # Re number variation of parameters
    (_, cm1, cm2, cm3) = params(Re, a.cmk_cvx, a.cmk_sgn, a.cmk_nterms, a.cmk_fits,)

    # CM calculation
    return cm1 + 2 * cm2 * cl + 3 * cm3 * cl^2
end
dcm_δ(a::XfoilFits, cl::Real, δ::Real, Re::Real) = dcm_α(a, cl, δ, Re)

"""
Negative value corresponds to the domain covered by the valid data.
It corresponds to dcl_α < 0.7 * clα (clα is the lift curve slope far away from stall)
"""
function clmaxrange(a::XfoilFits, δ::Real, Re::Real)
    logRe = log(Re)
    clstar = params(Re, a.clk_cvx[3], a.clk_sgn[3], a.clk_nterms[3], a.clk_fits[3],logRe=logRe)
    clmax  = params(Re, a.clk_cvx[5], a.clk_sgn[5], a.clk_nterms[5], a.clk_fits[5],logRe=logRe)
    return (clstar - (clstar-clmax)/sqrt(1-0.7^2))
end

# 2D parasitic drag model
DragModel2d(Re_ref, t_c) = (1+2*t_c) * 2.04 * (0.074/Re_ref^0.2)

# Thin airfoil
struct ThinAirfoil <: Airfoil end # singleton type
geom(::ThinAirfoil) = (; t_c   = 0.05, s_c = 1., x_c = 0.5)
liftcoefficient(a::ThinAirfoil, α::Real, δ::Real, Re::Real) = 6.3*α*π/180#2π*α*π/180#*(1 + 0.77*(a.t_c))*α*π/180
dcl_α(a::ThinAirfoil, α::Real, δ::Real, Re::Real) = 6.3#*(1 + 0.77*(a.t_c))
dcm_α(::ThinAirfoil, α::Real, δ::Real, Re::Real) = 0.
dcl_δ(a::ThinAirfoil, α::Real, δ::Real, Re::Real) = dcl_α(a, α, δ, Re) # full flying
dcm_δ(::ThinAirfoil, α::Real, δ::Real, Re::Real) = 0.
dragcoefficient(a::ThinAirfoil, α::Real, δ::Real, Re::Real) = DragModel2d(Re, 0.12)#DragModel2d(Re, a.t_c)
momentcoefficient(::ThinAirfoil, α::Real, δ::Real, Re::Real) = 0.
clmax(::ThinAirfoil, δ, Re) = 0.82
clmaxrange(::ThinAirfoil, δ, Re) = 0.77/1.5

# Simpler RG15
# cl only depends on alpha and delta, no Re dependence
# clmax has a re dependance
# Cm only depends on Re, not on cl 
struct SimplerXfoilFits{TCL,TCD,TCM} <: Airfoil
    a::XfoilFits{TCL,TCD,TCM}
    Re::Float64
    η::Float64 # effectiveness (reduced effectiveness can come from being in the slipstream of another body)
    function SimplerXfoilFits(a::XfoilFits{TCL,TCD,TCM}; Re::Float64=2e6, η::Float64=1.0) where {TCL,TCD,TCM}
        new{TCL,TCD,TCM}(a, Re, η)
    end
end
geom(a::SimplerXfoilFits) = geom(a.a)

# Input α and δ in degrees, but derivatives are in per rad.
liftcoefficient(a::SimplerXfoilFits, α::Real, δ::Real, Re::Real) = dcl_α(a.a, 0., 0., a.Re)*α*π/180 + liftcoefficient(a.a, 0., 0., a.Re)
dcl_α(a::SimplerXfoilFits, α::Real, δ::Real, Re::Real) = dcl_α(a.a, 0., 0., a.Re)*a.η
dcl_δ(a::SimplerXfoilFits, α::Real, δ::Real, Re::Real) = dcl_δ(a.a, 0., 0., a.Re)*a.η
clmax(a::SimplerXfoilFits, δ, Re) = clmax(a.a, δ, Re)*a.η
clmaxrange(a::SimplerXfoilFits, δ, Re) = clmaxrange(a.a, δ, Re)*a.η

dragcoefficient(a::SimplerXfoilFits, cl::Real, δ::Real, Re::Real) = dragcoefficient(a.a, cl, δ, Re)*a.η
dcd_cl(a::SimplerXfoilFits, cl::Real, δ::Real, Re::Real) = dcd_cl(a.a, cl, δ, Re)
dcd_δ(a::SimplerXfoilFits, cl::Real, δ::Real, Re::Real) = dcd_δ(a.a, cl, δ, Re)

momentcoefficient(a::SimplerXfoilFits, cl::Real, δ::Real, Re::Real) = momentcoefficient(a.a, 0., 0., a.Re)*a.η
dcm_cl(a::SimplerXfoilFits, cl::Real, δ::Real, Re::Real) = 0.
dcm_δ(a::SimplerXfoilFits, cl::Real, δ::Real, Re::Real) = 0.

SimplerRG15(;Re::Float64=2e6, η::Float64=1.0) = SimplerXfoilFits(RG15(), Re=Re, η=η)
SimplerNACA0009(;Re::Float64=2e6, η::Float64=1.0) = SimplerXfoilFits(NACA0009(), Re=Re, η=η)
SimplerFlatPlate(;Re::Float64=2e6, η::Float64=1.0) = SimplerXfoilFits(FlatPlate(), Re=Re, η=η)
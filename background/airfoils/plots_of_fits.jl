using Plots
using LaTeXStrings
using JLD2
using Printf
##
pltpath = "/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/plots_09232022"
mkpath(pltpath)
cllims = [-0.5, 1.5]
cdlims = [0., 0.05]
cmlims = [-0.05, 0.1]

function params(re, coef_cvx, coef_sgn, coef_nterms, coef_fits)
    return map(1:length(coef_fits)) do i
        k = coef_fits[i]
        nterms = coef_nterms[i]
        coef_sgn[i] * sum(exp.(k[j]*log(re) .+ k[nterms+j]) for j=1:nterms).^(coef_cvx[i])
    end
end

## Saturation for CL
x = LinRange(-5,1,100)
plot(x,t->-sqrt(t^2+1), label="", linewidth=2, 
        # title="Common Saturation Function " * (L" f: x \Rightarrow  -\sqrt{x^2+1}"),
        title="Common Saturation Function ",
        xlabel=L"x",
        ylabel=L"f(x)") 
savefig("$pltpath/simplesaturation.pdf")

## Example CL plot
cls = 1.
clmax = 0.9
cla = 5.
cl0 = 0.1
al = LinRange(-3, 13., 100)
clexample = cls .- sqrt.((cla*al*π/180 .+ (cl0 - cls)).^2 .+ (clmax - cls)^2)
plot(al,clexample, label="", 
        linewidth=2, 
        # thickness_scaling = 2,
        # title=L"c_l(\alpha) = c_l^* + \sqrt{(c_{l\alpha}\alpha + c_{l0} - c_l^*)^2 + (c_l^* - c_{l\sf MAX})^2}",
        title="Lift Coefficient Representation",
        xlabel=L"\alpha {\sf \; (deg)}",
        yaxis=(L"c_l(\alpha)",(-.2, 1.1)),
        ylims = (-.2, 1.2),
        )
savefig("$pltpath/typical_cl_zoomed.pdf")

ylims!(cllims...)
savefig("$pltpath/typical_cl.pdf")

## Divergence for CD
x = LinRange(0,1,100)
plot(x,c->1/(1-c) - 1. - c, label="", linewidth=2, 
        # title="Common Saturation Function " * (L" f: x \Rightarrow  -\sqrt{x^2+1}"),
        title="Possible Function to Represent Drag Divergence",
        xlabel=L"x",
        ylabel=L"f(x)")
savefig("$pltpath/simpledivergence.pdf")

## Example drag
cd0 = 0.005
cd1 = -0.001
cd2 = 0.01
cddiv = 1e-3
cl = LinRange(-.2, clmax, 100)
cdq = map(cl->cd0+cd1*cl+cd2*cl^2, cl)
cdd = map(cl->cddiv*(1/(1-cl/clmax)-1-cl/clmax), cl)
cdexample = cdq+cdd
plot(cl, cdq, linestyle=:dash, label="quadratic term",linewidth=2, )
plot!(cl, cdd, linestyle=:dash, label="divergence term",linewidth=2, )
plot!(cl,cdexample, label="quadratic+divergence", 
        linewidth=2, 
        # thickness_scaling = 1.2,
        title="Drag Coefficient Representation",
        ylabel=L"c_d",
        xaxis=(L"c_l",(-.2, 1.1)),
        ylims = (0, 0.025),
        legend=:topleft
        )
savefig("$pltpath/typical_cd_zoomed.pdf")

ylims!(cdlims...)
# xlims!(cllims...)
savefig("$pltpath/typical_cd.pdf")

##
airfoil = "naca0009"
symmetric=true
airfoilpath = "/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil"

## Available data
@load "$airfoilpath/xfoildata/$(airfoil)xfoil.jld2" c_l c_m c_d re α
@load "$airfoilpath/$(airfoil)_clfit.jld2" cldirect
(nal,nre)=size(cldirect)

legend = map(re) do r
    @sprintf "Re %6d" r
end

# All data
scatter(α,c_l[:,1:nre],label = reshape(legend,1,:), linewidth=2, 
        thickness_scaling=1.2, legend=:topleft, marker=:+,
        xlabel=L"α \; (deg)", ylabel=L"c_l",
        title="Lift Coefficient")
savefig("$airfoilpath/xfoildata/liftcoefficient_$(airfoil)data_zoomed.pdf")
cllims_all = ylims()

scatter(α,c_l[:,1:nre],label = reshape(legend,1,:), linewidth=2, 
        thickness_scaling=1.2, marker=:+,
        xlabel=L"α \; (deg)", ylabel=L"c_l", 
        title="Lift Coefficient",
        legend=:outerright,
        size=(700,400))
        ylims!(cllims...)
savefig("$airfoilpath/xfoildata/liftcoefficient_$(airfoil)data.pdf")

##
scatter(c_l[:,1:nre],c_d[:,1:nre],label = reshape(legend,1,:), linewidth=2, 
        thickness_scaling=1.2, legend=:topleft, marker=:+,
        xlabel=L"c_l", ylabel=L"c_d",
        title="Drag Coefficient")
savefig("$airfoilpath/xfoildata/dragcoefficient_$(airfoil)data_zoomed.pdf")


scatter(c_l[:,1:nre],c_d[:,1:nre],label = reshape(legend,1,:), linewidth=2, 
        thickness_scaling=1.2, marker=:+,
        xlabel=L"c_l", ylabel=L"c_d",
        title="Drag Coefficient",
        legend=:outerright,
        size=(700,400))
xlims!(cllims...)
ylims!(cdlims...)
savefig("$airfoilpath/xfoildata/dragcoefficient_$(airfoil)data.pdf")

##
scatter(α,-c_m[:,1:nre],label = reshape(legend,1,:), linewidth=2, 
        thickness_scaling=1.2, legend=:bottomleft, marker=:+,
        xlabel=L"α \; (deg)", ylabel=L"-c_{m,c/4}",
        title="Moment Coefficient At Quarter Chord",)
savefig("$airfoilpath/xfoildata/momentcoefficient_$(airfoil)data_zoomed.pdf")

scatter(α,-c_m[:,1:nre],label = reshape(legend,1,:), linewidth=2, 
        thickness_scaling=1.2, marker=:+,
        xlabel=L"α \; (deg)", ylabel=L"-c_{m,c/4}",
        title="Moment Coefficient At Quarter Chord",
        legend=:outerright,
        size=(700,400),
        ylims=cmlims,
        )
savefig("$airfoilpath/xfoildata/momentcoefficient_$(airfoil)data.pdf")

## Only data before stall
@load "$airfoilpath/$(airfoil)_cdfit.jld2" imax

P = plot( legend=:topleft,
xlabel=L"α \; (deg)", ylabel=L"c_l", 
title="Lift Coefficient")
for i=1:nre
    scatter!(α[1:imax[i]],c_l[1:imax[i],i],label = legend[i], linewidth=2, 
    thickness_scaling=1.2,marker=:+)
end
cllims_tight = ylims()
savefig("$airfoilpath/xfoildata/liftcoefficient_$(airfoil)data_prestall_zoomed.pdf")

P = plot( legend=:topleft,
xlabel=L"α \; (deg)", ylabel=L"c_l", 
title="Lift Coefficient")
for i=1:nre
    scatter!(α[1:imax[i]],c_l[1:imax[i],i],label = legend[i], linewidth=2, 
    thickness_scaling=1.2,marker=:+, legend=:outerright, size=(700,400))
end
ylims!(cllims...)
savefig("$airfoilpath/xfoildata/liftcoefficient_$(airfoil)data_prestall.pdf")
##
P = plot(  legend=:topleft,
xlabel=L"c_l", ylabel=L"c_d",
title="Drag Coefficient")
for i=1:nre
    scatter!(c_l[1:imax[i],i], c_d[1:imax[i],i],label = legend[i], linewidth=2, 
    thickness_scaling=1.2,marker=:+)
end
cdlims_tight = ylims()
savefig("$airfoilpath/xfoildata/dragcoefficient_$(airfoil)data_prestall_zoomed.pdf")

P = plot(xlabel=L"c_l", ylabel=L"c_d",
title="Drag Coefficient", ylims=cdlims, xlims=cllims, legend=:outerright, size=(700,400))
for i=1:nre
    scatter!(c_l[1:imax[i],i], c_d[1:imax[i],i],label = legend[i], linewidth=2, 
    thickness_scaling=1.2,marker=:+)
end
savefig("$airfoilpath/xfoildata/dragcoefficient_$(airfoil)data_prestall.pdf")
##
P = plot(legend=:bottomleft,
xlabel=L"c_l", ylabel=L"-c_{m,c/4}",
title="Moment Coefficient At Quarter Chord")
for i=1:nre
    scatter!(c_l[1:imax[i],i], -c_m[1:imax[i],i],label = legend[i], linewidth=2, 
        thickness_scaling=1.2,marker=:+)
end
cmlims_tight = ylims()
savefig("$airfoilpath/xfoildata/momentcoefficient_$(airfoil)data_prestall_zoomed.pdf")

P = plot(legend=:outerright, xlims=cllims, ylims=cmlims,size=(700,400),
xlabel=L"c_l", ylabel=L"-c_{m,c/4}",
title="Moment Coefficient At Quarter Chord")
for i=1:nre
    scatter!(c_l[1:imax[i],i], -c_m[1:imax[i],i],label = legend[i], linewidth=2, 
        thickness_scaling=1.2,marker=:+)
end
cmlims_tight = ylims()
savefig("$airfoilpath/xfoildata/momentcoefficient_$(airfoil)data_prestall.pdf")

##
# Lift curve slope
@load "$airfoilpath/$(airfoil)_cdfit.jld2" imax
@load "$airfoilpath/$(airfoil)_clfit.jld2" raw_clk clk_cvx clk_sgn clk_nterms clk_fits

P = plot(
    title = "Lift Curve Slope Predicted by Fits to " *L"c_l(\alpha)",
    ylabel=L"dc_l/d\alpha", xlabel=L"\alpha \; (deg)",
    legend=:bottomleft
)
for i=1:nre
    cli = c_l[:,i]
    (clα, αstall, clstar, α0, clmax, cl0) = params(re[i], clk_cvx, clk_sgn, clk_nterms, clk_fits)
    cla = map(α ->  clα * ((clstar - cl0) - clα * α * π/180) / sqrt(((clstar - cl0) - clα * α * π/180)^2 + (clstar-clmax)^2), α)
    plot!(α[1:imax[i]], cla[1:imax[i]], label=legend[i], linewidth=2, thickness_scaling=1.2)
end
savefig("$airfoilpath/xfoildata/clalpha.pdf")

## Direct Fitting
mkpath("$airfoilpath/cl_directfit")
for i=1:nre
    title = @sprintf "Re %6i" re[i]
    plot(α, c_l[:,i], label="",title=title,xlabel=L"\alpha \; ({\sf deg})", ylabel=L"C_L",marker=:+, linewidth=2, thickness_scaling=1.2)
    plot!(α, cldirect[:,i], label="", linewidth=2, thickness_scaling=1.2)
    savefig("$airfoilpath/cl_directfit/re$i.pdf")
end

## Coefts versus Re
@load "$airfoilpath/$(airfoil)_clfit.jld2" raw_clk clk_cvx clk_sgn clk_nterms clk_fits

# (clα, αstall, clstar, α0, clmax, cl0)
title = (L"c_{lα}", L"αstall", L"c_l^*", L"α0", L"c_{l_{\sf MAX}}", L"c_{l0}",)
Re = LinRange(re[1] , re[nre], 100)
logRe = log.(Re)
lims = (
    cla    = (3.5, 6.5),
    αstall = (8., 17.),
    clstar = (.6, 1.6),
    α0     = (-10., 2.),
    clmax  = (.6, 1.6),
    cl0    = (0., .5),
)

mkpath("$airfoilpath/cl_coefsfit")
for i=1:6
    nterms = clk_nterms[i]
    k = clk_fits[i]
    pred = clk_sgn[i] * sum(exp.(k[j]*logRe .+ k[nterms+j]) for j=1:nterms).^(clk_cvx[i])
    
    plot(re[1:nre], raw_clk[i], label="result of fixed Re data fitting",
        xlabel=L"Re", 
        ylabel=title[i],
        marker=:+, linewidth=2, thickness_scaling=1.2)
    plot!(Re, pred, label="power law interpolation",
        linewidth=2, thickness_scaling=1.2, 
        xguidefontsize=14, yguidefontsize=14,
    )
    savefig("$airfoilpath/cl_coefsfit/coef$(i)_zoomed.pdf")
        
    plot(re[1:nre], raw_clk[i], label="result of fixed Re data fitting",
        xlabel=L"Re", 
        ylabel=title[i],
        marker=:+, linewidth=2, thickness_scaling=1.2,
        ylims=lims[i],
        legend=:outertop,
        size=(600,450)
    )
    plot!(Re, pred, label="power law interpolation",
        linewidth=2, thickness_scaling=1.2, xguidefontsize=14, yguidefontsize=14,
        legend=:outertop,
        size=(600,450)
        )
    savefig("$airfoilpath/cl_coefsfit/coef$i.pdf")

    ## log version
    plot(re[1:nre], clk_sgn[i]*raw_clk[i], label="",
    xlabel=L"Re", 
    ylabel=(clk_sgn[i]>0) ? title[i] : "-"*title[i],
    marker=:+, linewidth=2, thickness_scaling=1.2,
    yaxis=:log10,
    xaxis=:log10, 
    )
    plot!(Re, clk_sgn[i]*pred, label="",
        linewidth=2, thickness_scaling=1.2, xguidefontsize=14, yguidefontsize=14)
    savefig("$airfoilpath/cl_coefsfit/log_coef$i.pdf")
end

## cl predictions at all Re numbers
mkpath("$airfoilpath/cl_end2endprediction")
for ire=1:nre
    title = @sprintf "Re %6i" re[ire]
    cli = c_l[:,ire]
    (clα, αstall, clstar, α0, clmax, cl0) = params(re[ire], clk_cvx, clk_sgn, clk_nterms, clk_fits)
    cl = map(α -> clstar - sqrt(((clstar - cl0) - clα * α * π/180)^2 + (clstar-clmax)^2), α)
    plot(α, cli, label="xfoil data", marker=:+, 
        ylabel=L"c_l", xlabel=L"\alpha \; (deg)", title=title,
        linewidth=2, thickness_scaling=1.2, legend=:bottomright)
    plot!(α, cl, label="prediction",
        linewidth=2, thickness_scaling=1.2,
        ylims=cllims_all)
    savefig("$airfoilpath/cl_end2endprediction/re$(ire)_zoomed.pdf")


    cl = map(α -> clstar - sqrt(((clstar - cl0) - clα * α * π/180)^2 + (clstar-clmax)^2), α)
    plot(α, cli, label="xfoil data", marker=:+, 
        ylabel=L"c_l", xlabel=L"\alpha \; (deg)", title=title,
        linewidth=2, thickness_scaling=1.2, legend=:bottomright)
    plot!(α, cl, label="prediction",
        linewidth=2, thickness_scaling=1.2,
        ylims=cllims)
    savefig("$airfoilpath/cl_end2endprediction/re$ire.pdf")
end

## Drag
# Coefts versus Re
@load "$airfoilpath/$(airfoil)_cdfit.jld2" raw_cdk cdk_cvx cdk_sgn cdk_nterms cdk_fits imax

# (cd0, cd1, cd2, cde)
title = (L"c_{d0}", L"c_{d1}", L"c_{d2}", L"c_{d\sf DIV}")
logRe = log.(Re)

lims=(
    cd0 = (0., 0.05),
    cd1 = (-0.05, 0.0),
    cd2 = (0., 0.05),
    cd3 = (0., 0.05),
)

mkpath("$airfoilpath/cd_coefsfit")
for i=1:4
    plot(re[1:nre], raw_cdk[i],label="result of fixed Re data fitting",
        xlabel="Re", 
        ylabel=title[i],
        marker=:+, linewidth=2, thickness_scaling=1.2)
    nterms = cdk_nterms[i]
    k = cdk_fits[i]
    pred = cdk_sgn[i] * sum(exp.(k[j]*logRe .+ k[nterms+j]) for j=1:nterms).^(cdk_cvx[i])
    plot!(Re, pred,label="power law interpolation",
        linewidth=2, thickness_scaling=1.2,
        xguidefontsize=14, yguidefontsize=14)
    savefig("$airfoilpath/cd_coefsfit/coefs$(i)_zoomed.pdf")
    
    plot(re[1:nre], raw_cdk[i],label="result of fixed Re data fitting",
        xlabel="Re", 
        ylabel=title[i],
        marker=:+, linewidth=2, thickness_scaling=1.2,
        ylims=lims[i],
        legend=:outertop,
        size=(600,450))
    nterms = cdk_nterms[i]
    k = cdk_fits[i]
    pred = cdk_sgn[i] * sum(exp.(k[j]*logRe .+ k[nterms+j]) for j=1:nterms).^(cdk_cvx[i])
    plot!(Re, pred,label="power law interpolation",
        linewidth=2, thickness_scaling=1.2,
        xguidefontsize=14, yguidefontsize=14)
    savefig("$airfoilpath/cd_coefsfit/coefs$i.pdf")
end

## end2end prediction
mkpath("$airfoilpath/cd_end2endprediction")
ftr(c) = 1/(1-c) - 1. - c
for ire=1:nre
    title = @sprintf "Re %6i" re[ire]
    cli = c_l[1:imax[ire],ire]
    cdi = c_d[1:imax[ire],ire]
    (cd0, cd1, cd2, cde) = params(re[ire], cdk_cvx, cdk_sgn, cdk_nterms, cdk_fits)
    clmax = params(re[ire], clk_cvx, clk_sgn, clk_nterms, clk_fits)[5]
    cd = map(cl -> cd0 + cd1 * cl + cd2 * cl^2 + cde * ftr(cl / clmax), cli)
    plot(cli, cdi, label="xfoil data", marker=:+, 
        ylabel=L"C_D", xlabel=L"c_L", title=title,
        linewidth=2, thickness_scaling=1.2, legend=:topleft)
    plot!(cli, cd, label="prediction",
        linewidth=2, thickness_scaling=1.2,
        xlims=cllims_tight,
        ylims=cdlims_tight)
    savefig("$airfoilpath/cd_end2endprediction/re$(ire)_zoomed.pdf")

    plot(cli, cdi, label="xfoil data", marker=:+, 
        ylabel=L"C_D", xlabel=L"c_L", title=title,
        linewidth=2, thickness_scaling=1.2, legend=:topleft)
    plot!(cli, cd, label="prediction",
        linewidth=2, thickness_scaling=1.2,
        xlims=cllims,
        ylims=cdlims)
    savefig("$airfoilpath/cd_end2endprediction/re$ire.pdf")
end

## Moment
@load "$airfoilpath/$(airfoil)_cdfit.jld2" imax
@load "$airfoilpath/$(airfoil)_cmfit.jld2" raw_cmk cmk_cvx cmk_sgn cmk_nterms cmk_fits

title = (L"c_{m0}", L"c_{m1}", L"c_{m2}", L"c_{m3}")
logRe = log.(Re)

mkpath("$airfoilpath/cm_coefsfit")
for i=1:4
    plot(re[1:nre], raw_cmk[i], label="result of fixed Re data fitting",
        # title=title[i],
        xlabel="Re", 
        ylabel=title[i],
        marker=:+, linewidth=2, thickness_scaling=1.2,
        xguidefontsize=14, yguidefontsize=14)
    nterms = cmk_nterms[i]
    k = cmk_fits[i]
    pred = cmk_sgn[i] * sum(exp.(k[j]*logRe .+ k[nterms+j]) for j=1:nterms).^(cmk_cvx[i])
    plot!(Re, pred, label="power law interpolation",
        linewidth=2, thickness_scaling=1.2,
        ylims=-[cmlims[2],cmlims[1]],
        legend=:outertop,
        size=(600,450))
    savefig("$airfoilpath/cm_coefsfit/re$i.pdf")

    ## log version
    plot(re[1:nre], cmk_sgn[i]*raw_cmk[i], label="",
        xlabel="Re", 
        ylabel=(cmk_sgn[i]>0) ? title[i] : ("-"*title[i]),
        marker=:+, linewidth=2, thickness_scaling=1.2,
    yaxis=:log10,
    xaxis=:log10,
    xguidefontsize=14, yguidefontsize=14
    )
    plot!(Re, cmk_sgn[i]*pred, label="",
        linewidth=2, thickness_scaling=1.2, xguidefontsize=14, yguidefontsize=14)
    savefig("$airfoilpath/cm_coefsfit/log_re$i.pdf")
end

## end2end prediction
@load "$airfoilpath/$(airfoil)_cdfit.jld2" imax
@load "$airfoilpath/$(airfoil)_cmfit.jld2" raw_cmk cmk_cvx cmk_sgn cmk_nterms cmk_fits
@load "$airfoilpath/$(airfoil)_clfit.jld2" raw_clk clk_cvx clk_sgn clk_nterms clk_fits
mkpath("$airfoilpath/cm_end2endprediction")

function liftcoefficent(α, clstar, clα, cl0, clmax)
    f(α) = (clstar - sqrt((clstar-(cl0 + clα * α * π/180))^2 + (clstar-clmax)^2))
    if symmetric
        α0 = 1/clα * (clstar - sqrt((2*clstar-clmax)*clmax)-cl0)*180/π # very nearly 0 but including for smoothness
        return (α > α0) ? f(α) : f(0.)-f(α0-α)
    else
        return f(α)
    end
end

for ire=1:nre
    title = @sprintf "Re %6i" re[ire]
    cli = c_l[1:imax[ire],ire]
    cmi = c_m[1:imax[ire],ire]
    (cm0, cm1, cm2, cm3) = params(re[ire], cmk_cvx, cmk_sgn, cmk_nterms, cmk_fits)
    (clα, αstall, clstar, α0, clmax, cl0) = params(re[ire], clk_cvx, clk_sgn, clk_nterms, clk_fits)
    clpred(α) = liftcoefficent(α, clstar, clα, cl0, clmax)
    α = -LinRange(-12.,12.,100)
    cm  = map(α -> cm0 + cm1 * clpred(α) + cm2 * clpred(α)^2 + cm3 * clpred(α)^3, α)
    cl = map(clpred, α)
    plot(cli,-cmi, label="xfoil data", marker=:+, 
        ylabel=L"-c_{m,c/4}", xlabel=L"c_l", title=title,
        linewidth=2, thickness_scaling=1.2, legend=:bottomleft,
        xlims=cllims,
        ylims=cmlims
        )
    plot!(cl, -cm, label="prediction",
        linewidth=2, thickness_scaling=1.2,
        xlims=cllims,
        ylims=cmlims
        )
    savefig("$airfoilpath/cm_end2endprediction/re$ire.pdf")

    plot(cli,-cmi, label="xfoil data", marker=:+, 
        ylabel=L"-c_{m,c/4}", xlabel=L"c_l", title=title,
        linewidth=2, thickness_scaling=1.2, legend=:bottomleft,
        xlims=cllims_tight,
        ylims=cmlims_tight
        )
    plot!(cl, -cm, label="prediction",
        linewidth=2, thickness_scaling=1.2,
        xlims=cllims_tight,
        ylims=cmlims_tight
        )
    savefig("$airfoilpath/cm_end2endprediction/re$(ire)_zoomed.pdf")
end

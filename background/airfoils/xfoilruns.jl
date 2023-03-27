using Xfoil, Plots, Printf
using JLD2
using LinearAlgebra
using Statistics

# read airfoil coordinates from a file
airfoil = "rg15"
symmetric = true

function xfoilrun(airfoil::String, savename="$(airfoil)xfoil")
    # read coordinates
    x, y = open("$(pwd())/UIUC_data/$airfoil/$airfoil.dat", "r") do f
        x = Float64[]
        y = Float64[]
        for line in eachline(f)
            entries = split(chomp(line))
            try
                push!(x, parse(Float64, entries[1]))
                push!(y, parse(Float64, entries[2]))
            catch
                continue
            end
        end
        x, y
    end

    # load airfoil coordinates into XFOIL
    Xfoil.set_coordinates(x,y)

    # repanel using XFOIL's `PANE` command
    xr, yr = Xfoil.pane()
    plot(xr, yr, label="", framestyle=:none, aspect_ratio=1.0, show=true, color=:black, linewidth=3)
    savefig("$(pwd())/UIUC_data/$(airfoil)_geom.pdf")

    # prealloate arrays
    nal = 110
    nre = 13
    α   = LinRange(-3., 12, nal) # range of angle of attacks, in degrees
    re  = 10. .^ LinRange(4.6, 5.5,nre)

    c_l       = zeros(nal, nre)
    c_d       = zeros(nal, nre)
    c_dp      = zeros(nal, nre)
    c_m       = zeros(nal, nre)
    xsep_u    = zeros(nal, nre)
    xsep_l    = zeros(nal, nre)
    converged = zeros(Bool, nal, nre)

    # inviscid calculations
    c_li = zeros(nal)
    c_mi = zeros(nal)
    for i=eachindex(α)
        c_li[i], c_mi[i] = Xfoil.solve_alpha(α[i])
    end

    # viscous calculations
    for j=1:nre
        reinit = true
        for i = 1:nal
            c_l[i,j], c_d[i,j], c_dp[i,j], c_m[i,j], converged[i,j] = Xfoil.solve_alpha(α[i], re[j]; iter=100, reinit=reinit, ncrit=5, xtrip=(0.05,0.05))
            # xsep_u[i,j], xsep_l[i,j] = Xfoil.get_xsep(xtol=1e-4)
        
            reinit = false
            if converged[i,j] < 1
                @printf "%i\t%i\t%i\n" i j converged[i,j]
                reinit = true
                if (i>5) && sum(converged[i-3:i,j])==0
                    break
                end
            end
        end
    end

    # fill-in unconverged runs
    for j=1:nre
        for i=1:nal
            if !converged[i,j]
                if 1==i
                    c_l[1,j] = 2*c_l[2,j] - c_l[3,j]
                    c_d[1,j] = 2*c_d[2,j] - c_d[3,j]
                    c_dp[1,j] = 2*c_dp[2,j] - c_dp[3,j]
                    c_m[1,j] = 2*c_m[2,j] - c_m[3,j]
                elseif i==nal
                    c_l[nal,j] = 2*c_l[nal-1,j] - c_l[nal-2,j]
                    c_d[nal,j] = 2*c_d[nal-1,j] - c_d[nal-2,j]
                    c_dp[nal,j] = 2*c_dp[nal-1,j] - c_dp[nal-2,j]
                    c_m[nal,j] = 2*c_m[nal-1,j] - c_m[nal-2,j]
                else
                    c_l[i,j] = (c_l[i+1,j]+c_l[i-1,j])/2
                    c_d[i,j] = (c_d[i+1,j]+c_d[i-1,j])/2
                    c_dp[i,j] = (c_dp[i+1,j]+c_dp[i-1,j])/2
                    c_m[i,j] = (c_m[i+1,j]+c_m[i-1,j])/2
                end
            end
        end
    end

    # Save
    @save "/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil/xfoildata/$savename.jld2" c_l c_d α re c_dp converged c_m xsep_u xsep_l
end
##
xfoilrun(airfoil, "$(airfoil)xfoil_oct20")
##
# a   = load("/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil/xfoildata/$(airfoil)xfoil.jld2")
a   = load("/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil/xfoildata/$(airfoil)xfoil_oct20.jld2")

nre = 12
c_l = a["c_l"][:,1:nre]
c_d = a["c_d"][:,1:nre]
c_m = a["c_m"][:,1:nre]
α   = a["α"]
re  = a["re"][1:nre]
nal = length(α)
al  = α*π/180;

## Direct Fit
clk =  zeros(4, nre)
for r in 1:nre
    X     = [2*c_l[:,r] al.^2 -al -ones(nal)];
    Nx    = map(norm, eachcol(X))
    cl2   = c_l[:,r].^2
    Ncl2  = norm(cl2)
    map((x,nx)->(x ./= nx),eachcol(X),Nx)
    cl2 ./= Ncl2
    clk[:, r] = (X'X) \ (X'cl2)
    @. clk[:, r] *= Ncl2/Nx
end
cldirect  = (clk[1,:] .- sqrt.(clk[2,:] .* (al').^2 - clk[3,:] .* (al') .- clk[4,:] .+ clk[1,:].^2))';
clαdirect = - ((clk[2,:] .* (al') .- clk[3,:]/2) ./ sqrt.(clk[2,:] .* (al').^2 - clk[3,:] .* (al') .- clk[4,:] .+ clk[1,:].^2))';

## Coefficients vs Re
clα         = sqrt.(clk[2,:])
αstall      = clk[3,:] / 2 ./ clk[2,:] * 180 / π
clstar      = clk[1,:]
α0          = αstall - clk[1,:] ./ clα * 180 / π
clmax       = clk[1,:] .- sqrt.(clk[2,:] .* (αstall*π/180).^2 - clk[3,:] .* (αstall*π/180) .- clk[4,:] .+ clk[1,:].^2)
cl0         = - clα .* α0*π/180
raw_clk   = (clα, αstall, clstar, α0, clmax, cl0)
# α0     = 1/clα * (clstar - sqrt((2*clstar-clmax)*clmax)-cl0)*180/π
##
using LsqFit
function fit_1dposynomial(xdata, ydata; nterms=2, p0=(2*rand(2*nterms).-1)/1000)
    ndata = length(xdata)
    model!(F,x,p) = map!(x->log(sum(exp(p[i]*x + p[nterms+i]) for i=1:nterms)),F, x)
    function jacobian!(J,x,p)
        foreach(i->map!(x->exp(p[i]*x + p[nterms+i]),view(J,1:ndata,nterms+i), x), 1:nterms) # all biases
        foreach(i->map!((Ji, x)-> Ji * x, view(J,1:ndata,i), x, view(J,1:ndata,nterms+i)), 1:nterms) # all biases
        for Jj = eachrow(J)
            sJ = sum(view(Jj,nterms+1:2*nterms))
            @. Jj /= sJ 
        end
    end
    # curve_fit(model!, jacobian!, xdata, ydata, p0, inplace=true)
    curve_fit(model!, xdata, ydata, p0, inplace=true)
end

if airfoil=="rg15"
    # (clα, αstall, clstar, α0, clmax, cl0)
    clk_cvx       = (-1, 1, 1, 1,-1,-1) # -1 for concave in loglog plot, 1 for convex
    clk_sgn       = ( 1, 1, 1,-1, 1, 1) # 0 for positive,                +1 for negative
    clk_nterms    = ( 2, 2, 2, 2, 1, 2) # number of terms in posynomial regression
elseif airfoil=="naca0009"
    # (clα, αstall, clstar, α0, clmax, cl0)
    clk_cvx       = (-1, 1, 1, 1, 1, 1) # -1 for concave in loglog plot, 1 for convex
    clk_sgn       = ( 1, 1, 1,-1, 1, 1) # 0 for positive,                +1 for negative
    clk_nterms    = ( 2, 2, 2, 2, 2, 2) # number of terms in posynomial regression
end
logRe          = log.(re[1:nre])
fits           = ntuple(i->fit_1dposynomial(logRe, clk_cvx[i]*log.(clk_sgn[i]*raw_clk[i]), nterms=clk_nterms[i]), 6)
clk_fits     = ntuple(i->fits[i].param, 6)
clk_fits_res = ntuple(i->fits[i].resid, 6)
@save "/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil/$(airfoil)_clfit.jld2" cldirect raw_clk clk_cvx clk_sgn clk_nterms clk_fits clk_fits_res

## DRAG

# Will not consider data if clα is less than 70% of clα at 0 lift.
imax = map(zip(eachcol(cldirect),eachcol(clαdirect))) do (c,cα)
    i0lift = findfirst(t->(t>0.), c)
    findfirst(t->(t<cα[i0lift]*0.7), cα)
end

##
ftr(c) = 1/(1-c) - 1. - c
cd0 = zeros(nre)
cd1 = zeros(nre)
cd2 = zeros(nre)
cde = zeros(nre)
cddirect = [zeros(i) for i=imax]
for j=(1:nre)
    c_d1 = c_d[1:imax[j],j]
    c_l1 = c_l[1:imax[j],j]
    clm  = clmax[j]
    x = [ones(imax[j]) (c_l1).^2 ftr.(c_l1/clm) c_l1 c_d1]
    
    # Fit with total least squares
    U, S, V = svd(x)
    n = size(x,2)-1
    VXY = V[1:n, 1+n:end]
    VYY = V[1+n:end, 1+n:end]
    B   = -VXY / VYY

    cd0[j] = B[1]
    cd1[j] = B[4]
    cd2[j] = B[2]
    cde[j] = B[3]
    cddirect[j] .= x[:,1:n]*B
end
raw_cdk = (cd0, cd1, cd2,cde);

P = plot()
for j=eachindex(imax)
    scatter!(P, c_l[1:imax[j],j], c_d[1:imax[j],j], label="", marker=:+, color=j)
    plot!(P, c_l[1:imax[j],j], cddirect[j], label="", color=j)
end
P
##
if airfoil == "rg15"
    cdk_cvx       = ( 1, 1, 1, 1) # -1 for concave in loglog plot, 1 for convex
    cdk_sgn       = ( 1,-1, 1, 1) # 1 for positive, +1 for negative
    cdk_nterms    = ( 2, 1, 1, 1) # number of terms in posynomial regression
elseif airfoil == "naca0009"
    cdk_cvx       = ( 1, 1, 1, 1) # -1 for concave in loglog plot, 1 for convex
    cdk_sgn       = ( 1,-1, 1, 1) # 1 for positive, +1 for negative
    cdk_nterms    = ( 2, 1, 2, 1) # number of terms in posynomial regression
end
logRe          = log.(re[1:nre])
fits           = ntuple(i->fit_1dposynomial(logRe, cdk_cvx[i]*log.(cdk_sgn[i]*raw_cdk[i]), nterms=cdk_nterms[i]), 4)
cdk_fits     = ntuple(i->fits[i].param, 4)
cdk_fits_res = ntuple(i->fits[i].resid, 4)
@save "/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil/$(airfoil)_cdfit.jld2" cddirect raw_cdk cdk_cvx cdk_sgn cdk_nterms cdk_fits cdk_fits_res imax

##
function params(re, coef_cvx, coef_sgn, coef_nterms, coef_fits)
    return map(1:length(coef_fits)) do i
        k = coef_fits[i]
        nterms = coef_nterms[i]
        coef_sgn[i] * sum(exp.(k[j]*log(re) .+ k[nterms+j]) for j=1:nterms).^(coef_cvx[i])
    end
end

function liftcoefficent(α, clstar, clα, cl0, clmax)
    f(α) = (clstar - sqrt((clstar-(cl0 + clα * α * π/180))^2 + (clstar-clmax)^2))
    if symmetric
        α0 = 1/clα * (clstar - sqrt((2*clstar-clmax)*clmax)-cl0)*180/π # very nearly 0 but including for smoothness
        return (α > α0) ? f(α) : f(0.)-f(α0-α)
    else
        return f(α)
    end
end


P=plot()
cm0 = -ones(nre) * eps()
cm1 =  ones(nre) * eps()
cm2 =- ones(nre) * eps()
cm3 =  ones(nre) * eps()
cmdirect = [zeros(i) for i=imax]
e = 1e-10
for i=1:nre
    (clα, αstall, clstar, α0, clmax, cl0) = params(re[i], clk_cvx, clk_sgn, clk_nterms, clk_fits)
    clpred(α) = liftcoefficent(α, clstar, clα, cl0, clmax)
    cl  = map(clpred, α)
    if symmetric
        x = [cl[1:imax[i]] cl[1:imax[i]].^3]
        k = (x'x + e*I) \ (x'*(c_m[1:imax[i],i]))
        cm1[i], cm3[i] = k
    else
        x = [ones(imax[i]) cl[1:imax[i]] cl[1:imax[i]].^2 cl[1:imax[i]].^3]
        k = (x'x + e*I) \ (x'*(c_m[1:imax[i],i]))
        cm0[i], cm1[i], cm2[i], cm3[i] = k
    end
    cmdirect[i] .= x*k
    plot!(P, α[1:imax[i]], -c_m[1:imax[i],i], label="", marker=:+, alpha=.5)
    plot!(P, α[1:imax[i]], -x*k, label="")
end
raw_cmk = (cm0, cm1, cm2, cm3)
P
##
if airfoil == "naca0009"
    cmk_sgn    = (-1, 1,-1, 1) # 1 for positive, +1 for negative
    cmk_cvx    = (-1, 1,-1,-1) # -1 for concave in loglog plot, 1 for convex
    cmk_nterms = ( 1, 2, 1, 2) # number of terms in posynomial regression
elseif airfoil == "rg15"
    cmk_sgn    = (-1, 1,-1, 1) # 1 for positive, +1 for negative
    cmk_cvx    = (-1, 1,-1, 1) # -1 for concave in loglog plot, 1 for convex
    cmk_nterms = ( 2, 1, 1, 1) # number of terms in posynomial regression
end
logRe        = log.(re[1:nre])
fits         = ntuple(i->fit_1dposynomial(logRe, cmk_cvx[i]*log.(cmk_sgn[i]*raw_cmk[i]), nterms=cmk_nterms[i]), 4)
cmk_fits     = ntuple(i->fits[i].param, 4)
cmk_fits_res = ntuple(i->fits[i].resid, 4)
@save "/Users/jeandebecdelievre/.julia/dev/AeroDesign/UIUC_data/$airfoil/$(airfoil)_cmfit.jld2" cmdirect raw_cmk cmk_cvx cmk_sgn cmk_nterms cmk_fits cmk_fits_res imax

## NACA
# f(x) = 0.025*((x<0.5) ? x/(x+0.01) : (1-x)/(.05+(1-x)))
f(x) = 0.025*( x/(x+0.01) + (1-x)/(.05+(1-x))-1.)
n = 100
x = LinRange(0., 1., 100)
y = map(f,x)
y .-= y[1]
y = [y; 0.]
x = [x; x[end]+x[2]] / (x[end]+x[2])
y = [reverse(y); -y[2:end]]
x = [reverse(x); x[2:end]]
# x = [reverse(x); x[2:end]]

scatter(x,y,aspectratio=1)


##
# load airfoil coordinates into XFOIL
Xfoil.set_coordinates(x,y)

# repanel using XFOIL's `PANE` command
xr, yr = Xfoil.pane()
scatter(xr, yr, label="", framestyle=:none, aspect_ratio=1.0, show=true, color=:black, linewidth=3)
##
# prealloate arrays
nal = 110
nre = 13
α   = LinRange(-3., 12, nal) # range of angle of attacks, in degrees
re  = 10. .^ LinRange(3.6, 4.5,nre)
c_l       = zeros(nal, nre)
c_d       = zeros(nal, nre)
c_dp      = zeros(nal, nre)
c_m       = zeros(nal, nre)
xsep_u    = zeros(nal, nre)
xsep_l    = zeros(nal, nre)
converged = zeros(Bool, nal, nre);

Xfoil.solve_alpha(0., 100000; iter=100, reinit=false)#, ncrit=5, xtrip=(0.01,0.01))

##
# viscous calculations
for j=1:nre
    reinit = true
    for i = 1:nal
        c_l[i,j], c_d[i,j], c_dp[i,j], c_m[i,j], converged[i,j] = Xfoil.solve_alpha(α[i], re[j]; iter=200, reinit=reinit, ncrit=5, xtrip=(0.01,0.01))
        # xsep_u[i,j], xsep_l[i,j] = Xfoil.get_xsep(xtol=1e-4)
    
        reinit = false
        if converged[i,j] < 1
            @printf "%i\t%i\t%i\n" i j converged[i,j]
            reinit = true
        end
    end
end

# fill-in unconverged runs
for j=1:nre
    for i=1:nal
        if !converged[i,j]
            if 1==i
                c_l[1,j] = 2*c_l[2,j] - c_l[3,j]
                c_d[1,j] = 2*c_d[2,j] - c_d[3,j]
                c_dp[1,j] = 2*c_dp[2,j] - c_dp[3,j]
                c_m[1,j] = 2*c_m[2,j] - c_m[3,j]
            elseif i==nal
                c_l[nal,j] = 2*c_l[nal-1,j] - c_l[nal-2,j]
                c_d[nal,j] = 2*c_d[nal-1,j] - c_d[nal-2,j]
                c_dp[nal,j] = 2*c_dp[nal-1,j] - c_dp[nal-2,j]
                c_m[nal,j] = 2*c_m[nal-1,j] - c_m[nal-2,j]
            else
                c_l[i,j] = (c_l[i+1,j]+c_l[i-1,j])/2
                c_d[i,j] = (c_d[i+1,j]+c_d[i-1,j])/2
                c_dp[i,j] = (c_dp[i+1,j]+c_dp[i-1,j])/2
                c_m[i,j] = (c_m[i+1,j]+c_m[i-1,j])/2
            end
        end
    end
end

##
plot(c_l, c_d)
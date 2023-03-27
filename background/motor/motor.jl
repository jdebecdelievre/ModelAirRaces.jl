using Plots
using DataFrames
using LinearAlgebra
using CSV
using LaTeXStrings
using JLD2
##
motor = CSV.read(joinpath(@__DIR__,"motordata.csv"),DataFrame,delim=',', header=1)
scatter(motor.Kv, motor.I0, motor.R0)
n = length(motor.mass)
##
Kv = motor.Kv * 2π/60 # Kv in rad/s/V
pmax_Kv = motor.Pmax ./ Kv
mass = motor.mass / 1000

a = SA[n sum(pmax_Kv); sum(pmax_Kv) sum(t->t^2,pmax_Kv)]
m = a \ SA[sum(mass); sum(p*m for (p,m)=zip(pmax_Kv,mass))]
scatter(pmax_Kv, mass, xlabel=L"P_{\max}/K_v \; (WVs)", ylabel=L"{\sf mass } \; (kg)", label="T-Motor AT Series", legend=:bottomright)
plot!(LinRange(minimum(pmax_Kv), maximum(pmax_Kv),50), t->(m[1]+m[2]*t), label="Fit "*L"{\sf mass}(P_{\max}/K_v) = {m_0+m_1P_{\max}/K_v}", ylims=(0.,0.2), xlims=(0.,9.))
title!("Mass vs. Max. Power over Motor Constant \n for T-Motor AT Series ")
savefig(joinpath(@__DIR__,"motormass.pdf"))
##
R0 = motor.R0 / 1000
a = SA[n sum(R0); sum(R0) sum(t->t^2,R0)]
i0 = a \ SA[sum(1/i for i=motor.I0), sum(r/i for (i,r)=zip(motor.I0, R0))]

# scatter(R0, 1.0 ./motor.I0, xlabel=L"R_0 \; (\Omega)", ylabel=L"1/I_0 \; (A^{-1})", label="T-Motor AT Series")
# plot!(LinRange(minimum(R0), maximum(R0),50), t->i0[1]+i0[2]*t, label="Linear Fit", legend=:bottomright)
##
scatter(R0, motor.I0, xlabel=L"R_0 \; (\Omega)", ylabel=L"I_0 \; (A)", label="T-Motor AT Series")
plot!(LinRange(minimum(R0), maximum(R0),50), t->1/(i0[1]+i0[2]*t), label="Fit "*L"I(R_0) = \frac{1}{i_0+i_1R_0}", ylims=(0.,3.), xlims=(0.,0.26))
title!("Internal Resistance vs. Idle Current \n for T-Motor AT Series ")
savefig(joinpath(@__DIR__,"motori0.pdf"))
##
plot!(LinRange(minimum(R0), maximum(R0),50), t->(0.0467 / t)^(1/1.892))
##
save(joinpath(@__DIR__,"motor.jld2"),"i0",i0,"m",m)
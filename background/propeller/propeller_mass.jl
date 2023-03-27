"""
Use mass of thin electric propellers to come up with a mass model
"""

using JLD2
using CSV
using DataFrames
using Dates
using Plots
using LinearAlgebra
using StaticArrays
##
dat = CSV.read(joinpath(@__DIR__,"APC_data/performance.csv"), DataFrame)
cond = (d-> (endswith(d,"E") && !(startswith(d,"B")))&& !endswith(d,"WE"))
thinelec = map(cond,dat."Product Name")
dat = dat[thinelec,:]
##
mass = map(d->parse(Float64, d[1:end-3]), (dat."Weight")) * 28.35 / 1000
diam = map(d->parse(Float64, d), (dat."Diameter (INCHES)")) * 0.0254
pitch = map(d->parse(Float64, d),(dat."Pitch (INCHES)")) * 0.0254
##
dmax = 15* 0.0254
dmin = 6* 0.0254
dm = map(d->isless(d,dmax)&&isless(dmin,d),diam)
mass = mass[dm]
diam = diam[dm]
pitch = pitch[dm]
n = length(mass)
##
# d1 = diam
# d3 = diam.^3
# XX = SA[n sum(d1) sum(d3); sum(d1) d1⋅d1 d1⋅d3; sum(d3) d1⋅d3 d3⋅d3]
# Xy = SA[sum(mass),d1⋅mass,d3⋅mass]
# k = XX \ Xy
# ##
# scatter(diam,mass,xlabel="Diameter (m)",ylabel="Mass (kg)",marker_z=pitch, zlabel="Pitch (m)", label="APC propellers", legend=:bottomright)
# plot!(LinRange(dmin,dmax,50), d->k[1]+d*k[2]+d^3*k[3], label="linear fit")
# 
##
d = diam
XX = SA[n sum(d); sum(d) d⋅d]
Xy = SA[sum(mass),d⋅mass]
k = XX \ Xy
##
scatter(diam,mass,xlabel="Diameter (m)",ylabel="Mass (kg)",
        marker_z=pitch./diam, colorbar_title="Pitch (m)", label="APC propellers", 
        legend=:bottomright, title="Mass Variations of APC Thin Electric Propellers")
plot!(LinRange(dmin,dmax,50), d->k[1]+d*k[2], label="linear fit")
# savefig(joinpath(@__DIR__,"massmodel.pdf"))
# save(joinpath(@__DIR__,"propmassmodel.jld2"),"bounds",(dmin,dmax), "k", k)
##
d = load(joinpath(@__DIR__,"propmassmodel.jld2"))
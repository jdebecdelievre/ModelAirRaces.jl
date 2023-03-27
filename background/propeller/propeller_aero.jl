"""
Use mass of thin electric propellers to come up with a model of ct and cp versus J, diameter, pitch
"""

using JLD2
using CSV
using DataFrames
using Dates
using Plots
using LinearAlgebra
using StaticArrays

##
datadir = joinpath(@__DIR__,"UIUC_data")
files = readdir(datadir)
a = DataFrame[]
for f = files
    if occursin("static",f) || startswith(f,"README") || endswith(f,"geom.txt")
        continue
    end
    
    ff = split(f, '_')
    dxp = split(ff[2], 'x')
    diam = parse(Float64,dxp[1])
    pitch = parse(Float64,dxp[2])
    if diam > 15
        continue
    end

    df = CSV.read(joinpath(datadir,f), DataFrame, delim="   ", skipto=2, header=["J","CT","CP","eta"])
    df.diam .= diam * 0.0254
    df.pitch .= pitch * 0.0254
    df.rpm .= parse(Float64,ff[end][1:end-4])
    push!(a, df)
end
dat = vcat(a...);

##
n = length(dat.J)
x1 = ones(n)
x2 = dat.J
x3 = dat.J.^2
x4 = (dat.pitch ./ dat.diam)
x5 = (dat.pitch ./ dat.diam).^2
x6 = (dat.pitch ./ dat.diam).*dat.J
xt = [x1 x2 x3 x4 x5 x6]; ct = xt'xt\(xt'dat.CT);
##
scatter(dat.J, dat.CT, marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch")
scatter!(dat.J, xt*ct, marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch", marker=:+)
##
n = length(dat.J)
x1 = ones(n)
x2 = dat.J
x3 = dat.J.^2
x4 = (dat.pitch ./ dat.diam)
x5 = (dat.pitch ./ dat.diam).^2
x6 = (dat.pitch ./ dat.diam).*dat.J
xp = [x1 x2 x3 x4 x5 x6]; cp = xp'xp\(xp'dat.CP);
##
scatter(dat.J, dat.CP, marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch/diameter")
scatter!(dat.J, xp*cp, marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch/diameter", marker=:+)
##
save(joinpath(@__DIR__, "propaeromodel.jld2"),"cp",cp,"ct",ct,"pitch_diam_bounds",(0.5,1.0))
##
scatter(dat.J, dat.CT .* dat.J ./ dat.CP, marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch/diameter")
scatter!(dat.J, xt*ct .* dat.J ./ (xp*cp), marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch/diameter", marker=:+, ylims=(0.,1.))


##
scatter(dat.J, dat.CT .* dat.J ./ dat.CP, marker_z=dat.pitch./dat.diam, label="", colorbar_title="pitch", ylims=(0.,1.))


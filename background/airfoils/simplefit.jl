## RG15 Data at Re = 60400
import Plots
using Printf
const Pl=Plots
using JLD2
#=
Data from UIUC attached files
=#
Re =60000
# Re = (59900+60300+60100)/3
α = [[-6.50, -4.96, -3.43, -1.93, -0.39, 1.15, 2.72, 4.26, 5.80, 7.31, 8.84],
     [-6.53,-4.98,-3.46,-1.92,-0.41, 1.12, 2.70, 4.23, 5.76, 7.27, 8.80],
     [-8.56, -7.01, -5.46, -3.97, -2.42, -0.89, 0.66, 2.25, 3.79, 5.30, 6.80]]

cl =[[-0.346, -0.235, -0.129, -0.026, 0.129, 0.305, 0.573, 0.705, 0.844, 0.950, 1.027],
[-0.316, -0.143, -0.025, 0.134, 0.286, 0.415, 0.660, 0.796, 0.892, 0.990, 1.061],
[-0.309, -0.109, 0.050, 0.114, 0.251, 0.410, 0.590, 0.915, 1.085, 1.146, 1.205]]

cd =[[0.0330, 0.0227, 0.0145, 0.0145, 0.0158, 0.0215, 0.0211, 0.0210, 0.0214, 0.0257, 0.0371],
[0.0300, 0.0232, 0.0156, 0.0164, 0.0190, 0.0242, 0.0232, 0.0251, 0.0293, 0.0365, 0.0460],
[0.0510, 0.0290, 0.0227, 0.0215, 0.0238, 0.0252, 0.0285, 0.0296, 0.0252, 0.0315, 0.0417]]

cm = [-0.058, -0.11, -0.13];
data = [(; Re, α, cl, cd, cm)]

Re= 100000
# Re = (100100+99800+99600)/3
α = [[-6.48, -4.97, -3.45, -1.93, -0.37, 1.18, 2.72, 4.26, 5.79, 7.32, 8.83],
     [-6.48, -4.98, -3.46, -1.93, -0.37, 1.20, 2.73, 4.26, 5.78, 7.30, 8.82],
     [-8.55, -6.99, -5.48, -3.95, -2.41, -0.87, 0.66, 2.25, 3.78, 5.30, 6.79]]

cl =[[-0.366, -0.250, -0.129, -0.025, 0.170, 0.435, 0.587, 0.728, 0.867, 0.977, 1.047],
     [-0.252, -0.094, 0.023, 0.161, 0.358, 0.633, 0.807, 0.921, 1.023, 1.090, 1.156],
     [-0.279, -0.091, 0.038, 0.137, 0.286, 0.458, 0.664, 0.929, 1.038, 1.117, 1.170]]

cd =[[0.0304, 0.0217, 0.0150, 0.0128, 0.0164, 0.0157, 0.0148, 0.0158, 0.0189, 0.0236, 0.0362],
     [0.0276, 0.0224, 0.0166, 0.0148, 0.0180, 0.0177, 0.0144, 0.0178, 0.0232, 0.0312, 0.0454],
     [0.0539, 0.0298, 0.0198, 0.0204, 0.0194, 0.0201, 0.0253, 0.0177, 0.0215, 0.0266, 0.03617]]

push!(data, (; Re, α, cl, cd, cm))

Re = 200000
# Re = (200000+199400+199200)/3
α = [[-6.49, -4.99, -3.45, -1.90, -0.36, 1.19, 2.72, 4.26, 5.80, 7.34, 8.86],
     [-6.44, -4.96, -3.37, -0.61, -0.34, 1.19, 2.76, 4.26, 5.79, 7.31],
     [-8.50, -6.89, -5.53, -3.90, -2.32, -0.83, 0.81, 2.23, 3.74, 5.29, 6.81]]

cl =[[-0.387, -0.267, -0.124, 0.069, 0.244, 0.395, 0.553, 0.724, 0.867, 0.978, 1.059],
     [-0.091, 0.033, 0.162, 0.464, 0.498, 0.670, 0.818, 0.947, 1.051, 1.128],
     [-0.090, 0.029, 0.143, 0.277, 0.425, 0.581, 0.837, 0.954, 1.059, 1.141, 1.208]]

cd =[[0.0291, 0.0180, 0.0124, 0.0096, 0.0088, 0.0090, 0.0100, 0.0124, 0.0154, 0.0218, 0.0315],
     [0.0183, 0.0154, 0.0121, 0.0111, 0.0103, 0.0101, 0.0120, 0.0155, 0.0206, 0.0286],
     [0.0266, 0.0206, 0.0175, 0.0172, 0.0177, 0.0182, 0.0104, 0.0134, 0.0181, 0.0251, 0.0333]]

push!(data, (; Re, α, cl, cd, cm))

Re =300000
# Re = (300000+299100+298700)/3
α =[[-6.57, -4.98, -3.56, -1.91, -0.37, 1.19, 2.75, 4.27, 5.78, 7.30, 8.85],
    [-6.47, -4.97, -3.42, -2.04, -0.40, 1.19, 2.68, 4.23, 5.76, 7.32, 8.81],
    [-7.48, -5.93, -4.41, -2.89, -1.33, 0.23, 1.71, 3.30, 4.67, 6.29, 7.89]]

cl =[[-0.390, -0.246, -0.084, 0.090, 0.227, 0.429, 0.595, 0.744, 0.880, 0.991, 1.077],
     [-0.094, 0.042, 0.191, 0.330, 0.489, 0.652, 0.808, 0.933, 1.030, 1.110, 1.177],
     [0.022, 0.177, 0.301, 0.407, 0.646, 0.810, 0.930, 1.044, 1.122, 1.198, 1.262]]

cd =[[0.0290, 0.0153, 0.0113, 0.0078, 0.0069, 0.0077, 0.0089, 0.0110, 0.0143, 0.0205, 0.0293],
     [0.0171, 0.0124, 0.0094, 0.0086, 0.0099, 0.0105, 0.0115, 0.0142, 0.0190, 0.0261, 0.0368],
     [0.0173, 0.0137, 0.0146, 0.0158, 0.0131, 0.0104, 0.0117, 0.0154, 0.0206, 0.0282, 0.0394]]

push!(data, (; Re, α, cl, cd, cm))
##
P = plot()
for D=data
    plot!(P,D.cl[1], D.cd[1], label="")
end
P
## Cm versus δ
function simplefit(dat)
    (; cl, cd, α, Re) = dat
    δ =  [0, 5, 10.]
    δ_ = [δ*pi/180  [1., 1., 1]]
    cmδ, cm0 = (δ_'*δ_) \ (δ_'*cm) # least squares
    println("cmδ = $cmδ, cm0 = $cm0 at Re=$Re.")
    str = @sprintf "Fit of Pitching Moment versus Flap Deflection \n at Re=%6d cm0 = %2.4f, cmδ = %2.4f" Re cm0 cmδ
    
    p1=Pl.scatter(δ,cm, marker=:+, label="")
    Pl.plot!(δ,cm0 .+ cmδ * δ*pi/180, label="")
    Pl.xlabel!("Flap deflection (deg)")
    Pl.ylabel!("cm")
    Pl.title!(str,titlefontsize=12)
    Pl.savefig("UIUC_data/rg15_re$(Re)_cm.pdf")

    ## Cl versus α and δ
    nδ = 3
    δ = [0, 5, 10.]
    nα = length.(α)

    δα = hcat([fill(δ[1], nα[1]); fill(δ[2], nα[2]); fill(δ[3], nα[3])]*pi/180,
        vcat(α...)*pi/180,
        ones(sum(nα)))

    clδ, clα, cl0 = (δα'*δα) \ (δα'*vcat(cl...)) # least squares
    α0 = cl0 / clα
    println("clδ = $clδ, clα = $clα, cl0 = $cl0, α0 = $α0 at Re=$Re.")

    p2 = Pl.plot()
    colors =[:blue, :red, :green]
    for i=1:nδ
        Pl.scatter!(α[i],cl[i], marker=:+, label="", color=colors[i])
        Pl.plot!(α[i],cl0 .+ clδ * δ[i]*pi/180 .+ clα * α[i]*pi/180, label="δ=$(δ[i])",color=colors[i],legend=:bottomright)
     end
     Pl.xlabel!("Angle of attack (deg)")
     Pl.ylabel!("cl")
     str = @sprintf "Fit of Lift Coefficient Versus AoA and \n Flap Deflection at Re=%6d \n cl0 = %2.4f, clα = %2.4f, clδ = %2.4f" Re cl0 clα clδ
     Pl.title!(str,titlefontsize=12)
    Pl.savefig("UIUC_data/rg15_re$(Re)_cl.pdf")

    ## cd versus cl
    cl_ = hcat(vcat(cl...).^2, vcat(cl...), ones(sum(nα)))
    cd2, cd1, cd0 = (cl_'*cl_) \ (cl_'*vcat(cd...)) # least squares
    println("cd2 = $cd2, cd1 = $cd1, cd0 = $cd0 at Re=$Re.")
    str = @sprintf "Fit of Drag Versus Lift Coefficients \n at Re=%6d cd2 = %2.4f, cd1 = %2.4f, cd0 = %2.4f" Re cd2 cd1 cd0
    p3 = Pl.plot()
    colors =[:blue, :red, :green]
    for i=1:nδ
        Pl.scatter!(cl[i], cd[i], marker=:+, label="", color=colors[i])
        Pl.plot!(cl[i],cd0 .+ cd1 * cl[i] .+ cd2 * cl[i].^2, label="δ=$(δ[i])",color=colors[i],legend=:bottomright)
     end
     Pl.xlabel!("cl")
     Pl.ylabel!("cd")
     Pl.title!(str,titlefontsize=12)
    Pl.savefig("UIUC_data/rg15_re$(Re)_cd.pdf")
    return (; cmδ, cm0, cd2, cd1, cd0, clδ, clα, cl0, α0, Re)
end
rg15fit = simplefit.(data)
##
function refit(coefs, key)
     c = getindex.(coefs, key)
     re = getindex.(coefs, :Re)
     nre = length(re)
     dre = hcat(re, ones(nre))
     c1, c0 = (dre'*dre)\(dre'*c)
     # @show (c'*c)\(c'*dre)
     Pl.scatter(re, c, marker=:+, label="")
     Pl.plot!(re,c1.*re.+c0, label="")
     Pl.ylabel!("$key")
     Pl.xlabel!("Re")
     str = @sprintf "Fit of %s versus Reynold's Number \n %s=%.3ef Re + %.3ef" key key c1 c0
     Pl.title!(str,titlefontsize=12)
     Pl.savefig("UIUC_data/rg15_$(key)_re.pdf")
     return c, c0, c1
end
coefs = (:cd0, :cd1, :cd2, :clα, :clδ, :cl0, :cm0, :cmδ)
rg15fits = NamedTuple{coefs}(map(c->refit(rg15fit, c),coefs))
@save "UIUC_data/rg15.jld2" rg15fits


##
# Results
# cmδ = -0.4125296124941925, cm0 = -0.06333333333333337 at Re=60400.
# clδ = 1.947700976189669, clα = 5.602383974802183, cl0 = 0.20386090626164127, α0 = 0.03638824243010575 at Re=60400.
# cd2 = 0.028300487439960303, cd1 = -0.018621740679992036, cd0 = 0.02239894904987368 at Re=60400.

##

# DRG data

αd = [-6.48, -4.97, -3.45, -1.93, -0.37, 1.18, 2.72, 4.26, 5.79, 7.32, 8.83]
cld =[-0.366, -0.250, -0.129, -0.025, 0.170, 0.435, 0.587, 0.728, 0.867, 0.977, 1.047]
cd =[0.0304, 0.0217, 0.0150, 0.0128, 0.0164, 0.0157, 0.0148, 0.0158, 0.0189, 0.0236, 0.0362]
plot(αd, cld, marker=:o,label="")
xlabel!("α (deg)")
ylabel!("cl")
title!("RG15 Lift Coefficient at Re 100,000")
savefig("UIUC_data/available_data_comp_rg15_cl_100000.pdf")

plot(cld, cd, marker=:o,label="")
xlabel!("cl")
ylabel!("cd")
title!("RG15 Drag Coefficient at Re 100,000")
savefig("UIUC_data/available_data_comp_rg15_cd_100000.pdf")

## LFT data
D = [-6.57  -0.426  -0.0580
    -5.01  -0.262  -0.0580
    -4.03  -0.185  -0.0580
    -3.05  -0.112  -0.0580
    -2.07  -0.049  -0.0580
    -1.09   0.052  -0.0580
    0.08   0.252  -0.0580
    1.06   0.387  -0.0580
    2.04   0.483  -0.0580
    3.01   0.570  -0.0580
    4.26   0.675  -0.0580
    5.18   0.762  -0.0580
    6.14   0.843  -0.0580
    7.14   0.921  -0.0580
    8.09   0.967  -0.0580
    9.06   1.016  -0.0580
    10.27   1.062  -0.0580
    11.12   1.055  -0.0580
    12.13   1.013  -0.0580
    13.05   0.991  -0.0580
    14.27   0.943  -0.0580
    15.21   0.951  -0.0580
    16.29   0.941  -0.0580
    17.25   0.944  -0.0580
    18.20   0.977  -0.0580
    16.40   0.969  -0.0580
    15.56   0.947  -0.0580
    14.56   0.973  -0.0580
    13.48   0.942  -0.0580
    12.48   0.994  -0.0580
    11.49   1.078  -0.0580
    10.53   1.072  -0.0580
    9.68   1.036  -0.0580
    8.48   0.976  -0.0580
    7.55   0.932  -0.0580
    6.50   0.865  -0.0580
    5.43   0.771  -0.0580
    4.59   0.698  -0.0580
    3.36   0.594  -0.0580
    2.33   0.510  -0.0580
    1.34   0.417  -0.0580
    0.32   0.289  -0.0580
    -0.69   0.117  -0.0580
    -1.65  -0.013  -0.0580
    -2.62  -0.084  -0.0580
    -3.82  -0.179  -0.0580
    -4.86  -0.260  -0.0580
    -5.90  -0.343  -0.0580
    -6.74  -0.409  -0.0580]
α = D[:,1]; cl = D[:,2]; cm = D[:,3];
plot(α, cl, marker=:o,label="")
xlabel!("α (deg)")
ylabel!("cl")
title!("RG15 Lift Coefficient at Re 100,000")
savefig("UIUC_data/available_data_comp_rg15_clcm_100000.pdf")

## Together
plot(α, cl, marker=:o,label="")
plot!(αd, cld, marker=:o,label="")
xlabel!("α (deg)")
ylabel!("cl")
title!("RG15 Lift Coefficient at Re 100,000")
savefig("UIUC_data/available_data_comp_rg15_both_100000.pdf")
##
P = plot()
for d=data
    plot!(d.cl[3], d.cd[3], label="$(d.Re)", marker=:+, legend=:bottomleft)
end
P

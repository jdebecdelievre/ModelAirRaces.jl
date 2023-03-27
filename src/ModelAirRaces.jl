module ModelAirRaces

using StaticArrays
using ForwardDiff
const FD = ForwardDiff
using JLD2
using LinearAlgebra
using FastGaussQuadrature
using OptimUtils
using LatinHypercubeSampling

const envt = (
      rho_atm         = 1.225,           # kg/m^3
      g               = 9.807,           # N/kg
      nu              = 1.46E-5,         # m^2/s
      mu              = 1.789E-5,        # kg/m/s
    # battery_energy  = 117000.0,
      marathonlength  = 21.0975*1e3,
      defaultaltitude = 50.,
      batteryenergy   = 61e3, #2.2*3*3.7 *3600 * 0.7, # Joules (Ah*nCells*Vcell*3600)
      nturns          = 40,
      sidelength      = 21.0975e3/40/3,
      turnradius      = 21.0975e3/40/3/cosd(30)/2,

    length_propulsion = 0.065 # motor length with margin for wires and mount

)
const deriv = (:_, :α, :p, :q, :r)
const r2d = 180/pi
const d2r = pi/180
export d2r, r2d

# Airfoils fit
include("airframe/airfoilproperties.jl")
export RG15, NACA0009, FlatPlate, Airfoil, liftcoefficient, momentcoefficient, dragcoefficient, striptheory, geom
export SimplerRG15, SimplerNACA0009, SimplerFlatPlate

# Weissinger Lifting Line
include("airframe/weissingerresults.jl")
include("airframe/weissinger.jl")
export envt, WeissingerResults, WeissingerResultsAD, weissinger!, deriv, AbstractWeissingerResults, get_array, get_results, Output, pack!, primalize
export DragModel2d, weiss2vlm, wing_tail_yspacing!, wing_tail_yspacing, elementwise_cl, quad_wing_tail_yspacing!

# Strutural analysis routines
include("airframe/structuralanalysis.jl")
export cantilevered_beam_bending!, mass_estimation, updatecrosssection!, cantilevered_beam_torsion!, CircularCrossSection, RectangularCrossSection, npanels, inertia_estimation
export StructMesh, MaterialProperties, showplane

# Airframe subspace (combines aero and structures routines)
include("airframe/fixedfuselage.jl")
include("airframe/airframeanalysis.jl")
include("airframe/airframesubspace.jl")
export mass_estimation, aeroanalysis, bendinganalysis, fixedfuselage

# Propulsion subspace
include("propulsion/propulsionanalysis.jl")
include("propulsion/propulsionsubspace.jl")
export LinearMotorQuadraticPropeller, FixedEfficiencyProp, propulsion_analysis

# Trajectory subspace
include("trajectory/simpletrajanalysis.jl")
include("trajectory/simpletrajsubspace.jl")
export flattraj, flattraj_constraints, showsimpletraj, Airplane, airplane_fit!, TrajElement

end

# ModelAirRaces

Model Air Races aircraft design problem. All the modeling details can be found in the following PhD thesis:

*Aerospace vehicle design with Bayesian collaborative optimization*, Jean de Becdelièvre, 2023. Stanford University.

The package is currently not registered, so you will need to clone the repo in your `.julia/dev/` repo, open a REPL,
access the pacakge manager with `]`, and run:
`dev ModelAirRaces`.

Similarly, you will need to clone the OptimUtils package from `https://github.com/jdebecdelievre/OptimUtils.jl`
and run `dev OptimUtils`.

To run the optimization example using IPOPT,
run the `idf.jl` file under the `examples` directory.
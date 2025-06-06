```@meta
CurrentModule = AMOCStochasticBoxModel
```

# AMOCStochasticBoxModel

Documentation for [AMOCStochasticBoxModel](https://github.com/UCL/AMOCStochasticBoxModel.jl).

Julia implementation of _Atlantic meridional overturning circulation_ (AMOC) stochastic box model described in 
[Soons, Grafke & Dijkstra (2024)](https://doi.org/10.1175/JPO-D-23-0234.1).
This is a stochastic extension of the five compartment AMOC box model described in 
[Wood et al. (2019)](https://doi.org/10.1007/s00382-019-04956-1).

## Usage example

To simulate the model using the default parameters, matching those in Soons, Grafke & Dijkstra (2024), 
and plot a visualisation of the solution, the following snippet can be run

```@example
using AMOCStochasticBoxModel
using Plots # hide

parameters = AMOCStochasticBoxModelParameters()
solution = simulate_model(; parameters)
plot_solution(solution)
savefig("example-solution.svg"); nothing # hide
```
This produces the plot below with 

- three-dimensional projections of non-dimensional state trajectories in top panel,
- simulated box salinities over time in second panel, 
- simulated AMOC strength over time in third panel,
- and simulated stochastic freshwater forcing over time in bottom panel.

![Example module solution](example-solution.svg)

The AMOC strength can be seen to be stochastically switching between low and high levels,
corresponding to bistable states of the model, 
with the stochastic freshwater forcing causing the model state to transition between the
bistable states.

## API reference

```@index
```

```@autodocs
Modules = [AMOCStochasticBoxModel]
```

## References

1. Soons, J., Grafke, T., & Dijkstra, H. A. (2024). Optimal transition paths for AMOC
   collapse and recovery in a stochastic box model. _Journal of Physical Oceanography_,
   54(12), 2537-2552.
2. Wood, R. A., Rodríguez, J. M., Smith, R. S., Jackson, L. C., & Hawkins, E. (2019). 
   Observable, low-order dynamical controls on thresholds of the Atlantic meridional
   overturning circulation. _Climate Dynamics_, 53, 6A815-6834.

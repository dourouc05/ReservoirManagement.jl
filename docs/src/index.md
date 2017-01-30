# ReservoirManagement Documentation

ReservoirManagement.jl is a [Julia](http://julialang.org/) package that can be used to manage reservoirs (such as dams) using mathematical optimisation tools. Its design should allow easy extensibility for new use cases and kinds of management. For now, it focuses on determining rule curves, based on inflow to the reservoir and purpose needs.

These rule curves are determined based on mathematical optimisation models to derive the minimum level in the reservoir that is allowable. The solutions are largely impacted by the uncertainty in the inflow, which is taken into account in two ways:

  * stochastic programming, with scenarios
  * robust programming, with an uncertainty set

The number of available scenarios is often quite limited, as precise measures have not always been performed for a long period of time. To this end, *scenario generation* techniques can be used: instead of using the shear scenarios, some algorithms are applied on them to generate a larger set of scenarios; these techniques break some intra-year correlations, and their results are often more extreme scenarios that what is available in the data set. Robust programming is another way to work around this lack of data: they are only used to derive some confidence intervals for the inflow, using statistical modelling; if the confidence level is very high, then the intervals are automatically more conservative.

Once the rule curves are determined by optimisation, this package proposes different ways to evaluate these solutions and to plot them. This part is fully integrated within a report-generation facility that also allows comparing the solutions when some parameters vary.

## Installation

This package relies on mathematical optimisation, and thus requires a solver to be installed. For most features, a simple LP solver is enough (in rare cases, a MIP solver is required). Any solver that implements the [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) interface works ([see JuliaOpt's website for more information](http://www.juliaopt.org/)).

**Caution:** For now, Gurobi is hardcoded into the source code, occurrences of `GurobiSolver()` must be replaced by the corresponding solver object. This will be addressed in a near release.

For some statistical operations, this package uses unpublished Julia packages, which thus need to be installed separately (namely, [ExtremeValueDistributions.jl](https://github.com/sammorris81/ExtremeValueDistributions.jl/tree/develop) for fitting extreme value distributions). Execute those commands from a Julia shell:

```Julia
Pkg.clone("https://github.com/sammorris81/MetropolisUpdaters.jl.git")
Pkg.clone("https://github.com/sammorris81/DataTransformations.jl.git")
Pkg.clone("https://github.com/sammorris81/ExtremeValueDistributions.jl.git")
Pkg.add("Optim")

Pkg.checkout("ExtremeValueDistributions", "develop")
```

## Using this package

The first step is to import it, along with SIUnits to provide the units for physical parameters.

```Julia
using ReservoirManagement
using SIUnits
using SIUnits.ShortUnits
```

Then, the first step is to build a `Reservoir` data structure that contains all the relevant parameters about the reservoir, its purposes and its tributaries. After, this object can be used to call the optimisation solvers, or the integrated report generator.

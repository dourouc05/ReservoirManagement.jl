# ReservoirManagement

[![Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive) [![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](http://opensource.org/licenses/MIT)

ReservoirManagement.jl is a [Julia](http://julialang.org/) package that can be used to manage reservoirs (such as dams) using mathematical optimisation tools. Its design should allow easy extensibility for new use cases and kinds of management. For now, it focuses on determining rule curves, based on inflow to the reservoir and purpose needs. 

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

## Usage

```Julia
using ReservoirManagement
using SIUnits
using SIUnits.ShortUnits
```

The first step is to define the topology of the reservoir: the tributaries. This package defines two kinds of rivers: those that directly flow into the reservoir (the most common case), represented by `NaturalRiver`, and diverted rivers (whose contribution can be decided up to some level), `DivertedRiver`. 

The main information stored for each river is a series of inflow scenarios, each scenario being represented by an array of time series (whose type is [`TimeArray`](https://github.com/JuliaStats/TimeSeries.jl/)). 

For example, the [Vesdre reservoir](https://en.wikipedia.org/wiki/Lake_Eupen) has two main tributaries (the rivers Vesdre and Getzbach), plus one diverted river (Helle, which has a minimum environmental flow that must remain in the river after the diversion): 

```Julia
Vesdre = NaturalRiver(name="Vesdre", scenarios=scenariosVesdre)
Getzbach = NaturalRiver(name="Getzbach", scenarios=scenariosGetzbach)
Helle = DivertedRiver(name="Vesdre", scenarios=scenariosHelle, environmental_flow=0.5m^3/s, maximum_flow=15.0m^3/s)
```

Once this topology is modelled, a reservoir can be defined. The main parameters are: 

  * the purposes of the reservoir (i.e. its need for water); for this example, mostly drinking water
  * the outputs from the reservoir (i.e. how much water it can release); for this example, mostly the bottom outlets
  * the minimum and maximum allowable levels

```Julia
purpose = DeterministicPurposes(drinkingWater=0.5m^3/s)
out = ConstantDamOutputs(bottomOutlets=100m^3/s)
VesdreReservoir = Reservoir(name="Vesdre", capacity=(2500000.m^3, 25000000.m^3),
                            purposes=purpose, outputs=out, 
                            rivers_in=[Vesdre, Getzbach], rivers_diverted=[Helle])
```

The setup phase is done: with this `Reservoir` data structure, this package can provide optimisation results with the `report()` function. For now, it focuses on the determination of minimum rule curves (that span over `ruleDuration`) that guarantee the purposes for some duration (`guarantee`). If provided, it can compare its results to the current rule curve (`currentRuleCurve`). The results are mainly a series of plots (enabled by `savePlots=true` and saved into `saveFolder`), but also the intermediate results (enabled by `saveIntermediate=true` and also saved into `saveFolder`)

```Julia
report(VesdreReservoir, ruleDuration=Year(1), guarantee=Year(2), 
       currentRuleCurve=current_target,
       saveFolder="…/outputs/", savePlots=true, saveIntermediate=true)
```

**Remark:** The values in this example are not the true ones (albeit they are realistic). 

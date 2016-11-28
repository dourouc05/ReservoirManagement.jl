["""

For time series: `meta` corresponds to a function in `Base.Dates` (https://github.com/JuliaLang/julia/blob/master/base/Dates.jl#L20).

"""]


### For the article:
###     CODE     PLOTS     PLOT DATA     ARTICLE
###     X        X         X             X           Various CI levels (95..99: steps of 1 or .5)
###                                                  Analysis by month, week: differences? Go for two-three months, seasons, days?
###                                                  Mixing with seasons instead of years.
###                                                  Safe solver: evolution of the objective function when considering more and more scenarios. A bit like https://github.com/JuliaSmoothOptimizers/BenchmarkProfiles.jl? 
###                                                  Rework the engine so that the previous is possible: work by configurable blocks of scenarios (between one block for all scenarios, as is done currently, to one scenario per block); may then add parallelism in the code.
###    X         X         X            X            Basic solver: compare solutions with scenarios, indicate which scenarios are support (solid vs dashed)
###    X         X         X            X            Basic solver: compare solutions with scenarios, as large as an A4 page (legend on the side); completely rework the legend (two lines, a gradient)
### High priority:
###     TODO REFACTOR [REPORT]: rework evaluation data structures (evaluations and globalEvaluations) to allow easy sorting on the evaluator algorithm (and no ugly loop over loop with test on the evaluator to check that it is output at the right time). I.e. have them being Dict(evaluator -> Dict(parameters -> evaluation)).
###     TODO REFACTOR [REPORT]: enable configuration of output folders, especially for plots. For example, give one "main" output folder, allow to create a subfolder for the current report() call (include many parameters in the folder names, such as time steps?). Allow (or force?) creating subfolders for the different kinds of plots?                Try to be generic: specify the directory structure by an abstract data structure, with a specific interface to allow the user to implement their own folder structure (getFolder(:plots), getFolder(:plots, :confidenceIntervals), etc.).            Allow to add automatically the date within the outermost folder.                   Or rather use hash maps instead of full-blown data structures? (Or a main data structure containing mainly a hashmap, so that dispatch is still possible, without having to specify many parameters, i.e. the parameters variable and the folder chooser?)
###     TODO REFACTOR [REPORT]: should be runnable in multiple passes: one to compute things, serialise the results; one that prints a textual report on the console (and in a file); one that prepares the figures. All three could be run in one shot, without reading the results from disk, but the postprocessing parts should also work when reading the data from a serialised save.
### Short-term:
###     TODO REFACTOR: Migration to Julia 0.5.
###     TODO DOC [SOLVERS]: describe generic options for the solvers once and for all (doc).
###     TODO DESIGN [SOLVER_DS]: complete the null pattern for solver data structures (those zero objects should behave as real objects, but empty). I.e. follow documentation to provide basic implementation of each object, so that ony those fields are present (and contain nothing).
###     TODO REFACTOR [UTILS] [UPSTREAM]: teach TimeSeries.collapse() how to work with Base.Dates.Period. (Only used in computeCI_sub.)
###     TODO DESIGN [SOLVER_EVALUATION]: the minLevel/maxLevel variables for the evaluation solvers do not need to be InverseSolution, do they? Generalise to anything that has a `getLevel` method? (Currently doing strange things for hydropower evaluation.)
###     TODO REFACTOR [SOLVERS]: how to use another solver than Gurobi?
###     TODO MODELLING [SOLVER_HYDROPOWER]: hydropower: what about pumping? Modify Reservoir constructor to allow no river, but still pumping-storage (e.g.: Plate Taille, one of the lake of Eau d'Heure).     Use the same data structure as for hydropower, with convenience constructors to keep compatibility? This would be a PumpingStoragePurpose, with variable efficiency depending on whether it is pumping or producing hydropower. If this set of parameters is null (nullable sub-type for each part?), then not present (e.g. pumping).
###     TODO REFACTOR [EVERYWHERE]: Refactor periods handling? No more multiplying, dividing a bit everywhere, please... Function to go from a period to another? What about sub-day values?      Extend period2seconds, e.g.
###     TODO REFACTOR [SOLVER_EVALUATION]: Get rid of minLevel/maxLevel arguments and replace them by purposes? Would remove quite a few quirks (reduce the number of arguments, remove quirks in purposeshortage). Makes sense: a constraint the reservoir must meet.                   Requirement, even if the purpose way is not chosen: helpers to implement the rule curves.
###     TODO REFACTOR [SOLVER_DATASTRUCTURES]: Could get rid of the feasibility field in ReservoirSolution? Some evaluations do not need it (mc_feasibility and purposeshortage). Or just say it should stay to true in those cases?
###     TODO REFACTOR [UNITS HANDLING]: Switch to something more developed than SIUnits, like http://ajkeller34.github.io/Unitful.jl/. Requires migrating to Julia 0.5.                        The new units package should allow to refactor some timing issues (i.e. converting Hour(1) to either "1 hour" or "1 h", or even to "hour" and "h", which are often used to generate names).
### Medium-term:
###     TODO DEBUG [STATS]: Bring Adam2015 back to working state.
###     TODO FUNCTIONALITY [SOLVERS]: Outputting model files should only take a prefix as argument (force it to be a folder?): the solver decides the filename and the extension. Allows a single solver to output multiple models! Required for evaluate_purposeshortage             Nice functions: 1/ take a solver ID (as Symbol) and an optional supplement to create a file name (extension as positional parameter); 2/ takes a model, a solver ID, an optional supplement, but no extension, and makes the LP file on its own.
###     TODO FUNCTIONALITY [PURPOSES]: Purposes whose needs depend on time and time step. (E.g., discharge required for fish spawn only at some periods of the year.)
###     TODO FUNCTIONALITY [PURPOSES]: Other purposes: bounds and max deltas for reservoir level within the reservoir (e.g. recreational use); variations within a day IFF THE OPTIMISATION TIME STEP ALLOWS.
###     TODO DESIGN [SOLVERS] [WIP]: how to refactor the models? Many constraints and other variable definitions are basically copied from a file to another!
###     TODO DESIGN [SOLVER_DATASTRUCTURES]: remove variant descriptions and replace by a list of symbols describing the chain to this point? (could be easier to generate strings for plots)
###     TODO DESIGN [EVERYWHERE]: replace verbosity by Lumberjack.jl or Logging.jl, or a better journaling system (akin to Ipopt's Journalist http://web.mit.edu/ipopt_v3.11.6/Ipopt-doxydoc-3.11.6/classIpopt_1_1Journalist.html). Would solve indenting problems too in the outputs. Maybe write a package for this.        Or just use simple logger, with a function wrapping around info()? Redirect solver's output: https://github.com/JuliaOpt/JuMP.jl/issues/325#issuecomment-64974955        What about sub-solver calls? When some output stream is enabled for the main solver, should the sub-solvers also enable it?           How to deal with many iterations for a solver? A parameter to let output every so often? A discretised list of possible frequencies? (What if the user asks for a total number of iterations that is much higher than initially thought of?) Much better: what duration between the messages?      Use some reservoir level parameter when calling the log() function? Informational messages will have a high values; called solvers will have some increment wrt to the calling solver, so fewer messages are displayed for them. Have beginFunction() or so to automatically increase the base figure, so that each solver only deals with its own messages? Exactly like CPLEX or Gurobi: when solving a MIP, almost no feedback on the root node (except if long); but if only a LP of the same size, then there will be much more output from the simplex.              Similar to https://github.com/timholy/ProgressMeter.jl partly?
###     TODO DESIGN [SOLVER_DATASTRUCTURES]: define operator [] on reservoirs to select scenarios and time steps?
###     TODO FUNCTIONALITY [PLOTTING]: refactor the plotting facilities around the notion of recipe? https://juliaplots.github.io/recipes/
###     TODO REFACTOR [REPORT]: allow the user to add things externally, and avoid always doing the Cartesian product of all inputs. Idea: the main report() function takes a vector of analyses to run (may be generated by another function that does the Cartesian product).
###     TODO FUNCTIONALITY [SOLVERS]: add parameters for hydropower, so that the maximum discharge is limited by the hydraulic head just like the other outputs.
###     TODO FUNCTIONALITY [SOLVER_SHORTAGE]: new index, the PDSI (Palmer drought severity index). Def: http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450(1984)023%3C1100%3ATPDSIL%3E2.0.CO%3B2, sec 2. Original definition: https://www.ncdc.noaa.gov/temp-and-precip/drought/docs/palmer.pdf. Relevant for our purposes?
### Long-term:
###     TODO FUNCTIONALITY [PLOTTING]: develop a Web interface based on Escher.jl https://github.com/shashi/Escher.jl (or R's Shiny, or something else).
###     TODO DESIGN [EVERYWHERE]: think about units everywhere.
###     TODO MODELLING [SOLVER_HELPERS]: uncertainty for the purpose needs?
###     TODO FUNCTIONALITY [DATASTRUCTURES]: implement basins, purposes over a basin rather than only one dam.
###     TODO FUNCTIONALITY [STATS]: implement other methods for confidence intervals, e.g. based on Monte-Carlo approach with resampling.
###     TODO FUNCTIONALITY [SOLVERS]: think about goal programming to order the purposes (see RTC-tools https://oss.deltares.nl/web/rtc-tools). Impact on the way to solve problems, but also to evaluate solutions (e.g. shortage index).
###     TODO REFACTOR [REPORT]: write a package to write arrays to the console (a bit like PHP's http://symfony.com/doc/current/components/console/helpers/table.html), that can also output other formats (like CSV or HTML/CALS tables, or even LaTeX tables? LyX can import CVS "spreadsheets"); different output formats for the table (compact, nicer, etc.). This will allow report() to have a better output. Also include a caption for the table!          To decide the output format, reimplement show/display/writemime/...? Will those functions stand the test of time (https://github.com/JuliaLang/julia/issues/14052)?

using TimeSeries
using Base.Dates
using DataStructures
using SIUnits
using SIUnits.ShortUnits

using JuMP
using Gurobi

using StatsBase
using HypothesisTests
using Distributions
using ExtremeValueDistributions # Strange dependency (not published in METADATA), try to reimplement somewhere
# (and make it fully compatible with Julia 0.4, without warnings).

using DataFrames
using Colors
using Plots
using StatPlots

# The generalised distributions are defined in both Distributions and ExtremeValueDistributions, but those are identical.
GeneralizedExtremeValue = Distributions.GeneralizedExtremeValue
GeneralizedPareto = Distributions.GeneralizedPareto

# Conflict between standard Dates and SIUnits.
Second = Dates.Second

# To allow overriding those functions.
import Base: copy, isempty, start, next, done, eltype, length

# Order of inclusion is important!
include("utils.jl")
include("reservoir/_ds.jl")
include("read.jl")
include("stats/_stats.jl")
include("scenariosgeneration.jl")
include("solvers/_solvers.jl")
include("plots_dataframe.jl")
include("plots_plots.jl")
include("report.jl")

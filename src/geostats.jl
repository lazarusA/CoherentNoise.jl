using Meshes: Point
using Meshes: GridTopology
using Meshes: CartesianGrid
using Meshes: domain, topology
using Meshes: centroid, coordinates
using GeoStatsBase

import GeoStatsBase: preprocess, solvesingle

"""
    CNS(var₁=>param₁, var₂=>param₂, ...)

Corehent noise simulation.

## Parameters

* `sampler` - Coherent noise sampler

### References

* MISSING REFERENCE
"""
@simsolver CNS begin
  @param sampler
end

function preprocess(problem::SimulationProblem, solver::CNS)
  # retrieve domain of simulation
  pdomain = domain(problem)

  # sanity check
  @assert topology(pdomain) isa GridTopology "simulation only defined for domains with grid topology"

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine algorithm for variable
      sampler = varparams.sampler

      preproc[var] = (sampler=sampler,)
    end
  end

  preproc
end

function solvesingle(problem::SimulationProblem, covars::NamedTuple, ::CNS, preproc)
  # retrieve domain of simulation
  pdomain = domain(problem)

  # define corresponding Cartesian grid with
  # corners at (0,0,...) and (1,1,...)
  dims   = size(topology(pdomain))
  ndim   = length(dims)
  start  = Point(ntuple(i->0.0, ndim))
  finish = Point(ntuple(i->1.0, ndim))
  grid   = CartesianGrid(start, finish, dims=dims)

  # loop over simulation variables
  varreal = map(covars.names) do var
    # unpack preprocessed parameters
    sampler = preproc[var].sampler

    # perform simulation
    vals = map(1:prod(dims)) do i
      coords = coordinates(centroid(grid, i))
      sample(sampler, coords...)
    end

    # save result
    var => vals
  end

  Dict(varreal)
end
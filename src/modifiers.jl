using ..Common: AbstractSampler, RandomState, Seed, random_state, lerp, cubic_interpolate, curve3
using ..Fractals: FBM
using ..Patterns: Constant
import ..Common: sample

"""
Supertype for all `N`-dimensional modifier samplers.
"""
abstract type ModifierSampler{N} <: AbstractSampler{N} end

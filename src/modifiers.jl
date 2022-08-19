using ..Common: AbstractSampler, RandomState, Seed, random_state, lerp, cubic_interpolate, curve3
using ..Fractals: FBMFractal
using ..Patterns: constant_1d
import ..Common: sample

"""
Supertype for all modifier samplers.
"""
abstract type ModifierSampler{N} <: AbstractSampler{N} end

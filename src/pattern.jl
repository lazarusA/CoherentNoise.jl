using ..Common: AbstractSampler, RandomState, Seed
import ..Common: sample

"""
Supertype for all pattern samplers.
"""
abstract type PatternSampler{N} <: AbstractSampler{N} end

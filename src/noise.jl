using ..Common: AbstractSampler

"""
Supertype for all `N`-dimensional noise samplers.
"""
abstract type NoiseSampler{N} <: AbstractSampler{N} end

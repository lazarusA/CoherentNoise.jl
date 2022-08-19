using ..Common: AbstractSampler

"""
Supertype for all noise algorithm samplers.
"""
abstract type NoiseSampler{N} <: AbstractSampler{N} end

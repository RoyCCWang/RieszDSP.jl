module RieszDSP

using LinearAlgebra
using FFTW

import Combinatorics

# constant values.
function twopi(::Type{Float32})::Float32
    return 6.2831855f0 #convert(T, 2*π)
end

function twopi(::Type{Float64})::Float64
    return 6.283185307179586 #convert(T, 2*π)
end

function twopi(::Type{T})::T where T <: AbstractFloat
    return convert(T, 2*π)
end

function onepi(::Type{Float32})::Float32
    return 3.1415927f0
end

function onepi(::Type{Float64})::Float64
    return 3.141592653589793
end

function onepi(::Type{T})::T where T <: AbstractFloat
    return convert(T, π)
end

include("steerable/wavelet.jl")
include("steerable/riesz_filters.jl")
include("steerable/monogenicanalysis.jl")
include("steerable/riesz_wavelet.jl")
include("steerable/front_end.jl")

include("sym_tensor.jl")
include("filters/dft_filters.jl")

export  getdefaultscale,
    # getRTfilters, gethigherorderRTfilters,
    # RieszAnalysis, RieszSynthesis,
    # RieszAnalysisLimited, RieszSynthesisLimited,
    # directionalHilbert, monogenicanalysis, instantfreq,
    # getprefilters,
    #waveletanalysis, waveletsynthesis,
    #getmonogenicwaveletanalysis
    rieszwaveletanalysis, #convert𝓡ψtoψ𝓡,
    rieszwaveletsynthesis #convert𝓡ψtoψ𝓡vectorfield,

end # module RieszFilters

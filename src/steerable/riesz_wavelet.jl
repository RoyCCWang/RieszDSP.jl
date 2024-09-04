# functions for combining Riesz transform with isotropic, redundant wavelet
#   transform. Current implementation does Riesz transform first,
#   then redundant wavelet.
#
# Current implementation recreates the wavelet and Riesz filters as needed,
#   as opposed to storing them.

#abstract type OperationMode end

#struct RieszWavelet <: OperationMode end


"""
    rieszwaveletanalysis(y::Array{T}, N_scales::Int) where T <: AbstractFloat
    
Return types:
- `Vector{Vector{Array{T,D}}}`
- `Array{T<:AbstractFloat,D}`

returns WRY, residual.
WRY[d][s][n], where:
- d is dimension index, up to D.
- s is scale index, up to N_scales.
- n is sampling position index, for a D-dim array.

`residual` is the portion of `y` that doesn't under go the transform. Keep this if you want to do synthesis.

# Example: forward Riesz transform (rieszwaveletanalysis) and its inverse transform (rieszwaveletsynthesis)
```julia
import RieszDSP as RZ
using LinearAlgebra

# generate input data.
T = Float64
y = randn(T, 45, 512, 3)

# specification for the number of wavelet subbands.
N_scales = round(Int, log2( maximum(size(y))))

# forward transform.
WRY, residual = RZ.rieszwaveletanalysis(y, N_scales)
# residual is the portion that does not under go the Riesz-wavelet transform.

# inverse transform.
yr = RZ.rieszwaveletsynthesis(WRY, residual)

# this is a perfect reconstruction (up to numerical precision), if `residual` is kept.
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()
```
"""
function rieszwaveletanalysis(y::Array{T}, N_scales::Int) where T <: AbstractFloat

    A, r, _ = rieszwaveletanalysis(y, N_scales, 1)
    return A, r
end



"""
rieszwaveletanalysis(
    y::Array{T,D},
    N_scales::Integer,
    order::Integer,
    ) where {T <: AbstractFloat,D}
    

returns WRY, residual, a_array.

WRY[d][s][n], where:
- d is dimension index, up to D.
- s is scale index, up to N_scales.
- n is sampling position index, for a D-dim array.

`residual` is the portion of `y` that doesn't under go the transform. Keep this if you want to do synthesis.

a_array: a 1-D array of D-tuples. Each tuple sums to `order`. The `d`-th entry of the tuple takes on an integer between 0 and `order`, and reflects the number of times the Riesz transform is applied in the `d`-th coordinate axis. In other words, `a_array` is the index set for the multi-indices of a `D`-dimensional order-`L` symmetric tensor.

# Example: higher-order Riesz-wavelet transform, round-trip.
```julia
import RieszDSP as RZ
using LinearAlgebra

# generate input data.
T = Float64
y = randn(T, 45, 512, 3)

# specification for the number of wavelet subbands.
N_scales = round(Int, log2( maximum(size(y))))


# Specify the Riesz transform order.
order = 7

WRY, residual, a_array = RZ.rieszwaveletanalysis(y, N_scales, order)

# each entry in a_array specifies the number of times the Riesz transform is iterated per dimension.
# Since our in put `y` is a 3-D array, each entry of `a_array` is a 3-tuple that sums to the order we specified.
println("The Riesz order states are:", a_array)

yr, a_array_rec = RZ.rieszwaveletsynthesis(WRY, residual, order)
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()

# The returned Riesz order states from the analysis and synthesis directions of the transform should be the same. The following should be zero.
@show norm(a_array - a_array_rec)
```
"""
function rieszwaveletanalysis(
    y::Array{T,D},
    N_scales::Integer,
    order::Integer,
    ) where {T <: AbstractFloat,D}

    LP, HP = getprefilters(T, Val(D), size(y))

    # pre-filter.
    Y = real.(ifft(fft(y).*LP)) # bandlimited version of y.
    residual = real.(ifft(fft(y).*HP))

    # Riesz transform.
    H, a_array = gethigherorderRTfilters(T, Val(D), size(Y), order)
    RY = RieszAnalysisLimited(Y,H)

    # wavelet transform.
    WRY = collect( waveletanalysis(RY[j], N_scales) for j in eachindex(RY) )

    return WRY, residual, a_array
end

# # convert data structures
# function convert𝓡ψtoψ𝓡(WRY::Vector{Vector{Array{T,D}}})::Vector{Vector{Array{T,D}}} where {T,D}
#     @assert !isempty(WRY)

#     N_scales = length(WRY[1])
#     ψRY = Vector{Vector{Array{T,D}}}(undef, N_scales)

#     for s = 1:N_scales
#         ψRY[s] = Vector{Array{T,D}}(undef, D)

#         for d = 1:D
#             ψRY[s][d] = WRY[d][s]
#         end
#     end

#     return ψRY
# end

# # convert data structures
# function convert𝓡ψtoψ𝓡vectorfield(WRY::Vector{Vector{Array{T,D}}})::Vector{Array{Vector{T},D}} where {T,D}
#     @assert !isempty(WRY)

#     N_scales = length(WRY[1])
#     ψRY = Array{Array{Vector{T},D}}(undef, N_scales)

#     for s = 1:N_scales
#         ψRY[s] = Array{Vector{T}}(undef, size(WRY[1][s]))

#         for i = 1:length(WRY[1][s])
#             ψRY[s][i] = Vector{T}(undef, D)

#             for d = 1:D
#                 ψRY[s][i][d] = WRY[d][s][i]
#             end
#         end
#     end

#     return ψRY
# end

"""
rieszwaveletsynthesis(
    WRY::Vector{Vector{Array{T,D}}},
    residual::Array{T,D},
    )::Array{T,D} where {T <: AbstractFloat, D}
    
See rieszwaveletanalysis for a description of the inputs.

returns the original input, y of type Array{T,D}.

# Example: first-order Riesz-wavelet transform, round-trip.
```julia
import RieszDSP as RZ
using LinearAlgebra

# generate input data.
T = Float64
y = randn(T, 45, 512, 3)

# specification for the number of wavelet subbands.
N_scales = round(Int, log2( maximum(size(y))))

# forward transform.
WRY, residual = RZ.rieszwaveletanalysis(y, N_scales)
# residual is the portion that does not under go the Riesz-wavelet transform.

# inverse transform.
yr = RZ.rieszwaveletsynthesis(WRY, residual)

# this is a perfect reconstruction (up to numerical precision), if `residual` is kept.
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()
```
"""
function rieszwaveletsynthesis(
    WRY::Vector{Vector{Array{T,D}}},
    residual::Array{T,D},
    )::Array{T,D} where {T <: AbstractFloat, D}

    A, _ = rieszwaveletsynthesis(WRY, residual, 1)
    return A
end


"""
rieszwaveletsynthesis(
    WRY::Vector{Vector{Array{T,D}}},
    residual::Array{T,D},
    )::Array{T,D} where {T <: AbstractFloat, D}

See rieszwaveletanalysis for a description of the inputs.

returns the original input, y of type Array{T,D}.
# Example: higher-order Riesz-wavelet transform, round-trip.
```julia
import RieszDSP as RZ
using LinearAlgebra

# generate input data.
T = Float64
y = randn(T, 45, 512, 3)

# specification for the number of wavelet subbands.
N_scales = round(Int, log2( maximum(size(y))))


# Specify the Riesz transform order.
order = 7

WRY, residual, a_array = RZ.rieszwaveletanalysis(y, N_scales, order)

# each entry in a_array specifies the number of times the Riesz transform is iterated per dimension.
# Since our in put `y` is a 3-D array, each entry of `a_array` is a 3-tuple that sums to the order we specified.
println("The Riesz order states are:", a_array)

yr, a_array_rec = RZ.rieszwaveletsynthesis(WRY, residual, order)
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()

# The returned Riesz order states from the analysis and synthesis directions of the transform should be the same. The following should be zero.
@show norm(a_array - a_array_rec)
```
"""
function rieszwaveletsynthesis(
    WRY::Vector{Vector{Array{T,D}}},
    residual::Array{T,D},
    order::Integer,
    ) where {T <: AbstractFloat, D}

    RY_rec = collect( waveletsynthesis(WRY[j]) for j in eachindex(WRY) )

    sz_Y = size(RY_rec[begin])
    @assert sz_Y == size(residual)

    H, a_array = gethigherorderRTfilters(T, Val(D), sz_Y, order)
    Y_rec = RieszSynthesisLimited(RY_rec, H)

    LP, HP = getprefilters(T, Val(D), size(Y_rec))
    A_rec = real.(ifft(fft(Y_rec).*LP + fft(residual).*HP))

    return A_rec, a_array
end

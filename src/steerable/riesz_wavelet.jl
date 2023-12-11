# functions for combining Riesz transform with isotropic, redundant wavelet
#   transform. Current implementation does Riesz transform first,
#   then redundant wavelet.
#
# Current implementation recreates the wavelet and Riesz filters as needed,
#   as opposed to storing them.

#abstract type OperationMode end

#struct RieszWavelet <: OperationMode end


"""
rieszwaveletanalysis(
    y::Array{T,D},
    N_scales::Int,
    )::Tuple{Vector{Vector{Array{T,D}}},Array{T<:AbstractFloat,D}}

returns WRY, residual.
WRY[d][s][n], where:
- d is dimension index, up to D.
- s is scale index, up to N_scales.
- n is sampling position index, for a D-dim array.
"""
function rieszwaveletanalysis(
    y::Array{T,D},
    N_scales::Int,
    ) where {T <: AbstractFloat,D}

    A, r, _ = rieszwaveletanalysis(y, N_scales, 1)
    return A, r
end

function rieszwaveletanalysis(
    y::Array{T,D},
    N_scales::Integer,
    order::Integer,
    ) where {T <: AbstractFloat,D}

    LP, HP = getprefilters(T, Val(D), size(y))
    Y = real.(ifft(fft(y).*LP)) # bandlimited version of y.
    residual = real.(ifft(fft(y).*HP))

    H, a_array = gethigherorderRTfilters(T, Val(D), size(Y), order)
    RY = RieszAnalysisLimited(Y,H)

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

function rieszwaveletsynthesis(
    WRY::Vector{Vector{Array{T,D}}},
    residual::Array{T,D},
    )::Array{T,D} where {T <: AbstractFloat, D}

    A, _ = rieszwaveletsynthesis(WRY, residual, 1)
    return A
end

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

# Only Simoncelli lowpass/highpass pairs are used.
# The wavelet analysis assumes an input signal f: ℝᴰ → ℝ.
# Since Simoncelli filters have only real components, and we're doing
#   linear filtering, the wavelet coefficients are also real-valued.

function getprefilters(
    ::Type{T},
    ::Val{D},
    sz_Y::Tuple,
    ) where {T,D}

    return getSimoncellifilters(T, Val(D), sz_Y, convert(T, 2))
end

# isotropic band reductions across scales.
# uses flters that are strictly real and symmetric.
function waveletanalysis(
    Y::Array{T,D},
    N_scales::Int,
    )::Vector{Array{T,D}} where {T,D}

    Ŷ = fft(Y)
    Ŷ_LP = copy(Ŷ)
    ψY = Vector{Array{T,D}}(undef, N_scales)

    for s = 1:(N_scales-1)
        LP, HP = getSimoncellifilters(T, Val(D), size(Y), convert(T, 1/2^(s-1)))

        ψY[s] = real.(ifft(Ŷ_LP.*HP))
        Ŷ_LP = Ŷ_LP.*LP
    end
    ψY[end] = real.(ifft(Ŷ_LP))

    return ψY
end


function waveletsynthesis(ψY::Vector{Array{T,D}})::Array{T,D} where {T,D}

    N_scales = length(ψY)
    ψŶ = fft(ψY[end]) # this is fft(ψY[N_scales])

    for s = (N_scales-1):-1:1
        LP, HP = getSimoncellifilters(T, Val(D), size(ψY[s]), convert(T, 1.0/2^(s-1)))

        ψŶ = ψŶ.*LP + fft(ψY[s]).*HP
    end
    Y = real.(ifft(ψŶ))

    return Y
end

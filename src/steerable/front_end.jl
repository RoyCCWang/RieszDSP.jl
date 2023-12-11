function getmonogenicwaveletanalysis(
    y::Array{T,D},
    N_scales = round(Int, log2( maximum(size(y)) )),
    scale_select = cld(N_scales,2),
    ) where {T <: AbstractFloat, D}

    # get Riesz-wavelet transform.
    WRY,residual = rieszwaveletanalysis(y, N_scales)

    # get wavelet bands. Need this for the real part of the monogenic signal.
    LP,HP = getprefilters(T, Val(D), size(y))
    Y = real.(ifft(fft(y).*LP)) # bandlimited version of y.
    residual = real.(ifft(fft(y).*HP))
    ψY = waveletanalysis(Y, N_scales)

    # convert data structure.
    ψRY = convert𝓡ψtoψ𝓡(WRY)

    HᵤWRY_ds = directionalHilbert(ψRY[scale_select])
    Aᵤ, ϕᵤ = monogenicanalysis(HᵤWRY_ds, ψY[scale_select])
    ζᵤ, ∇ϕᵤ = instantfreq(ϕᵤ, ψRY[scale_select], ψY[scale_select])

    ∇ϕᵤ_norm = collect( norm(∇ϕᵤ[i]) for i in eachindex(∇ϕᵤ) )

    return WRY, ψY, HᵤWRY_ds, Aᵤ, ϕᵤ, ζᵤ, ∇ϕᵤ, ∇ϕᵤ_norm
end

function Rieszreconstructiondemo(A::Array{T,D}) where {T <: AbstractFloat, D}

    LP,HP = getprefilters(T, Val(D), size(A))
    Y = real.(ifft(fft(A).*LP))
    residual = real.(ifft(fft(A).*HP))
    println("Demo for Riesz transform reconstruction.")
    println("Y is isotropically bandlimited version of A.")

    H = getRTfilters(T, D, size(Y))
    B = RieszAnalysisLimited(Y,H)
    Yr = RieszSynthesisLimited(B,H)
    println("Discard imaginary parts: discrepancy between Y and Yr: ", sum(abs.(Y-Yr)) )

    H = getRTfilters(T, D, size(A))
    B = RieszAnalysisLimited(A,H)
    Ar = RieszSynthesisLimited(B,H)
    println("Discard imaginary parts: discrepancy between A and Ar: ", sum(abs.(A-Ar)), ". This should not be zero for a non-bandlimited A." )

    H = getRTfilters(T, D, size(A))
    B = RieszAnalysis(A,H)
    Ar = RieszSynthesis(B,H)
    println("Discard nothing: discrepancy between A and Ar: ", sum(abs.(A-Ar)) )
    println()

    return residual
end

# for one frequency band.
function filterpairreconstructiondemo(A::Array{T,D}; N_tests = 100) where {T <: AbstractFloat, D}
    #(h,w,d) = size(A)

    total_discrepancy = zero(T)
    for n = 1:N_tests
        s = rand(T)
        LP, HP = getSimoncellifilters(T, Val(D), size(A), convert(T, 1/2^(s-1)))
        Y = real.(ifft(fft(A).*LP))
        residual = real.(ifft(fft(A).*HP))
        total_discrepancy += sum(abs.(ifft(fft(Y).*LP+fft(residual).*HP)-A))
    end
    println("Demo for Lowpass and highpass reconstruction.")
    println("Total reconstruction discrepancies of ", N_tests, " trials: ")
    println(total_discrepancy)
    println()
end

# for multiple frequency bands, which together makes a redundant wavelet analysis.
function redundantwaveletreconstructiondemo(A::Array{T,D}) where {T <: AbstractFloat, D}
    #(h,w,d) = size(A)

    LP,HP = getprefilters(T, Val(D), size(A))
    Y = real.(ifft(fft(A).*LP))
    residual = real.(ifft(fft(A).*HP))

    levels = log2( maximum(size(Y)) )
    N_scales = round(Int, levels)
    ψY = waveletanalysis(Y, N_scales)
    Yr = waveletsynthesis(ψY)

    Ar = real.(ifft(fft(Yr).*LP + fft(residual).*HP))
    discrepancy = sum(abs.(A-Ar))
    println("Demo for Redundant wavelet reconstruction.")
    println("discrepancy: ", discrepancy)
    println()

    return discrepancy
end

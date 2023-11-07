
Random.seed!(25)

T = Float64
#T = Float32 # there is a noticeable discrepancy with Float32.

y = randn(T, 45, 512, 3)

N_scales = round(Int, log2( maximum(size(y))))

println("Demo for Riesz-wavelet reconstruction.")

𝓡ψY, residual = RZ.rieszwaveletanalysis(y, N_scales)

yr = RZ.rieszwaveletsynthesis(𝓡ψY, residual)
println("discrepancy between y and yr: ", sum(abs.(y-yr)) )
println()

nothing
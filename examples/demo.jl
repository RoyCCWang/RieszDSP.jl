
Random.seed!(25)

T = Float64
#T = Float32 # there is a noticeable discrepancy with Float32.

y = randn(T, 45, 512, 3)

N_scales = round(Int, log2( maximum(size(y))))

println("Demo for Riesz-wavelet reconstruction.")

WRY, residual = RZ.rieszwaveletanalysis(y, N_scales)

yr = RZ.rieszwaveletsynthesis(WRY, residual)
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()

nothing
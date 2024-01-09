# RieszDSP.jl
This is an implementation of the higher-order Riesz-wavelet transform for grid-sampled signals; i.e.  the input is a multi-dimensional array like a grayscale image or video. The Riesz-wavelet is a perfect reconstruction transform, provided that a residual signal that the Riesz-wavelet transform does not operate on is saved.

A portion of the code here is based on the [Generalized Riesz-Wavelet Toolbox for Matlab](https://bigwww.epfl.ch/demo/steerable-wavelets/) authored by Nicolas Chenouard, Dimitri Van De Ville and Michael Unser. I made modifications so that the wavelet subbands are not critically sampled.

# Documentation
The [documentation](https://royccwang.github.io/RieszDSP.jl/) has further details and an image decomposition demo.

# Quick demo: round-trip discrepancy
Apply the forward and inverse transform and compare against the original input, `y`.

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
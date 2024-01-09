
Random.seed!(25)

T = Float64
#T = Float32 # there is a noticeable discrepancy with Float32.

y = randn(T, 45, 512, 3)

# specification for the number of wavelet subbands.
N_scales = round(Int, log2( maximum(size(y))))

# # First-order version.
println("Demo for Riesz-wavelet reconstruction.")

WRY, residual = RZ.rieszwaveletanalysis(y, N_scales)
#residual is the portion that does not under go the Riesz-wavelet transform.

# inverse transform.
yr = RZ.rieszwaveletsynthesis(WRY, residual)
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()

# # Higher-order version.

# Specify the Riesz transform order.
order = 7

println("Demo for $order-order Riesz-wavelet reconstruction.")

# forward transform.
WRY, residual, a_array = RZ.rieszwaveletanalysis(y, N_scales, order)

# each entry in a_array specifies the number of times the Riesz transform is iterated per dimension.
# Since our in put `y` is a 3-D array, each entry of `a_array` is a 3-tuple that sums to the order we specified.
println("The Riesz order states are:", a_array)

# inverse transform.
yr, a_array_rec = RZ.rieszwaveletsynthesis(WRY, residual, order)
println("relative discrepancy between y and yr: ", norm(y-yr)/norm(y) )
println()

# The returned Riesz order states from the analysis and synthesis directions of the transform should be the same. The following should be zero.
@show norm(a_array - a_array_rec)

nothing
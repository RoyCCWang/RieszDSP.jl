
function computefreqgrid(dummy::Array{T,D})::Vector{Vector{T}} where {T <: AbstractFloat ,D}

    sz_Y = size(dummy)

    ω_grid = Vector{Vector{T}}(undef, D)
    for d = 1:D
        N = sz_Y[d]
        step = convert(T, twopi(T)/N)
        half = div(N,2) # same thing as floor(N/2)

        if iseven(N)
            positive_freqs = collect(LinRange(zero(T), onepi(T), half+1))
            pop!(positive_freqs)

            negative_freqs = collect(LinRange(-onepi(T), -step, half))

            ω_grid[d] = vcat(positive_freqs, negative_freqs)
        else
            half_step = step/2

            positive_freqs = collect(LinRange(zero(T), onepi(T)-half_step, half+1))

            negative_freqs = collect(LinRange(-onepi(T)+half_step, zero(T), half+1))
            pop!(negative_freqs)

            ω_grid[d] = vcat(positive_freqs,negative_freqs)
        end
    end

    return ω_grid
end


# we want to move the band from [π/4,π/2] to freq_bound_factor*[π/4,π/2].
function getSimoncellifilters(
    dummy::Array{T,D},
    freq_bound_factor::T,
    )::Tuple{Array{Complex{T},D},Array{Complex{T},D}} where {T <: AbstractFloat, D}

    sz_Y = size(dummy)
    ω_grid = computefreqgrid(dummy)

    ###lowpass frequency rsp.
    LP = zeros(Complex{T}, sz_Y)
    #ω_norm = zeros(T,sz_Y)

    π½ = convert(T, onepi(T)/2)
    π¼ = convert(T, onepi(T)/4)

    multiindex_LUT = CartesianIndices(size(LP))

    # at non-DC frequencies.
    for linear_i in eachindex(LP)

        # compute the corresponding multi-linear index
        i = multiindex_LUT[linear_i]

        tuple_i = Tuple(i)
        if !mapreduce( x -> (x==1 ? true : false), &, tuple_i ) # skip DC.

            # get norm.
            ω_norm = zero(T)
            for dd = 1:D
                ω_norm += ω_grid[dd][i[dd]]^2
            end
            ω_norm = sqrt(ω_norm)

            # correct for the specified frequency band.
            # division by freq_bound_factor is the same as ω_grid[d] = ω_grid[d]/freq_bound_factor
            ω_norm = ω_norm/freq_bound_factor

            # store filter to array.
            if ω_norm <= π¼
                LP[linear_i] = Complex(one(T), zero(T))
            elseif π¼ < ω_norm <= π½
                LP[linear_i] = Complex( cos(π½ * log2(2*ω_norm/π½)) , zero(T))
            end

        end

    end

    # at DC frequency.
    LP[begin] = Complex( one(T), zero(T))

    ### high-pass frequency rsp. Designed to allow perfect reconstruction.
    HP = sqrt.(1 .- abs2.(LP))

    return LP, HP
end
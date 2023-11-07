# SymTensorTools

# type unstable due to inability to dispatch on D.
function getsubscriptarray(L::Int, D::Int)

    α_array = Vector{Vector{Int}}(undef,binomial(L+D-1,L))
    mapsymtensor!(α_array, L, D)

    return α_array
end

function mapsymtensor!(out::Vector{Vector{Int}}, L, D)::Nothing
    @assert length(out) == binomial(D+L-1,L)
    ν = zeros(Int,D)
    ν[end] = L
    ν_linear = 1
    for ν_linear = 1:length(out)
        # update the last index.
        ν[end] = L - sum(ν) + ν[end]

        # compute the outer product
        out[ν_linear] = copy(ν)

        # update indices. Start from the 2nd last index, work our way towards the first.
        updatesubscript!(ν,length(ν)-1, ν[end])
        ν_linear += 1
    end

    return nothing
end

# devised from
# N = 4 case: ν1 = 0:L, ν2 = 0:L-ν1, ν3 = 0:L-ν1-ν2, ν4 = L-ν1-ν2-ν3
# ν1 = 0:L, ν2 = 0:L-ν1, ν3 = 0:L-ν1-ν2, ν4 = L-ν1-ν2-ν3
#collect((ν1,ν2,ν3,ν4) for ν1 = 0:L for ν2 = 0:L-ν1 for ν3 = 0:L-ν1-ν2 for ν4 = L-ν1-ν2-ν3)
function updatesubscript!(ν,d, limit)::Nothing
    if d < 1
        return nothing
    end

    limit += ν[d]
    ν[d] += 1

    if ν[d] > limit
        ν[d] = 0
        updatesubscript!(ν,d-1,limit)
    end
    return nothing
end
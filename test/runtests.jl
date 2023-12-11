
import RieszDSP as RZ

import Random
Random.seed!(25)

using LinearAlgebra
using Test


@testset "order 1, round-trip" begin

    N_tests = 10
    Ts = [Float64; Float32]
    eps_tol_multiplier = 100

    for T in Ts
        rel_tol = eps(T)*eps_tol_multiplier

        for _ = 1:N_tests

            y = randn(T, 45, 512, 3)

            N_scales = round(Int, log2( maximum(size(y))))

            WRY, residual = RZ.rieszwaveletanalysis(y, N_scales)
            yr = RZ.rieszwaveletsynthesis(WRY, residual)

            @test norm(yr-y)/norm(y) < rel_tol

        end
    end
end

@testset "orders 2 to 10, round-trip" begin

    N_tests = 10
    Ts = [Float64; Float32]
    eps_tol_multiplier = 100
    orders = [2; 3; 4; 5; 6; 7; 8; 9; 10;]

    for T in Ts
        rel_tol = eps(T)*eps_tol_multiplier

        for order in orders
            for _ = 1:N_tests

                y = randn(T, 45, 72, 3)

                N_scales = round(Int, log2( maximum(size(y))))

                WRY, res, a_array = RZ.rieszwaveletanalysis(y, N_scales, order)
                y_rec, a_array2 = RZ.rieszwaveletsynthesis(WRY, res, order)

                @test norm(y_rec-y)/norm(y) < rel_tol

            end
        end
    end
end
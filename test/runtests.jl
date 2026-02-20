using Hydrogen
using AtomicLevels
using LinearAlgebra
using WignerSymbols
using PrettyTables
using Test

include("burgess.jl")

function test_approx_eq(a, b; on_fail::Union{Nothing,Function}=nothing, kwargs...)
    size(a) == size(b) || throw(DimensionMismatch("Cannot compare objects of sizes $(size(a)) and $(size(b))"))
    if !isapprox(a, b; kwargs...)
        @error "Approximate equality failed:"
        na = norm(a)
        nb = norm(b)
        Δ = norm(a-b)
        relΔ = Δ/max(na,nb)
        pretty_table(["|a|" na
                      "|b|" nb
                      "Abs. Δ" Δ
                      "Rel. Δ" relΔ],
                     show_header=false,
                     alignment=[:r,:l],tf=tf_borderless)
        isnothing(on_fail) || on_fail()
    end

    @test isapprox(a, b; kwargs...)
end

# The Lyman and Balmer series radial dipole moments are given by
# Eq. (63.4)

lyman²(n) = 2^8*n^7*(n-1)^(2n-5)/((n+1)^(2n+5))

balmer_s²(n) = 2^17*n^7*(n^2-1)*(n-2)^(2n-6)/((n+2)^(2n+6))

function balmer_p²(n,ℓ)
    if ℓ == 0
        2^15*n^9*(n-2)^(2n-6)/(3*(n+2)^(2n+6))
    elseif ℓ == 2
        2^19*n^9*(n^2-1)*(n-2)^(2n-7)/(3*(n+2)^(2n+7))
    end
end

function lyman_balmer²(o, n, ℓ)
    # Eq. (BS63.5)
    n == o.n && return (3/2*n*√(n^2-max(ℓ,o.ℓ)^2))^2
    n = float(n)
    if o == o"1s"
        lyman²(n)
    elseif o == o"2s"
        balmer_s²(n)
    elseif o == o"2p"
        balmer_p²(n, ℓ)
    end
end

@testset "Hydrogen.jl" begin
    @testset "Factorial ratios" begin
        for a = 1:6
            for b = 1:6
                @test Hydrogen.factorial_ratio(a,b) == factorial(a)/factorial(b)
            end
        end
    end

    @testset "Radial dipole moments" begin
        # We cannot compare with Table 13 of Bethe and Salpeter 1977,
        # since it provides so few figures and some of the values
        # appear to be wrong. Instead, we compare with the specific
        # formulas for the Lyman and Balmer series.

        for (j,(o,ℓ′)) in enumerate([(o"1s",1), (o"2s",1), (o"2p",0), (o"2p",2)])
            for n′ = 1+ℓ′:8
                test_approx_eq(radial_dipole_moment(o, Orbital(n′, ℓ′))^2,
                               lyman_balmer²(o, n′, ℓ′))
            end
        end
    end

    @testset "Photoionization dipole radial integrals" begin
        for n = 1:4
            channels = [Orbital(n, ℓ) => Orbital(:k, ℓ′)
                        for ℓ ∈ 0:n-1 for ℓ′ ∈ (ℓ == 0 ? (1,) : (ℓ-1,ℓ+1))]
            ref = [Burgess[ch] for ch in channels]
            local κ = ref[1].k
            g = Hydrogen.burgess_g(n, κ)
            # We test each number separately, since their magnitudes are quite
            # different.
            for (i,ref) in enumerate(ref)
                for (j,κ) in enumerate(κ)
                    test_approx_eq(g[i,j], ref.absg[j], rtol=3e-5,
                                   on_fail=() -> begin
                                       @info "Channel $(channels[i])" κ j
                                   end)
                end
            end
        end
    end

    @testset "Orbitals" begin
        # Error conditions:
        @test_throws ArgumentError non_relativistic_orbital(0, 0)
        @test_throws ArgumentError non_relativistic_orbital(-1, 0)
        @test_throws ArgumentError non_relativistic_orbital(1, 1)
        @test_throws ArgumentError non_relativistic_orbital(2, -1)
        @test_throws MethodError non_relativistic_orbital(2.5, -1)
        # Check values:
        P10 = non_relativistic_orbital(1, 0)
        @test P10(0) ≈ 0.0
    end
end

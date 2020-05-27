using Hydrogen
using AtomicLevels
using Test

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
                @test radial_dipole_moment(o, Orbital(n′, ℓ′))^2 ≈
                    lyman_balmer²(o, n′, ℓ′)
            end
        end
    end
end

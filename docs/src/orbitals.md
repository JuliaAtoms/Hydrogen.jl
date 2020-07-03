# Orbitals

[`non_relativistic_orbital`](@ref) is can be used to calculate the non-relativistic reduced
radial orbitals:

```@example
using Hydrogen, Plots; pyplot()
function plot_nros(n; rmax=70)
    rs = range(0, rmax; length=500)
    p = plot(legend=false, xaxis="r (a.u.)", yaxis="\$P_{n\\ell}(r)\$  for  \$n=$n\$")
    for ℓ in 0:n-1
        P = non_relativistic_orbital(n, ℓ)
        plot!(rs, P.(rs))
    end
    return p
end
plot([plot_nros(n) for n=1:4]..., layout=(4,1), size=(800, 1200))
```

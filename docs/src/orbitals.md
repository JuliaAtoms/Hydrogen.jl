# Orbitals

[`non_relativistic_orbital`](@ref) is can be used to calculate the non-relativistic reduced
radial orbitals:

```@setup
using Hydrogen, CairoMakie, LaTeXStrings
function plot_nros!(ax, n; rmax=70)
    rs = range(0, sqrt(rmax); length=500).^2

    for ℓ in 0:n-1
        P = non_relativistic_orbital(n, ℓ)
        lines!(rs, P.(rs))
    end
end
fig = Figure(size=(800,1200))
for n = 1:4
    last_plot = n == 4
    ax = Axis(fig[n,1],
              xlabel=last_plot ? L"$r$ [Bohr]" : "",
              xscale=sqrt,
              xticksvisible=last_plot, xticklabelsvisible=last_plot,
              ylabel=L"$P_{n\ell}(r)$ [Bohr$^{-1/2}$]",
              title=L"n=%$(n)")
    plot_nros!(ax, n)
end
save("radial-orbitals.svg", fig)
```

![](radial-orbitals.svg)

@doc raw"""
    factorial_ratio(a,b)

Compute

```math
\frac{a!}{b!}
```

via simple looping. This is efficient if ``\abs{a-b}`` is not too
large.
"""
function factorial_ratio(a,b)
    a < b && return inv(factorial_ratio(b,a))
    d = a-b
    isinteger(d) || throw(DomainError("factorial_ratio only valid for integer distances"))
    a == b && return one(a)

    f = a
    while a > b+1
        a -= 1
        f *= a
    end
    f
end

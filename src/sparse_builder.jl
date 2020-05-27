function sparse_builder(fun::Function, ::Type{T}, M, N;
                        issymmetric=false, ishermitian=false) where T
    I = Int[]
    J = Int[]
    V = T[]

    set! = if issymmetric || ishermitian && T <: Real
        (i,j,v) -> begin
            iszero(v) && return
            push!(I, i)
            push!(J, j)
            push!(V, v)

            push!(I, j)
            push!(J, i)
            push!(V, v)
        end
    elseif ishermitian
        (i,j,v) -> begin
            iszero(v) && return
            push!(I, i)
            push!(J, j)
            push!(V, v)

            push!(I, j)
            push!(J, i)
            push!(V, conj(v))
        end
    else
        (i,j,v) -> begin
            iszero(v) && return
            push!(I, i)
            push!(J, j)
            push!(V, v)
        end
    end

    fun(set!)

    sparse(I, J, V, M, N)
end

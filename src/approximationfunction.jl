export ApproximationFunction
"""
    ApproximationFunction

Type used to store a function defined on elements of basis and returning a Number.
"""
struct ApproximationFunction{BT<:MB.AbstractPolynomialBasis}
    func::Function
    vars::Vector{<:MP.AbstractVariable}
    basis_type::Type{BT}
end

function (f::ApproximationFunction)(p::AbstractPolynomialLike)
    basis = maxdegree_basis(f.basis_type, f.vars, maxdegree(p))
    coefs, basis = MB.change_basis(p, basis)
    return dot(coefs, f.func.(basis))
end

# Linear operations on ApproximationFunctions
function compatible(f1::ApproximationFunction, f2::ApproximationFunction)
    return f1.basis_type == f2.basis_type && f1.vars == f2.vars
end

function Base.sum(v::Vector{<:ApproximationFunction})
    f = first(v)
    for g in v[2:end]
        @assert compatible(f, g) "Only functions acting on same basis can be added."
    end
    return ApproximationFunction(
        x -> sum(g.func(x) for g in v),
        f.vars,
        f.basis_type,
    )
end

function Base.:+(f1::ApproximationFunction, f2::ApproximationFunction)
    return sum([f1, f2])
end
function Base.:*(a::Number, f::ApproximationFunction)
    return ApproximationFunction(x -> a * f.func(x), f.vars, f.basis_type)
end
Base.:*(f::ApproximationFunction, a::Number) = a * f
Base.:/(f::ApproximationFunction, a::Number) = inv(a) * f
Base.:-(f::ApproximationFunction) = (-1) * f
Base.:-(f1::ApproximationFunction, f2::ApproximationFunction) = f1 + (-f2)

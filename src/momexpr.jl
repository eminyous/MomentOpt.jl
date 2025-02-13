abstract type AbstractMomentExpressionLike <: AbstractGMPScalar end

function cover_variables(c::Number, e::GMPVariableRef) end
function cover_variables(c::MP.AbstractPolynomialLike, e::GMPVariableRef)
    @assert variables(c) ⊆ variables(e) "$e does not act on $c."
end

export MomentExpr
"""
    MomentExpr{T, C}

Type representing a linear combination of Moments.
"""
mutable struct MomentExpr{T<:Union{Number,MP.AbstractPolynomialLike}} <:
               AbstractMomentExpressionLike
    coefs::Vector{T}
    meass::Vector{GMPVariableRef}
    function MomentExpr(
        coefs::Vector{T},
        meass::Vector{GMPVariableRef},
    ) where {T}
        same_length(coefs, meass)
        for (c, m) in zip(coefs, meass)
            cover_variables(c, m)
        end
        return new{T}(coefs, meass)
    end
end

MP.coefficients(me::MomentExpr) = me.coefs
measures(me::MomentExpr) = me.meass

Base.zero(::MomentExpr{T}) where {T} = MomentExpr(T[], GMPVariableRef[])
Base.zero(::Type{MomentExpr{T}}) where {T} = MomentExpr(T[], GMPVariableRef[])

Base.iszero(me::MomentExpr) = length(me) == 0
function Base.:(==)(me1::MomentExpr, me2::MomentExpr)
    return all(coefficients(me1) .== coefficients(me2)) &&
           all(measures(me1) == measures(me2))
end

Base.eltype(::MomentExpr{T}) where {T} = Tuple{T,GMPVariableRef}
Base.length(me::MomentExpr) = length(coefficients(me))
function Base.iterate(me::MomentExpr)
    return iszero(me) ? nothing :
           ((first(coefficients(me)), first(measures(me))), 2)
end
function Base.iterate(me::MomentExpr, i)
    return length(me) < i ? nothing :
           ((coefficients(me)[i], measures(me)[i]), i + 1)
end

function momexp_by_measexp(me::MomentExpr)
    C = []
    M = []
    for (c, m) in me
        idx = findfirst(x -> x == c, C)
        if idx isa Nothing
            push!(C, c)
            push!(M, m)
        else
            M[idx] += m
        end
    end
    return C, M
end

function JuMP.function_string(mode, me::MomentExpr)
    str = ""
    if length(me) == 0
        str = "⟨0, 0⟩"
    else
        C, V = momexp_by_measexp(me)
        for (i, (c, v)) in enumerate(zip(C, V))
            if i == 1
                s = ""
            else
                s = " + "
            end
            str = str * s * "⟨" * sprint(show, c) * ", " * sprint(show, v) * "⟩"
        end
    end
    return str
end

function Base.getindex(me::MomentExpr{T}, mu::GMPVariableRef) where {T}
    idx = findall(x -> x == mu, measures(me))
    if isempty(idx)
        return zero(T)
    else
        return sum(coefficients(me)[idx])
    end
end

function Base.promote_rule(
    ::Type{MomentExpr{T}},
    ::Type{MomentExpr{S}},
) where {T,S}
    return MomentExpr{promote_type(T, S)}
end
function Base.convert(::Type{MomentExpr{T}}, m::MomentExpr) where {T}
    return MomentExpr(convert.(T, coefficients(m)), measures(m))
end

function Base.:*(a::Number, e::MomentExpr)
    return MomentExpr(a .* coefficients(e), measures(e))
end

function Base.sum(mev::Vector{MomentExpr{T}}) where {T}
    if T <: Number
        coefs = T[]
    else
        coefs = polynomial_type(T)[]
    end
    vars = GMPVariableRef[]
    for me in mev
        for (c, m) in me
            index = findfirst(x -> x == m, vars)
            if index isa Nothing
                push!(coefs, c)
                push!(vars, m)
            else
                coefs[index] += c
            end
        end
    end
    index = findall(x -> !iszero(x), coefs)
    return MomentExpr(coefs[index], vars[index])
end

function MomentExpr(
    cv::Vector{<:Union{Number,MP.AbstractPolynomialLike}},
    mv::Vector{<:MeasExpr},
)
    return sum([
        MomentExpr(c .* coefficients(m), measures(m)) for (c, m) in zip(cv, mv)
    ])
end
function MomentExpr(
    cv::Union{Number,MP.AbstractPolynomialLike},
    mv::Union{GMPVariableRef,MeasExpr},
)
    return MomentExpr([cv], [mv])
end

degree(e::MomentExpr) = maximum(MP.maxdegree.(coefficients(e)))

export Mom
const Mom = MomentExpr

mutable struct AffMomentExpr{S,T} <: AbstractMomentExpressionLike
    expr::MomentExpr{S}
    cons::T
end

expr(ame::AffMomentExpr) = ame.expr
JuMP.constant(ame::AffMomentExpr) = ame.cons

expr(ame::MomentExpr) = ame
JuMP.constant(ame::MomentExpr) = 0

MP.coefficients(ame::AffMomentExpr) = coefficients(expr(ame))
measures(ame::AffMomentExpr) = measures(expr(ame))

function Base.:(==)(me1::AffMomentExpr, me2::AffMomentExpr)
    return coefficients(me1) == coefficients(me2) &&
           measures(me1) == measures(me2)
end

function JuMP.function_string(mode, ame::AffMomentExpr)
    if constant(ame) == 0
        cs = ""
    else
        cs = " + " * string(JuMP.constant(ame))
    end
    return function_string(mode, expr(ame)) * cs
end

Base.:+(me::MomentExpr, meas::T) where {T<:Number} = AffMomentExpr(me, meas)
function Base.:+(me::AffMomentExpr, meas::T) where {T<:Number}
    return AffMomentExpr(expr(me), constant(me) + meas)
end
function Base.sum(amev::Vector{<:AffMomentExpr})
    return AffMomentExpr(sum(expr.(amev)), sum(constant.(amev)))
end
function Base.:-(me::Union{AffMomentExpr,MomentExpr}, meas::T) where {T<:Number}
    return me + (-meas)
end
function Base.:*(a::Number, me::AffMomentExpr)
    return AffMomentExpr(a * expr(me), a * constant(me))
end

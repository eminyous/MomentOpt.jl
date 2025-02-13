function SemialgebraicSets.inequalities(::SAS.AbstractAlgebraicSet)
    return []
end

function polynomials(s::SAS.AbstractSemialgebraicSet)
    return [equalities(s)..., inequalities(s)...]
end

function MultivariatePolynomials.maxdegree(s::SAS.AbstractSemialgebraicSet)
    pols = polynomials(s)
    return isempty(pols) ? 0 : maximum(maxdegree.(pols))
end

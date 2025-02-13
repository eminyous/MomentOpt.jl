"""
Approximating a discontinuous soltution to Burgers equation with GMP.jl
For theoretical background see : http://aimsciences.org/article/doi/10.3934/mcrf.2019032

Let f(y)= 1/4y^2, y_0(x) = 1 if x<0 and 0 if x>0. We want to find a solution y(t,x) for 

∂_t y(t,x) + ∂_x f(y(t,x)) = 0
y(0,x) = y_0(x)

on T×X = [0,1]×[-1/2,1/2]. In this example we can and will impose that y(-1/2,t)=1. In
addition, based on the initial condition, we can derive that y(t,x)∈[0,1] on T×X.

This problem can be equivalently defined as

find μ,μ_T,μ_R
<μ_T, φ(t,x)y> + <μ_R, φ(t,x)f(y)> - <μ,∂_tφ(t,x) y + ∂_xφ(t,x) f(y)> = ∫φ(0,y_0(x) dx + ∫φ(t,f(1)) dt
<μ,φ(t,x)> = ∫φ(t,x) dt dx
<μ_T,φ(x)> = ∫φ(x) dx
<μ_R,φ(t)> = ∫φ(t) dt     	
μ  ∈ M(T×X×Y)
μ_T∈ M({1}×X×Y))
μ_R∈ M(T×{-1/2}×Y))

where the intgral with respect to x is from -1/2 to 1/2 and the integral with respect to t
is from 0 to 1. 


From the moments of the measure μ the graph of the solution y(t,x) can be extracted via 
an approach described in : http://www.optimization-online.org/DB_HTML/2019/03/7146.html

"""

using DynamicPolynomials
using MomentOpt
using MosekTools

# Choose relaxation order
order = 6

flux(y) = 1 / 4 * y^2
T = [0, 1]
X = [-1 / 2, 1 / 2]
Y = [0, 1]

# Define polynomial variables
@polyvar t x y

gmp = GMPModel()
# Add measures

@variable gmp μ Meas(
    [t, x, y],
    support = @set(
        (t - T[1]) * (T[2] - t) >= 0 &&
        (x - X[1]) * (X[2] - x) >= 0 &&
        (y - Y[1]) * (Y[2] - y) >= 0
    )
)
@variable gmp μT Meas(
    [t, x, y],
    support = @set(
        t == 1 && (x - X[1]) * (X[2] - x) >= 0 && (y - Y[1]) * (Y[2] - y) >= 0
    )
)
@variable gmp μR Meas(
    [t, x, y],
    support = @set(
        (t - T[1]) * (T[2] - t) >= 0 &&
        x == 1 / 2 &&
        (y - Y[1]) * (Y[2] - y) >= 0
    )
)

# As T×X×Y is compact we can use the monomials as testfunctions φ. 
# We define the values for the right hand side of the constraint of the dynamic
rhs(i, j) = 0^i * (-(-1 / 2)^(j + 1)) / (j + 1) + (-1 / 2)^j / ((i + 1) * 4)

# DynamicPolynomials.jl provides the possibility to define a monomial vector.
mons = monomials([t, x], 0:2*order-2)
# The monomial vector is not of a polynomial type, we need to convert it first.
pons = polynomial.(mons)
wpde =
    vec(differentiate(pons, [t])) * y + vec(differentiate(pons, [x])) * flux(y)

# The fiels mons.Z provides the exponents of the monomials of mons (and pons)
@constraint gmp [i = 1:length(pons)] Mom(pons[i] * y, μT) +
                                     Mom(pons[i] * flux(y), μR) -
                                     Mom(wpde[i], μ) == rhs(mons.Z[i]...)

# Right hand sides for the marginal constraints
lebt(i) = (T[2]^(i + 1) - T[1]^(i + 1)) / (i + 1)
lebx(j) = (X[2]^(j + 1) - X[1]^(j + 1)) / (j + 1)
lebtx(i, j) = lebt(i) * lebx(j)
monstx = monomials([t, x], 0:2*order)
monsx = monomials(x, 0:2*order)
monst = monomials(t, 0:2*order)

@constraint gmp Mom.(polynomial.(monstx), μ) .==
                [lebtx(monstx.Z[i]...) for i in 1:length(monstx.Z)]
@constraint gmp Mom.(polynomial.(monsx), μT) .==
                [lebx(monsx.Z[i]...) for i in 1:length(monsx.Z)]
@constraint gmp Mom.(polynomial.(monst), μR) .==
                [lebt(monst.Z[i]...) for i in 1:length(monst.Z)]

# The GMP formulation of the PDE does not have any objective function as the measures are uniquely determined by the
# problem. However, when relaxing the problem the measures will not be unique anymore which is why we add an objective function.
# Note that a GMPModel always needs an objective to be relaxed.
# A common heuristic for an objective function is to minimize the trace of the moment matrix:

mmons = polynomial.(monomials([t, x, y], 0:order))
trace = mmons' * mmons

@objective gmp Min Mom(trace, μ) + Mom(trace, μT) + Mom(trace, μR)

set_approximation_mode(gmp, DUAL_STRENGTHEN_MODE())
set_optimizer(gmp, Mosek.Optimizer)
optimize!(gmp)

println(
    """
    The approximated function can be plotted, e.g., using the following commands:

    using Plots
    pyplot()

    T = range(0 , stop = 1, length = 10)
    X = range(-1/2, stop = 1/2, length = 10)
    f = graph(μ, [t,x]; regpar = 1e-6)
    plot(T, X, (tt, xx) -> f([tt, xx]) , st=:surface,camera=(-40,30))
    """,
)

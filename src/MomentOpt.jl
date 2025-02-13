module MomentOpt

using Reexport

import MathOptInterface as MOI

import MutableArithmetics as MA
using MutableArithmetics

using LinearAlgebra

import MultivariatePolynomials as MP
using MultivariatePolynomials
include("MPext/MPextra.jl")

Reexport.@reexport using JuMP

import SemialgebraicSets as SAS
Reexport.@reexport using SemialgebraicSets
include("SASext/SASextra.jl")

using SumOfSquares
export Sparsity

import MultivariateMoments as MM
Reexport.@reexport using MultivariateMoments


import MultivariateBases as MB
Reexport.@reexport using MultivariateBases

include("MBext/MBextra.jl")

abstract type AbstractGMPModel <: JuMP.AbstractModel end

include("approximationfunction.jl")
include("approximationscheme.jl")

include("objects.jl")
include("variables.jl")
include("defaultmeasures.jl")

include("measexpr.jl")
include("momexpr.jl")

abstract type GMPSubstitution <: JuMP.AbstractJuMPScalar end
# define MomentSubstitution

include("constraints.jl")
#
include("gmpmodel.jl")
include("MMext/MMextra.jl")
include("approximate.jl")
include("gmppostproc.jl")

end# module

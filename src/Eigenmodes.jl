#
#  Eigenmodes.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#

module Eigenmodes

# Julia packages
using LinearAlgebra
using SparseArrays
using FromFile

# Local modules
@from "Laplacians.jl" using Laplacians

function eigenmodesLf(R, A, B)
    Lf = makeLf(R, A, B)
    decomposition = (eigen(Matrix(Lf))).vectors
    return decomposition
end

function eigenmodesLc(R, A, B)
    Lc = makeLc(R, A, B)
    decomposition = (eigen(Matrix(Lc))).vectors
    return decomposition
end

function eigenmodesLt(R, A, B)
    Lₜ = makeLt(R, A, B)
    decomposition = (eigen(Matrix(Lₜ))).vectors
    return decomposition
end

function eigenmodesLv(R, A, B)
    Lᵥ = makeLv(R, A, B)
    decomposition = (eigen(Matrix(Lᵥ))).vectors
end 

export eigenmodesLf
export eigenmodesLc
export eigenmodesLt
export eigenmodesLv

end #end module 
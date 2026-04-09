#
#  Constants.jl
#  DiscreteCalculus
#
#

module Constants

using LinearAlgebra
using StaticArrays

const ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])

const ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])

export ϵᵢ
export ϵₖ

end
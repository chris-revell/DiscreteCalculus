#
#  CovarianDerivatives.jl
#  DiscreteCalculus
#
#
# Differential operators over primary and dual networks

module CovarianDerivatives

using DrWatson
using LinearAlgebra
using SparseArrays
using StaticArrays
using FromFile

@from "GeometryFunctions.jl" using GeometryFunctions
@from "TopologyFunctions.jl" using TopologyFunctions

function make𝐃c(R, A, B, 𝐯)
    𝐧ᵢⱼ = findCellOutwardNormals(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [outerProd(𝐧ᵢⱼ[i,j], 𝐯[j])/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

function make𝐃v(R, A, B, 𝐕ⱼ)
    𝐍ⱼₖ = findTriangleOutwardNormals(R, A, B)
    E = findCellLinkVertexTriangleAreas(R, A, B)
    tmp = [outerProd(𝐍ⱼₖ[j,k], 𝐕ⱼ[j])/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

export make𝐃c
export make𝐃v

end
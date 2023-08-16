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
using GeometryBasics
using DrWatson

# Local modules
@from "GeometryFunctions.jl" using GeometryFunctions
@from "Laplacians.jl" using Laplacians

function eigenmodesLt(nCells, nEdges, nVerts, R, A, Aᵀ, B, boundaryVertices, cellCentresOfMass, edgeMidpoints)    
    T = makeCellLinks(nCells, nEdges, B, cellCentresOfMass, edgeMidpoints)
    edgeTrapezia = makeEdgeTrapezia(nEdges, R, A, B, cellCentresOfMass)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(nCells, nVerts, R, A, B, C, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lₜ = makeLt(A,Aᵀ,T,linkTriangleAreas,trapeziumAreas)
    decomposition = (eigen(Matrix(Lₜ))).vectors
    return decomposition
end

function eigenmodesLf(nCells, nEdges, R, A, B, Bᵀ, C, cellAreas, edgeLengths, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    edgeTrapezia = makeEdgeTrapezia(nEdges, R, A, B, cellCentresOfMass)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    Lf = makeLf(nCells,B,Bᵀ,cellAreas,edgeLengths,trapeziumAreas)
    decomposition = (eigen(Matrix(Lf))).vectors
    return decomposition
end

function eigenmodesLv(nCells, nEdges, nVerts, R, A, Aᵀ, B, C, boundaryVertices, cellCentresOfMass, edgeMidpoints, edgeLengths)
    T = makeCellLinks(nCells, nEdges, B, cellCentresOfMass, edgeMidpoints)
    edgeTrapezia = makeEdgeTrapezia(nEdges, R, A, B, cellCentresOfMass)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(nCells, nVerts, R, A, B, C, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lᵥ = makeLv(A,Aᵀ,edgeLengths,linkTriangleAreas,trapeziumAreas)
    decomposition = (eigen(Matrix(Lᵥ))).vectors
end 

export eigenmodesLt, eigenmodesLf, eigenmodesLv

end #end module 
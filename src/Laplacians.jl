#
#  Laplacians.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#

module Laplacians

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using FromFile
using GeometryBasics

# Local modules
@from "SpatialData.jl" using SpatialData
@from "GeometryFunctions.jl" using GeometryFunctions

function makeLf(R, A, B)
    nCells = size(B,1)
    Bᵀ = sparse(transpose(B))
    cellAreas = findCellAreas(R, A, B)
    edgeLengths = findEdgeLengths(R, A) 
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    H = Diagonal(cellAreas)
    boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
    diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths.^2)./(2.0.*trapeziumAreas)))[:,1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
    Tₑ = Diagonal(diagonalComponent)
    invH = inv(H)
    Lf = invH*B*Tₑ*Bᵀ
    dropzeros!(Lf)    
    return Lf
end

function makeLc(R, A, B)
    nCells = size(B,1)
    Bᵀ = sparse(transpose(B))
    cellAreas = findCellAreas(R, A, B)
    T = makeCellLinks(R,A,B)
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    onesVec = ones(1,nCells)
    boundaryEdges = abs.(onesVec*B)
    boundaryEdgesFactor = abs.(boundaryEdges.-1)# =1 for internal vertices, =0 for boundary vertices
    H = Diagonal(cellAreas)
    Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
    invTₗ = inv(Tₗ)
    boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1,:])
    Lc = (H\B)*boundaryEdgesFactorMat*invTₗ*Bᵀ
    dropzeros!(Lc)
    return Lc
end

function makeLv(R, A, B)
    Aᵀ = sparse(transpose(A))
    edgeLengths = findEdgeLengths(R, A)
    linkTriangles = makeLinkTriangles(R, A, B)
    linkTriangleAreas = abs.(area.(linkTriangles))
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    E = Diagonal(linkTriangleAreas)
    Tₑ = Diagonal((edgeLengths.^2)./(2.0.*trapeziumAreas))
    Lᵥ = (E\Aᵀ)*(Tₑ\A)
    dropzeros!(Lᵥ)
    return Lᵥ
end

function makeLt(R, A, B)
    Aᵀ = sparse(transpose(A))
    T = makeCellLinks(R,A,B)
    linkTriangles = makeLinkTriangles(R, A, B)
    linkTriangleAreas = abs.(area.(linkTriangles))
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    E = Diagonal(linkTriangleAreas)
    Tₗ = Diagonal(((norm.(T)).^2)./(2.0.*trapeziumAreas))
    Lₜ = (E\Aᵀ)*Tₗ*A
    dropzeros!(Lₜ)
    return Lₜ
end

export makeLf
export makeLc
export makeLv
export makeLt

end #end module 
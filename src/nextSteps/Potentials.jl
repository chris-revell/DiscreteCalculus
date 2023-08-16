#
#  Potentials.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#

module Potentials

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using FromFile
using GeometryBasics

# Local modules
@from "GeometryFunctions.jl" using GeometryFunctions
@from "Laplacians.jl" using Laplacians

function psicPotential(nCells, nEdges, nVerts, R, A, B, C, F, boundaryVertices, cellCentresOfMass, cellAreas, edgeMidpoints, edgeTangents, ϵᵢ)
    T = makeCellLinks(nCells, nEdges, B, cellCentresOfMass, edgeMidpoints)
    edgeTrapezia = makeEdgeTrapezia(nEdges, R, A, B, cellCentresOfMass)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(nCells, nVerts, R, A, B, C, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lf = makeLf(nCells,B,Bᵀ,cellAreas,edgeLengths)
    cellDivs = -1.0.*calculateCellDivs(nCells, R, B, C, F, cellCentresOfMass, edgeMidpoints, edgeTangents, cellAreas, ϵᵢ)
    onesVec = ones(nCells)
    H = Diagonal(cellAreas)
    eigenvectors = (eigen(Matrix(Lf))).vectors
    eigenvalues = (eigen(Matrix(Lf))).values
    ḡ = ((onesVec'*H*cellDivs)/(onesVec'*H*ones(nCells))).*onesVec
    ğ = cellDivs.-ḡ
    ψ̆ = zeros(nCells)
    spectrum = Float64[]
    for k=2:nCells
        numerator = eigenvectors[:,k]'*H*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*H*eigenvectors[:,k])
        ψ̆ .-= (numerator/denominator).*eigenvectors[:,k]
        push!(spectrum,(numerator/denominator))
    end
    return ψ̆, spectrum
end

function psivPotential(nCells, nEdges, nVerts, R, A, Aᵀ, B, C, F, ϵᵢ, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    T = makeCellLinks(nCells, nEdges, B, cellCentresOfMass, edgeMidpoints)
    edgeTrapezia = makeEdgeTrapezia(nEdges, R, A, B, cellCentresOfMass)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(nCells, nVerts, R, A, B, C, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    linkTriangleAreas = abs.(area.(linkTriangles))
    q = makeSpokes(nVerts, nCells, R, C, cellCentresOfMass)
    Lₜ = makeLt(A,Aᵀ,T,linkTriangleAreas,trapeziumAreas)
    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values
    vertexDivs = -1.0.*calculateVertexDivs(nVerts, R, C, cellCentresOfMass, F, ϵᵢ, q, linkTriangleAreas)
    onesVec = ones(nVerts)
    E = Diagonal(linkTriangleAreas)
    ḡ = ((onesVec'*E*vertexDivs)/(onesVec'*E*ones(nVerts))).*onesVec
    ğ = vertexDivs.-ḡ
    ψ̆ = zeros(nVerts)
    spectrum = Float64[]
    for k=2:nVerts
        numerator = -eigenvectors[:,k]'*E*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(spectrum,(numerator/denominator))
    end
    return ψ̆, spectrum
end

function capitalPsivPotential(nCells, nEdges, nVerts, R, A, B, C, F, ϵᵢ, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    T = makeCellLinks(nCells, nEdges, B, cellCentresOfMass, edgeMidpoints)
    edgeTrapezia = makeEdgeTrapezia(nEdges, R, A, B, cellCentresOfMass)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    linkTriangles = makeLinkTriangles(nCells, nVerts, R, A, B, C, boundaryVertices, cellCentresOfMass, edgeMidpoints)
    linkTriangleAreas = abs.(area.(linkTriangles))
    Lₜ = makeLt(A,Aᵀ,T,linkTriangleAreas,trapeziumAreas)
    eigenvectors = (eigen(Matrix(Lₜ))).vectors
    eigenvalues = (eigen(Matrix(Lₜ))).values
    q = makeSpokes(nVerts, nCells, R, C, cellCentresOfMass)
    vertexCurls = calculateVertexCurls(nVerts, R, C, cellCentresOfMass, F, ϵᵢ, q, linkTriangleAreas)
    onesVec = ones(nVerts)
    E = Diagonal(linkTriangleAreas)
    ḡ = ((onesVec'*E*vertexCurls)/(onesVec'*E*ones(nVerts))).*onesVec
    ğ = vertexCurls.-ḡ
    ψ̆ = zeros(nVerts)
    spectrum = Float64[]
    for k=2:nVerts
        numerator = -eigenvectors[:,k]'*E*ğ
        denominator = eigenvalues[k]*(eigenvectors[:,k]'*E*eigenvectors[:,k])
        ψ̆ .+= (numerator/denominator).*eigenvectors[:,k]
        push!(spectrum,(numerator/denominator))
    end
    return ψ̆, spectrum
end

export psicPotential, psivPotential, capitalPsivPotential

end #end module 
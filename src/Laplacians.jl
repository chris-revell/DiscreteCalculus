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
using StaticArrays

# Local modules
@from "GeometryFunctions.jl" using GeometryFunctions
@from "TopologyFunctions.jl" using TopologyFunctions

function geometricLf(R, A, B)
    nCells = size(B, 1)
    Bᵀ = Transpose(B)
    cellAreas = findCellAreas(R, A, B)
    edgeLengths = findEdgeLengths(R, A)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    onesVec = ones(1, nCells)
    boundaryEdges = abs.(onesVec * B)
    H = Diagonal(cellAreas)
    boundaryEdgesFactor = abs.(boundaryEdges .- 1)# =1 for internal vertices, =0 for boundary vertices
    diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths .^ 2)./(2.0 .* trapeziumAreas)))[:, 1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
    Tₑ = Diagonal(diagonalComponent)
    invH = inv(H)
    Lf = invH * B * Tₑ * Bᵀ
    dropzeros!(Lf)
    return Lf
end

function geometricLc(R, A, B)
    nCells = size(B, 1)
    Bᵀ = Transpose(B)
    cellAreas = findCellAreas(R, A, B)
    T = findCellLinks(R, A, B)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    onesVec = ones(1, nCells)
    boundaryEdges = abs.(onesVec * B)
    boundaryEdgesFactor = abs.(boundaryEdges .- 1)# =1 for internal vertices, =0 for boundary vertices
    H = Diagonal(cellAreas)
    Tₗ = Diagonal(((norm.(T)) .^ 2) ./ (2.0 .* trapeziumAreas))
    invTₗ = inv(Tₗ)
    boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1, :])
    Lc = (H \ B) * boundaryEdgesFactorMat * invTₗ * Bᵀ
    dropzeros!(Lc)
    return Lc
end

function geometricLv(R, A, B)
    Aᵀ = Transpose(A)
    edgeLengths = findEdgeLengths(R, A)
    linkTriangles = findCellLinkTriangles(R, A, B)
    linkTriangleAreas = abs.(area.(linkTriangles))
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    E = Diagonal(linkTriangleAreas)
    Tₑ = Diagonal((edgeLengths .^ 2) ./ (2.0 .* trapeziumAreas))
    Lᵥ = (E \ Aᵀ) * (Tₑ \ A)
    dropzeros!(Lᵥ)
    return Lᵥ
end

function geometricLt(R, A, B)
    Aᵀ = Transpose(A)
    T = findCellLinks(R, A, B)
    linkTriangles = findCellLinkTriangles(R, A, B)
    linkTriangleAreas = abs.(area.(linkTriangles))
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    E = Diagonal(linkTriangleAreas)
    Tₗ = Diagonal(((norm.(T)) .^ 2) ./ (2.0 .* trapeziumAreas))
    Lₜ = (E \ Aᵀ) * Tₗ * A
    dropzeros!(Lₜ)
    return Lₜ
end

function topologicalLf(A, B)
    nCells = size(B, 1)
    Bᵀ = Transpose(B)
    onesVec = ones(1, nCells)
    boundaryEdges = abs.(onesVec * B)
    H = Diagonal(cellAreas)
    boundaryEdgesFactor = abs.(boundaryEdges .- 1)# =1 for internal vertices, =0 for boundary vertices
    diagonalComponent = boundaryEdgesFactor'
    Tₑ = Diagonal(diagonalComponent)
    invH = inv(H)
    Lf = invH * B * Tₑ * Bᵀ
    dropzeros!(Lf)
    return Lf
end

function topologicalLc(A, B)
    nCells = size(B, 1)
    Bᵀ = Transpose(B)
    onesVec = ones(1, nCells)
    boundaryEdges = abs.(onesVec * B)
    boundaryEdgesFactor = abs.(boundaryEdges .- 1)      # =1 for internal vertices, =0 for boundary vertices
    boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1, :])
    Lc = B * boundaryEdgesFactorMat * Bᵀ
    Lc = B * Bᵀ
    dropzeros!(Lc)
    return Lc
end

function topologicalLv(A, B)
    Aᵀ = Transpose(A)
    Lᵥ = Aᵀ * A
    dropzeros!(Lᵥ)
    return Lᵥ
end

function topologicalLt(A, B)
    Aᵀ = Transpose(A)
    Lₜ = Aᵀ * A
    dropzeros!(Lₜ)
    return Lₜ
end


function edgeMidpointLDirichlet(R, A, B)
    nCells = size(B,1); nEdges = size(B,2); nVerts = size(A,2)
    cellAreas = findCellAreas(R, A, B)
    edgeMidpointLinks = findEdgeMidpointLinks(R, A, B)
    edgeTangents = findEdgeTangents(R, A)
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2
    vertexAreas = findVertexAreas(R, A, B)
    L = spzeros(Float64,(nCells,nCells))
    for k=1:nVerts
        for i in findall(x->x!=0, C[:,k])
            for i′ in findall(x->x!=0, C[:,k])
                L[i, i′] += (edgeMidpointLinks[i,k]⋅edgeMidpointLinks[i′,k])/(cellAreas[i]*vertexAreas[k])
            end
        end
    end
    dropzeros!(L)
    return L
end


function edgeMidpointLNeumann(R, A, B)
    nCells = size(B,1); nEdges = size(B,2); nVerts = size(A,2)
    Aᵀ = Transpose(A)
    Āᵀ = abs.(Aᵀ)
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2
    cellAreas = findCellAreas(R, A, B)
    edgeMidpointLinks = findEdgeMidpointLinks(R, A, B)
    boundaryVertices = Āᵀ * abs.(sum.(eachcol(B))) .÷ 2
    for k in findall(x->x!=0, boundaryVertices)
        edgeMidpointLinks[:,k] .= fill(SVector{2,Float64}(zeros(2)),nCells)
    end
    edgeTangents = findEdgeTangents(R, A)
    vertexAreas = findVertexAreas(R, A, B)
    L = spzeros(Float64,(nCells,nCells))
    for k=1:nVerts
        for i in findall(x->x!=0, C[:,k])
            for i′ in findall(x->x!=0, C[:,k])
                L[i, i′] += (edgeMidpointLinks[i,k]⋅edgeMidpointLinks[i′,k])/(cellAreas[i]*vertexAreas[k])
            end
        end
    end
    dropzeros!(L)
    return L
end

#   Lⱼⱼ′=∑ₖAⱼₖ(t̂ⱼ⋅t̂ⱼ′)Aⱼ′ₖ/Eₖ
function scalarEdgeL(R, A, B)
    nCells = size(B,1); nEdges = size(B,2); nVerts = size(A,2)
    edgeTangentsNormalised = normalize.(findEdgeTangents(R, A))
    L = spzeros(nEdges,nEdges)
    vertexAreas = findVertexAreas(R, A, B)
    for k=1:nVerts
        k_js = findall(x->x!=0,A[:,k])
        for j in k_js
            for j′ in k_js 
                L[j,j′] += (edgeTangentsNormalised[j]⋅edgeTangentsNormalised[j′])*A[j,k]*A[j′,k]/vertexAreas[k]
            end
        end
    end
    dropzeros!(L)
    return L
end

# Lᵢᵢ′ = ∑ₖB̄ᵢⱼAⱼₖ(t̂ⱼ⋅t̂ⱼ′)(Aⱼ′ₖ/Eₖ)B̄ᵢ′ⱼ′
function uniformCellTensionL(R, A, B)
    nCells = size(B,1); nEdges = size(B,2); nVerts = size(A,2)
    edgeTangentsNormalised = normalize.(findEdgeTangents(R, A))
    L = spzeros(nCells,nCells)
    B̄ = abs.(B)
    C = findC(A,B)
    vertexAreas = findVertexAreas(R, A, B)
    for k=1:nVerts
        k_js = findall(x->x!=0,A[:,k])
        for j in k_js
            for j′ in k_js 
                j_is =findall(x->x!=0,B[:,j])#∩findall(x->x!=0,B[:,j′])
                j′_i′s =findall(x->x!=0,B[:,j′])#∩findall(x->x!=0,B[:,j′])
                for i in j_is
                    for i′ in j′_i′s 
                        L[i,i′] += (edgeTangentsNormalised[j]⋅edgeTangentsNormalised[j′])*A[j,k]*A[j′,k]*B̄[i,j]*B̄[i′,j′]/vertexAreas[k]
                    end
                end
            end
        end
    end
    dropzeros!(L)
    return L
end



# function edgeMidpointLfunction(ϕᵢ, R, A, B)
    
#     nCells = size(B,1); nEdges = size(B,2); nVerts = size(A,2)

#     cellAreas = findCellAreas(R, A, B)
    
#     sᵢₖ = findEdgeMidpointLinks(R, A, B)

#     edgeTangents = findEdgeTangents(R, A)
    
#     Ā = abs.(A)
#     B̄ = abs.(B)
#     C = B̄ * Ā .÷ 2

#     αₖ=Float64[]
#     for k=1:nVerts
#         k_is = findall(x->x!=0, C[:,k])
#         α = 0.5*norm([sᵢₖ[(k_is[1], k)]...,0.0]×[sᵢₖ[(k_is[2],k)]...,0.0])
#         push!(αₖ,α)
#     end

#     Lϕ = zeros(nCells)

#     for k=1:nVerts
#         for i in findall(x->x!=0, C[:,k])
#             for i′ in findall(x->x!=0, C[:,k])
#                 Lϕ[i] += (sᵢₖ[i,k]⋅sᵢₖ[i′,k])*ϕᵢ[i′]/(cellAreas[i]*αₖ[k])
#             end
#         end
#     end
    
#     return Lϕ
# end

export geometricLf
export geometricLc
export geometricLv
export geometricLt
export topologicalLf
export topologicalLc
export topologicalLv
export topologicalLt
export edgeMidpointLDirichlet
export edgeMidpointLNeumann
export scalarEdgeL
export uniformCellTensionL

end #end module 

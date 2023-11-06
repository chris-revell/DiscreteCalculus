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
@from "GeometryFunctions.jl" using GeometryFunctions

function geometricLf(R, A, B)
    nCells = size(B, 1)
    Bᵀ = Transpose(B)
    cellAreas = findCellAreas(R, A, B)
    edgeLengths = findEdgeLengths(R, A)
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
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
    T = makeCellLinks(R, A, B)
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
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
    linkTriangles = makeLinkTriangles(R, A, B)
    linkTriangleAreas = abs.(area.(linkTriangles))
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    E = Diagonal(linkTriangleAreas)
    Tₑ = Diagonal((edgeLengths .^ 2) ./ (2.0 .* trapeziumAreas))
    Lᵥ = (E \ Aᵀ) * (Tₑ \ A)
    dropzeros!(Lᵥ)
    return Lᵥ
end

function geometricLt(R, A, B)
    Aᵀ = Transpose(A)
    T = makeCellLinks(R, A, B)
    linkTriangles = makeLinkTriangles(R, A, B)
    linkTriangleAreas = abs.(area.(linkTriangles))
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    E = Diagonal(linkTriangleAreas)
    Tₗ = Diagonal(((norm.(T)) .^ 2) ./ (2.0 .* trapeziumAreas))
    Lₜ = (E \ Aᵀ) * Tₗ * A
    dropzeros!(Lₜ)
    return Lₜ
end

function topologicalLf(R, A, B)
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

function topologicalLc(R, A, B)
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

function topologicalLv(R, A, B)
    Aᵀ = Transpose(A)
    Lᵥ = Aᵀ * A
    dropzeros!(Lᵥ)
    return Lᵥ
end

function topologicalLt(R, A, B)
    Aᵀ = Transpose(A)
    Lₜ = Aᵀ * A
    dropzeros!(Lₜ)
    return Lₜ
end

function edgeMidpointL(R, A, B)
    nCells = size(B,1); nEdges = size(B,2); nVerts = size(A,2)

    cellAreas = findCellAreas(R, A, B)
    
    sᵢₖ = findEdgeMidpointLinks(R, A, B)

    edgeTangents = findEdgeTangents(R, A)
    
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2

    αₖ=Float64[]
    for k=1:nVerts
        k_is = findall(x->x!=0, C[:,k])
        # if length(k_is) == 2
        #     edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
        #     α = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
        #     edgesSharedBy_i2_And_k = findall(x->x!=0, B[k_is[2],:])∩findall(x->x!=0, A[:,k])
        #     α += 0.5^3*norm([edgeTangents[edgesSharedBy_i2_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i2_And_k[2]]...,0.0])
        # else
            α = 0.5*norm([sᵢₖ[(k_is[1], k)]...,0.0]×[sᵢₖ[(k_is[2],k)]...,0.0])
        # end
        push!(αₖ,α)
    end

    L = spzeros(Float64,(nCells,nCells))
    for k=1:nVerts
        for i in findall(x->x!=0, C[:,k])
            for i′ in findall(x->x!=0, C[:,k])
                L[i, i′] += (sᵢₖ[i,k]⋅sᵢₖ[i′,k])/(cellAreas[i]*αₖ[k])
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
export edgeMidpointL

# export edgeMidpointLfunction

end #end module 
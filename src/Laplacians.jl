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
@from "InnerOuterProduct.jl" using InnerOuterProduct

# Scalar Laplacian Lf with metric components 
function geometricLf(R, A, B)
    I = size(B, 1)
    Bᵀ = Transpose(B)
    cellAreas = findCellAreas(R, A, B)
    edgeLengths = findEdgeLengths(R, A)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    onesVec = ones(1, I)
    boundaryEdges = abs.(onesVec * B)
    H = Diagonal(cellAreas)
    boundaryEdgesFactor = abs.(boundaryEdges .- 1)# =1 for internal edges, =0 for boundary edges
    diagonalComponent = (boundaryEdgesFactor'.*((edgeLengths .^ 2)./(2.0 .* trapeziumAreas)))[:, 1] # Multiply by boundaryEdgesFactor vector to set boundary vertex contributions to zero
    Tₑ = Diagonal(diagonalComponent)
    invH = inv(H)
    Lf = invH * B * Tₑ * Bᵀ
    dropzeros!(Lf)
    return Lf
end

# Scalar Laplacian Lc with metric components 
function geometricLc(R, A, B)
    I = size(B, 1)
    Bᵀ = Transpose(B)
    cellAreas = findCellAreas(R, A, B)
    T = findCellLinks(R, A, B)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    onesVec = ones(1, I)
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

# Scalar Laplacian Lv with metric components 
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

# Scalar Laplacian Lt with metric components 
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

# Purely topological scalar Laplacian Lf without metric components 
function topologicalLf(A, B)
    I = size(B, 1)
    Bᵀ = Transpose(B)
    onesVec = ones(1, I)
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

# Purely topological scalar Laplacian Lc without metric components 
function topologicalLc(A, B)
    I = size(B, 1)
    Bᵀ = Transpose(B)
    onesVec = ones(1, I)
    boundaryEdges = abs.(onesVec * B)
    boundaryEdgesFactor = abs.(boundaryEdges .- 1)      # =1 for internal vertices, =0 for boundary vertices
    boundaryEdgesFactorMat = Diagonal(boundaryEdgesFactor[1, :])
    Lc = B * boundaryEdgesFactorMat * Bᵀ
    Lc = B * Bᵀ
    dropzeros!(Lc)
    return Lc
end

# Purely topological scalar Laplacian Lv without metric components 
function topologicalLv(A, B)
    Aᵀ = Transpose(A)
    Lᵥ = Aᵀ * A
    dropzeros!(Lᵥ)
    return Lᵥ
end

# Purely topological scalar Laplacian Lt without metric components 
function topologicalLt(A, B)
    Aᵀ = Transpose(A)
    Lₜ = Aᵀ * A
    dropzeros!(Lₜ)
    return Lₜ
end

# Vector Laplacian
function edgeLaplacianPrimal(R, A, B)
    aᵢ = findCellAreas(R, A, B)
    Eⱼ = findCellLinkTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    Aᵀ = transpose(A)
    Bᵀ = transpose(B)
    H⁻¹ = spdiagm(1.0./aᵢ)
    E⁻¹ = spdiagm(1.0./Eⱼ)
    Tₑ = spdiagm((tⱼ.^2)./Fⱼ)
    Tₑ⁻¹ = spdiagm(Fⱼ./(tⱼ.^2))
    # boundaryEdgesFactor = abs.(findBoundaryEdges .- 1)# =1 for internal vertices, =0 for boundary vertices
    # diagonalComponent = boundaryEdgesFactor'
    # Tₑ = Diagonal(diagonalComponent)
    Lprimal = A*E⁻¹*Aᵀ*Tₑ⁻¹ + Tₑ*Bᵀ*H⁻¹*B
    dropzeros!(Lprimal)
    return Lprimal
end
   
# Vector Laplacian
function edgeLaplacianDual(R, A, B)
    aᵢ = findCellAreas(R, A, B)
    Eⱼ = findCellLinkTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    Aᵀ = transpose(A)
    Bᵀ = transpose(B)
    H⁻¹ = spdiagm(1.0./aᵢ)
    E⁻¹ = spdiagm(1.0./Eⱼ)
    Tₗ = spdiagm((Tⱼ.^2)./Fⱼ)
    Tₗ⁻¹ = spdiagm(Fⱼ./(Tⱼ.^2))
    Ldual = Bᵀ*H⁻¹*B*Tₗ⁻¹ + Tₗ*A*E⁻¹*Aᵀ
    dropzeros!(Ldual)
    return Ldual
end

# Proc Roy Soc A20; 
function cotanL(R, A, B)
    I = size(B,1)
    J = size(B,2)
    K = size(A,2)
    aᵢ = findCellAreas(R, A, B)
    𝐬ᵢₖ = findEdgeMidpointLinks(R, A, B)
    Dₖ = findEdgeMidpointLinkVertexAreas(R, A, B)
    tmp = [𝐬ᵢₖ[i,k]⋅𝐬ᵢₖ[i′,k]/(aᵢ[i]*Dₖ[k]) for i=1:I, i′=1:I, k=1:K]
    L = sparse(dropdims(sum(tmp, dims=3), dims=3))
    return L
end

function cotan𝐋(R, A, B)
    I = size(B,1)
    J = size(B,2)
    K = size(A,2)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    aᵢ = findCellAreas(R, A, B)
    𝐬ᵢₖ = findEdgeMidpointLinks(R, A, B)
    Dₖ = findEdgeMidpointLinkVertexAreas(R, A, B)
    tmp = [(outerProd(𝐬ᵢₖ[i,k], 𝐬ᵢₖ[i,k′])/Dₖ[k] + outerProd(ϵᵢ*𝐬ᵢₖ[i,k], ϵᵢ*𝐬ᵢₖ[i,k′])/Dₖ[k])/aᵢ[i] for i=1:I, k=1:K, k′=1:K]
    𝐋 = sparse(dropdims(sum(tmp, dims=1), dims=1))
    return 𝐋
end


#!!!!!!!!!!!!

#   Lⱼⱼ′=∑ₖAⱼₖ(t̂ⱼ⋅t̂ⱼ′)Aⱼ′ₖ/Eₖ
function scalarEdgeL(R, A, B)
    I = size(B,1); J = size(B,2); K = size(A,2)
    edgeTangentsNormalised = normalize.(findEdgeTangents(R, A))
    L = spzeros(J,J)
    vertexAreas = findEdgeMidpointLinkVertexAreas(R, A, B)
    for k=1:K
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

# Lᵢᵢ′ = ∑ⱼⱼ,ₖ(1/Lᵢ)B̄ᵢⱼAⱼₖ(t̂ⱼ⋅t̂ⱼ′)(Aⱼ′ₖB̄ᵢ′ⱼ′/Dₖ)
# Lᵢ = ∑ⱼB̄ᵢⱼAⱼₖ𝐭̂ⱼ
# Proc Roy Soc A20; 
function uniformCellTensionL(R, A, B)
    I = size(B,1); J = size(B,2); K = size(A,2)
    𝐭 = findEdgeTangents(R, A)
    𝐭̂ = normalize.(𝐭)
    B̄ = abs.(B)
    Dₖ = findEdgeMidpointLinkVertexAreas(R, A, B)
    Lᵢtmp = [B̄[i,j]*𝐭[j] for i=1:I, j=1:J] # A11
    Lᵢ = dropdims(sum(Lᵢtmp, dims=2), dims=2)

    Ltmp = spzeros(I,I)
    Ltmp .= [B̄[i,j]*A[j,k]*(𝐭̂[j]⋅𝐭̂[j′])*A[j′,k]*B̄[i′,j′]./(Lᵢ[i].*Dₖ[k]) for i=1:I, i′=1:I, j=1:J, j′=1:J, k=1:K]
    return dropdims(sum(tmp, dims=(3:5)), dims=(3:5))
end

#!!!!!!!!!!!!


export geometricLf
export geometricLc
export geometricLv
export geometricLt
export topologicalLf
export topologicalLc
export topologicalLv
export topologicalLt
export edgeLaplacianPrimal
export edgeLaplacianDual

export cotanL
export cotan𝐋
export edgeMidpointLNeumannOld

export scalarEdgeL
export uniformCellTensionL

end #end module 



# for k=1:K
    #     k_js = findall(x->x!=0,A[:,k])
    #     for j in k_js
    #         for j′ in k_js 
    #             j_is =findall(x->x!=0,B[:,j])#∩findall(x->x!=0,B[:,j′])
    #             j′_i′s =findall(x->x!=0,B[:,j′])#∩findall(x->x!=0,B[:,j′])
    #             for i in j_is
    #                 for i′ in j′_i′s 
    #                     L[i,i′] += (𝐭̂[j]⋅𝐭̂[j′])*A[j,k]*A[j′,k]*B̄[i,j]*B̄[i′,j′]/vertexAreas[k]
    #                 end
    #             end
    #         end
    #     end
    # end
    # dropzeros!(L)
    # return L



# function edgeMidpointLNeumannOld(R, A, B)
#     I = size(B,1); J = size(B,2); K = size(A,2)
#     Aᵀ = Transpose(A)
#     Āᵀ = abs.(Aᵀ)
#     Ā = abs.(A)
#     B̄ = abs.(B)
#     C = B̄ * Ā .÷ 2
#     cellAreas = findCellAreas(R, A, B)
#     edgeMidpointLinks = findEdgeMidpointLinks(R, A, B)
#     boundaryVertices = Āᵀ * abs.(sum.(eachcol(B))) .÷ 2
#     for k in findall(x->x!=0, boundaryVertices)
#         edgeMidpointLinks[:,k] .= fill(SVector{2,Float64}(zeros(2)),I)
#     end
#     edgeTangents = findEdgeTangents(R, A)
#     vertexAreas = findEdgeMidpointLinkVertexAreas(R, A, B)
#     L = spzeros(Float64,(I,I))
#     for k=1:K
#         for i in findall(x->x!=0, C[:,k])
#             for i′ in findall(x->x!=0, C[:,k])
#                 L[i, i′] += (edgeMidpointLinks[i,k]⋅edgeMidpointLinks[i′,k])/(cellAreas[i]*vertexAreas[k])
#             end
#         end
#     end
#     dropzeros!(L)
#     return L
# end




# function edgeMidpointLfunction(ϕᵢ, R, A, B)
    
#     I = size(B,1); J = size(B,2); K = size(A,2)

#     cellAreas = findCellAreas(R, A, B)
    
#     𝐬ᵢₖ = findEdgeMidpointLinks(R, A, B)

#     edgeTangents = findEdgeTangents(R, A)
    
#     Ā = abs.(A)
#     B̄ = abs.(B)
#     C = B̄ * Ā .÷ 2

#     αₖ=Float64[]
#     for k=1:K
#         k_is = findall(x->x!=0, C[:,k])
#         α = 0.5*norm([𝐬ᵢₖ[(k_is[1], k)]...,0.0]×[𝐬ᵢₖ[(k_is[2],k)]...,0.0])
#         push!(αₖ,α)
#     end

#     Lϕ = zeros(I)

#     for k=1:K
#         for i in findall(x->x!=0, C[:,k])
#             for i′ in findall(x->x!=0, C[:,k])
#                 Lϕ[i] += (𝐬ᵢₖ[i,k]⋅𝐬ᵢₖ[i′,k])*ϕᵢ[i′]/(cellAreas[i]*αₖ[k])
#             end
#         end
#     end
    
#     return Lϕ
# end
#
#  Laplacians.jl
#  DiscreteCalculus
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

# Laplacian Lf
# Jensen and Revell 2023 Eq.6b
# Scalar dual network laplacian acting on cells
function geometricLf(R, A, B)
    Bᵀ = Transpose(B)
    aᵢ = findCellAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    Tₑ = spdiagm((tⱼ .^ 2)./Fⱼ)
    H⁻¹ = spdiagm(1.0./aᵢ)
    Ib = spdiagm(1 .- findPeripheralEdges(B))
    L_ℱ = H⁻¹*B*Ib*Tₑ*Bᵀ
    dropzeros!(L_ℱ)
    return L_ℱ
end

function geometricLfHat(R, A, B)
    I = size(B,1)
    J = size(A,1)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ

    aᵢ = findCellAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    
    Bⁿ = B[:,jⁿ.==1]
    Bⁿᵀ = Transpose(Bⁿ)
    Bⁱ = B[:, jⁱ.==1]
    Bⁱᵀ = Transpose(Bⁱ)

    H⁻¹ = spdiagm(1.0./aᵢ)
    Tₑⁿ = spdiagm((tⱼ[jⁿ.==1] .^ 2)./Fⱼ[jⁿ.==1])
    Tₑⁱ = spdiagm((tⱼ[jⁱ.==1] .^ 2)./Fⱼ[jⁱ.==1])

    L̂_ℱ = H⁻¹*(Bⁿ*Tₑⁿ*Bⁿᵀ + Bⁱ*Tₑⁱ*Bⁱᵀ)
    dropzeros!(L̂_ℱ)

    return L̂_ℱ, collect(1:I)
end

# Laplacian Lc
# Jensen and Revell 2023 Eq.6b
# Scalar dual network laplacian acting on cells
function geometricLc(R, A, B)
    Bᵀ = Transpose(B)
    aᵢ = findCellAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    H⁻¹ = spdiagm(1.0./aᵢ)
    Tₗ⁻¹ = spdiagm(Fⱼ./(Tⱼ.^2))
    Ib = spdiagm(1 .- findPeripheralEdges(B))
    L_𝒞 = H⁻¹*B*Ib*Tₗ⁻¹*Bᵀ
    dropzeros!(L_𝒞)
    return L_𝒞
end

function geometricLcHat(R, A, B)
    I = size(B,1)
    J = size(A,1)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ

    aᵢ = findCellAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)

    Bⁿ = B[:,jⁿ.==1]
    Bⁿᵀ = Transpose(Bⁿ)
    Bⁱ = B[:, jⁱ.==1]
    Bⁱᵀ = Transpose(Bⁱ)

    H⁻¹ = spdiagm(1.0./aᵢ)
    Tₗⁿ⁻¹ = spdiagm(Fⱼ[jⁿ.==1]./(Tⱼ[jⁿ.==1].^2))
    Tₗⁱ⁻¹ = spdiagm(Fⱼ[jⁱ.==1]./(Tⱼ[jⁱ.==1].^2))
    
    L̂_𝒞 = H⁻¹*(Bⁿ*Tₗⁿ⁻¹*Bⁿᵀ + Bⁱ*Tₗⁱ⁻¹*Bⁱᵀ)
    dropzeros!(L̂_𝒞)
    
    return L̂_𝒞, collect(1:I)
end

# Laplacian Lv
# Jensen and Revell 2023 Eq.6a
# Scalar primal network laplacian acting on vertices
function geometricLv(R, A, B)
    Aᵀ = Transpose(A)
    tⱼ = findEdgeLengths(R, A)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    E⁻¹ = spdiagm(1.0./Eₖ)
    Tₑ⁻¹ = spdiagm(Fⱼ./(tⱼ .^ 2))
    L_𝒱 = E⁻¹*Aᵀ*Tₑ⁻¹*A
    dropzeros!(L_𝒱)
    return L_𝒱
end

function geometricLvHat(R, A, B)
    J = size(A,1)
    K = size(A,2)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ
    kᵖ = findPeripheralVertices(A,B)
    kⁱ = ones(K).-kᵖ

    tⱼ = findEdgeLengths(R, A)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)

    Aⁿⁱ = A[jⁿ.==1, kⁱ.==1]
    Aⁿⁱᵀ = Transpose(Aⁿⁱ)
    Aⁱⁱ = A[jⁱ.==1, kⁱ.==1]
    Aⁱⁱᵀ = Transpose(Aⁱⁱ)

    Eⁱ⁻¹ = spdiagm(1.0./Eₖ[kⁱ.==1])
    Tₑⁿ⁻¹ = spdiagm(Fⱼ[jⁿ.==1]./(tⱼ[jⁿ.==1] .^ 2))
    Tₑⁱ⁻¹ = spdiagm(Fⱼ[jⁱ.==1]./(tⱼ[jⁱ.==1] .^ 2))

    L̂_𝒱 = Eⁱ⁻¹*(Aⁿⁱᵀ*Tₑⁿ⁻¹*Aⁿⁱ + Aⁱⁱᵀ*Tₑⁱ⁻¹*Aⁱⁱ)
    dropzeros!(L̂_𝒱)
    return L̂_𝒱, collect(1:K)[kⁱ.==1]
end

# Laplacian Lt
# Jensen and Revell 2023 Eq.6a
# Scalar dual network laplacian acting on vertices
function geometricLt(R, A, B)
    Aᵀ = Transpose(A)
    Tⱼ = findCellLinkLengths(R, A, B)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = findEdgeQuadrilateralAreas(R, A, B)
    E⁻¹ = spdiagm(1.0./Eₖ)
    Tₗ = spdiagm((Tⱼ.^2)./Fⱼ)
    L_𝒯 = E⁻¹*Aᵀ*Tₗ*A
    dropzeros!(L_𝒯)
    return L_𝒯
end

function geometricLtHat(R, A, B)
    J = size(A,1)
    K = size(A,2)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ
    kᵖ = findPeripheralVertices(A,B)
    kⁱ = ones(K).-kᵖ

    Tⱼ = findCellLinkLengths(R, A, B)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = findEdgeQuadrilateralAreas(R, A, B)
    
    Aⁿⁱ = A[jⁿ.==1, kⁱ.==1]
    Aⁿⁱᵀ = Transpose(Aⁿⁱ)
    Aⁱⁱ = A[jⁱ.==1, kⁱ.==1]
    Aⁱⁱᵀ = Transpose(Aⁱⁱ)
    
    Eⁱ⁻¹ = spdiagm(1.0./Eₖ[kⁱ.==1])
    Tₗⁿ = spdiagm((Tⱼ[jⁿ.==1].^2)./Fⱼ[jⁿ.==1])
    Tₗⁱ = spdiagm((Tⱼ[jⁱ.==1].^2)./Fⱼ[jⁱ.==1])

    L̂_𝒯 = Eⁱ⁻¹*(Aⁿⁱᵀ*Tₗⁿ*Aⁿⁱ + Aⁱⁱᵀ*Tₗⁱ*Aⁱⁱ)
    dropzeros!(L̂_𝒯)
    
    return L̂_𝒯, collect(1:K)[kⁱ.==1]
end

# Purely topological scalar Laplacians Lf, Lc, Lv, Lt without metric components. Without metric component, Lf=Lc, Lv=Lt
topologicalLf(A, B) = dropzeros(B * Transpose(B))
topologicalLc(A, B) = topologicalLf(A, B)
topologicalLv(A, B) = dropzeros(Aᵀ * Transpose(A))
topologicalLt(A, B) = topologicalLv(A, B)    

# Vector Laplacian \mathsf{L}_{\mathcal{E}}
# Jensen and Revell 2026 Eq.2.48
function edgeLaplacianPrimal(R, A, B)
    aᵢ = findCellAreas(R, A, B)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    Aᵀ = transpose(A)
    Bᵀ = transpose(B)
    H⁻¹ = spdiagm(1.0./aᵢ)
    E⁻¹ = spdiagm(1.0./Eₖ)
    Tₑ = spdiagm((tⱼ.^2)./Fⱼ)
    Tₑ⁻¹ = spdiagm(Fⱼ./(tⱼ.^2))
    L_ℰ = A*E⁻¹*Aᵀ*Tₑ⁻¹ + Tₑ*Bᵀ*H⁻¹*B
    dropzeros!(L_ℰ)
    return L_ℰ
end

# Vector Laplacian \mathsf{L}_{\mathcal{E}}
# Jensen and Revell 2026 Eq.2.48
function edgeLaplacianPrimalHat(R, A, B)
    aᵢ = findCellAreas(R, A, B)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    H⁻¹ = spdiagm(1.0./aᵢ)
    
    jᵖ = findPeripheralEdges(B)
    kᵖ = findPeripheralVertices(A,B)
    kⁱ = 1 .-kᵖ
    Â = findÂ(A, B)
    Âᵀ = Transpose(Â)
    B̂ = findB̂(A, B)
    B̂ᵀ = Transpose(B̂)
    T̂ₑ = spdiagm((tⱼ[jᵖ.==0].^2)./Fⱼ[jᵖ.==0])
    T̂ₑ⁻¹ = spdiagm(Fⱼ[jᵖ.==0]./(tⱼ[jᵖ.==0].^2))
    Ê⁻¹ = spdiagm(1.0./Eₖ[kⁱ.==1])
    
    L̂_ℰ = Â*Ê⁻¹*Âᵀ*T̂ₑ⁻¹ + T̂ₑ*B̂ᵀ*H⁻¹*B̂
    dropzeros!(L̂_ℰ)
    return L̂_ℰ
end

# Vector Laplacian \mathsf{L}_{\mathcal{L}}
# Jensen and Revell 2026 Eq.2.54
function edgeLaplacianDual(R, A, B)
    aᵢ = findCellAreas(R, A, B)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    Aᵀ = transpose(A)
    Bᵀ = transpose(B)
    H⁻¹ = spdiagm(1.0./aᵢ)
    E⁻¹ = spdiagm(1.0./Eₖ)
    Tₗ = spdiagm((Tⱼ.^2)./Fⱼ)
    Tₗ⁻¹ = spdiagm(Fⱼ./(Tⱼ.^2))
    L_ℒ = Bᵀ*H⁻¹*B*Tₗ⁻¹ + Tₗ*A*E⁻¹*Aᵀ
    dropzeros!(L_ℒ)
    return L_ℒ
end

# Vector Laplacian \mathsf{L}_{\mathcal{L}}
# Jensen and Revell 2026 Eq.2.54
function edgeLaplacianDualHat(R, A, B)
    aᵢ = findCellAreas(R, A, B)
    Eₖ = findCellLinkVertexTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    H⁻¹ = spdiagm(1.0./aᵢ)

    jᵖ = findPeripheralEdges(B)
    kᵖ = findPeripheralVertices(A,B)
    kⁱ = 1 .-kᵖ
    Â = findÂ(A, B)
    Âᵀ = Transpose(Â)
    B̂ = findB̂(A, B)
    B̂ᵀ = Transpose(B̂)
    T̂ₗ = spdiagm((Tⱼ[jᵖ.==0].^2)./Fⱼ[jᵖ.==0]) 
    T̂ₗ⁻¹ = spdiagm(Fⱼ[jᵖ.==0]./(Tⱼ[jᵖ.==0].^2)) 
    Ê⁻¹ = spdiagm(1.0./Eₖ[kⁱ.==1])

    L̂_ℒ = B̂ᵀ*H⁻¹*B̂*T̂ₗ⁻¹ + T̂ₗ*Â*Ê⁻¹*Âᵀ
    dropzeros!(L̂_ℒ)
    return L̂_ℒ
end

# Cowley, N., Revell, C. K., Johns, E., Woolner, S. & Jensen, O. E. Spectral approaches to stress 
# relaxation in epithelial monolayers. Proceedings of the Royal Society A (2024) doi:10.1098/rspa.2024.0224.
# Eq A20
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
    aᵢ = findCellAreas(R, A, B)
    𝐬ᵢₖ = findEdgeMidpointLinks(R, A, B)
    Dₖ = findEdgeMidpointLinkVertexAreas(R, A, B)
    tmp = [(outerProd(𝐬ᵢₖ[i,k], 𝐬ᵢₖ[i,k′])/Dₖ[k] + outerProd(ϵᵢ*𝐬ᵢₖ[i,k], ϵᵢ*𝐬ᵢₖ[i,k′])/Dₖ[k])/aᵢ[i] for i=1:I, k=1:K, k′=1:K]
    𝐋 = sparse(dropdims(sum(tmp, dims=1), dims=1))
    return 𝐋
end

export geometricLf
export geometricLfHat
export geometricLc
export geometricLcHat
export geometricLv
export geometricLvHat
export geometricLt
export geometricLtHat
export topologicalLf
export topologicalLc
export topologicalLv
export topologicalLt
export edgeLaplacianPrimal
export edgeLaplacianPrimalHat
export edgeLaplacianDual
export edgeLaplacianDualHat

export cotanL
export cotan𝐋

end #end module 


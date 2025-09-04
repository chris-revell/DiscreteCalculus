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
    Bᵀ = Transpose(B)
    aᵢ = findCellAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    jⁱ = abs.(findPeripheralEdges(B).-1) # =1 for internal edges, =0 for boundary edges
    Tₑ = spdiagm((jⁱ.*((tⱼ .^ 2)./(Fⱼ)))) # Multiply by jⁱ vector to set boundary vertex contributions to zero
    H⁻¹ = spdiagm(1.0./aᵢ)
    L_ℱ = H⁻¹*B*Tₑ*Bᵀ
    dropzeros!(L_ℱ)
    return L_ℱ
end


function geometricLfHat(R, A, B)
    J = size(A,1)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ

    aᵢ = findCellAreas(R, A, B)
    tⱼ = findEdgeLengths(R, A)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    Tₑ = spdiagm((tⱼ .^ 2)./(Fⱼ))
    H⁻¹ = spdiagm(1.0./aᵢ)

    Bⁿ = copy(B)
    Bⁿ[:,jⁿ.==0].=0
    dropzeros!(Bⁿ)
    Bⁿᵀ = Transpose(Bⁿ)
    Bⁱ = copy(B)
    Bⁱ[:, jⁱ.==0].=0
    dropzeros!(Bⁱ)
    Bⁱᵀ = Transpose(Bⁱ)

    Tₑⁿ = spdiagm(jⁿ.*(tⱼ .^ 2)./(Fⱼ))
    Tₑⁱ = spdiagm(jⁱ.*(tⱼ .^ 2)./(Fⱼ))

    L̂_ℱ = H⁻¹*(Bⁿ*Tₑⁿ*Bⁿᵀ + Bⁱ*Tₑⁱ*Bⁱᵀ)

    dropzeros!(L̂_ℱ)
    return L̂_ℱ
end

function geometricLfHatReduced(R, A, B)
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

# Scalar Laplacian Lc with metric components 
function geometricLc(R, A, B)
    Bᵀ = Transpose(B)
    aᵢ = findCellAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    jⁱ = abs.(findPeripheralEdges(B).-1) # =1 for internal edges, =0 for boundary edges
    H⁻¹ = spdiagm(1.0./aᵢ)
    Tₗ⁻¹ = spdiagm(jⁱ.*Fⱼ./(Tⱼ.^2)) # Multiply by jⁱ vector to set boundary vertex contributions to zero
    L_𝒞 = H⁻¹*B*Tₗ⁻¹*Bᵀ
    dropzeros!(L_𝒞)
    return L_𝒞
end

function geometricLcHat(R, A, B)
    J = size(A,1)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ

    aᵢ = findCellAreas(R, A, B)
    Tⱼ = findCellLinkLengths(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    H⁻¹ = spdiagm(1.0./aᵢ)

    Bⁿ = copy(B)
    Bⁿ[:,jⁿ.==0].=0
    dropzeros!(Bⁿ)
    Bⁿᵀ = Transpose(Bⁿ)
    Bⁱ = copy(B)
    Bⁱ[:, jⁱ.==0].=0
    dropzeros!(Bⁱ)
    Bⁱᵀ = Transpose(Bⁱ)

    Tₗⁿ⁻¹ = spdiagm(jⁿ.*Fⱼ./(Tⱼ.^2))
    Tₗⁱ⁻¹ = spdiagm(jⁱ.*Fⱼ./(Tⱼ.^2))
    
    L̂_𝒞 = H⁻¹*(Bⁿ*Tₗⁿ⁻¹*Bⁿᵀ + Bⁱ*Tₗⁱ⁻¹*Bⁱᵀ)
    dropzeros!(L̂_𝒞)
    return L̂_𝒞
end

function geometricLcHatReduced(R, A, B)
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

# Scalar Laplacian Lv with metric components 
function geometricLv(R, A, B)
    Aᵀ = Transpose(A)
    tⱼ = findEdgeLengths(R, A)
    Eₖ = findCellLinkTriangleAreas(R, A, B)
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
    Eₖ = findCellLinkTriangleAreas(R, A, B)
    Fⱼ = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    E⁻¹ = spdiagm(1.0./Eₖ)
    
    Aⁿⁱ = copy(A)
    Aⁿⁱ[jⁿ.==0, :].=0
    Aⁿⁱ[:, kⁱ.==0].=0
    dropzeros!(Aⁿⁱ)
    Aⁿⁱᵀ = Transpose(Aⁿⁱ)
    Aⁱⁱ = copy(A)
    Aⁱⁱ[jⁱ.==0, :].=0
    Aⁱⁱ[:, kⁱ.==0].=0
    dropzeros!(Aⁱⁱ)
    Aⁱⁱᵀ = Transpose(Aⁱⁱ)

    Tₑⁿ⁻¹ = spdiagm(jⁿ.*Fⱼ./(tⱼ .^ 2))
    Tₑⁱ⁻¹ = spdiagm(jⁱ.*Fⱼ./(tⱼ .^ 2))

    L̂_𝒱 = E⁻¹*(Aⁿⁱᵀ*Tₑⁿ⁻¹*Aⁿⁱ + Aⁱⁱᵀ*Tₑⁱ⁻¹*Aⁱⁱ)
    dropzeros!(L̂_𝒱)
    return L̂_𝒱
end

function geometricLvHatReduced(R, A, B)
    J = size(A,1)
    K = size(A,2)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ
    kᵖ = findPeripheralVertices(A,B)
    kⁱ = ones(K).-kᵖ

    tⱼ = findEdgeLengths(R, A)
    Eₖ = findCellLinkTriangleAreas(R, A, B)
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

# Scalar Laplacian Lt with metric components 
function geometricLt(R, A, B)
    Aᵀ = Transpose(A)
    Tⱼ = findCellLinkLengths(R, A, B)
    Eₖ = findCellLinkTriangleAreas(R, A, B)
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
    Eₖ = findCellLinkTriangleAreas(R, A, B)
    Fⱼ = findEdgeQuadrilateralAreas(R, A, B)
    E⁻¹ = spdiagm(1.0./Eₖ)
    
    Aⁿⁱ = copy(A)
    Aⁿⁱ[jⁿ.==0, :].=0
    Aⁿⁱ[:, kⁱ.==0].=0
    dropzeros!(Aⁿⁱ)
    Aⁿⁱᵀ = Transpose(Aⁿⁱ)
    Aⁱⁱ = copy(A)
    Aⁱⁱ[jⁱ.==0, :].=0
    Aⁱⁱ[:, kⁱ.==0].=0
    dropzeros!(Aⁱⁱ)
    Aⁱⁱᵀ = Transpose(Aⁱⁱ)

    Tₗⁿ = spdiagm(jⁿ.*(Tⱼ.^2)./Fⱼ)
    Tₗⁱ = spdiagm(jⁱ.*(Tⱼ.^2)./Fⱼ)

    L̂_𝒯 = E⁻¹*(Aⁿⁱᵀ*Tₗⁿ*Aⁿⁱ + Aⁱⁱᵀ*Tₗⁱ*Aⁱⁱ)
    dropzeros!(L̂_𝒯)
    return L̂_𝒯
end

function geometricLtHatReduced(R, A, B)
    J = size(A,1)
    K = size(A,2)
    jⁿ = findNormalEdges(A,B)
    jᵖ = findPeripheralEdges(B)
    jⁱ = ones(J).-jⁿ.-jᵖ
    kᵖ = findPeripheralVertices(A,B)
    kⁱ = ones(K).-kᵖ

    Tⱼ = findCellLinkLengths(R, A, B)
    Eₖ = findCellLinkTriangleAreas(R, A, B)
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

# Check topological versions 
#!!!!!!!!!
# Purely topological scalar Laplacian Lf without metric components 
function topologicalLf(A, B)
    Bᵀ = Transpose(B)
    onesVec = ones(1, I)
    b = spdiagm(abs.(findPeripheralEdges(B) .- 1)) # =1 for internal vertices, =0 for boundary vertices
    Lf = B * b * Bᵀ
    dropzeros!(Lf)
    return Lf
end
# Purely topological scalar Laplacian Lc without metric components 
function topologicalLc(A, B)
    Bᵀ = Transpose(B)
    onesVec = ones(1, I)
    b = spdiagm(abs.(findPeripheralEdges(B) .- 1)) # =1 for internal vertices, =0 for boundary vertices
    Lf = B * b * Bᵀ
    dropzeros!(Lf)
    return Lf
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
#!!!!!!!!!

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
    # boundaryEdgesFactor = abs.(findPeripheralEdges .- 1)# =1 for internal vertices, =0 for boundary vertices
    # diagonalComponent = boundaryEdgesFactor'
    # Tₑ = spdiagm(diagonalComponent)
    L_ℰ = A*E⁻¹*Aᵀ*Tₑ⁻¹ + Tₑ*Bᵀ*H⁻¹*B
    dropzeros!(L_ℰ)
    return L_ℰ
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
    L_ℒ = Bᵀ*H⁻¹*B*Tₗ⁻¹ + Tₗ*A*E⁻¹*Aᵀ
    dropzeros!(L_ℒ)
    return L_ℒ
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
    # 𝐪ₖ = eachcol(spdiagm(ones(Int64, K)))
    # tmp = [outerProd(𝐪ₖ[k],𝐪ₖ[k′]).*(outerProd(𝐬ᵢₖ[i,k], 𝐬ᵢₖ[i,k′])/Dₖ[k] + outerProd(ϵᵢ*𝐬ᵢₖ[i,k], ϵᵢ*𝐬ᵢₖ[i,k′])/Dₖ[k])/aᵢ[i] for i=1:I, k=1:K, k′=1:K]
    # 𝐋 = sparse(sum(tmp, dims=1))
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
export geometricLfHat
export geometricLfHatReduced
export geometricLc
export geometricLcHat
export geometricLcHatReduced
export geometricLv
export geometricLvHat
export geometricLvHatReduced
export geometricLt
export geometricLtHat
export geometricLtHatReduced
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
#     aᵢ = findCellAreas(R, A, B)
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
#                 L[i, i′] += (edgeMidpointLinks[i,k]⋅edgeMidpointLinks[i′,k])/(aᵢ[i]*vertexAreas[k])
#             end
#         end
#     end
#     dropzeros!(L)
#     return L
# end




# function edgeMidpointLfunction(ϕᵢ, R, A, B)
    
#     I = size(B,1); J = size(B,2); K = size(A,2)

#     aᵢ = findCellAreas(R, A, B)
    
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
#                 Lϕ[i] += (𝐬ᵢₖ[i,k]⋅𝐬ᵢₖ[i′,k])*ϕᵢ[i′]/(aᵢ[i]*αₖ[k])
#             end
#         end
#     end
    
#     return Lϕ
# end
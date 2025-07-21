#
#  DifferentialOperators.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
# Differential operators over primary and dual networks

# Definitions: 
# {cocurlᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
# {divᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
# {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
# {cocurlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ          With boundary considerations (Jensen and Revell 2023 Eq 12): {cocurlᵛ 𝐛}ₖ = -∑ᵢⱼₖBᵢⱼAⱼₖ(ϵₖqᵢₖ)⋅𝐛ⱼ/Eₖ
# {curlᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
# {codivᶜ 𝐛}ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
# {codivᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
# {curlᵛ 𝐛}ₖ = -∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ               With boundary considerations (Jensen and Revell 2023 Eq 12): {curlᵛ 𝐛}ₖ = ∑ᵢⱼₖBᵢⱼAⱼₖqᵢₖ⋅𝐛ⱼ/Eₖ
# {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
# {cogradᶜ f}ⱼ = ∑ᵢϵₖ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
# {gradᵛ ϕ }ⱼ = = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
# {cogradᵛ ϕ}ⱼ = ∑ₖϵᵢ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
# {corotᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
# {rotᶜ f }ⱼ = -∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
# {corotᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ            With boundary considerations (Jensen and Revell 2023 Eq 12): {corotᵛ ϕ}ⱼ = -∑ᵢⱼₖBᵢⱼAⱼₖϵₖ𝐪ᵢₖϕₖ/Fⱼ
# {rotᵛ ϕ}ⱼ = -∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ                 With boundary considerations (Jensen and Revell 2023 Eq 12): {rotᵛ ϕ}ⱼ = ∑ᵢⱼₖBᵢⱼAⱼₖ𝐪ᵢₖϕₖ/Fⱼ

# Old name to new name conversions
# Old   => new:
# -divᶜ => cocurlᶜ
# -d̃ivᶜ => -divᶜ
# -d̃ivᵛ => -divᵛ
# -divᵛ => cocurlᵛ
# curlᶜ => -curlᶜ
# C̃URLᶜ => codivᶜ
# c̃urlᵛ => codivᵛ 
# CURLᵛ => -curlᵛ
# gradᶜ => gradᶜ
# CURLᶜ => -cogradᶜ 
# gradᵛ => gradᵛ
# curlᵛ => -cogradᵛ
# g̃radᶜ => corotᶜ
# c̃urlᶜ => -rotᶜ
# g̃radᵛ => corotᵛ
# C̃URLᵛ => -rotᵛ

# Component definitions
# aᵢ => Cell areas, for each cell i
# Aⱼₖ => Incidence matrix, , for each edge j and vertex k
# Bᵢⱼ => Incidence matrix, for each cell i and edge j
# 𝐭ⱼ => Edge tangents, for each edge j
# tⱼ => Edge lengths, for each edge j
# 𝐓ⱼ => Cell centre links, for each edge j
# Tⱼ => Cell centre lengths, for each edge j 
# Eₖ => Cell centre link triangle areas, for each vertex k
# Fⱼ => 2x Edge quadrilateral area, for each edge j 
# Kᵢₖ => Vertex kite area, for each cell-vertex pair i, k.
# 𝐜ⱼ => Edge midpoint for each edge j 
# 𝐫ₖ => Vertex position for each vertex k 

# ϕ some scalar field over vertices k 
# 𝐛 some vector field over edges j 
# f some scalar field over cells i 

module DifferentialOperators

using DrWatson
using SparseArrays
using StaticArrays
using LinearAlgebra
using FromFile 

@from "GeometryFunctions.jl" using GeometryFunctions
@from "TopologyFunctions.jl" using TopologyFunctions

# {cocurlᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
function cocurlᶜ(R, A, B, 𝐛)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    tmp = [B[i,j]*(ϵᵢ*𝐭[j])⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end
 
# {divᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
function divᶜ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    T = findCellLinkLengths(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [-B[i,j]*(F[j]/(T[j]^2))*𝐓[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end 
# With boundary component suppression 
function divᶜsuppress(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    T = findCellLinkLengths(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    a = findCellAreas(R, A, B)
    b = abs.(findBoundaryEdges(B).-1) # Indicator function 
    tmp = [-B[i,j]*(F[j]/(T[j]^2))*b[j]*𝐓[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end 

# {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
function divᵛ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    t = findEdgeLengths(R, A)
    𝐭 = findEdgeTangents(R, A)
    E = findCellLinkTriangleAreas(R, A, B)
    tmp = [-A[j,k]*(F[j]/(t[j]^2))*𝐭[j]⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end
# With boundary component suppression 
function divᵛsuppress(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    t = findEdgeLengths(R, A)
    𝐭 = findEdgeTangents(R, A)
    E = findCellLinkTriangleAreas(R, A, B)
    b = abs.(findBoundaryVertices(A, B) .-1) # Indicator function 
    tmp = [-A[j,k]*(F[j]/(t[j]^2))*b[k]*𝐭[j]⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# {cocurlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ
function cocurlᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    𝐓 = findCellLinks(R, A, B)
    tmp = [A[j,k]*(ϵₖ*𝐓[j])⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end
# With boundary considerations (Jensen and Revell 2023 Eq 12): {cocurlᵛ 𝐛}ₖ = -∑ᵢⱼₖBᵢⱼAⱼₖ(ϵₖqᵢₖ)⋅𝐛ⱼ/Eₖ
function cocurlᵛspokes(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    tmp = [-B[i,j]*A[j,k].*(ϵₖ*q[i,k])⋅𝐛[j]./E[k] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,2)), dims=(1,2))
end

# {curlᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
function curlᶜ(R, A, B, 𝐛)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    tmp = [-B[i,j]*𝐭[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# {codivᶜ 𝐛}ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
function codivᶜ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
        0.0 -1.0
        1.0 0.0
    ])
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [-B[i,j]*(F[j]/(T[j]^2))*(ϵₖ*𝐓[j])⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end
# With boundary component suppression 
function codivᶜsuppress(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
        0.0 -1.0
        1.0 0.0
    ])
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    a = findCellAreas(R, A, B)
    b = abs.(findBoundaryEdges(B).-1) # Indicator function 
    tmp = [-B[i,j]*(F[j]/(T[j]^2))*b[j]*(ϵₖ*𝐓[j])⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end


# {codivᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
function codivᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [-A[j,k]*(F[j]/(t[j]^2))*(ϵᵢ*𝐭[j])⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end
# With boundary component suppression 
function codivᵛsuppress(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    b = abs.(findBoundaryVertices(A, B) .-1) # Indicator function 
    tmp = [-A[j,k]*(F[j]/(t[j]^2))*b[k]*(ϵᵢ*𝐭[j])⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# {curlᵛ 𝐛}ₖ = -∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ
function curlᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    tmp = [-A[j,k]*(𝐓[j]⋅𝐛[j])/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end 
# With boundary considerations (Jensen and Revell 2023 Eq 12): {curlᵛ 𝐛}ₖ = ∑ᵢⱼₖBᵢⱼAⱼₖqᵢₖ⋅𝐛ⱼ/Eₖ
function curlᵛspokes(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    q = findSpokes(R, A, B)
    tmp = [B[i,j]*A[j,k]*(q[i,k]⋅𝐛[j])/E[k] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,2)), dims=(1,2))
end 

# {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
function gradᶜ(R, A, B, f)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    tmp = [B[i,j].*(𝐓[j]./(T[j]^2)).*f[i] for j=1:size(B,2), i=1:size(B,1)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# {cogradᶜ f}ⱼ = ∑ᵢϵₖ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
function cogradᶜ(R, A, B, f)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    tmp = [ϵₖ*(𝐓[j]./(T[j]^2)).*B[i,j]*f[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# {gradᵛ ϕ }ⱼ = = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
function gradᵛ(R, A, ϕ)    
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [A[j,k].*(𝐭[j]./(t[j]^2)).*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    # gradᵛ = [A[j,k].*(𝐭[j]./(t[j]^2)) for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
    return gradᵛ*ϕ
end

# {cogradᵛ ϕ}ⱼ = ∑ₖϵᵢ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
function cogradᵛ(R, A, B, ϕ)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [ϵᵢ*(𝐭[j]./(t[j]^2)).*A[j,k]*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# {corotᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
function corotᶜ(R, A, B, f)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    tmp = [B[i,j].*ϵᵢ*(𝐭[j]./F[j]).*f[i] for j=1:size(B,2), i=1:size(B,1)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# {rotᶜ f }ⱼ = -∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
function rotᶜ(R, A, B, f)
    𝐭 = findEdgeTangents(R, A)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tmp = [-B[i,j].*(𝐭[j]./F[j]).*f[i] for j=1:size(B,2), i=1:size(B,1)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# {corotᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ
function corotᵛ(R, A, B, ϕ)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    tmp = [A[j,k].*(ϵₖ*𝐓[j]).*ϕ[k]/F[j] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end
# With boundary considerations (Jensen and Revell 2023 Eq 12): {corotᵛ ϕ}ⱼ = -∑ᵢⱼₖBᵢⱼAⱼₖϵₖ𝐪ᵢₖϕₖ/Fⱼ
function corotᵛspokes(R, A, B, ϕ)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    tmp = [-B[i,j]*A[j,k].*ϵₖ*q[i,k].*ϕ[k]/F[j] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,3)), dims=(1,3))
end

# {rotᵛ ϕ}ⱼ = -∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ
function rotᵛ(R, A, B, ϕ)
    𝐓 = findCellLinks(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tmp = [-A[j,k].*𝐓[j].*ϕ[k]/F[j] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end 
# With boundary considerations (Jensen and Revell 2023 Eq 12): {rotᵛ ϕ}ⱼ = ∑ᵢⱼₖBᵢⱼAⱼₖ𝐪ᵢₖϕₖ/Fⱼ
function rotᵛspokes(R, A, B, ϕ)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    q = findSpokes(R, A, B)
    tmp = [B[i,j]*A[j,k].*q[i,k].*ϕ[k]/F[j] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,3)), dims=(1,3))
end 


export gradᵛ
export cogradᵛ
export corotᶜ
export rotᶜ
export gradᶜ
export cogradᶜ
export corotᵛ
export corotᵛspokes
export rotᵛ
export rotᵛspokes
export divᵛ
export divᵛsuppress
export codivᵛ
export codivᵛsuppress
export cocurlᶜ
export curlᶜ
export divᶜ
export divᶜsuppress
export codivᶜ
export codivᶜsuppress
export cocurlᵛ
export cocurlᵛspokes
export curlᵛ
export curlᵛspokes

end
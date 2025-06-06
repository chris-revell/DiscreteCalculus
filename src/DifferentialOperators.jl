#
#  DifferentialOperators.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
# Differential operators over primary and dual networks

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

# Primary operators:
# grad
# curl (associated with vectors parallel to edges and links)
# cog 
# cocurl (associated with vectors perpendicult to edges and links)

# Derived operators (adjoint under inner products associated with Hodge stars)
# -div 
# rot 
# -cod 
# corot 

module DifferentialOperators

using SparseArrays
using StaticArrays
using LinearAlgebra
using DiscreteCalculus

# Old: {gradᵛ ϕ }ⱼ = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
# Old => new: gradᵛ => gradᵛ
# New: {gradᵛ ϕ }ⱼ = = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
function gradᵛ(R, A)
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    gradᵛ = [A[j,k].*(𝐭[j]./(t[j]^2)) for j=1:size(A,1), k=1:size(A,2)]
    return gradᵛ
end
function gradᵛ(R, A, ϕ)    
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [A[j,k].*(𝐭[j]./(t[j]^2)).*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    # gradᵛ = [A[j,k].*(𝐭[j]./(t[j]^2)) for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
    return gradᵛ*ϕ
end

# Old: {curlᵛ ϕ}ⱼ = ∑ₖϵₖ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
# Old => new: curlᵛ => -cogᵛ
# New: {cogᵛ ϕ}ⱼ = -∑ₖϵₖ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
function cogᵛ(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    cogᵛ = [-ϵₖ*(𝐭[j]./(t[j]^2)).*A[j,k] for j=1:size(A,1), k=1:size(A,2)]
    return cogᵛ
end
function cogᵛ(R, A, B, ϕ)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [-ϵₖ*(𝐭[j]./(t[j]^2)).*A[j,k]*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {g̃radᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
# Old => new: g̃radᶜ => corotᶜ
# New: {corotᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
function corotᶜ(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    corotᶜ = [B[i,j].*ϵᵢ*(𝐭[j]./F[j]) for j=1:size(B,2), i=1:size(B,1)]
    return corotᶜ
end
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

# Old: {c̃urlᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
# Old => new: c̃urlᶜ => -rotᶜ
# New: {rotᶜ f }ⱼ = -∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
function rotᶜ(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    rotᶜ = [-B[i,j].*(𝐭[j]./F[j]) for j=1:size(B,2), i=1:size(B,1)]
    return rotᶜ
end
function rotᶜ(R, A, B, f)
    𝐭 = findEdgeTangents(R, A)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tmp = [-B[i,j].*(𝐭[j]./F[j]).*f[i] for j=1:size(B,2), i=1:size(B,1)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

#%%%%

# Old: {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
# Old => new: gradᶜ => gradᶜ
# New: {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
function gradᶜ(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    gradᶜ = [B[i,j].*(𝐓[j]./(T[j]^2)) for j=1:size(B,2), i=1:size(B,1)]
    return gradᶜ
end
function gradᶜ(R, A, B, f)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    tmp = [B[i,j].*(𝐓[j]./(T[j]^2)).*f[i] for j=1:size(B,2), i=1:size(B,1)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {CURLᶜ f}ⱼ = ∑ᵢϵᵢ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
# Old => new: CURLᶜ => -cogᶜ 
# New: {cogᶜ f}ⱼ = -∑ᵢϵᵢ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
function cogᶜ(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    cogᶜ = [-ϵᵢ*(𝐓[j]./(T[j]^2)).*B[i,j] for j=1:size(B,2), i=1:size(B,1)]
    return cogᶜ
end
function cogᶜ(R, A, B, f)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    tmp = [-ϵᵢ*(𝐓[j]./(T[j]^2)).*B[i,j]*f[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {g̃radᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ
# Old => new: g̃radᵛ => corotᵛ
# New: {corotᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ
# With boundary consideration: {corotᵛ ϕ}ⱼ = -∑ᵢⱼₖBᵢⱼAⱼₖϵₖ𝐪ᵢₖϕₖ/Fⱼ

# function corotᵛnaive(R, A, B)
#     F = spdiagm(1.0./(2.0.*findEdgeQuadrilateralAreas(R, A, B)))
#     ϵₖ = SMatrix{2, 2, Float64}([
#                 0.0 -1.0
#                 1.0 0.0
#             ])
#     ϵ = spdiagm(fill(ϵₖ, size(A,1)))
#     𝐓 = spdiagm(findCellLinks(R, A, B))
#     corotᵛ = ϵ*𝐓*F*A
#     return corotᵛ
# end
function corotᵛ(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    tmp = [-B[i,j]*A[j,k].*ϵₖ*q[i,k]./F[j] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    corotᵛ = dropdims(sum(tmp, dims=1), dims=1)
    return corotᵛ
end
# function corotᵛnaive(R, A, B, ϕ)
#     F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
#     ϵₖ = SMatrix{2, 2, Float64}([
#                 0.0 -1.0
#                 1.0 0.0
#             ])
#     𝐓 = findCellLinks(R, A, B)
#     tmp = [A[j,k].*ϵₖ*(𝐓[j]./F[j]).*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
#     return dropdims(sum(tmp, dims=2), dims=2)
# end
# With boundary considerations, Jensen and Revell 2023 Eq 12
function corotᵛ(R, A, B, ϕ)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    tmp = [-B[i,j]*A[j,k].*ϵₖ*q[i,k].*ϕ[k]/F[j] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,3)), dims=(1,3))
end

# Old: {C̃URLᵛ ϕ }ⱼ = ∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ
# Old => new: C̃URLᵛ => -rotᵛ
# New: {rotᵛ ϕ}ⱼ = -∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ
function rotᵛ(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    rotᵛ = [A[j,k].*𝐓[j]./F[j] for j=1:size(A,1), k=1:size(A,2)]
    return rotᵛ
end 
# function rotᵛ(R, A, B, ϕ)
#     F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
#     𝐓 = findCellLinks(R, A, B)
#     tmp = [A[j,k].*𝐓[j].*ϕ[k]/F[j] for j=1:size(A,1), k=1:size(A,2)]
#     return dropdims(sum(tmp, dims=2), dims=2)
# end 
# With boundary considerations, Jensen and Revell 2023 Eq 12
function rotᵛ(R, A, B, ϕ)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    q = findSpokes(R, A, B)
    tmp = [B[i,j]*A[j,k].*q[i,k].*ϕ[k]/F[j] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,3)), dims=(1,3))
end 

#%%%%

# Old: {d̃ivᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
# Old => new: -d̃ivᵛ => -divᵛ
# New: {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
function divᵛ(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    t = findEdgeLengths(R, A)
    𝐭 = findEdgeTangents(R, A)
    E = findCellLinkTriangleAreas(R, A, B)
    divᵛ = [-A[j,k]*(F[j]/(t[j]^2)).*transpose(𝐭[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return divᵛ
end
function divᵛ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    t = findEdgeLengths(R, A)
    𝐭 = findEdgeTangents(R, A)
    E = findCellLinkTriangleAreas(R, A, B)
    tmp = [-A[j,k]*(F[j]/(t[j]^2))*𝐭[j]⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {c̃urlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵₖ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
# Old => new: c̃urlᵛ => codᵛ 
# New: {codᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵₖ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
function codᵛ(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
        0.0 -1.0
        1.0 0.0
    ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    codᵛ = [A[j,k]*(F[j]/(t[j]^2)).*transpose(ϵₖ*𝐭[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return codᵛ
end
function codᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
        0.0 -1.0
        1.0 0.0
    ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [A[j,k]*(F[j]/(t[j]^2))*(ϵₖ*𝐭[j])⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old {divᶜ 𝐛}ᵢ = -∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
# Old => new: -divᶜ => cocurlᶜ
# New: {cocurlᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
function cocurlᶜ(R, A, B)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    cocurlᶜ = [B[i,j].*transpose(ϵᵢ*𝐭[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return cocurlᶜ
end
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

# Old: {curlᶜ 𝐛 }ᵢ = ∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
# Old => new: curlᶜ => -curlᶜ
# New: {curlᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
function curlᶜ(R, A, B)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    curlᶜ = [-B[i,j].*transpose(𝐭[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return curlᶜ
end
function curlᶜ(R, A, B, 𝐛)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    tmp = [-B[i,j]*𝐭[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

#%%%%

# Old: {d̃ivᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
# Old => new: -d̃ivᶜ => -divᶜ
# New: {divᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
function divᶜ(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    T = findCellLinkLengths(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    a = findCellAreas(R, A, B)
    divᶜ = [-B[i,j]*(F[j]/(T[j]^2)).*transpose(𝐓[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return divᶜ
end 
function divᶜ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    T = findCellLinkLengths(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [-B[i,j]*(F[j]/(T[j]^2))*𝐓[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end 

# Old: {C̃URLᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵᵢ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
# Old => new: C̃URLᶜ => codᶜ
# New: {codᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵᵢ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
function codᶜ(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    a = findCellAreas(R, A, B)
    codᶜ = [B[i,j]*(F[j]/(T[j]^2)).*transpose(ϵᵢ*𝐓[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return codᶜ
end
function codᶜ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [B[i,j]*(F[j]/(T[j]^2))*(ϵᵢ*𝐓[j])⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ
# Old => new: -divᵛ => cocurlᵛ
# New: {cocurlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ
function cocurlᵛ(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    𝐓 = findCellLinks(R, A, B)
    cocurlᵛ = [A[j,k].*transpose(ϵₖ*𝐓[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return cocurlᵛ
end
# function cocurlᵛ(R, A, B, 𝐛)
#     E = findCellLinkTriangleAreas(R, A, B)
#     ϵₖ = SMatrix{2, 2, Float64}([
#                 0.0 -1.0
#                 1.0 0.0
#             ])
#     ϵ = fill(ϵₖ, size(A,2))
#     𝐓 = findCellLinks(R, A, B)
#     tmp = [A[j,k]*(ϵ[k]*𝐓[j])⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
#     return dropdims(sum(tmp, dims=1), dims=1)
# end
# With boundary considerations, Jensen and Revell 2023 Eq 12
function cocurlᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    tmp = [-B[i,j]*A[j,k].*(ϵₖ*q[i,k])⋅𝐛[j]./E[k] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,2)), dims=(1,2))
end

# Old: {CURLᵛ 𝐛 }ₖ = ∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ
# Old => new: CURLᵛ => -curlᵛ
# New: {curlᵛ 𝐛}ₖ = -∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ
function curlᵛ(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    curlᵛ = [-A[j,k].*transpose(𝐓[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return curlᵛ
end 
# function curlᵛ(R, A, B, 𝐛)
#     E = findCellLinkTriangleAreas(R, A, B)
#     𝐓 = findCellLinks(R, A, B)
#     tmp = [-A[j,k]*(𝐓[j]⋅𝐛[j])/E[k] for j=1:size(A,1), k=1:size(A,2)]
#     return dropdims(sum(tmp, dims=1), dims=1)
# end 
# With boundary considerations, Jensen and Revell 2023 Eq 12
function curlᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    q = findSpokes(R, A, B)
    tmp = [B[i,j]*A[j,k]*(q[i,k]⋅𝐛[j])/E[k] for i=1:size(B,1), j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=(1,2)), dims=(1,2))
end 

export gradᵛ
export cogᵛ
export corotᶜ
export rotᶜ
export gradᶜ
export cogᶜ
export corotᵛ
export rotᵛ
export divᵛ
export codᵛ
export cocurlᶜ
export curlᶜ
export divᶜ
export codᶜ
export cocurlᵛ
export curlᵛ

end
#
#  DifferentialOperators.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
# A set of functions to derive objects that depend on system topology and spatial information 

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

using SparseArrays
using StaticArrays
using LinearAlgebra
using DiscreteCalculus

# Old: {gradᵛ ϕ }ⱼ = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
# Old => new: gradᵛ => gradᵛ
# New: {gradᵛ ϕ }ⱼ = = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
function gradᵛ(R, A)
    𝐭 = spdiagm(findEdgeTangents(R, A))
    t = spdiagm(1.0./findEdgeLengths(R, A))
    gradᵛ = 𝐭*t*t*A
    return gradᵛ
end
function gradᵛ(R, A, ϕ)    
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [A[j,k].*(𝐭[j]./(t[j]^2)).*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {curlᶜ 𝐛 }ᵢ = ∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
# Old => new: curlᶜ => -curlᶜ
# New: {curlᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
function curlᶜ(R, A, B)
    a = spdiagm(1.0./findCellAreas(R, A, B))
    𝐭 = spdiagm(transpose.(findEdgeTangents(R, A)))
    curlᶜ = -a*B*𝐭
    return curlᶜ
end
function curlᶜ(R, A, B, 𝐛)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    tmp = [-B[i,j]*𝐭[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
# Old => new: gradᶜ => gradᶜ
# New: {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
function gradᶜ(R, A, B)
    𝐓 = spdiagm(findCellCentreLinks(R, A, B))
    T = spdiagm(1.0./findCellCentreLinkLengths(R, A, B))
    Bᵀ = transpose(B)
    gradᶜ = 𝐓*T*T*Bᵀ
    return gradᶜ
end
function gradᶜ(R, A, B, f)
    𝐓 = findCellCentreLinks(R, A, B)
    T = findCellCentreLinkLengths(R, A, B)
    tmp = [B[i,j].*(𝐓[j]./(T[j]^2)).*f[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {CURLᵛ 𝐛 }ₖ = ∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ
# Old => new: CURLᵛ => -curlᵛ
# New: {curlᵛ 𝐛}ₖ = -∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ
function curlᵛ(R, A, B)
    E = spdiagm(1.0./findCellLinkTriangleAreas(R, A, B))
    𝐓 = spdiagm(transpose.(findCellCentreLinks(R, A, B)))
    Aᵀ = transpose(A)
    curlᵛ = -E*Aᵀ*𝐓
    return curlᵛ
end 
function curlᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    𝐓 = findCellCentreLinks(R, A, B)
    tmp = [-A[j,k]*(𝐓[j]⋅𝐛[j])/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end 

# Old: {C̃URLᵛ ϕ }ⱼ = ∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ
# Old => new: C̃URLᵛ => -rotᵛ
# New: {rotᵛ ϕ}ⱼ = -∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ
function rotᵛ(R, A, B)
    F = spdiagm(1.0./(2.0.*findEdgeQuadrilateralAreas(R, A, B)))
    𝐓 = spdiagm(findCellCentreLinks(R, A, B))
    rotᵛ = -𝐓*F*A
    return rotᵛ
end 
function rotᵛ(R, A, B, ϕ)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    𝐓 = findCellCentreLinks(R, A, B)
    tmp = [-A[j,k].*𝐓[j].*ϕ[k]/F[j] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end 

# Old: {d̃ivᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
# Old => new: -d̃ivᶜ => -divᶜ
# New: {divᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
function divᶜ(R, A, B)
    F = spdiagm(2.0.*findEdgeQuadrilateralAreas(R, A, B))
    T = spdiagm(1.0./findCellCentreLinkLengths(R, A, B))
    𝐓 = spdiagm(transpose.(findCellCentreLinks(R, A, B)))
    a = spdiagm(1.0./findCellAreas(R, A, B))
    divᶜ = -a*B*F*T*T*𝐓
    return divᶜ
end 
function divᶜ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    T = findCellCentreLinkLengths(R, A, B)
    𝐓 = findCellCentreLinks(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [-B[i,j]*(F[j]/(T[j]^2))*𝐓[j]⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end 

# Old: {c̃urlᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
# Old => new: c̃urlᶜ => -rotᶜ
# New: {rotᶜ f }ⱼ = -∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
function rotᶜ(R, A, B)
    𝐭 = spdiagm(findEdgeTangents(R, A))
    F = spdiagm(1.0./(2.0.*findEdgeQuadrilateralAreas(R, A, B)))
    Bᵀ = transpose(B)
    rotᶜ = -𝐭*F*Bᵀ
    return rotᶜ
end
function rotᶜ(R, A, B, f)
    𝐭 = findEdgeTangents(R, A)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    tmp = [-B[i,j].*(𝐭[j]./F[j])*f[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {d̃ivᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
# Old => new: -d̃ivᵛ => -divᵛ
# New: {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
function divᵛ(R, A, B)
    F = spdiagm(2.0.*findEdgeQuadrilateralAreas(R, A, B))
    t = spdiagm(1.0./findEdgeLengths(R, A))
    𝐭 = spdiagm(transpose.(findEdgeTangents(R, A)))
    E = spdiagm(1.0./findCellLinkTriangleAreas(R, A, B))
    Aᵀ = transpose(A)
    divᵛ = -E*Aᵀ*F*t*t*𝐭
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

# Old: {curlᵛ ϕ}ⱼ = ∑ₖϵₖ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
# Old => new: curlᵛ => -cogᵛ
# New: {cogᵛ ϕ}ⱼ = -∑ₖϵₖ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
function cogᵛ(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    ϵ = spdiagm(fill(ϵₖ, size(A,1))) # Assume that ϵₖ does not vary with k and therefore we can fill this matrix with size JxJ
    𝐭 = spdiagm(findEdgeTangents(R, A))
    t = spdiagm(1.0./findEdgeLengths(R, A))
    cogᵛ = -ϵ*𝐭*t*t*A
    return cogᵛ
end
function cogᵛ(R, A, B, ϕ)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    ϵ = fill(ϵₖ, size(A,2))
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    tmp = [-ϵ[k]*(𝐭[j]./(t[j]^2)).*A[j,k]*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old {divᶜ 𝐛}ᵢ = -∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
# Old => new: -divᶜ => cocurlᶜ
# New: {cocurlᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
function cocurlᶜ(R, A, B)
    a = spdiagm(1.0./findCellAreas(R, A, B))
    𝐭 = spdiagm(transpose.(findEdgeTangents(R, A)))
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    ϵ = spdiagm(fill(transpose(ϵᵢ), size(B,2)))
    divᶜᵢ = a*B*𝐭*ϵ 
    return divᶜᵢ
end
function cocurlᶜ(R, A, B, 𝐛)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    ϵ = fill(ϵᵢ, size(B,1))
    tmp = [B[i,j]*(ϵ[i]*𝐭[j])⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {CURLᶜ f}ⱼ = ∑ᵢϵᵢ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
# Old => new: CURLᶜ => -cogᶜ 
# New: {cogᶜ f}ⱼ = -∑ᵢϵᵢ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
function cogᶜ(R, A, B)
    𝐓 = spdiagm(findCellCentreLinks(R, A, B))
    T = spdiagm(1.0./findCellCentreLinkLengths(R, A, B))
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    ϵ = spdiagm(fill(ϵᵢ, size(B,2)))
    Bᵀ = transpose(B)
    cogᶜ = -ϵ*𝐓*T*T*Bᵀ
    return cogᶜ
end
function cogᶜ(R, A, B, f)
    𝐓 = findCellCentreLinks(R, A, B)
    T = findCellCentreLinkLengths(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    ϵ = fill(ϵᵢ, size(B,1))
    tmp = [-ϵ[i]*(𝐓[j]./(T[j]^2)).*B[i,j]*f[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ
# Old => new: -divᵛ => cocurlᵛ
# New: {cocurlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ
function cocurlᵛ(R, A, B)
    E = spdiagm(1.0./findCellLinkTriangleAreas(R, A, B))
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    ϵ = spdiagm(fill(transpose(ϵₖ), size(A,1)))
    𝐓 = spdiagm(transpose.(findCellCentreLinks(R, A, B)))
    Aᵀ = transpose(A)
    cocurlᵛ = E*Aᵀ*𝐓*ϵ
    return cocurlᵛ
end
function cocurlᵛ(R, A, B, 𝐛)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    ϵ = fill(ϵₖ, size(A,2))
    𝐓 = findCellCentreLinks(R, A, B)
    tmp = [A[j,k]*(ϵ[k]*𝐓[j])⋅𝐛[j]/E[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {g̃radᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ
# Old => new: g̃radᵛ => corotᵛ
# New: {corotᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ
function corotᵛ(R, A, B)
    F = spdiagm(1.0./(2.0.*findEdgeQuadrilateralAreas(R, A, B)))
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    ϵ = spdiagm(fill(ϵₖ, size(A,1)))
    𝐓 = spdiagm(findCellCentreLinks(R, A, B))
    corotᵛ = ϵ*𝐓*F*A
    return corotᵛ
end
function corotᵛ(R, A, B, ϕ)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    𝐓 = findCellCentreLinks(R, A, B)
    tmp = [A[j,k].*ϵₖ*(𝐓[j]./F[j]).*ϕ[k] for j=1:size(A,1), k=1:size(A,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {C̃URLᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵᵢ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
# Old => new: C̃URLᶜ => codᶜ
# New: {codᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵᵢ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
function codᶜ(R, A, B)
    F = spdiagm(2.0.*findEdgeQuadrilateralAreas(R, A, B))
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    ϵ = spdiagm(fill(transpose(ϵᵢ), size(B,2)))
    𝐓 = spdiagm(transpose.(findCellCentreLinks(R, A, B)))
    T = spdiagm(1.0./findCellCentreLinkLengths(R, A, B))
    a = spdiagm(1.0./findCellAreas(R, A, B))
    codᶜ = a*B*F*T*T*𝐓*ϵ
    return codᶜ
end
function codᶜ(R, A, B, 𝐛)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐓 = findCellCentreLinks(R, A, B)
    T = findCellCentreLinkLengths(R, A, B)
    a = findCellAreas(R, A, B)
    tmp = [B[i,j]*(F[j]/(T[j]^2))*(ϵᵢ*𝐓[j])⋅𝐛[j]/a[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=2), dims=2)
end

# Old: {g̃radᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
# Old => new: g̃radᶜ => corotᶜ
# New: {corotᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
function corotᶜ(R, A, B)
    F = spdiagm(1.0./(2.0.*findEdgeQuadrilateralAreas(R, A, B)))
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    ϵ = spdiagm(fill(ϵᵢ, size(B,2)))
    𝐭 = spdiagm(findEdgeTangents(R, A))
    Bᵀ = transpose(B)
    corotᶜ = ϵ*𝐭*F*Bᵀ 
    return corotᶜ
end
function corotᶜ(R, A, B, f)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐭 = findEdgeTangents(R, A)
    tmp = [B[i,j].*ϵᵢ*(𝐭[j]./F[j]).*f[i] for i=1:size(B,1), j=1:size(B,2)]
    return dropdims(sum(tmp, dims=1), dims=1)
end

# Old: {c̃urlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵₖ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
# Old => new: c̃urlᵛ => codᵛ 
# New: {codᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵₖ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
function codᵛ(R, A, B)
    E = spdiagm(1.0./findCellLinkTriangleAreas(R, A, B))
    F = spdiagm(2.0.*findEdgeQuadrilateralAreas(R, A, B))
    ϵₖ = SMatrix{2, 2, Float64}([
        0.0 -1.0
        1.0 0.0
    ])
    ϵ = spdiagm(fill(transpose(ϵₖ), size(B,2)))
    𝐭 = spdiagm(transpose.(findEdgeTangents(R, A)))
    t = spdiagm(1.0./findEdgeLengths(R, A))
    Aᵀ = transpose(A)
    codᵛ = E*Aᵀ*F*t*t*𝐭*ϵ
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

export gradᵛ
# export gradᵛ2
export curlᶜ
# export curlᶜ2
export curlᶜ3
export gradᶜ
# export gradᶜ2
export curlᵛ
# export curlᵛ2
export rotᵛ
# export rotᵛ2
export divᶜ
# export divᶜ2
export rotᶜ
# export rotᶜ2
export divᵛ
# export divᵛ2
export cogᵛ
# export cogᵛ2
export cocurlᶜ
# export cocurlᶜ2
export cogᶜ
# export cogᶜ2
export cocurlᵛ
# export cocurlᵛ2
export corotᵛ
# export corotᵛ2
export codᶜ
# export codᶜ2
export corotᶜ
# export corotᶜ2
export codᵛ
# export codᵛ2

end
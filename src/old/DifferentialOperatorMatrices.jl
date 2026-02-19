#
#  DifferentialOperatorMatrices.jl
#  DiscreteCalculus
#
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

module DifferentialOperatorMatrices

using SparseArrays
using StaticArrays
using LinearAlgebra
using DiscreteCalculus

# {cocurlᶜ 𝐛}ᵢ = ∑ⱼBᵢⱼ(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/aᵢ
function cocurlᶜMat(R, A, B)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    cocurlᶜ = [B[i,j].*transpose(ϵᵢ*𝐭[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return cocurlᶜ
end
 
# {divᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)𝐓ⱼ⋅𝐛ⱼ/aᵢ
function divᶜMat(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    T = findCellLinkLengths(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    a = findCellAreas(R, A, B)
    divᶜ = [-B[i,j]*(F[j]/(T[j]^2)).*transpose(𝐓[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return divᶜ
end 

# {divᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)𝐭ⱼ⋅𝐛ⱼ/Eₖ
function divᵛMat(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    t = findEdgeLengths(R, A)
    𝐭 = findEdgeTangents(R, A)
    E = findCellLinkTriangleAreas(R, A, B)
    divᵛ = [-A[j,k]*(F[j]/(t[j]^2)).*transpose(𝐭[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return divᵛ
end

# {cocurlᵛ 𝐛}ₖ = ∑ⱼAⱼₖ(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/Eₖ
function cocurlᵛMat(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    𝐓 = findCellLinks(R, A, B)
    cocurlᵛ = [A[j,k].*transpose(ϵₖ*𝐓[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return cocurlᵛ
end
# With boundary considerations (Jensen and Revell 2023 Eq 12): {cocurlᵛ 𝐛}ₖ = -∑ᵢⱼₖBᵢⱼAⱼₖ(ϵₖqᵢₖ)⋅𝐛ⱼ/Eₖ
function cocurlᵛboundaryMat(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    cocurlᵛboundary = [-B[i,j]*A[j,k].*transpose(ϵₖ*q[i,k])./E[k] for k=1:size(A,2), j=1:size(A,1), i=1:size(B,1)]
    return dropdims(sum(cocurlᵛboundary, dims=3), dims=3)
end

# {curlᶜ 𝐛 }ᵢ = -∑ⱼBᵢⱼ𝐭ⱼ⋅𝐛ⱼ/aᵢ
function curlᶜMat(R, A, B)
    a = findCellAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    curlᶜ = [-B[i,j].*𝐭[j]'./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return curlᶜ
end

# {codivᶜ 𝐛}ᵢ = -∑ⱼBᵢⱼ(Fⱼ/Tⱼ²)(ϵₖ𝐓ⱼ)⋅𝐛ⱼ/aᵢ
function codivᶜMat(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
        0.0 -1.0
        1.0 0.0
    ])
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    a = findCellAreas(R, A, B)
    codivᶜ = [-B[i,j]*(F[j]/(T[j]^2)).*transpose(ϵₖ*𝐓[j])./a[i] for i=1:size(B,1), j=1:size(B,2)]
    return codivᶜ
end

# {codivᵛ 𝐛}ₖ = -∑ⱼAⱼₖ(Fⱼ/tⱼ²)(ϵᵢ𝐭ⱼ)⋅𝐛ⱼ/Eₖ
function codivᵛMat(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    codivᵛ = [-A[j,k]*(F[j]/(t[j]^2)).*transpose(ϵᵢ*𝐭[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return codivᵛ
end

# {curlᵛ 𝐛}ₖ = -∑ⱼAⱼₖ𝐓ⱼ⋅𝐛ⱼ/Eₖ
function curlᵛMat(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    curlᵛ = [-A[j,k].*transpose(𝐓[j])./E[k] for k=1:size(A,2), j=1:size(A,1)]
    return curlᵛ
end 
# With boundary considerations (Jensen and Revell 2023 Eq 12): {curlᵛ 𝐛}ₖ = ∑ᵢⱼₖBᵢⱼAⱼₖqᵢₖ⋅𝐛ⱼ/Eₖ
function curlᵛboundaryMat(R, A, B)
    E = findCellLinkTriangleAreas(R, A, B)
    q = findSpokes(R, A, B)
    curlᵛboundary = [B[i,j]*A[j,k].*transpose(q[i,k])./E[k] for k=1:size(A,2), j=1:size(A,1), i=1:size(B,1)]
    return dropdims(sum(curlᵛboundary, dims=3), dims=3)
end 

# {gradᶜ f }ⱼ = ∑ᵢBᵢⱼ(𝐓ⱼ/Tⱼ²)fᵢ
function gradᶜMat(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    gradᶜ = [B[i,j].*(𝐓[j]./(T[j]^2)) for j=1:size(B,2), i=1:size(B,1)]
    return gradᶜ
end

# {cogradᶜ f}ⱼ = ∑ᵢϵₖ(𝐓ⱼ/Tⱼ²)Bᵢⱼfᵢ
function cogradᶜMat(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    cogradᶜ = [ϵₖ*(𝐓[j]./(T[j]^2)).*B[i,j] for j=1:size(B,2), i=1:size(B,1)]
    return cogradᶜ
end

# {gradᵛ ϕ }ⱼ = = ∑ₖAⱼₖ(𝐭ⱼ/tⱼ²)ϕₖ
function gradᵛMat(R, A)    
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    gradᵛ = [A[j,k].*(𝐭[j]./(t[j]^2)) for j=1:size(A,1), k=1:size(A,2)]
    return gradᵛ
end

# {cogradᵛ ϕ}ⱼ = ∑ₖϵᵢ(𝐭ⱼ/tⱼ²)Aⱼₖϕₖ
function cogradᵛMat(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
                0.0 1.0
                -1.0 0.0
            ])
    𝐭 = findEdgeTangents(R, A)
    t = findEdgeLengths(R, A)
    cogradᵛ = [ϵᵢ*(𝐭[j]./(t[j]^2)).*A[j,k] for j=1:size(A,1), k=1:size(A,2)]
    return cogradᵛ
end

# {corotᶜ f}ⱼ = ∑ᵢBᵢⱼϵᵢ(𝐭ⱼ/Fⱼ)fᵢ
function corotᶜMat(R, A, B)
    ϵᵢ = SMatrix{2, 2, Float64}([
        0.0 1.0
        -1.0 0.0
    ])
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    corotᶜ = [B[i,j].*ϵᵢ*(𝐭[j]./F[j]) for j=1:size(B,2), i=1:size(B,1)]
    return corotᶜ
end

# {rotᶜ f }ⱼ = -∑ᵢBᵢⱼ(𝐭ⱼ/Fⱼ)fᵢ
function rotᶜMat(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    rotᶜ = [-B[i,j].*(𝐭[j]./F[j]) for j=1:size(B,2), i=1:size(B,1)]
    return rotᶜ
end

# {corotᵛ ϕ}ⱼ = ∑ₖAⱼₖϵₖ(𝐓ⱼ/Fⱼ)ϕₖ
function corotᵛMat(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    T = findCellLinkLengths(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    corotᵛ = [A[j,k].*(ϵₖ*𝐓[j])./F[j] for k=1:size(A,2), j=1:size(A,1)]
    return corotᵛ
end
# With boundary considerations (Jensen and Revell 2023 Eq 12): {corotᵛ ϕ}ⱼ = -∑ᵢⱼₖBᵢⱼAⱼₖϵₖ𝐪ᵢₖϕₖ/Fⱼ
function corotᵛboundaryMat(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    ϵₖ = SMatrix{2, 2, Float64}([
                0.0 -1.0
                1.0 0.0
            ])
    q = findSpokes(R, A, B)
    corotᵛboundary = [-B[i,j]*A[j,k].*ϵₖ*q[i,k]./F[j] for j=1:size(A,1), k=1:size(A,2), i=1:size(B,1)]
    return dropdims(sum(corotᵛboundary, dims=3), dims=3)
end

# {rotᵛ ϕ}ⱼ = -∑ₖAⱼₖ𝐓ⱼϕₖ/Fⱼ
function rotᵛMat(R, A, B)
    𝐓 = findCellLinks(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    rotᵛ = [-A[j,k].*𝐓[j]./F[j] for j=1:size(A,1), k=1:size(A,2)]
    return rotᵛ
end 
# With boundary considerations (Jensen and Revell 2023 Eq 12): {rotᵛ ϕ}ⱼ = ∑ᵢⱼₖBᵢⱼAⱼₖ𝐪ᵢₖϕₖ/Fⱼ
function rotᵛboundaryMat(R, A, B)
    F = 2.0.*findEdgeQuadrilateralAreas(R, A, B)
    q = findSpokes(R, A, B)
    rotᵛboundary = [B[i,j]*A[j,k].*q[i,k]./F[j] for j=1:size(A,1), k=1:size(A,2), i=1:size(B,1)]
    return dropdims(sum(rotᵛboundary, dims=3), dims=3)
end 


export gradᵛMat
export cogradᵛMat
export corotᶜMat
export rotᶜMat
export gradᶜMat
export cogradᶜMat
export corotᵛMat
export corotᵛboundaryMat
export rotᵛMat
export rotᵛboundaryMat
export divᵛMat
export codivᵛMat
export cocurlᶜMat
export curlᶜMat
export divᶜMat
export codivᶜMat
export cocurlᵛMat
export cocurlᵛboundaryMat
export curlᵛMat
export curlᵛboundaryMat

end
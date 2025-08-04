#
#  GeometryFunctions.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
# A set of functions to derive objects that depend on both system topology and vertex positions
# K => Number of vertices in network
# J => Number of edges in network
# I => Number of cells in network

# 𝐭ⱼ    Aⱼₖ𝐫ₖ       findEdgeTangents
# tⱼ    |𝐭ⱼ|        findEdgeLengths
# 𝐜ⱼ    Āⱼₖ𝐫ₖ/2     findEdgeMidpoints
# 𝐑ᵢ    Cᵢₖ𝐫ₖ/Zᵢ    indCellCentresOfMass
# lᵢ    B̄ᵢⱼtⱼ       findCellPerimeterLengths
# 𝐚ᵢ                findCellPolygons
# aᵢ                findCellAreas
# 𝐓ⱼ    ∑ᵢBᵢⱼ(𝐑ᵢ-𝐜ᵖⱼ), cᵖⱼ = ∑ᵢBᵢⱼ𝐜ⱼ    findCellLinks
# 𝐂ⱼ                findCellLinkMidpoints
# Cⱼ                findCellLinkLengths
# 𝐄ₖ findCellLinkTriangles
# Eₖ findCellLinkTriangleAreas
# 𝐅ⱼ                findEdgeQuadrilaterals
# Fⱼ/2              findEdgeQuadrilateralAreas
# 𝐪ᵢₖ               findSpokes
#                   findEdgeMidpointCellPolygons
# 𝐬ᵢₖ findEdgeMidpointLinks
# Dₖ findEdgeMidpointLinkVertexAreas
#  findEdgeLinkMidpoints

module GeometryFunctions

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using StaticArrays
using GeometryBasics
using FromFile
using CircularArrays

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "TopologyFunctions.jl" using TopologyFunctions

# 𝐭ⱼ = Aⱼₖ𝐫ₖ
# Returns a vector of SVectors corresponding to the tangent vectors of each edge in the network, with length equal to the length of the edge
findEdgeTangents(R, A) = A*R

# tⱼ = |𝐭ⱼ| 
# Returns a vector of floats corresponding to the lengths of each edge in the network
findEdgeLengths(R, A) = norm.(A*R)

# 𝐜ⱼ = Āⱼₖ𝐫ₖ/2
# Returns a vector of SVectors corresponding to the midpoint locations of each edge in the network
findEdgeMidpoints(R, A) = 0.5.*abs.(A)*R

# 𝐑ᵢ = Cᵢₖ𝐫ₖ/Zᵢ
# Returns a vector of SVectors corresponding to the centre of mass locations of each cell in the network, assuming that cells have mass only in their vertices, and each vertex has the same mass. 
findCellCentresOfMass(R, A, B) = findC(A, B)*R./findCellEdgeCount(B)

# lᵢ = B̄ᵢⱼtⱼ 
# Returns a vector of floats corresponding to the perimeter lengths of each cell in the network
findCellPerimeterLengths(R, A, B) = abs.(B)*norm.(A*R)

# Returns a vector of polygons for each cell, where each polygon is a vector of Point{2,Float64} objects from GeometryBasics.jl. This construction is primarily useful for plotting or for area calculations.
function findCellPolygons(R, A, B)
    cellPolygons = Vector{Point{2,Float64}}[]
    for i = 1:size(B, 1)
        orderedVertices, orderedEdges = orderAroundCell(A, B, i)
        push!(cellPolygons, Point{2, Float64}.(R[orderedVertices]))
    end
    return cellPolygons
end

# aᵢ
# Returns a vector of floats corresponding to the areas of each cell in the network, using the GeometryBasics area() function for simplicity
function findCellAreas(R, A, B)
    cellPolygons = findCellPolygons(R, A, B)
    return abs.(area.(cellPolygons))
end

# 𝐓ⱼ = ∑ᵢBᵢⱼ(𝐑ᵢ-𝐜ᵖⱼ), cᵖⱼ = ∑ᵢBᵢⱼ𝐜ⱼ
# Returns a vector of SVectors, indexed by edge j, corresponding to the vector separating adjacent cell centres, or cell centres and edge midpoints at the periphery
function findCellLinks(R, A, B)
    B̄ = abs.(B)
    𝐑ᵢ = findCellCentresOfMass(R, A, B)
    𝐜ⱼ = findEdgeMidpoints(R, A)
    𝐜ᵖ = dropdims(sum([B̄[i,j].*𝐜ⱼ[j] for i=1:size(B,1), j=1:size(B,2)], dims=1), dims=1)
    𝐓ⱼ = dropdims(sum([B[i,j].*(𝐑ᵢ[i] .- 𝐜ᵖ[j]) for i=1:size(B,1), j=1:size(B,2)], dims=1), dims=1)
    return 𝐓ⱼ
end

# 𝐂ⱼ
# Returns a vector of SVectors, indexed by edge j, corresponding to midpoint position of cell links 𝐓ⱼ
function findCellLinkMidpoints(R, A, B)
    𝐓ⱼ = findCellLinks(R, A, B)
    𝐑ᵢ = findCellCentresOfMass(R, A, B)
    𝐂ⱼ = []
    for j=1:size(B,2)
        i = findall(x->x!=0, B[:,j])[1]
        push!(𝐂ⱼ, 𝐑ᵢ[i].-0.5.*B[i,j].*𝐓ⱼ[j])
    end
    return 𝐂ⱼ
end

# Tⱼ
# Returns a vector of floats, indexed by edge j, corresponding to the lengths of cell links 𝐓ⱼ
function findCellLinkLengths(R, A, B)
    T = findCellLinks(R, A, B)
    return norm.(T)
end

# Conventional variable name 𝐄ₖ
# Function that returns a vector with length equal to the number of vertices, 
# where each components is a vector of Point objects forming the triangle formed around the corresonding vertex 
# by the lines connecting adjacent cell centroids. In the case of peripheral vertices, the function instead forms 
# a polygon around the vertex from adjacent cell centroids and peripheral edge midpoints.
function findCellLinkTriangles(R, A, B)
    boundaryVertices = findBoundaryVertices(A, B)
    boundaryEdges = findBoundaryEdges(B)
    C = findC(A, B)
    𝐑ᵢ = findCellCentresOfMass(R, A, B)
    𝐜ⱼ = findEdgeMidpoints(R, A)
    linkTriangles = Vector{Point{2,Float64}}[]
    for k = 1:size(A,2)
        if boundaryVertices[k] == 0
            # If this vertex is not at the system boundary, link triangle is easily formed from the positions of surrounding cells
            # Note that there is no need to work out an order for these cells since any ordering of triangle vertices has a valid orientation
            k_is = findall(x -> x != 0, C[:, k])
            push!(linkTriangles, Point{2,Float64}.(𝐑ᵢ[k_is]))
        else
            # If this vertex is at the system boundary, we must form a more complex kite from surrounding cell centres and midpoints of surrounding boundary edges
            k_js = [j for j in findall(x -> x != 0, A[:, k]) if boundaryEdges[j]==1]
            k_is = [findall(x->x!=0, B[:, j])[1] for j in k_js]
            kiteVertices = [R[k], 𝐜ⱼ[k_js[1]], 𝐑ᵢ[k_is]..., 𝐜ⱼ[k_js[2]]]
            push!(linkTriangles, Point{2,Float64}.(kiteVertices))
        end
    end
    return linkTriangles
end

# Eₖ
function findCellLinkTriangleAreas(R, A, B)
    Eₖ = findCellLinkTriangles(R, A, B)
    return abs.(area.(Eₖ))
end

# Function that returns a vector of floats with length J where each element corresponds to an edge and gives the area of the quadrilateral 
# formed by the vertices at each end of the edge and the centroids of the cells adjacent to the edge. In the case of a peripheral edge this is instead 
# a triangle formed with only one cell centroid 
# This area is equivalent to 0.5*Fⱼ in the conventional naming convention 
function findEdgeQuadrilaterals(R, A, B)
    J = size(B, 2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    edgeQuadrilaterals = Vector{Point{2,Float64}}[]
    for j = 1:J
        edgeCells = findall(x -> x != 0, B[:, j])
        edgeVertices = findall(x -> x != 0, A[j, :])
        if length(edgeCells) > 1
            push!(edgeQuadrilaterals, Point{2,Float64}.([R[edgeVertices[1]], cellCentresOfMass[edgeCells[1]], R[edgeVertices[2]], cellCentresOfMass[edgeCells[2]]]))
        else
            push!(edgeQuadrilaterals, Point{2,Float64}.([R[edgeVertices[1]], cellCentresOfMass[edgeCells[1]], R[edgeVertices[2]]]))
        end
    end
    return edgeQuadrilaterals
end

# 0.5*Fⱼ
function findEdgeQuadrilateralAreas(R, A, B)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    return abs.(area.(edgeQuadrilaterals))
end

# 𝐪ᵢₖ = C̄ᵢₖ*(𝐫ₖ - 𝐑ᵢ)
function findSpokes(R, A, B)
    K = size(A, 2)
    I = size(B, 1)
    C = abs.(B) * abs.(A) .÷ 2
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    𝐪ᵢₖ = Matrix{SVector{2,Float64}}(undef, I, K)
    for i = 1:I
        for k = 1:K
            𝐪ᵢₖ[i, k] = abs(C[i, k]) .* (R[k] .- cellCentresOfMass[i])
        end
    end
    return 𝐪ᵢₖ
end


function findEdgeMidpointCellPolygons(R, A, B)
    I = size(B, 1)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeMidpointPolygons = Vector{Point{2,Float64}}[]
    for i = 1:I
        orderedVertices, orderedEdges = orderAroundCell(A, B, i)
        push!(edgeMidpointPolygons, Point{2,Float64}.(edgeMidpoints[orderedEdges]))
    end
    return edgeMidpointPolygons
end

# 𝐬ᵢₖ = ∑ⱼ(1/2)BᵢⱼtⱼĀⱼₖ
function findEdgeMidpointLinks(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    tmp = [0.5*B[i,j]*𝐭[j]*abs(A[j,k]) for i=1:size(B,1), j=1:size(B,2), k=1:size(A,2)]
    𝐬ᵢⱼ = sparse(dropdims(sum(tmp, dims=2), dims=2))
    return 𝐬ᵢⱼ
end

#!!!!!!!!!!!!!!!
# Function that returns a vector of floats corresponding to the areas surrounding each vertex in the network.
# For internal vertices this area is bounded by the lines connecting adjacent edge midpoints. For peripheral vertices, 
# it is a the area of the quadrilateral formed by adjacent edge midpoint links and 2 adjacent peripheral edges.
# Conventional variable name Dₖ
function findEdgeMidpointLinkVertexAreas(R, A, B)
    K = size(A,2)
    edgeTangents = findEdgeTangents(R, A)
    𝐬ᵢₖ = findEdgeMidpointLinks(R, A, B)
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2 # (NB Integer division)
    dropzeros!(C)
    vertexAreas = zeros(K)
    for k=1:K
        k_is = findall(x->x!=0, C[:,k]) # Cells i surrounding vertex k
        if length(k_is) == 1
            # If peripheral vertex with only one adjacent cell 
            k_js = findall(x->x!=0, A[:,k]) # Edges j around vertex k
            vertexAreas[k] = 0.5^3*norm([edgeTangents[k_js[1]]...,0.0]×[edgeTangents[k_js[2]]...,0.0]) # Triangle area from cross product of adjacent edge tangents
        elseif length(k_is) == 2
            # If peripheral vertex with 2 adjacent cells 
            edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
            edgesSharedBy_i2_And_k = findall(x->x!=0, B[k_is[2],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] += 0.5^3*norm([edgeTangents[edgesSharedBy_i2_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i2_And_k[2]]...,0.0])
        else
            # If internal vertex with 3 adjacent cells 
            vertexAreas[k] = 0.5*norm([𝐬ᵢₖ[k_is[1], k]...,0.0]×[𝐬ᵢₖ[k_is[2],k]...,0.0])
        end
    end
    return vertexAreas
end

#!!!!!!!!!!!!!!!
function findEdgeLinkMidpoints(R, A, B, ϵᵢ)
    J = size(B, 2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    edgeTangents = findEdgeTangents(R, A)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    T = findCellLinks(R, A, B)
    # Rotation matrix around vertices is the opposite of that around cells
    ϵₖ = -1 * ϵᵢ
    onesVec = ones(1, I)
    boundaryEdges = abs.(onesVec * B)
    cᵖ = boundaryEdges' .* edgeMidpoints
    intersections = SVector{2,Float64}[]
    for j = 1:J
        if boundaryEdges[j] == 0
            k = findall(x -> x < 0, A[j, :])[1]
            i = findall(x -> x < 0, B[:, j])[1]
            mⱼ = R[k] .+ ((cellCentresOfMass[i] .- R[k]) ⋅ (ϵₖ * T[j])) / (2.0 * trapeziumAreas[j]) .* edgeTangents[j]
            push!(intersections, mⱼ)
        else
            push!(intersections, cᵖ[j])
        end
    end
    return intersections
end

export findEdgeTangents
export findEdgeLengths
export findEdgeMidpoints
export findCellCentresOfMass
export findCellPerimeterLengths
export findCellPolygons
export findCellAreas
export findCellLinks
export findCellLinkMidpoints
export findCellLinkLengths
export findCellLinkTriangles
export findCellLinkTriangleAreas
export findEdgeQuadrilaterals
export findEdgeQuadrilateralAreas
export findSpokes
export findEdgeMidpointCellPolygons
export findEdgeMidpointLinks
export findEdgeMidpointLinkVertexAreas
export findEdgeLinkMidpoints

end 
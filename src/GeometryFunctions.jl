#
#  GeometryFunctions.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
# A set of functions to derive objects that depend on system topology and spatial information 
# nVerts => Number of vertices in network
# nEdges => Number of edges in network
# nCells => Number of cells in network

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

# Function that returns a vector of SVectors corresponding to the tangent vectors of each edge in the network, with length equal to the length of the edge
# Conventional variable name 𝓉ⱼ
# 𝐭 = AR
findEdgeTangents(R, A) = A*R

# Function that returns a vector of floats corresponding to the lengths of each edge in the network
# Conventional variable name ???
# t = |AR| 
findEdgeLengths(R, A) = norm.(A*R)

# Function that returns a vector of SVectors corresponding to the midpoint locations of each edge in the network
# Conventional variable name 𝒸ⱼ
# = ĀR/2
findEdgeMidpoints(R, A) = 0.5.*abs.(A)*R

# Function that returns a vector of SVectors corresponding to the centre of mass locations of each cell in the network,
# assuming that cells have mass only in their vertices, and each vertex has the same mass. 
# Conventional variable name ???
# = CR/n
findCellCentresOfMass(R, A, B) = findC(A, B)*R./findCellEdgeCount(B)

# Function that returns a vector of floats corresponding to the perimeter lengths of each cell in the network
# Conventional variable name ???
# = B̄t 
findCellPerimeterLengths(R, A, B) = abs.(B)*norm.(A*R)

# Function that returns a vector of floats corresponding to the areas surrounding each vertex in the network.
# For internal vertices this area is bounded by the lines connecting adjacent edge midpoints. For peripheral vertices, 
# it is a the area of the quadrilateral formed by adjacent edge midpoint links and 2 adjacent peripheral edges.
# Conventional variable name Dₖ
function findVertexAreas(R, A, B)
    nVerts = size(A,2)
    edgeTangents = findEdgeTangents(R, A)
    edgeMidpointLinks = findEdgeMidpointLinks(R, A, B)
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2 # (NB Integer division)
    dropzeros!(C)
    vertexAreas = zeros(nVerts)
    for k=1:nVerts
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
            vertexAreas[k] = 0.5*norm([edgeMidpointLinks[k_is[1], k]...,0.0]×[edgeMidpointLinks[k_is[2],k]...,0.0])
        end
    end

    return vertexAreas
end

# Function that returns a vector of polygons for each cell, where each polygon is a vector of Point{2,Float64} objects from GeometryBasics.jl.
# This construction is primarily useful for plotting or for area calculations.
function findCellPolygons(R, A, B)
    cellPolygons = Vector{Point{2,Float64}}[]
    for i = 1:size(B, 1)
        orderedVertices, orderedEdges = orderAroundCell(A, B, i)
        push!(cellPolygons, Point{2, Float64}.(R[orderedVertices]))
    end
    return cellPolygons
end

# Function that returns a vector of floats corresponding to the areas of each cell in the network, using the GeometryBasics area() function for simplicity
# Conventional variable name Aᵢ
function findCellAreas(R, A, B)
    cellPolygons = findCellPolygons(R, A, B)
    return abs.(area.(cellPolygons))
end

# Function that returns a vector of SVectors with nEdges components, where each component 
# corresponds to the vector connecting cell centres of mass across the corresponding edge. 
# Conventional variable name Tⱼ
# function findCellLinks(R, A, B)
#     nCells = size(B, 1)
#     nEdges = size(B, 2)
#     cellCentresOfMass = findCellCentresOfMass(R, A, B)
#     edgeMidpoints = findEdgeMidpoints(R, A)
#     boundaryEdges = findBoundaryEdges(B)
#     cᵖ = boundaryEdges' .* edgeMidpoints
#     T = SVector{2,Float64}[]
#     for j = 1:nEdges
#         Tⱼ = @SVector zeros(2)
#         for i = 1:nCells
#             Tⱼ = Tⱼ + B[i, j] * (cellCentresOfMass[i] .- cᵖ[j])
#         end
#         push!(T, Tⱼ)
#     end
#     return T
# end

# 𝐓ⱼ = ∑ᵢBᵢⱼ(𝐑ᵢ-𝐜ᵖⱼ), cᵖⱼ = ∑ᵢBᵢⱼ𝐜ⱼ
function findCellLinks(R, A, B)
    B̄ = abs.(B)
    𝐑ᵢ = findCellCentresOfMass(R, A, B)
    𝐜ⱼ = findEdgeMidpoints(R, A)
    𝐜ᵖ = dropdims(sum([B̄[i,j].*𝐜ⱼ[j] for i=1:size(B,1), j=1:size(B,2)], dims=1), dims=1)
    𝐓ⱼ = dropdims(sum([B[i,j].*(𝐑ᵢ[i] .- 𝐜ᵖ[j]) for i=1:size(B,1), j=1:size(B,2)], dims=1), dims=1)
    return 𝐓ⱼ
end

function findCellLinkLengths(R, A, B)
    T = findCellLinks(R, A, B)
    return norm.(T)
end

# Function that returns a vector with length equal to the number of vertices, 
# where each components is a vector of Point objects forming the triangle formed around the corresonding vertex 
# by the lines connecting adjacent cell centroids. In the case of peripheral vertices, the function instead forms 
# a polygon around the vertex from adjacent cell centroids and peripheral edge midpoints.
# Conventional variable name Eₖ
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

function findCellLinkTriangleAreas(R, A, B)
    Eₖ = findCellLinkTriangles(R, A, B)
    return abs.(area.(Eₖ))
end

# Function that returns a vector of floats with length nEdges where each element corresponds to an edge and gives the area of the quadrilateral 
# formed by the vertices at each end of the edge and the centroids of the cells adjacent to the edge. In the case of a peripheral edge this is instead 
# a triangle formed with only one cell centroid 
# This area is equivalent to 0.5*Fⱼ in the conventional naming convention 
function findEdgeQuadrilaterals(R, A, B)
    nEdges = size(B, 2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    edgeQuadrilaterals = Vector{Point{2,Float64}}[]
    for j = 1:nEdges
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

function findEdgeQuadrilateralAreas(R, A, B)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    return abs.(area.(edgeQuadrilaterals))
end

function findSpokes(R, A, B)
    nVerts = size(A, 2)
    nCells = size(B, 1)
    C = abs.(B) * abs.(A) .÷ 2
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    q = Matrix{SVector{2,Float64}}(undef, nCells, nVerts)
    for i = 1:nCells
        for k = 1:nVerts
            q[i, k] = abs(C[i, k]) .* (R[k] .- cellCentresOfMass[i])
        end
    end
    return q
end





function findEdgeMidpointPolygons(R, A, B)
    nCells = size(B, 1)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeMidpointPolygons = Vector{Point{2,Float64}}[]
    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, i)
        push!(edgeMidpointPolygons, Point{2,Float64}.(edgeMidpoints[orderedEdges]))
    end
    return edgeMidpointPolygons
end

# sᵢₖ = ∑ⱼ(1/2)BᵢⱼtⱼĀⱼₖ
function findEdgeMidpointLinks(R, A, B)
    spzeros(SVector{2,Float64}, size(B,1), size(A,2))
    Ā = abs.(A)
    C = findC(A, B)
    edgeTangents = findEdgeTangents(R, A)
    cellVertexOrders  = fill(CircularVector(Int64[]), size(B,1)),
    for i = 1:size(B,1)
        cellVertexOrders[i], cellEdgeOrders[i] = orderAroundCell(A, B, i)
    end
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1], nzC[2])
    for (i, k) in ikPairs
        for j in cellEdgeOrders[i]
            edgeMidpointLinks[i, k] = edgeMidpointLinks[i, k] .+ 0.5 .* B[i, j] .* edgeTangents[j] .* Ā[j, k]
        end
    end
    return edgeMidpointLinks
end

function findEdgeLinkMidpoints(R, A, B, ϵᵢ)
    nEdges = size(B, 2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    edgeTangents = findEdgeTangents(R, A)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    trapeziumAreas = abs.(area.(edgeQuadrilaterals))
    T = findCellLinks(R, A, B)
    # Rotation matrix around vertices is the opposite of that around cells
    ϵₖ = -1 * ϵᵢ
    onesVec = ones(1, nCells)
    boundaryEdges = abs.(onesVec * B)
    cᵖ = boundaryEdges' .* edgeMidpoints
    intersections = SVector{2,Float64}[]
    for j = 1:nEdges
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


export findEdgeMidpointPolygons
export findEdgeLinkMidpoints
export findSpokes
export findEdgeMidpointLinks
export findCellCentresOfMass
export findEdgeTangents
export findEdgeMidpoints
export findCellPolygons
export findCellLinks
export findCellLinkLengths
export findCellPerimeterLengths
export findCellAreas
export findEdgeLengths
export findVertexAreas
export findCellLinkTriangles
export findCellLinkTriangleAreas
export findEdgeQuadrilaterals
export findEdgeQuadrilateralAreas


end #end module 






# function findCellAreas(R, A, B)
#     nCells = size(B,1)
#     Bᵀ = sparse(Transpose(B)) # Have to convert transpose type to sparse type here because nzrange won't operate on a transpose type object
#     edgeTangents = A*R
#     edgeMidpoints = 0.5.*abs.(A)*R
#     # Calculate oriented cell areas    
#     cellOrientedAreas = fill(SMatrix{2,2}(zeros(2,2)),nCells)
#     cellAreas = zeros(nCells)
#     for i=1:nCells
#         for j in nzrange(Bᵀ,i)
#             cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
#         end
#         cellAreas[i] = cellOrientedAreas[i][1,2]
#     end
#     return cellAreas 
# end
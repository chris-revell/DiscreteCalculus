#
#  GeometryFunctions.jl
#  DiscreteCalculus
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
# 𝐄ₖ                findCellLinkVertexTriangles
# Eₖ                findCellLinkVertexTriangleAreas
# 𝐅ⱼ                findEdgeQuadrilaterals
# Fⱼ/2              findEdgeQuadrilateralAreas
# 𝐪ᵢₖ               findSpokes
#                   findEdgeMidpointCellPolygons
# 𝐬ᵢₖ               findEdgeMidpointLinks
# Dₖ                findEdgeMidpointLinkVertexAreas
#                   findEdgeLinkIntersections

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
@from "Constants.jl" using Constants

# 𝐭ⱼ = Aⱼₖ𝐫ₖ
# Returns a vector of SVectors corresponding to the tangent vectors of each 
# edge in the network, with length equal to the length of the edge
findEdgeTangents(R, A) = A*R
findEdgeTangents!(R, A, 𝐭ⱼ) = 𝐭ⱼ .= A*R

# tⱼ = |𝐭ⱼ| 
# Returns a vector of floats corresponding to the lengths of each edge in the network
findEdgeLengths(R, A) = norm.(A*R)
findEdgeLengths!(R, A, tⱼ) = tⱼ .= norm.(A*R)

# 𝐜ⱼ = Āⱼₖ𝐑ₖ/2
# Returns a vector of SVectors corresponding to the midpoint locations of each edge in the network
findEdgeMidpoints(R, A) = 0.5.*abs.(A)*R
findEdgeMidpoints!(R, A, 𝐜ⱼ) = 𝐜ⱼ .= 0.5.*abs.(A)*R

# 𝐧ᵢⱼ = -Bᵢⱼϵᵢ𝐭ⱼ
# Returns a sparse matrix of normal vectors to the surface of each cell at each of that 
# cell's edges, with empty components where B[i,j] = 0
function findCellOutwardNormals(R, A, B)
    𝐭ⱼ = findEdgeTangents(R, A)
    𝐧ᵢⱼ = [-B[i,j]*ϵᵢ*𝐭ⱼ[j] for i=1:size(B,1), j=1:size(B,2)]
    return sparse(𝐧ᵢⱼ)
end

# 𝐓ⱼ = ∑ᵢBᵢⱼ(𝐑ᵢ-𝐜ᵖⱼ), cᵖⱼ = ∑ᵢBᵢⱼ𝐜ⱼ
# Returns a vector of SVectors, indexed by edge j, corresponding to the vector separating adjacent 
# cell centres, or cell centres and edge midpoints at the periphery
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

# Tⱼ = |𝐓ⱼ|
# Returns a vector of floats, indexed by edge j, corresponding to the lengths of cell links 𝐓ⱼ
function findCellLinkLengths(R, A, B)
    T = findCellLinks(R, A, B)
    return norm.(T)
end

# 𝐍ⱼₖ = -Aⱼₖϵₖ𝐓ⱼ
# Returns a sparse matrix of normal vectors to the surface of each triangle formed by cell-cell links around a vertex,
# indexed by edge j and vertex k
function findCellLinkTriangleOutwardNormals(R, A, B)
    𝐓ⱼ = findCellLinks(R, A, B)
    𝐍ⱼₖ = [-A[j,k]*ϵₖ*𝐓ⱼ[j] for j=1:size(A,1), k=1:size(A,2)]
    return sparse(𝐍ⱼₖ)
end

# 𝐫ᵢ = Cᵢₖ𝐑ₖ/Zᵢ
# Returns a vector of SVectors corresponding to the centre of mass locations of each cell in 
# the network, assuming that cells have mass only in their vertices, and each vertex has the same mass. 
findCellCentresOfMass(R, A, B) = findC(A, B)*R./findCellEdgeCount(B)
findCellCentresOfMass!(R, A, B, 𝐫ᵢ) = 𝐫ᵢ .= findC(A, B)*R./findCellEdgeCount(B)

# lᵢ = B̄ᵢⱼtⱼ 
# Returns a vector of floats corresponding to the perimeter lengths of each cell in the network.
findCellPerimeterLengths(R, A, B) = abs.(B)*findEdgeLengths(R, A)

# Returns a vector of polygons for each cell, where each polygon is a vector of 
# Point{2,Float64} objects from GeometryBasics.jl. This construction is primarily 
# useful for plotting or for area calculations.
function findCellPolygons(R, A, B)
    cellPolygons = Vector{Point{2,Float64}}[]
    for i = 1:size(B, 1)
        orderedVertices, orderedEdges = orderAroundCell(A, B, i)
        push!(cellPolygons, Point{2, Float64}.(R[orderedVertices]))
    end
    return cellPolygons
end

# aᵢ
# Returns a vector of floats corresponding to the areas of each cell in the network, 
# using the GeometryBasics area() function for simplicity
function findCellAreas(R, A, B)
    cellPolygons = findCellPolygons(R, A, B)
    return abs.(area.(cellPolygons))
end

# 𝐄ₖ
# Function that returns a vector with length equal to the number of vertices, K, containing polygons 
# defined by the cell links connecting adjacent cells to each vertex. 
# In the case of peripheral vertices, the function instead forms a polygon around 
# the vertex comprising adjacent cell centroids and peripheral edge midpoints.
function findCellLinkVertexTriangles(R, A, B)
    boundaryVertices = findPeripheralVertices(A, B)
    peripheralEdges = findPeripheralEdges(B)
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
            k_js = [j for j in findall(x -> x != 0, A[:, k]) if peripheralEdges[j]==1]
            k_is = [findall(x->x!=0, B[:, j])[1] for j in k_js]
            kiteVertices = [R[k], 𝐜ⱼ[k_js[1]], 𝐑ᵢ[k_is]..., 𝐜ⱼ[k_js[2]]]
            push!(linkTriangles, Point{2,Float64}.(kiteVertices))
        end
    end
    return linkTriangles
end

# Eₖ
# Function that returns areas of cell link vertex triangles. 
function findCellLinkVertexTriangleAreas(R, A, B)
    Eₖ = findCellLinkVertexTriangles(R, A, B)
    return abs.(area.(Eₖ))
end

# 𝐛ⱼ 
# Returns a vector of SVectors corresponding to the intersection of edge 𝐭ⱼ with cell link 𝐓ⱼ
function findEdgeLinkIntersections(R, A, B)
    J = size(B,2)
    𝐭ⱼ = findEdgeTangents(R, A)
    𝐓ⱼ = findCellLinks(R, A, B)
    peripheralEdges = findPeripheralEdges(B).==1
    Rᵢ = findCellCentresOfMass(R, A, B)
    𝐛ⱼ = fill(SVector{2,Float64}(zeros(2)), J)
    for j = 1:J
        # Find vertices and cells surrounding edge j
        j_ks = findall(x -> x != 0, A[j, :])
        j_is = findall(x -> x != 0, B[:, j])
        # mⱼ = R[k] .+ ((cellCentresOfMass[i] .- R[k]) ⋅ (ϵₖ * T[j])) / (2.0 * trapeziumAreas[j]) .* edgeTangents[j]

        # Form Line object for edge
        𝐭Line = Line(Point{2,Float64}.(R[j_ks])...)
        # Form Line object for cell link. 
        # Note that by starting the line at Rᵢ[j_is[1]] and then 
        # ending it at Rᵢ[j_is[1]].-2.0*B[j_is[1],j].*𝐓ⱼ[j]], the 
        # line is guaranteed to intersect with 𝐭Line, even if 𝐭Line is a peripheral edge, 
        # in which case we can't construct the line with a second cell position, and using the edge midpoint will make the line 
        # too short to find an intersection 
        𝐓Line = Line(Point{2,Float64}.([Rᵢ[j_is[1]], Rᵢ[j_is[1]].-2.0*B[j_is[1],j].*𝐓ⱼ[j]])...)
        @show intersects(𝐭Line, 𝐓Line)[1]
        𝐛ⱼ[j] = SVector{2,Float64}(intersects(𝐭Line, 𝐓Line)[2])
    end
    return 𝐛ⱼ
end

# Function that returns a vector polygons corresponding to the quadrilateral formed by vertices 
# and cell centres adjacent to each edge. Note that this is a triangle at the monolayer periphery.
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
# Returns a vector of scalar areas for the quadrilaterals formed by vertices and cell centres adjacent to each edge.
function findEdgeQuadrilateralAreas(R, A, B)
    edgeQuadrilaterals = findEdgeQuadrilaterals(R, A, B)
    return abs.(area.(edgeQuadrilaterals))
end

# 𝐪ᵢₖ = C̄ᵢₖ*(𝐑ₖ - 𝐫ᵢ)
# Returns a vector of vectors corresponding to the spokes connecting cell centres and the corresonding vertices,
# indexed by cell and vertex pairs.
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

# 𝐬ᵢₖ = ∑ⱼ(1/2)BᵢⱼtⱼĀⱼₖ
# Returns a sparse matrix of vectors connecting midpoints of adjacent edges, 
# indexed by the corresponding cell i and vertex k. 
function findEdgeMidpointLinks(R, A, B)
    𝐭 = findEdgeTangents(R, A)
    𝐬ᵢₖ = spzeros(SVector{2,Float64}, size(B,1), size(A,2))
    C = findC(A, B)
    rows = rowvals(C)
    I, K = size(C)
    for k = 1:K
        for ind in nzrange(C, k)
            i = rows[ind]
            ik_js = findall(x->x!=0, B[i,:])∩findall(x->x!=0, A[:,k]) # Edges j shared by cell i and vertex k
            for j in ik_js
                𝐬ᵢₖ[i,k] += 0.5*B[i,j]*𝐭[j]*abs(A[j,k])
            end
        end
    end
    return 𝐬ᵢₖ
end

# Returns a vector of polygons corresponding to the ordered edge midpoints around each cell
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

# 𝐧ᵢₖ = -ϵᵢ𝐬ᵢₖ
# For each triangle formed by edge midpoint links 𝐬ᵢₖ around a vertex k, returns a sparse matrix of 
# vectors normal to the triangle, indexed by corresonding cells i and vertices k. 
function findEdgeMidpointLinkTriangleOutwardNormals(R, A, B)
    𝐬ᵢₖ = findEdgeMidpointLinks(R, A, B)
    𝐧ᵢₖ = [-ϵᵢ*𝐬ᵢₖ[i,k] for i=1:size(B,1), k=1:size(A,2)]
    return sparse(𝐧ᵢₖ)
end

# Dₖ
# Function that returns a vector of floats corresponding to the areas surrounding each vertex in the network.
# For internal vertices this area is bounded by the lines connecting adjacent edge midpoints. For peripheral vertices, 
# it is a the area of the quadrilateral formed by adjacent edge midpoint links and 2 adjacent peripheral edges.
function findEdgeMidpointLinkVertexAreas(R, A, B)
    K = size(A,2)
    𝐜ⱼ = findEdgeMidpoints(R, A)
    C = findC(A, B)
    Dₖ = zeros(K)
    for k=1:K
        k_js = findall(x->x!=0, A[:,k]) # Edges j around vertex k
        k_is = findall(x->x!=0, C[:,k]) # Cells i surrounding vertex k
        if length(k_is) == 1
            # If peripheral vertex with only one adjacent cell 
            polygon = Point{2,Float64}.([R[k], 𝐜ⱼ[k_js[1]], 𝐜ⱼ[k_js[2]]])
            Dₖ[k] = area(polygon)
        elseif length(k_is) == 2
            # If peripheral vertex with 2 adjacent cells 
            sharedEdge = (findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, B[k_is[2],:]))[1]
            notSharedEdges = setdiff(k_js, sharedEdge)
            polygon = Point{2,Float64}.([R[k], 𝐜ⱼ[notSharedEdges[1]], 𝐜ⱼ[sharedEdge], 𝐜ⱼ[notSharedEdges[2]]])
            Dₖ[k] = area(polygon)
        else
            # If internal vertex with 3 adjacent cells 
            polygon = Point{2,Float64}.(𝐜ⱼ[k_js])
            Dₖ[k] = area(polygon)
        end
    end
    return abs.(Dₖ)
end

export findEdgeTangents
export findEdgeTangents!
export findEdgeLengths
export findEdgeLengths!
export findEdgeMidpoints
export findEdgeMidpoints!
export findCellOutwardNormals
export findCellLinks
export findCellLinkMidpoints
export findCellLinkLengths
export findCellLinkTriangleOutwardNormals
export findCellCentresOfMass
export findCellCentresOfMass!
export findCellPerimeterLengths
export findCellPolygons
export findCellAreas
export findCellLinkVertexTriangles
export findCellLinkVertexTriangleAreas
export findEdgeLinkIntersections
export findEdgeQuadrilaterals
export findEdgeQuadrilateralAreas
export findSpokes
export findEdgeMidpointLinks
export findEdgeMidpointCellPolygons
export findEdgeMidpointLinkTriangleOutwardNormals
export findEdgeMidpointLinkVertexAreas

end 
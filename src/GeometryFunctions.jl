#
#  GeometryFunctions.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
#

module GeometryFunctions

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using StaticArrays
using GeometryBasics
# using Random
using FromFile
# using Colors

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell
@from "SpatialData.jl" using SpatialData

# function getRandomColor(seed)
#     Random.seed!(seed)
#     rand(RGB{})
# end

function makeCellPolygons(R, A, B)
    nCells = size(B,1)
    cellPolygons = Vector{Point2f}[]
    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, indexCell)
        push!(cellPolygons, Point2f.(R[orderedVertices]))
    end
    return cellPolygons
end

function makeCellLinks(R,A,B)
    nCells = size(B,1)
    nEdges = size(B,2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B) 
    edgeMidpoints = findEdgeMidpoints(R, A)
    onesVec = ones(1, nCells)
    boundaryEdges = abs.(onesVec * B)
    cᵖ = boundaryEdges' .* edgeMidpoints
    T = SVector{2,Float64}[]
    for j = 1:nEdges
        Tⱼ = @SVector zeros(2)
        for i = 1:nCells
            Tⱼ = Tⱼ + B[i, j] * (cellCentresOfMass[i] .- cᵖ[j])
        end
        push!(T, Tⱼ)
    end
    return T
end

function makeLinkTriangles(R, A, B)
    nCells = size(B,1)
    nVerts = size(A,2)
    Aᵀ = sparse(transpose(A))
    Āᵀ = abs.(Aᵀ)
    C = abs.(B) * abs.(A) .÷ 2
    boundaryVertices = Āᵀ * abs.(sum.(eachcol(B))) .÷ 2
    cellCentresOfMass = findCellCentresOfMass(R, A, B) 
    edgeMidpoints = findEdgeMidpoints(R, A)
    onesVec = ones(1, nCells)
    linkTriangles = Vector{Point2f}[]
    boundaryEdges = abs.(onesVec * B)
    for k = 1:nVerts
        if boundaryVertices[k] == 0
            # If this vertex is not at the system boundary, link triangle is easily formed from the positions of surrounding cells
            vertexCells = findall(x -> x != 0, C[:, k])
            push!(linkTriangles, Point2f.(cellCentresOfMass[vertexCells]))
        else
            # If this vertex is at the system boundary, we must form a more complex kite from surrounding cell centres and midpoints of surrounding boundary edges
            vertexCells = findall(x -> x != 0, C[:, k])
            vertexEdges = findall(x -> x != 0, A[:, k])
            boundaryVertexEdges = intersect(vertexEdges, findall(x -> x != 0, boundaryEdges[1, :]))
            kiteVertices = [edgeMidpoints[boundaryVertexEdges]; cellCentresOfMass[vertexCells]]
            push!(kiteVertices, R[k])
            com = sum(kiteVertices) ./ length(kiteVertices)
            angles = Float64[]
            for p = 1:length(kiteVertices)
                angle = atan((kiteVertices[p] .- com)...)
                push!(angles, angle)
            end
            kiteVertices .= kiteVertices[sortperm(angles)]
            push!(linkTriangles, Point2f.(kiteVertices))
        end
    end
    return linkTriangles
end

function makeEdgeTrapezia(R, A, B)
    nEdges = size(B,2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B) 
    edgeTrapezia = Vector{Point2f}[]
    for j = 1:nEdges
        edgeCells = findall(x -> x != 0, B[:, j])
        edgeVertices = findall(x -> x != 0, A[j, :])
        if length(edgeCells) > 1
            push!(edgeTrapezia, Point2f.([R[edgeVertices[1]], cellCentresOfMass[edgeCells[1]], R[edgeVertices[2]], cellCentresOfMass[edgeCells[2]]]))
        else
            push!(edgeTrapezia, Point2f.([R[edgeVertices[1]], cellCentresOfMass[edgeCells[1]], R[edgeVertices[2]]]))
        end
    end
    return edgeTrapezia
end

function makeEdgeMidpointPolygons(R, A, B)
    nCells = size(B,1)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeMidpointPolygons = Vector{Point2f}[]
    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, indexCell)
        push!(edgeMidpointPolygons, Point2f.(edgeMidpoints[orderedEdges]))
    end
    return edgeMidpointPolygons
end

function makeCellVerticesDict(A, B)
    nCells = size(B,1)
    cellVerticesDict = Dict()
    for i = 1:nCells
        cellVertices, cellEdges = orderAroundCell(A, B, indexCell)
        # Store sorted cell vertices for this cell
        cellVerticesDict[i] = cellVertices
    end
    return cellVerticesDict
end

function findEdgeLinkMidpoints(R, A, B, ϵᵢ)
    nEdges = size(B,2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B) 
    edgeTangents = findEdgeTangents(R, A)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    T = makeCellLinks(R,A,B)
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

function makeSpokes(R, A, B)
    nVerts = size(A,2)
    nCells = size(B,1)
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

export makeCellPolygons
export makeCellLinks
export makeLinkTriangles
export makeEdgeTrapezia
export makeEdgeMidpointPolygons
export makeCellVerticesDict
export findEdgeLinkMidpoints
export makeSpokes

end #end module 
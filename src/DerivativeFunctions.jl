#
#  DerivativeFunctions.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 16/08/2023.
#
#

module DerivativeFunctions

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
@from "GeometryFunctions.jl" using GeometryFunctions

# {curlᶜb}ᵢ
# Calculate curl on each cell
function calculateCellCurls(nCells, R, C, cellCentresOfMass, edgeMidpoints, edgeTangents, cellAreas)
    cellCurls = Float64[]
    for c = 1:nCells
        cellVertices = findall(x -> x != 0, C[c, :])
        vertexAngles = zeros(size(cellVertices))
        for (k, v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v] .- cellCentresOfMass[c])...)
        end
        m = minimum(vertexAngles)
        vertexAngles .-= m
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        cellEdges = findall(x -> x != 0, B[c, :])
        edgeAngles = zeros(size(cellEdges))
        for (k, e) in enumerate(cellEdges)
            edgeAngles[k] = atan((edgeMidpoints[e] .- cellCentresOfMass[c])...)
        end
        edgeAngles .+= (2π - m)
        edgeAngles .= edgeAngles .% (2π)
        cellEdges .= cellEdges[sortperm(edgeAngles)]
        h = @SVector [0.0, 0.0]
        curlSum = 0
        for (i, e) in enumerate(cellEdges)
            h = h + ϵᵢ * F[cellVertices[i], c]
            curlSum += B[c, e] * (h ⋅ edgeTangents[e]) / cellAreas[c]
        end
        push!(cellCurls, curlSum)
    end
    return cellCurls
end

# {divᶜb}ᵢ
# Calculate div on each cell
function calculateCellDivs(nCells, R, B, C, F, cellCentresOfMass, edgeMidpoints, edgeTangents, cellAreas, ϵᵢ)
    cellDivs = Float64[]
    for c = 1:nCells
        cellVertices = findall(x -> x != 0, C[c, :])
        vertexAngles = zeros(size(cellVertices))
        for (k, v) in enumerate(cellVertices)
            vertexAngles[k] = atan((R[v] .- cellCentresOfMass[c])...)
        end
        m = minimum(vertexAngles)
        vertexAngles .-= m
        cellVertices .= cellVertices[sortperm(vertexAngles)]
        cellEdges = findall(x -> x != 0, B[c, :])
        edgeAngles = zeros(size(cellEdges))
        for (k, e) in enumerate(cellEdges)
            edgeAngles[k] = atan((edgeMidpoints[e] .- cellCentresOfMass[c])...)
        end
        edgeAngles .+= (2π - m)
        edgeAngles .= edgeAngles .% (2π)
        cellEdges .= cellEdges[sortperm(edgeAngles)]
        h = @SVector [0.0, 0.0]
        divSum = 0
        for (i, e) in enumerate(cellEdges)
            h = h + ϵᵢ * F[cellVertices[i], c]
            divSum -= B[c, e] * (h ⋅ (ϵᵢ * edgeTangents[e])) / cellAreas[c]
        end
        # divSum *= (-0.5)
        push!(cellDivs, divSum)
    end
    return cellDivs
end

# {divᵛb}ₖ
# Calculate div at each vertex
function calculateVertexDivs(nVerts, R, C, cellCentresOfMass, F, ϵᵢ, q, linkTriangleAreas)
    vertexDivs = Float64[]
    for k = 1:nVerts
        divSum = 0
        vertexCells = findall(x -> x != 0, C[:, k])
        cellAngles = zeros(length(vertexCells))
        for i = 1:length(cellAngles)
            cellAngles[i] = atan((cellCentresOfMass[vertexCells[i]] .- R[k])...)
        end
        vertexCells .= vertexCells[sortperm(cellAngles, rev=true)]
        for i in vertexCells
            divSum += ((ϵᵢ * q[i, k]) ⋅ (ϵᵢ * F[k, i])) / linkTriangleAreas[k]
        end
        push!(vertexDivs, divSum)
    end
    return vertexDivs
end

# {CURLᵛb}ₖ
# Calculate curl at each vertex
function calculateVertexCurls(nVerts, R, C, cellCentresOfMass, F, ϵᵢ, q, linkTriangleAreas)
    vertexCurls = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k = 1:nVerts
        curlSum = 0
        vertexCells = findall(x -> x != 0, C[:, k])
        cellAngles = zeros(length(vertexCells))
        for i = 1:length(cellAngles)
            cellAngles[i] = atan((cellCentresOfMass[vertexCells[i]] .- R[k])...)
        end
        vertexCells .= vertexCells[sortperm(cellAngles, rev=true)]
        for i in vertexCells
            curlSum += (q[i, k] ⋅ (ϵᵢ * F[k, i])) / linkTriangleAreas[k]
        end
        push!(vertexCurls, curlSum)
    end
    return vertexCurls
end

function calculateVertexMidpointCurls(nVerts, nCells, nEdges, A, B, ϵᵢ, intersections, linkTriangleAreas, q)
    vertexMidpointCurls = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k = 1:nVerts
        curlSum = 0
        for i = 1:nCells
            for j = 1:nEdges
                curlSum -= B[i, j] * A[j, k] * (q[i, k] ⋅ intersections[j]) / linkTriangleAreas[k]
            end
        end
        push!(vertexMidpointCurls, curlSum)
    end
    return vertexMidpointCurls
end

function calculateVertexMidpointDivs(nVerts, nCells, nEdges, A, B, ϵᵢ, intersections, linkTriangleAreas, q)
    vertexMidpointDivs = Float64[]
    # Working around a given vertex, an h force space point from a cell is mapped to the next edge anticlockwise from the cell
    for k = 1:nVerts
        divSum = 0
        for i = 1:nCells
            for j = 1:nEdges
                divSum -= B[i, j] * A[j, k] * ((ϵᵢ * q[i, k]) ⋅ intersections[j]) / linkTriangleAreas[k]
            end
        end
        push!(vertexMidpointDivs, divSum)
    end
    return vertexMidpointDivs
end

function calculateCellMidpointDivs(nCells, B, cellAreas, edgeTangents, intersections, q)
    cellMidpointDivs = Float64[]
    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, indexCell)
        divSum = 0
        for j in orderedEdges
            divSum -= B[i, j] * (intersections[j] ⋅ (ϵᵢ * edgeTangents[j])) / cellAreas[i]
        end
        push!(cellMidpointDivs, divSum)
    end
    return cellMidpointDivs
end

function calculateCellMidpointCurls(nCells, B, cellAreas, edgeTangents, intersections, q)
    cellMidpointCurls = Float64[]
    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, indexCell)
        curlSum = 0
        for j in orderedEdges
            curlSum += B[i, j] * (intersections[j] ⋅ edgeTangents[j]) / cellAreas[i]
        end
        push!(cellMidpointCurls, curlSum)
    end
    return cellMidpointCurls
end

export calculateCellCurls
export calculateCellDivs
export calculateVertexDivs
export calculateVertexCurls
export makeCellVerticesDict
export calculateVertexMidpointCurls
export calculateVertexMidpointDivs
export calculateCellMidpointDivs
export calculateCellMidpointCurls

end #end module 
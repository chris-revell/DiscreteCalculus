#
#  HNetwork.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 01/05/2025.
#
#

module HNetwork

# Julia packages
using StaticArrays
using LinearAlgebra
using SparseArrays
using CircularArrays
using InvertedIndices
using FromFile

@from "TopologyFunctions.jl" using TopologyFunctions
@from "GeometryFunctions.jl" using GeometryFunctions
@from "OrderAroundCell.jl" using OrderAroundCell

function hAroundCell!(ii, B, h, ϵ, F, currentNeighbourShell, traversedCells, traversedEdges, cellEdgeOrders, cellVertexOrders)
    # Find an edge belonging to cell ii and shared by any other cell in traversedCells
    if length(traversedCells)>0
        startEdge = (findnz(B[traversedCells, :])[2] ∩ findnz(B[ii, :])[1])[1]
    else 
        startEdge = cellEdgeOrders[ii][1]
    end
    # Find the index of this edge within the ordering of edges around cell ii 
    startInd = findall(x->x==startEdge, cellEdgeOrders[ii])[1]
    # Remember in clockwise ordering, cellEdgeOrders[i][1] precedes cellVertexOrders[1]
    for vertexInd = startInd:(startInd+length(cellVertexOrders[ii])-1)
        h[cellEdgeOrders[ii][vertexInd+1]] = h[cellEdgeOrders[ii][vertexInd]] .+ ϵ*F[cellVertexOrders[ii][vertexInd], ii]
        push!(traversedEdges, cellEdgeOrders[ii][vertexInd]) 
    end
    return nothing 
end

function hNetwork(R, A, B, F)
    I = size(B,1)
    J = size(B,2)
    K = size(A,2)
    peripheralEdges = findPeripheralEdges(B)
    peripheralCells = findnz(B[:, peripheralEdges.==1])[1]
    cellVertexOrders  = fill(CircularVector(Int64[]), I)
    cellEdgeOrders    = fill(CircularVector(Int64[]), I)
    for i = 1:I
        cellVertexOrders[i], cellEdgeOrders[i] = orderAroundCell(A, B, i)
    end
    ϵ = SMatrix{2, 2, Float64}([
            0.0 1.0
            -1.0 0.0
        ])
    Ā = abs.(A)
    B̄ = abs.(B)

    # Ensure we don't start with a boundary cell
    # startCell = rand(collect(1:I)[Not(peripheralCells)])
    startCell = 1 
    traversedCells = Int64[]
    traversedEdges = Int64[]
    cellNeighbourMatrix = B*transpose(B)
    h = fill(SVector{2, Float64}(zeros(2)), J)

    hAroundCell!(startCell, B, h, ϵ, F, Int64[], traversedCells, traversedEdges, cellEdgeOrders, cellVertexOrders)
    push!(traversedCells, startCell)
    while length(traversedCells) < I
        currentNeighbourShell = setdiff(findnz(cellNeighbourMatrix[:, traversedCells])[1], traversedCells)
        for ii in currentNeighbourShell
            hAroundCell!(ii, B, h, ϵ, F, currentNeighbourShell, traversedCells, traversedEdges, cellEdgeOrders, cellVertexOrders)
            push!(traversedCells, ii)
        end
    end

    for j=2:J
        h[j] = h[j] - h[1]
    end
    h[1] = @SVector zeros(2)

    return h
end

export hNetwork

end
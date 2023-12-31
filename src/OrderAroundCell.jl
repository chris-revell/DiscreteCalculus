#
#  OrderAroundCell.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#
# Function to return a *circular* array of vertex labels in the correct clockwise order around a cell, 
# and a *circular* array of edge labels in the correct clockwise order around a cell,
# noting that the first vertex in orderedVertices is on the clockwise end of the first edge in orderedEdges
# ie orderedEdges[1] is the most anticlockwise-ward of the full set of edges and vertices

module OrderAroundCell

using LinearAlgebra
using SparseArrays
using CircularArrays

function orderAroundCell(A, B, indexCell)

    # Find all edges and vertices for cell indexCell
    cellEdges = findall(j -> j != 0, B[indexCell, :])

    orderedVertices = Int64[] # Ordered list of vertices around cell indexCell in clockwise direction 
    orderedEdges = Int64[]    # Ordered list of edges around cell indexCell in clockwise direction

    # Pick a random vertex to start ordering from 
    nextVertex = [findall(j -> j != 0, A[cellEdges[1], :])[1]]

    for _ = 1:length(cellEdges)
        allNeighbouringEdges = findall(j -> j != 0, A[:, nextVertex[end]]) # Find all neighbouring edges around vertex (could be up to 3)                
        iNeighbourEdges = allNeighbouringEdges ∩ findall(k -> k != 0, B[indexCell, :])    # Find the intersection of all neighbours with edges of cell indexCell to give 2 relevant edges
        # Testing both edges in iNeighbourEdges
        # Downstream clockwise if B[indexCell,j]>0 and A[j,k]<0 or B[indexCell,j]<0 and A[j,k]>0
        if B[indexCell, iNeighbourEdges[1]] > 0 && A[iNeighbourEdges[1], nextVertex[end]] < 0
            # iNeighbourEdges[1] is downstream clockwise
            push!(orderedEdges, iNeighbourEdges[1])
            push!(orderedVertices, findall(k -> k > 0, A[iNeighbourEdges[1], :])[1]) # Since A[iNeighbourEdges[1],nextVertex[end]]<0, the next vertex downstream must have A[iNeighbourEdges[1],k]>0
        elseif B[indexCell, iNeighbourEdges[1]] < 0 && A[iNeighbourEdges[1], nextVertex[end]] > 0
            # iNeighbourEdges[1] is downstream clockwise
            push!(orderedEdges, iNeighbourEdges[1])
            push!(orderedVertices, findall(k -> k < 0, A[iNeighbourEdges[1], :])[1]) # Since A[iNeighbourEdges[1],nextVertex[end]]>0, the next vertex downstream must have A[iNeighbourEdges[1],k]<0
        else
            # iNeighbourEdges[2] is downstream clockwise
            push!(orderedEdges, iNeighbourEdges[2])
            # Find the other vertex surrounding iNeighbourEdges[2] that isn't orderedVertices[end]
            if A[iNeighbourEdges[2], nextVertex[end]] > 0
                push!(orderedVertices, findall(k -> k < 0, A[iNeighbourEdges[2], :])[1])
            else
                push!(orderedVertices, findall(k -> k > 0, A[iNeighbourEdges[2], :])[1])
            end
        end
        nextVertex[1] = orderedVertices[end]
    end

    # Convert to circular arrays 
    orderedVertices = CircularArray(orderedVertices)
    orderedEdges = CircularArray(orderedEdges)

    return orderedVertices, orderedEdges

end

export orderAroundCell

end
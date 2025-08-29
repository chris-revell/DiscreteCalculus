#
#  TopologyFunctions.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#
# A set of functions to derive objects that depend only on system topology

module TopologyFunctions

# Julia packages
using LinearAlgebra
using SparseArrays
using CircularArrays
using FromFile 
using DrWatson

@from "OrderAroundCell.jl" using OrderAroundCell

# Non-mutating functions 
findĀ(A) = abs.(A)
findB̄(B) = abs.(B)
findC(A, B) = dropzeros(abs.(B) * abs.(A) .÷ 2)
findAᵀ(A) = Transpose(A)
findĀᵀ(A) = abs.(Transpose(A))
findBᵀ(B) = Transpose(B)
findB̄ᵀ(B) = abs.(Transpose(B))
findCellEdgeCount(B) = sum.(eachrow(abs.(B))) # Zᵢ
findBoundaryVertices(A,B) = abs.(Transpose(A)) * abs.(sum.(eachcol(B))) .÷ 2
findBoundaryEdges(B) = abs.([sum(x) for x in eachcol(B)]) # qⱼᵇ = Bᵀ𝟙ᵢ
function findBoundaryCells(B)
    boundaryEdges = findBoundaryEdges(B)
    boundaryCellIndices = findnz(B[:, boundaryEdges.==1])[1]
    boundaryCells = zeros(Int64, size(B,1))
    boundaryCells[boundaryCellIndices] .= 1
    return boundaryCells
end
function findPerpendicularEdges(A, B)
    kₚ = findBoundaryVertices(A, B).==1
    jₚ = findBoundaryEdges(B).==1
    tmp = findall(x->x!=0, A[:,kₚ])
    tmp2 = unique(getindex.(tmp,1))
    jₙinds = setdiff(tmp2, findall(x->x, jₚ))
    jₙ = zeros(size(A,1))
    jₙ[jₙinds] .= 1
end

# Mutating versions 
findĀ!(A, Ā) = Ā.=abs.(A)
findB̄!(B, B̄) = B̄.=abs.(B)
findC!(Ā, B̄, C) = C.=dropzeros(B̄*Ā.÷2)
findAᵀ!(A, Aᵀ) = Aᵀ.=Transpose(A)
findĀᵀ!(A, Āᵀ) = Āᵀ.=abs.(Transpose(A))
findBᵀ!(B, Bᵀ) = Bᵀ.=Transpose(B)
findB̄ᵀ!(B, B̄ᵀ) = B̄ᵀ.=abs.(Transpose(B))
findCellEdgeCount!(B, Zᵢ) = Zᵢ.=sum.(eachrow(abs.(B))) # Zᵢ
findBoundaryVertices!(A, B, bₖ) = bₖ.=abs.(Transpose(A)) * abs.(sum.(eachcol(B))) .÷ 2
findBoundaryEdges!(B, bⱼ) = bⱼ.=abs.([sum(x) for x in eachcol(B)])

function findBoundaryCells!(B, bᵢ)
    boundaryEdges = findBoundaryEdges(B)
    boundaryCellIndices = findnz(B[:, boundaryEdges.==1])[1]
    bᵢ .= zeros(Int64, size(B,1))
    bᵢ[boundaryCellIndices] .= 1
    return nothing
end

# Function to check for nonzero values in B*A
function senseCheck(A, B; marker="")
    test = B*A
    dropzeros!(test)
    if length(findnz(test)[1]) > 0
        throw("Non-zero values in BA: $(marker)")
    else
        return nothing
    end
end

export topologyMatrices
export findĀ
export findB̄
export findC
export findAᵀ
export findĀᵀ
export findBᵀ
export findB̄ᵀ
export findCellEdgeCount
export findBoundaryVertices
export findBoundaryEdges
export findBoundaryCells
export findPerpendicularEdges

export findĀ!
export findB̄!
export findC!
export findAᵀ!
export findĀᵀ!
export findBᵀ!
export findB̄ᵀ!
export findCellEdgeCount!
export findBoundaryVertices!
export findBoundaryEdges!
export findBoundaryCells!

export senseCheck

end


# function topologyMatrices(A, B)

#     # Find adjacency matrices from incidence matrices
#     Ā = abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)Ā :==>= a:C=>bs.:(A=>(A) :  =>   #: A=> All: =>: 1 :components co=>components conv:erted to +1 (In =>erted to +1 (In ot:her words, create adjacency matrix Ā from incidence matrix A=>her words, create adjacency matrix Ā from incidence matrix A)
#     B̄ = abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)
#     # C adjacency matrix. Rows => cells; Columns => vertices.
#     C = B̄ * Ā .÷ 2 # (NB Integer division)
#     dropzeros!(C)
#     # Update transpose matrices
#     Aᵀ = Transpose(A)
#     Āᵀ = abs.(Aᵀ)
#     Bᵀ = Transpose(B)
#     B̄ᵀ = abs.(Bᵀ)

  
#     # Calculate additional topology data
#     # Number of edges around each cell found by summing columns of B̄
#     cellEdgeCount = sum.(eachrow(B̄))

#     # Find boundary vertices
#     # Summing each column of B finds boundary edges (for all other edges, cell orientations on either side cancel);
#     # multiplying by Aᵀ gives nonzero values only where a vertex (row) has nonzero values at columns (edges) corresponding to nonzero values in the list of boundary edges.
#     # Note that the abs is needed in case the direction of boundary edges cancel at a vertex
#     boundaryVertices = Āᵀ * abs.(sum.(eachcol(B))) .÷ 2

#     # Find list of edges at system periphery
#     boundaryEdges = abs.([sum(x) for x in eachcol(B)])

#     cellVertexOrders  = fill(CircularVector(Int64[]), size(B, 1))
#     cellEdgeOrders    = fill(CircularVector(Int64[]), size(B, 1))
#     for i = 1:size(B,1)
#         cellVertexOrders[i], cellEdgeOrders[i] = orderAroundCell(matrices, i)
#     end

#     return Dict(:Ā=>Ā, :B̄=>B̄, :C=>C, :Aᵀ=>Aᵀ, :Āᵀ=>Āᵀ, :Bᵀ=>Bᵀ, :B̄ᵀ=>B̄ᵀ, :cellEdgeCount=>cellEdgeCount, :boundaryVertices=>boundaryVertices, :boundaryEdges=>boundaryEdges, :cellVertexOrders=>cellVertexOrders, :cellEdgeOrders=>cellEdgeOrders)
# end
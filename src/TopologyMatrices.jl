#
#  TopologyMatrices.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#

module TopologyMatrices

# Julia packages
using LinearAlgebra
using SparseArrays
using FromFile
using DrWatson

# Local modules
@from "SenseCheck.jl" using SenseCheck

function topologyMatrices(A, B)

    # Find adjacency matrices from incidence matrices
    Ā = abs.(A)    # All -1 components converted to +1 (In other words, create adjacency matrix Ā from incidence matrix A)
    B̄ = abs.(B)    # All -1 components converted to +1 (In other words, create adjacency matrix B̄ from incidence matrix B)

    # C adjacency matrix. Rows => cells; Columns => vertices.
    C = B̄ * Ā .÷ 2 # (NB Integer division)

    # Update transpose matrices
    Aᵀ = sparse(transpose(A))
    Āᵀ = abs.(Aᵀ)
    Bᵀ = sparse(transpose(B))
    B̄ᵀ = abs.(Bᵀ)

    dropzeros!(A)
    dropzeros!(B)
    dropzeros!(C)
    dropzeros!(Ā)
    dropzeros!(B̄)
    dropzeros!(C)
    dropzeros!(Aᵀ)
    dropzeros!(Āᵀ)
    dropzeros!(Bᵀ)
    dropzeros!(B̄ᵀ)

    # Calculate additional topology data
    # Number of edges around each cell found by summing columns of B̄
    cellEdgeCount = sum.(eachrow(B̄))

    # Find boundary vertices
    # Summing each column of B finds boundary edges (for all other edges, cell orientations on either side cancel);
    # multiplying by Aᵀ gives nonzero values only where a vertex (row) has nonzero values at columns (edges) corresponding to nonzero values in the list of boundary edges.
    # Note that the abs is needed in case the direction of boundary edges cancel at a vertex
    boundaryVertices = Āᵀ * abs.(sum.(eachcol(B))) .÷ 2

    # Find list of edges at system periphery
    boundaryEdges = abs.([sum(x) for x in eachcol(B)])

    return Ā, B̄, C, Aᵀ, Āᵀ, Bᵀ, B̄ᵀ, cellEdgeCount, boundaryVertices, boundaryEdges

end

export topologyMatrices

end

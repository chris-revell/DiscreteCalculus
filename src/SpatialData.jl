#
#  SpatialData.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
#
#
# Function to calculate spatial data including tangents, lengths, midpoints, tensions etc from incidence and vertex position matrices.

module SpatialData

# Julia packages
using LinearAlgebra
using StaticArrays
using SparseArrays

function findCellCentresOfMass(R, A, B) 
    C = abs.(B) * abs.(A) .÷ 2
    cellEdgeCount = sum.(eachrow(abs.(B)))
    return C*R./cellEdgeCount
end 

findEdgeTangents(R, A) = A*R

findEdgeLengths(edgeTangents) = norm.(edgeTangents)
findEdgeLengths(R, A) = norm.(A*R)

findEdgeMidpoints(R, A) = 0.5.*abs.(A)*R # 0.5.*Ā*R
    
function findCellPerimeterLengths(R, A, B) 
    edgeLengths = norm.(A*R)
    return abs.(B)*edgeLengths # B̄*edgeLengths
end 

function findCellAreas(R, A, B)
    nCells = size(B,1)
    Bᵀ = sparse(transpose(B))
    edgeTangents = A*R
    edgeMidpoints = 0.5.*abs.(A)*R
    # Calculate oriented cell areas    
    cellOrientedAreas = fill(SMatrix{2,2}(zeros(2,2)),nCells)
    cellAreas = zeros(nCells)
    for i=1:nCells
        for j in nzrange(Bᵀ,i)
            cellOrientedAreas[i] += B[i,rowvals(Bᵀ)[j]].*edgeTangents[rowvals(Bᵀ)[j]]*edgeMidpoints[rowvals(Bᵀ)[j]]'            
        end
        cellAreas[i] = cellOrientedAreas[i][1,2]
    end
    return cellAreas 
end

export findCellCentresOfMass
export findEdgeTangents
export findEdgeLengths
export findEdgeLengths
export findEdgeMidpoints
export findCellPerimeterLengths
export findCellAreas

end

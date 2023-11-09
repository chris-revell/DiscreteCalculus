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
using FromFile

# Local modules
@from "OrderAroundCell.jl" using OrderAroundCell

function findCellCentresOfMass(R, A, B) 
    C = abs.(B) * abs.(A) .÷ 2
    cellEdgeCount = sum.(eachrow(abs.(B)))
    return C*R./cellEdgeCount
end 

findEdgeTangents(R, A) = A*R

# findEdgeLengths(edgeTangents) = norm.(edgeTangents)
findEdgeLengths(R, A) = norm.(A*R)

findEdgeMidpoints(R, A) = 0.5.*abs.(A)*R # 0.5.*Ā*R
    
function findCellPerimeterLengths(R, A, B) 
    edgeLengths = norm.(A*R)
    return abs.(B)*edgeLengths # B̄*edgeLengths
end 

function findCellAreas(R, A, B)
    nCells = size(B,1)
    Bᵀ = sparse(Transpose(B)) # Have to convert transpose type to sparse type here because nzrange won't operate on a transpose type object
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
        k_is = findall(x->x!=0, C[:,k])
        if length(k_is) == 2
            edgesSharedBy_i1_And_k = findall(x->x!=0, B[k_is[1],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] = 0.5^3*norm([edgeTangents[edgesSharedBy_i1_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i1_And_k[2]]...,0.0])
            edgesSharedBy_i2_And_k = findall(x->x!=0, B[k_is[2],:])∩findall(x->x!=0, A[:,k])
            vertexAreas[k] += 0.5^3*norm([edgeTangents[edgesSharedBy_i2_And_k[1]]...,0.0]×[edgeTangents[edgesSharedBy_i2_And_k[2]]...,0.0])
        else
            vertexAreas[k] = 0.5*norm([edgeMidpointLinks[k_is[1], k]...,0.0]×[edgeMidpointLinks[k_is[2],k]...,0.0])
        end
    end
    return vertexAreas
end

function makeCellPolygons(R, A, B)
    nCells = size(B, 1)
    cellPolygons = Vector{Point2f}[]
    for indexCell = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, indexCell)
        push!(cellPolygons, Point2f.(R[orderedVertices]))
    end
    return cellPolygons
end

function makeCellLinks(R, A, B)
    nCells = size(B, 1)
    nEdges = size(B, 2)
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
    nCells = size(B, 1)
    nVerts = size(A, 2)
    Aᵀ = Transpose(A)
    Āᵀ = Transpose(abs.(A))
    C = abs.(B) * abs.(A) .÷ 2
    boundaryVertices = Āᵀ * abs.(sum.(eachcol(B))) .÷ 2
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    edgeMidpoints = findEdgeMidpoints(R, A)
    onesVec = ones(Int64, 1, nCells)
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
    nEdges = size(B, 2)
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
    nCells = size(B, 1)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeMidpointPolygons = Vector{Point2f}[]
    for i = 1:nCells
        orderedVertices, orderedEdges = orderAroundCell(A, B, indexCell)
        push!(edgeMidpointPolygons, Point2f.(edgeMidpoints[orderedEdges]))
    end
    return edgeMidpointPolygons
end

function makeCellVerticesDict(A, B)
    nCells = size(B, 1)
    cellVerticesDict = Dict()
    for i = 1:nCells
        cellVertices, cellEdges = orderAroundCell(A, B, indexCell)
        # Store sorted cell vertices for this cell
        cellVerticesDict[i] = cellVertices
    end
    return cellVerticesDict
end

function findEdgeLinkMidpoints(R, A, B, ϵᵢ)
    nEdges = size(B, 2)
    cellCentresOfMass = findCellCentresOfMass(R, A, B)
    edgeTangents = findEdgeTangents(R, A)
    edgeMidpoints = findEdgeMidpoints(R, A)
    edgeTrapezia = makeEdgeTrapezia(R, A, B)
    trapeziumAreas = abs.(area.(edgeTrapezia))
    T = makeCellLinks(R, A, B)
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

# sᵢₖ = ∑ⱼ(1/2)BᵢⱼtⱼĀⱼₖ
function findEdgeMidpointLinks(R, A, B)
    nCells = size(B,1)
    nEdges = size(B,2)
    nVerts = size(A,2)
    edgeTangents = A*R
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1],nzC[2])
    edgeMidpointLinks = fill(SVector{2,Float64}(zeros(2)), (nCells,nVerts))
    nzC = findnz(C)
    ikPairs = tuple.(nzC[1],nzC[2])
    for (i,k) in ikPairs
        for j=1:nEdges
            edgeMidpointLinks[i,k] = edgeMidpointLinks[i,k] .+ 0.5.*B[i,j]*edgeTangents[j]*Ā[j,k]
        end
    end
    return edgeMidpointLinks
end

export findCellCentresOfMass
export findEdgeTangents
export findEdgeLengths
export findEdgeLengths
export findEdgeMidpoints
export findCellPerimeterLengths
export findCellAreas
export findVertexAreas
export makeCellPolygons
export makeCellLinks
export makeLinkTriangles
export makeEdgeTrapezia
export makeEdgeMidpointPolygons
export makeCellVerticesDict
export findEdgeLinkMidpoints
export makeSpokes
export findEdgeMidpointLinks

end #end module 
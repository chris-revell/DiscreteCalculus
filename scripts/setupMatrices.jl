using DrWatson
using DiscreteCalculus
using DelaunayTriangulation
using CairoMakie
using GeometryBasics
using SparseArrays
using GR: delaunay
using StaticArrays
using UnPack
using JLD2
using DifferentialEquations


R = [SVector(x, 0.0) for x=1:9]
for j=1:4
    for i=1:9-j
        push!(R,SVector(i+0.5*j, j*sqrt(1-0.5^2)))
        push!(R,SVector(i+0.5*j, -j*sqrt(1-0.5^2)))
    end
end
xs = [x[1] for x in R]
ys = [x[2] for x in R]

# Delaunay triangulation of centroid locations using function from GR
n, tri = delaunay(xs, ys)

pairs1 = [r[1:2] for r in eachrow(tri)]
pairs2 = [r[2:3] for r in eachrow(tri)]
pairs3 = [r[1:2:3] for r in eachrow(tri)]

pairs = unique([sort.(pairs1)..., sort.(pairs2)..., sort.(pairs3)...])

nCells = n
nEdges = length(pairs)
nVerts = length(R)

A = spzeros(Int64, nEdges, nVerts)
for (i,p) in enumerate(pairs)
    A[i,p[1]] = 1
    A[i,p[2]] = -1
end

B = spzeros(Int64, nCells, nEdges)
for i=1:nCells
    for s=1:3
        p = circshift(tri[i,:],s)
        # tri gives vertices ordered clockwise around triangle
        jEdge = findall(x->x[1]!=0&&x[2]!=0, eachrow(A[:,p[1:2]]))[1]
        if A[jEdge, p[1]] == 1
            B[i, jEdge] = 1
        else
            B[i, jEdge] = -1
        end
    end
end

edgeFlip = 120

A[edgeFlip,6] = A[edgeFlip,5]
A[edgeFlip,5] = 0
A[edgeFlip,17] = A[edgeFlip,19]
A[edgeFlip,19] = 0

B[96,112] = B[76,112]
B[76,112] = 0
B[76,81] = B[96,81]
B[96,81] = 0

dropzeros!(A)
dropzeros!(B)

function model!(du, u, p, t)
    A, B, ϵ = p
    Ā = abs.(A)
    B̄ = abs.(B)
    edgeTangents = findEdgeTangents(u, A)
    edgeLengths = findEdgeLengths(u, A)
    cellPerimeterLengths = findCellPerimeterLengths(u, A, B) 
    cellAreas = findCellAreas(u, A, B)
    cellTensions   = (3.0 .- cellPerimeterLengths)
    cellPressures  = cellAreas .- sqrt(0.75)
    for k=1:nVerts
        for j in nzrange(A,k)
            for i in nzrange(B,rowvals(A)[j])
                # Force components from cell pressure perpendicular to edge tangents 
                du[k] += 0.5*cellPressures[rowvals(B)[i]]*B[rowvals(B)[i],rowvals(A)[j]]*Ā[rowvals(A)[j],k].*(ϵ*edgeTangents[rowvals(A)[j]])
                # Force components from cell membrane tension parallel to edge tangents 
                du[k] += cellTensions[rowvals(B)[i]]*B̄[rowvals(B)[i],rowvals(A)[j]]*A[rowvals(A)[j],k].*edgeTangents[rowvals(A)[j]]./edgeLengths[rowvals(A)[j]]                
            end
        end
    end
end

ϵ = @SMatrix [  # Clockwise rotation matrix setting orientation of cell faces
0.0 1.0
-1.0 0.0
]
prob = ODEProblem(model!,R,(0.0,0.1),(A,B,ϵ))
sol = solve(prob,Tsit5())


cellPolygons = makeCellPolygons(sol.u[end],A,B)

# fig = Figure(resolution=(1000,1000))
# ax = Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(ax); hidespines!(ax)
# for i=1:nCells
#     poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=2)
# end
# display(fig)

R = sol.u[end]

matrices = @dict A B R
save("defectMatrices.jld2",matrices)




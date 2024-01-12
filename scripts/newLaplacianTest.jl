using DrWatson
using DiscreteCalculus
using DelaunayTriangulation
using CairoMakie
using GeometryBasics
using SparseArrays
using GR: delaunay
using StaticArrays
using UnPack
using LinearAlgebra
using Printf
using StaticVectors

@unpack matrices, R = load(projectdir("scripts","referenceMatrices3.jld2"))
@unpack A, B = matrices

nCells = size(B,1)
nEdges = size(B,2)
nVerts = size(A,2)

ϵᵢ = @SMatrix [
0.0 1.0
-1.0 0.0
]

sᵢₖ = findEdgeMidpointLinks(R, A, B)

cellPolygons = makeCellPolygons(R,A,B)
edgeMidpoints = findEdgeMidpoints(R, A)

# fig = Figure(resolution=(1000,1000))
# ax = Axis(fig[1,1],aspect=DataAspect())
# hidedecorations!(ax); hidespines!(ax)
# for i=1:nCells
#     poly!(ax,cellPolygons[i],color=(:white,0.0),strokecolor=(:black,1.0),strokewidth=2)
# end
# for i=1:nCells
#     orderedVertices, orderedEdges = orderAroundCell(A, B, i)
#     for kk=1:length(orderedVertices)
#         lines!(ax, [Point2(edgeMidpoints[orderedEdges[kk]]...), Point2( (edgeMidpoints[orderedEdges[kk]].+sᵢₖ[(i,orderedVertices[kk])])... )], linestyle=:dot, color=:black)
#     end
# end


ϕᵢ = zeros(Float64, nCells)
cellPositions = findCellCentresOfMass(R, A, B) 
ϕᵢ = norm.(cellPositions)

function gradₖ(ϕᵢ, R, A, B, ϵᵢ)
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2 # (NB Integer division)
    dropzeros!(C)
    sᵢₖ = findEdgeMidpointLinks(R, A, B)
    αₖ = findVertexAreas(R, A, B)
    edgeTangents = findEdgeTangents(R, A)
    gradϕ = fill(SVector{2,Float64}(zeros(2)),nVerts)
    for k=1:nVerts
        for i=1:nCells
            for j=1:nEdges
                gradϕ[k] += SVector{2,Float64}(0.5*ϕᵢ[i].*B[i,j].*(ϵᵢ*edgeTangents[j]).*Ā[j,k])./αₖ[k]
            end
        end
    end
    return gradϕ
end

gradϕ = gradₖ(ϕᵢ, R, A, B, ϵᵢ)

matrices = topologyMatrices(A,B)
@unpack Āᵀ, boundaryVertices = matrices

for k in findall(x->x!=0, boundaryVertices)
    gradϕ[k] = SVector{2,Float64}(zeros(2))
end

fig2 = Figure(resolution=(1000,1000))
ax2 = Axis(fig2[1,1], aspect=DataAspect())
hidedecorations!(ax2); hidespines!(ax2)
clims = (minimum(ϕᵢ),maximum(ϕᵢ))
for i=1:nCells
    poly!(ax2,cellPolygons[i],color=ϕᵢ[i],colorrange=clims,strokewidth=2,strokecolor=(:black,0.25))
end
arrows!(ax2, [Point2(b...) for b in R], [Vec2(0.1.*a...) for a in gradϕ])
xlims!(minimum(first.(R)),maximum(first.(R)))
ylims!(minimum(last.(R)),maximum(last.(R)))
display(fig2)
save(projectdir("scripts","plots","gradϕ_Neumann3.png"),fig2)


function divᵢ(v, R, A, B, ϵᵢ)
    Aᵢ = findCellAreas(R, A, B)
    Ā = abs.(A)
    B̄ = abs.(B)
    C = B̄ * Ā .÷ 2 # (NB Integer division)
    dropzeros!(C)
    divv = zeros(nCells)
    for i=1:nCells
        for k in findall(x->x!=0, C[i,:])
            divv[i] += (ϵᵢ*sᵢₖ[i,k])⋅v[k]./Aᵢ[i]
        end
    end
    return divv
end
# Confirm div(R).==2.0
divv = divᵢ(R, R, A, B, ϵᵢ)


Lₐ = edgeMidpointLDirichlet(R, A, B)

test1 = L*ϕᵢ

# fig3 = Figure(resolution=(1000,1000))
# ax1 = Axis(fig3[1,1],aspect=DataAspect())
# hidedecorations!(ax1); hidespines!(ax1)
# lims1 = (minimum(test1),maximum(test1))
# for i=1:nCells
#     poly!(ax1,cellPolygons[i],color=test1[i],colorrange=lims1,colormap=:inferno,strokecolor=(:black,1.0),strokewidth=1)
# end
# save("Lϕ.png",fig3)
# display(fig3)


eigenvalues, eigenmodes = eigen(Matrix(Lₐ))

# fig4 = Figure(resolution=(1000,1000))
# ax4 = Axis(fig4[1,1], aspect=DataAspect())
# hidedecorations!(ax4); hidespines!(ax4)
# for mode=1:size(eigenmodes,2)
#     empty!(ax4)
#     lims = (-maximum(abs.(eigenmodes[:,mode])),maximum(abs.(eigenmodes[:,mode])))
#     for i=1:nCells
#         poly!(ax4,cellPolygons[i],color=eigenmodes[i,mode],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.25))
#     end
#     save("eigenmodes/mode$(@sprintf("%03d", mode)).png",fig4)
#     # display(fig4)
# end


Aᵢ = findCellAreas(R, A, B)

for mode1=1:size(eigenmodes,2)-1
    for mode2=mode1+1:size(eigenmodes,2)
        @show sum(Aᵢ.*eigenmodes[:,mode1].*eigenmodes[:,mode2])
    end
end


# fig5 = Figure()
# ax5 = Axis(fig5[1,1])
# scatter!(ax5, collect(1:length(eigenvalues)), eigenvalues)
# ax5.xlabel ="Eigenmode number"
# ax5.ylabel ="Eigenvalue"
# save("eigenvalueSpectrum.png",fig5)
# display(fig5)


# fig6 = 
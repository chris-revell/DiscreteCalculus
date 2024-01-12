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

@unpack matrices, R = load(projectdir("scripts","referenceMatrices2.jld2"))
@unpack A, B = matrices

nCells = size(B,1)
nEdges = size(B,2)
nVerts = size(A,2)

L = uniformCellTensionL(R, A, B)

eigenvalues, eigenmodes = eigen(Matrix(L))


# Confirm eigenmodes are orthogonal
eigenModeDotProducts = Float64[]
for j=1:nCells-1
    for j′=j+1:nCells
        push!(eigenModeDotProducts,eigenmodes[:,j]⋅eigenmodes[:,j′])
    end
end
@show maximum(eigenModeDotProducts)

cellPolygons = makeCellPolygons(R,A,B)
linkTriangles = makeLinkTriangles(R, A, B)
edgeTrapezia = makeEdgeTrapezia(R, A, B)

set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(2000,2000))

modesToUse = [collect(1:10)...,collect(nCells-9:nCells)...]
axDict = Dict()

for eigenvectorIndex in modesToUse
    absLim = maximum(abs.(eigenmodes[:,eigenvectorIndex]))
    lims = (-absLim,absLim)
    ax = Axis(fig,aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    for i=1:nCells
        poly!(ax,cellPolygons[i],color=eigenmodes[i,eigenvectorIndex],colorrange=lims,colormap=:bwr,strokecolor=(:black,1.0),strokewidth=1)
    end
    axDict[eigenvectorIndex] = ax
end

for y=1:2
    for x=1:5
        fig[y,x] = axDict[(y-1)*5+x]
        Label(fig[y,x,Bottom()],
            "$((y-1)*5+x)",
            fontsize = 32,
        )
    end
end
for y=1:2
    for x=1:5
        fig[y+2,x] = axDict[(y-1)*5+x+nCells-10]
        Label(fig[y+2,x,Bottom()],
            "$((y-1)*5+x+90)",
            fontsize = 32,
        )
    end
end

resize_to_layout!(fig)
# display(fig)
save(projectdir("scripts","plots","uniformCellTensionLEigenmodeTableau2.png"),fig)

fig2 = Figure()
ax2 = Axis(fig2[1,1])
scatter!(ax2, collect(1:length(eigenvalues)), eigenvalues)
ax2.xlabel ="Eigenmode number"
ax2.ylabel ="Eigenvalue"
save(projectdir("scripts","plots","uniformCellTensionLEigenmodeSpectrum2.png"),fig2)
# display(fig2)

# fig4 = Figure(size=(1000,1000))
# ax4 = Axis(fig4[1,1], aspect=DataAspect())
# hidedecorations!(ax4); hidespines!(ax4)
# isdir("test/plots/scalarEdgeNonlocalTensionLeigenmodes") ? nothing : mkdir("test/plots/scalarEdgeNonlocalTensionLeigenmodes")
# for mode=1:size(eigenmodes,2)
#     empty!(ax4)
#     lims = (-maximum(abs.(eigenmodes[:,mode])),maximum(abs.(eigenmodes[:,mode])))
#     for i=1:nCells
#         poly!(ax4,cellPolygons[i],color=[eigenmodes[i,mode]],colorrange=lims,colormap=:bwr,strokewidth=2,strokecolor=(:black,0.25))
#     end
#     save("test/plots/scalarEdgeNonlocalTensionLeigenmodes/scalarEdgeNonlocalTensionL_mode$(@sprintf("%03d", mode)).png",fig4)
#     # display(fig4)
# end
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

cellPolygons = makeCellPolygons(R,A,B)

L = edgeMidpointLNeumann(R, A, B)

eigenvalues, eigenmodes = eigen(Matrix(L))

set_theme!(figure_padding=1, backgroundcolor=(:white,1.0), font="Helvetica")
fig = Figure(size=(2000,2000))

modesToUse = [collect(1:10)...,collect(91:100)...]
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
        fig[y+2,x] = axDict[(y-1)*5+x+90]
        Label(fig[y+2,x,Bottom()],
            "$((y-1)*5+x+90)",
            fontsize = 32,
        )
    end
end

resize_to_layout!(fig)
# display(fig)
save(projectdir("scripts","plots","edgeMidpointLaplacianNeumannEigenmodeTableau3.png"),fig)

fig2 = Figure()
ax2 = Axis(fig2[1,1])
scatter!(ax2, collect(1:length(eigenvalues)), eigenvalues)
ax2.xlabel ="Eigenmode number"
ax2.ylabel ="Eigenvalue"
save(projectdir("scripts","plots","edgeMidpointLaplacianNeumannEigenvalueSpectrum3.png"),fig2)
# display(fig2)
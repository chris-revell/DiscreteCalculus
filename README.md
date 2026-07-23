# DiscreteCalculus

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> DiscreteCalculus

It is authored by chris-revell <christopher.revell@manchester.ac.uk>.

## Background

This package provides a user-friendly interface for operations described in [1]. Within `src/` are several modules that define functions that are then exported from `DiscreteCalculus.jl`. The functions are intended to operate on a discrete network of the form described in [1-5], defined by 3 matrices: a matrix `R` of vertex positions, signed incidence matrix `A` mapping edges to vertices, and signed incidence matrix `B` mapping cells to edges. This paradigm is used for exported data from the VertexModel.jl package [6]. 

Almost all of the functions exported by this package take all or some of `R`, `A`, and `B` as arguments. Many of the functions have both a mutating version, denoted by `!` added to the end of the function name, and non-mutating version with the same name but without the `!`. Mutating functions take an additional argument in which to store the result of the operation. 

Functions are divided into the following basic categories.

### Incidence matrix operations

These functions, listed below, act on the incidence matrices to obtain unsigned versions, denoted with an overbar, transposed versions, denoted with $^T$, or truncated versions, denoted by a hat.

```
findĀ
findB̄
findC
findAᵀ
findĀᵀ
findBᵀ
findB̄ᵀ
findÂ
findB̂
findĀ!
findB̄!
findC!
findAᵀ!
findĀᵀ!
findBᵀ!
findB̄ᵀ!
findÂ!
findB̂!
```

### Topology functions

This set of functions return properties of the network that depend only on the topology, given by the `A` and `B`.

```
findCellEdgeCount
findPeripheralVertices
findPeripheralEdges
findPeripheralCells
findNormalEdges
findCellEdgeCount!
findPeripheralVertices!
findPeripheralEdges!
findPeripheralCells!
findNormalEdges!
```

### Force potential network

```
hNetwork
```

### Geometry functions

This set of functions return properties of the network that depend on the positions of vertices, `R`.

```
findCellCentresOfMass
findEdgeTangents
findEdgeLengths
findEdgeMidpoints
findCellPerimeterLengths
findCellAreas
findEdgeMidpointLinkVertexAreas
findCellPolygons
findCellLinks
findCellLinkMidpoints
findCellLinkLengths
findCellLinkVertexTriangles
findCellLinkVertexTriangleAreas
findEdgeQuadrilaterals
findEdgeQuadrilateralAreas
findEdgeMidpointCellPolygons
findEdgeLinkIntersections
findSpokes
findEdgeMidpointLinks
findCellOutwardNormals
findCellLinkTriangleOutwardNormals
findEdgeMidpointLinkTriangleOutwardNormals
```

### Differential operators

```
gradᵛ
cogradᵛ
corotᶜ
rotᶜ
gradᶜ
cogradᶜ
corotᵛ
corotᵛspokes
rotᵛ
rotᵛspokes
divᵛ
divᵛsuppress
codivᵛ
codivᵛsuppress
cocurlᶜ
curlᶜ
divᶜ
divᶜsuppress
codivᶜ
codivᶜsuppress
cocurlᵛ
cocurlᵛspokes
curlᵛ
curlᵛspokes
𝐃c
𝐃v
𝐆c
grad_A
div_A
curl_A
rot_A
grad_L
div_L
```

### Laplacians

```
geometricLf
geometricLfHat
geometricLc
geometricLcHat
geometricLv
geometricLvHat
geometricLt
geometricLtHat
topologicalLf
topologicalLc
topologicalLv
topologicalLt
edgeLaplacianPrimal
edgeLaplacianPrimalHat
edgeLaplacianDual
edgeLaplacianDualHat
```

### Other exported functions and constants

```
senseCheck
orderAroundCell
cellNeighbourOrder
innerProd 
outerProd
ϵᵢ
ϵₖ
```


1. Jensen, O. E. & Revell, C. K. Harmonic fields and the mechanical response of a cellular monolayer to ablation. J. Math. Biol. 92, 61 (2026).

2. Jensen, O. E. & Revell, C. K. Couple stresses and discrete potentials in the vertex model of cellular monolayers. Biomech Model Mechanobiol 22, 1465–1486 (2023).

3. Nestor-Bergmann, A., Goddard, G., Woolner, S. & Jensen, O. E. Relating cell shape and mechanical stress in a spatially disordered epithelium using a vertex-based model. Mathematical Medicine and Biology: A Journal of the IMA 35, i1–i27 (2018).

4. Jensen, O. E., Johns, E. & Woolner, S. Force networks, torque balance and Airy stress in the planar vertex model of a confluent epithelium. Proc. R. Soc. A. 476, 20190716 (2020).

5. Cowley, N., Revell, C. K., Johns, E., Woolner, S. & Jensen, O. E. Spectral approaches to stress relaxation in epithelial monolayers. Proceedings of the Royal Society A https://doi.org/10.1098/rspa.2024.0224 (2024) doi:10.1098/rspa.2024.0224.

6. Revell, C. K. & Jensen, O. E. VertexModel.jl. (2025).










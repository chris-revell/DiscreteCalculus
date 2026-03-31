#
#  DiscreteCalculus.jl
#  DiscreteCalculus
#
#

module DiscreteCalculus

# Julia packages
using DrWatson
using LinearAlgebra
using SparseArrays
using FromFile

# Local modules
@from "GeometryFunctions.jl" using GeometryFunctions
@from "Laplacians.jl" using Laplacians
@from "OrderAroundCell.jl" using OrderAroundCell
@from "TopologyFunctions.jl" using TopologyFunctions
@from "HNetwork.jl" using HNetwork
@from "DifferentialOperators.jl" using DifferentialOperators
@from "InnerOuterProduct.jl" using InnerOuterProduct
@from "CovariantDerivatives.jl" using CovariantDerivatives

export findĀ
export findB̄
export findC
export findAᵀ
export findĀᵀ
export findBᵀ
export findB̄ᵀ
export findÂ
export findB̂
export findCellEdgeCount
export findPeripheralVertices
export findPeripheralEdges
export findPeripheralCells
export findNormalEdges

export findĀ!
export findB̄!
export findC!
export findAᵀ!
export findĀᵀ!
export findBᵀ!
export findB̄ᵀ!
export findÂ!
export findB̂!
export findCellEdgeCount!
export findPeripheralVertices!
export findPeripheralEdges!
export findPeripheralCells!
export findNormalEdges!

export senseCheck

export findCellCentresOfMass
export findEdgeTangents
export findEdgeLengths
export findEdgeLengths
export findEdgeMidpoints
export findCellPerimeterLengths
export findCellAreas
export findEdgeMidpointLinkVertexAreas
export findCellPolygons
export findCellLinks
export findCellLinkMidpoints
export findCellLinkLengths
export findCellLinkVertexTriangles
export findCellLinkVertexTriangleAreas
export findEdgeQuadrilaterals
export findEdgeQuadrilateralAreas
export findEdgeMidpointCellPolygons
export findEdgeLinkIntersections
export findSpokes
export findEdgeMidpointLinks
export findCellOutwardNormals
export findCellLinkTriangleOutwardNormals
export findEdgeMidpointLinkTriangleOutwardNormals

export orderAroundCell
export cellNeighbourOrder

export hNetwork

export geometricLf
export geometricLfHat
export geometricLfHatReduced
export geometricLc
export geometricLcHat
export geometricLcHatReduced
export geometricLv
export geometricLvHat
export geometricLvHatReduced
export geometricLt
export geometricLtHat
export geometricLtHatReduced
export topologicalLf
export topologicalLc
export topologicalLv
export topologicalLt
export edgeLaplacianPrimal
export edgeLaplacianPrimalHat
export edgeLaplacianDual
export edgeLaplacianDualHat

export edgeMidpointL
export edgeMidpointLNeumann
export scalarEdgeL
export uniformCellTensionL

export gradᵛ
export cogradᵛ
export corotᶜ
export rotᶜ
export gradᶜ
export cogradᶜ
export corotᵛ
export corotᵛspokes
export rotᵛ
export rotᵛspokes
export divᵛ
export divᵛsuppress
export codivᵛ
export codivᵛsuppress
export cocurlᶜ
export curlᶜ
export divᶜ
export divᶜsuppress
export codivᶜ
export codivᶜsuppress
export cocurlᵛ
export cocurlᵛspokes
export curlᵛ
export curlᵛspokes

export 𝐃c
export 𝐃v
export 𝐆v

export grad_A
export div_A
export curl_A
export rot_A
export grad_L
export div_L

export innerProd 
export outerProd

end 
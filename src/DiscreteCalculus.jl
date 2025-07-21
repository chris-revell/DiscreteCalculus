#
#  DiscreteCalculus.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 15/08/2023.
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
export senseCheck

export findCellCentresOfMass
export findEdgeTangents
export findEdgeLengths
export findEdgeLengths
export findEdgeMidpoints
export findCellPerimeterLengths
export findCellAreas
export findVertexAreas
export findCellPolygons
export findCellLinks
export findCellLinkMidpoints
export findCellLinkLengths
export findCellLinkTriangles
export findCellLinkTriangleAreas
export findEdgeQuadrilaterals
export findEdgeQuadrilateralAreas
export findEdgeMidpointCellPolygons
export findEdgeLinkMidpoints
export findSpokes
export findEdgeMidpointLinks

export orderAroundCell
export cellNeighbourOrder

export hNetwork

export geometricLf
export geometricLc
export geometricLv
export geometricLt
export topologicalLf
export topologicalLc
export topologicalLv
export topologicalLt
export edgeLaplacianPrimal
export edgeLaplacianDual

export edgeMidpointLDirichlet
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

export innerProd 
export outerProd

end 
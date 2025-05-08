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
export findCellCentreLinks
export findCellCentreLinkLengths
export findCellLinkTriangles
export findCellLinkTriangleAreas
export findEdgeQuadrilaterals
export findEdgeQuadrilateralAreas
export findEdgeMidpointPolygons
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
export edgeMidpointLDirichlet
export edgeMidpointLNeumann
export scalarEdgeL
export uniformCellTensionL

export gradᵛ
export curlᶜ
export gradᶜ
export curlᵛ
export rotᵛ
export divᶜ
export rotᶜ
export divᵛ
export cogᵛ
export cocurlᶜ
export cogᶜ
export cocurlᵛ
export corotᵛ
export codᶜ
export corotᶜ
export codᵛ 

end 
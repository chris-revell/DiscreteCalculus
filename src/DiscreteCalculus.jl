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

export topologyMatrices
export getRandomColor
export makeCellPolygons
export makeCellLinks
export makeLinkTriangles
export makeEdgeTrapezia
export makeEdgeMidpointPolygons
export calculateCellCurls
export calculateCellDivs
export calculateVertexDivs
export calculateVertexCurls
export makeCellVerticesDict
export findEdgeLinkMidpoints
export makeSpokes
export calculateVertexMidpointCurls
export calculateVertexMidpointDivs
export calculateCellMidpointDivs
export calculateCellMidpointCurls
export geometricLf
export geometricLc
export geometricLv
export geometricLt
export orderAroundCell
export psicPotential
export psivPotential
export capitalPsivPotential
export senseCheck
export findCellCentresOfMass
export findEdgeTangents
export findEdgeLengths
export findEdgeLengths
export findEdgeMidpoints
export findEdgeMidpointLinks
export findCellPerimeterLengths
export findCellAreas
export findVertexAreas
export topologyMatrices
export topologicalLf
export topologicalLc
export topologicalLv
export topologicalLt
export edgeMidpointLDirichlet
export edgeMidpointLNeumann
export scalarEdgeL
export uniformCellTensionL
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


# export edgeMidpointLfunction

end #end module 
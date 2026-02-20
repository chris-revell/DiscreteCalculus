#
#  TopologyFunctions.jl
#  DiscreteCalculus
#
#
# A set of functions to derive objects that depend only on system topology

module TopologyFunctions

# Julia packages
using LinearAlgebra
using SparseArrays
using CircularArrays
using FromFile 
using DrWatson

@from "OrderAroundCell.jl" using OrderAroundCell

# Non-mutating functions 
findĀ(A) = abs.(A)
findB̄(B) = abs.(B)
findC(A, B) = dropzeros(abs.(B) * abs.(A) .÷ 2)
findAᵀ(A) = Transpose(A)
findĀᵀ(A) = abs.(Transpose(A))
findBᵀ(B) = Transpose(B)
findB̄ᵀ(B) = abs.(Transpose(B))
findCellEdgeCount(B) = sum.(eachrow(abs.(B))) # Zᵢ
findPeripheralVertices(A,B) = abs.(Transpose(A)) * abs.(sum.(eachcol(B))) .÷ 2
findPeripheralEdges(B) = abs.([sum(x) for x in eachcol(B)]) # qⱼᵇ = Bᵀ𝟙ᵢ
function findPeripheralCells(B)
    jᵖ = findPeripheralEdges(B)
    peripheralCellIndices = findnz(B[:, jᵖ.==1])[1]
    iᵖ = zeros(Int64, size(B,1))
    iᵖ[peripheralCellIndices] .= 1
    return iᵖ
end
function findNormalEdges(A, B)
    kᵖ = findPeripheralVertices(A, B).==1
    jᵖ = findPeripheralEdges(B).==1
    tmp = findall(x->x!=0, A[:,kᵖ])
    tmp2 = unique(getindex.(tmp,1))
    jⁿinds = setdiff(tmp2, findall(x->x, jᵖ))
    jⁿ = zeros(Int64, size(A,1))
    jⁿ[jⁿinds] .= 1
    return jⁿ
end
function findÂ(A, B)
    jᵖ = findPeripheralEdges(B)
    kᵖ = findPeripheralVertices(A,B)
    Â = A[jᵖ.==0, kᵖ.==0]
    return Â
end
function findB̂(A, B)
    jᵖ = findPeripheralEdges(B)
    B̂ = B[:, jᵖ.==0]
    return B̂
end
    

# Mutating versions 
findĀ!(A, Ā) = Ā.=abs.(A)
findB̄!(B, B̄) = B̄.=abs.(B)
findC!(Ā, B̄, C) = C.=dropzeros(B̄*Ā.÷2)
findAᵀ!(A, Aᵀ) = Aᵀ.=Transpose(A)
findĀᵀ!(A, Āᵀ) = Āᵀ.=abs.(Transpose(A))
findBᵀ!(B, Bᵀ) = Bᵀ.=Transpose(B)
findB̄ᵀ!(B, B̄ᵀ) = B̄ᵀ.=abs.(Transpose(B))
findCellEdgeCount!(B, Zᵢ) = Zᵢ.=sum.(eachrow(abs.(B))) # Zᵢ
findPeripheralVertices!(A, B, bₖ) = bₖ.=abs.(Transpose(A)) * abs.(sum.(eachcol(B))) .÷ 2
findPeripheralEdges!(B, bⱼ) = bⱼ.=abs.([sum(x) for x in eachcol(B)])
function findPeripheralCells!(B, iᵖ)
    jᵖ = findPeripheralEdges(B)
    peripheralCellIndices = findnz(B[:, jᵖ.==1])[1]
    iᵖ .= zeros(Int64, size(B,1))
    iᵖ[peripheralCellIndices] .= 1
    return nothing
end
function findNormalEdges!(A, B, jⁿ)
    kᵖ = findPeripheralVertices(A, B).==1
    jᵖ = findPeripheralEdges(B).==1
    tmp = findall(x->x!=0, A[:,kᵖ])
    tmp2 = unique(getindex.(tmp,1))
    jⁿinds = setdiff(tmp2, findall(x->x, jᵖ))
    jⁿ .= zeros(Int64, size(A,1))
    jⁿ[jⁿinds] .= 1
    return nothing
end
function findÂ!(A, B, Â)
    jᵖ = findPeripheralEdges(B)
    kᵖ = findPeripheralVertices(A,B)
    Â .= A[jᵖ.==0, kᵖ.==0]
    return Â
end
function findB̂!(A, B, B̂)
    jᵖ = findPeripheralEdges(B)
    B̂ .= B[:, jᵖ.==0]
    return B̂
end

# Function to check for nonzero values in B*A
function senseCheck(A, B; marker="")
    test = B*A
    dropzeros!(test)
    if length(findnz(test)[1]) > 0
        throw("Non-zero values in BA: $(marker)")
    else
        return nothing
    end
end

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

end
#
#  InnerOuterProduct.jl
#  DiscreteCalculus
#
#  Created by Christopher Revell on 21/07/2025.
#
#

module InnerOuterProduct

using LinearAlgebra
using SparseArrays

innerProd(a, M, b) = transpose(a) * M * b # Inner product of a and b under metric M

function outerProd(a, b)
    return a*transpose(b)
end
function outerProd(a, M, b)
    return a*M*transpose(b)
end

export innerProd
export outerProd 

end
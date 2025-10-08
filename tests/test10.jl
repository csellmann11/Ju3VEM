


const VecOrNum{1} = Float64
const VecOrNum{U} = SVector{U,Float64}

struct Test{U} <: AbstractVector{VecOrNum{U}}

end



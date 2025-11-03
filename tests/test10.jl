@inline function binomial_ct(n::Int, k::Int)
    k = min(k, n - k)
    k == 0 ? 1 :
    k == 1 ? n :
    k == 2 ? n * (n - 1) ÷ 2 :
    k == 3 ? n * (n - 1) * (n - 2) ÷ 6 :
    k == 4 ? n * (n - 1) * (n - 2) * (n - 3) ÷ 24 :
    k == 5 ? n * (n - 1) * (n - 2) * (n - 3) * (n - 4) ÷ 120 :
    _binomial_general(n, k)
end

@inline function _binomial_general(n::Int, k::Int)
    result = 1
    for i in 1:k
        result = result * (n - i + 1) ÷ i
    end
    result
end

@inline function polynomial_order(::Val{L}, ::Val{D}) where {L, D}
    p = 0
    while p < 20  # reasonable upper bound
        binomial_ct(p + D, D) == L && return p
        p += 1
    end
    error("Order not found")
end

@code_llvm polynomial_order(Val(6), Val(2))
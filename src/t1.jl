using StaticArrays

"""
    MonomialIndexer{N,MaxDeg,OffsetSize,BinomSize}

Precomputed tables for ultra-fast monomial exponent to position mapping.
N is the number of variables, MaxDeg is the maximum total degree to precompute.
"""
struct MonomialIndexer{N,MaxDeg,OffsetSize,BinomSize}
    degree_offsets::SVector{OffsetSize,Int}  # Cumulative count of monomials up to each degree
    binomial_table::SMatrix{BinomSize,BinomSize,Int}  # Precomputed binomials
end

function MonomialIndexer(::Val{N}, ::Val{MaxDeg}) where {N,MaxDeg}
    OffsetSize = MaxDeg + 1
    BinomSize = MaxDeg + N
    
    # Precompute degree offsets
    offsets = zeros(Int, OffsetSize)
    offsets[1] = 0  # No monomials before degree 0
    for d in 1:MaxDeg
        # Number of monomials of degree d-1 in N variables
        offsets[d+1] = offsets[d] + binomial(N+d-2, d-1)
    end
    
    # Precompute binomial coefficients
    binom_table = zeros(Int, BinomSize, BinomSize)
    for n in 0:(BinomSize-1)
        binom_table[n+1,1] = 1  # C(n,0) = 1
        binom_table[n+1,n+1] = 1  # C(n,n) = 1
        for k in 1:min(n-1, BinomSize-1)
            binom_table[n+1,k+1] = binom_table[n,k] + binom_table[n,k+1]
        end
    end
    
    MonomialIndexer{N,MaxDeg,OffsetSize,BinomSize}(
        SVector{OffsetSize}(offsets),
        SMatrix{BinomSize,BinomSize}(binom_table)
    )
end

"""
    exp_to_position(indexer::MonomialIndexer{N}, exp::SVector{N,Int})

Convert monomial exponent vector to position in graded lexicographic ordering.
Returns 1-based index.

# Example
```julia
indexer = MonomialIndexer(Val(3), Val(10))  # 3 variables, max degree 10
pos = exp_to_position(indexer, SVector(2,1,0))  # x²y
```
"""
@inline function exp_to_position(indexer::MonomialIndexer{N,MaxDeg}, exp::SVector{N,Int}) where {N,MaxDeg}
    # Calculate total degree
    deg = sum(exp)
    
    # Early return for zero monomial
    deg == 0 && return 1
    
    # Get offset for all monomials of smaller degree
    @inbounds offset = deg <= MaxDeg ? indexer.degree_offsets[deg+1] : slow_degree_offset(N, deg)
    
    # Calculate position within degree using lexicographic ordering
    pos_in_degree = lex_position_in_degree(indexer, exp, deg)
    
    return offset + pos_in_degree
end

"""
    lex_position_in_degree(indexer, exp, deg)

Calculate position within monomials of the same degree using lexicographic ordering.
"""
@inline function lex_position_in_degree(indexer::MonomialIndexer{N,MaxDeg,OS,BS}, exp::SVector{N,Int}, deg::Int) where {N,MaxDeg,OS,BS}
    position = 1
    remaining_deg = deg
    
    @inbounds for i in 1:N-1
        # Count monomials that come before us with x_i having smaller exponent
        for e_i in 0:(exp[i]-1)
            new_remaining = remaining_deg - e_i
            if new_remaining >= 0
                # Number of ways to distribute new_remaining among variables i+1 to N
                position += get_binomial(indexer, N-i+new_remaining-1, new_remaining)
            end
        end
        remaining_deg -= exp[i]
    end
    
    return position
end

"""
    get_binomial(indexer, n, k)

Fast binomial coefficient lookup with fallback for large values.
"""
@inline function get_binomial(indexer::MonomialIndexer{N,MaxDeg,OS,BS}, n::Int, k::Int) where {N,MaxDeg,OS,BS}
    if n < BS && k < BS && n >= 0 && k >= 0
        @inbounds return indexer.binomial_table[n+1, k+1]
    else
        return binomial(n, k)  # Fallback for large values
    end
end

"""
    slow_degree_offset(N, deg)

Fallback for degrees beyond precomputed range.
"""
function slow_degree_offset(N::Int, deg::Int)
    offset = 0
    for d in 0:(deg-1)
        offset += binomial(N+d-1, d)
    end
    return offset
end

# Alternative: Even faster version for small, fixed N
"""
    MonomialIndexer2D{MaxDeg,OffsetSize}

Specialized ultra-fast indexer for 2D case.
"""
struct MonomialIndexer2D{MaxDeg,OffsetSize}
    degree_offsets::SVector{OffsetSize,Int}
end

function MonomialIndexer2D(::Val{MaxDeg}) where MaxDeg
    OffsetSize = MaxDeg + 1
    offsets = zeros(Int, OffsetSize)
    offsets[1] = 0
    for d in 1:MaxDeg
        offsets[d+1] = offsets[d] + d  # In 2D, degree d has d+1 monomials
    end
    MonomialIndexer2D{MaxDeg,OffsetSize}(SVector{OffsetSize}(offsets))
end

@inline function exp_to_position(indexer::MonomialIndexer2D{MaxDeg}, exp::SVector{2,Int}) where MaxDeg
    deg = exp[1] + exp[2]
    deg == 0 && return 1
    
    # In 2D, position within degree d is simply exp[2] + 1
    # (because monomials are ordered as x^d, x^(d-1)y, ..., y^d)
    @inbounds offset = indexer.degree_offsets[deg+1]
    return offset + exp[2] + 1
end

"""
    MonomialIndexer3D{MaxDeg,OffsetSize,TriangleSize}

Specialized ultra-fast indexer for 3D case.
"""
struct MonomialIndexer3D{MaxDeg,OffsetSize,TriangleSize}
    degree_offsets::SVector{OffsetSize,Int}
    triangle_numbers::SVector{TriangleSize,Int}  # Precomputed triangular numbers
end

function MonomialIndexer3D(::Val{MaxDeg}) where MaxDeg
    OffsetSize = MaxDeg + 1
    TriangleSize = MaxDeg + 1
    
    offsets = zeros(Int, OffsetSize)
    triangles = zeros(Int, TriangleSize)
    
    offsets[1] = 0
    for d in 0:MaxDeg
        triangles[d+1] = (d+1)*(d+2)÷2
        if d > 0
            offsets[d+1] = offsets[d] + triangles[d]
        end
    end
    
    MonomialIndexer3D{MaxDeg,OffsetSize,TriangleSize}(
        SVector{OffsetSize}(offsets),
        SVector{TriangleSize}(triangles)
    )
end

@inline function exp_to_position(indexer::MonomialIndexer3D{MaxDeg}, exp::SVector{3,Int}) where MaxDeg
    deg = exp[1] + exp[2] + exp[3]
    deg == 0 && return 1
    
    @inbounds offset = indexer.degree_offsets[deg+1]
    
    # Position within degree: count monomials with x exponent less than exp[1]
    pos = 1
    for i in 0:(exp[1]-1)
        remaining = deg - i
        @inbounds pos += indexer.triangle_numbers[remaining+1]
    end
    
    # Then add position for given x exponent (which is exp[3] + 1)
    pos += exp[3] + 1
    
    return offset + pos
end

# Example usage and benchmarking
function demo()
    # Generic N-dimensional version
    indexer = MonomialIndexer(Val(3), Val(20))  # 3 variables, max degree 20
    
    # Test cases
    println("Testing 3D monomials:")
    println("1 (constant): ", exp_to_position(indexer, SVector(0,0,0)))  # Should be 1
    println("x: ", exp_to_position(indexer, SVector(1,0,0)))  # Should be 2
    println("y: ", exp_to_position(indexer, SVector(0,1,0)))  # Should be 3
    println("z: ", exp_to_position(indexer, SVector(0,0,1)))  # Should be 4
    println("x²: ", exp_to_position(indexer, SVector(2,0,0)))  # Should be 5
    println("xy: ", exp_to_position(indexer, SVector(1,1,0)))  # Should be 6
    println("y²: ", exp_to_position(indexer, SVector(0,2,0)))  # Should be 7
    println("xz: ", exp_to_position(indexer, SVector(1,0,1)))  # Should be 8
    println("yz: ", exp_to_position(indexer, SVector(0,1,1)))  # Should be 9
    println("z²: ", exp_to_position(indexer, SVector(0,0,2)))  # Should be 10
    
    # Specialized 2D version (even faster)
    indexer2d = MonomialIndexer2D(Val(20))
    println("\nTesting 2D monomials:")
    println("1: ", exp_to_position(indexer2d, SVector(0,0)))  # 1
    println("x: ", exp_to_position(indexer2d, SVector(1,0)))  # 2
    println("y: ", exp_to_position(indexer2d, SVector(0,1)))  # 3
    println("x²: ", exp_to_position(indexer2d, SVector(2,0)))  # 4
    println("xy: ", exp_to_position(indexer2d, SVector(1,1)))  # 5
    println("y²: ", exp_to_position(indexer2d, SVector(0,2)))  # 6
    
    # Benchmark (comment out if BenchmarkTools not installed)
    
    exp = SVector(2, 0, 0)  # x⁵y³z²
    println("\nBenchmarking for x⁵y³z²:")
    @btime exp_to_position($indexer, $exp)
    
    exp2d = SVector(7, 3)  # x⁷y³
    println("\nBenchmarking 2D specialized for x⁷y³:")
    @btime exp_to_position($indexer2d, $exp2d)
end

# Alternative: Simple version using regular Arrays (if StaticArrays cause issues)
"""
    SimpleMonomialIndexer{N}

Simple version using regular Arrays for maximum compatibility.
"""
struct SimpleMonomialIndexer{N}
    max_deg::Int
    degree_offsets::Vector{Int}
    binomial_table::Matrix{Int}
    
    function SimpleMonomialIndexer{N}(max_deg::Int) where N
        # Precompute degree offsets
        offsets = zeros(Int, max_deg + 1)
        offsets[1] = 0
        for d in 1:max_deg
            offsets[d+1] = offsets[d] + binomial(N+d-2, d-1)
        end
        
        # Precompute binomial coefficients
        binom_size = max_deg + N
        binom_table = zeros(Int, binom_size, binom_size)
        for n in 0:(binom_size-1)
            binom_table[n+1,1] = 1
            binom_table[n+1,n+1] = 1
            for k in 1:min(n-1, binom_size-1)
                binom_table[n+1,k+1] = binom_table[n,k] + binom_table[n,k+1]
            end
        end
        
        new{N}(max_deg, offsets, binom_table)
    end
end

@inline function exp_to_position(indexer::SimpleMonomialIndexer{N}, exp::SVector{N,Int}) where N
    deg = sum(exp)
    deg == 0 && return 1
    
    # Get offset
    offset = deg <= indexer.max_deg ? indexer.degree_offsets[deg+1] : slow_degree_offset(N, deg)
    
    # Position within degree
    position = 1
    remaining_deg = deg
    
    @inbounds for i in 1:N-1
        for e_i in 0:(exp[i]-1)
            new_remaining = remaining_deg - e_i
            if new_remaining >= 0
                n, k = N-i+new_remaining-1, new_remaining
                if n < size(indexer.binomial_table, 1) && k < size(indexer.binomial_table, 2)
                    position += indexer.binomial_table[n+1, k+1]
                else
                    position += binomial(n, k)
                end
            end
        end
        remaining_deg -= exp[i]
    end
    
    return offset + position
end
demo()
# Run demo() to see it in action
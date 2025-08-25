using BenchmarkTools
using Random
using OrderedCollections

"""
    CompactIntDict

A storage-efficient dictionary for sparse integer keys and dense integer values.
Uses sorted parallel arrays with binary search.
"""
mutable struct CompactIntDict
    keys::Vector{Int32}    
    vals::Vector{UInt8}    
    
    CompactIntDict() = new(Int32[], UInt8[])
end

# Get value
function Base.getindex(sd::CompactIntDict, key::Int)
    idx = searchsortedfirst(sd.keys, Int32(key))
    (idx > length(sd.keys) || sd.keys[idx] != key) && throw(KeyError(key))
    return Int(sd.vals[idx])
end

# Set value
function Base.setindex!(sd::CompactIntDict, value::Int, key::Int)
    0 <= value <= 255 || throw(ArgumentError("Value must be in range 0-255"))
    
    idx = searchsortedfirst(sd.keys, Int32(key))
    
    if idx <= length(sd.keys) && sd.keys[idx] == key
        # Update existing
        sd.vals[idx] = UInt8(value)
    else
        # Insert new
        insert!(sd.keys, idx, Int32(key))
        insert!(sd.vals, idx, UInt8(value))
    end
    
    return value
end

# Check if key exists
Base.haskey(sd::CompactIntDict, key::Int) = 
    (idx = searchsortedfirst(sd.keys, Int32(key)); 
     idx <= length(sd.keys) && sd.keys[idx] == key)

# Get with default
Base.get(sd::CompactIntDict, key::Int, default) = 
    haskey(sd, key) ? sd[key] : default

# Delete key
function Base.delete!(sd::CompactIntDict, key::Int)
    idx = searchsortedfirst(sd.keys, Int32(key))
    if idx <= length(sd.keys) && sd.keys[idx] == key
        deleteat!(sd.keys, idx)
        deleteat!(sd.vals, idx)
    end
    return sd
end

# Length and iteration
Base.length(sd::CompactIntDict) = length(sd.keys)
Base.isempty(sd::CompactIntDict) = isempty(sd.keys)

function Base.iterate(sd::CompactIntDict, state=1)
    state > length(sd.keys) && return nothing
    return (Int(sd.keys[state]) => Int(sd.vals[state]), state + 1)
end

# Show
function Base.show(io::IO, sd::CompactIntDict)
    print(io, "CompactIntDict(")
    for i in 1:length(sd.keys)
        i > 1 && print(io, ", ")
        print(io, Int(sd.keys[i]), " => ", Int(sd.vals[i]))
    end
    print(io, ")")
end

"""
    PackedIntDict

Ultra-compact: packs key-value pairs into single integers.
22 bits for key (max ~4M), 10 bits for value (max 1023).
"""
struct PackedIntDict
    data::Vector{UInt32}  
    
    PackedIntDict() = new(UInt32[])
end

# Pack/unpack helpers
pack_kv(key::Int, val::Int) = UInt32((key << 10) | val)
unpack_key(packed::UInt32) = Int(packed >> 10)
unpack_val(packed::UInt32) = Int(packed & 0x3FF)

function Base.getindex(pd::PackedIntDict, key::Int)
    key > 0x3FFFFF && throw(ArgumentError("Key too large (max 4194303)"))
    
    # Search by key (high 22 bits)
    idx = searchsortedfirst(pd.data, pack_kv(key, 0))
    (idx > length(pd.data) || unpack_key(pd.data[idx]) != key) && throw(KeyError(key))
    return unpack_val(pd.data[idx])
end

function Base.setindex!(pd::PackedIntDict, value::Int, key::Int)
    key > 0x3FFFFF && throw(ArgumentError("Key too large (max 4194303)"))
    value > 0x3FF && throw(ArgumentError("Value too large (max 1023)"))
    
    packed = pack_kv(key, value)
    idx = searchsortedfirst(pd.data, packed, by = x -> x >> 10)
    
    if idx <= length(pd.data) && unpack_key(pd.data[idx]) == key
        pd.data[idx] = packed  # Update
    else
        insert!(pd.data, idx, packed)  # Insert
    end
    
    return value
end

Base.haskey(pd::PackedIntDict, key::Int) = 
    (idx = searchsortedfirst(pd.data, pack_kv(key, 0)); 
     idx <= length(pd.data) && unpack_key(pd.data[idx]) == key)

Base.length(pd::PackedIntDict) = length(pd.data)

function Base.iterate(pd::PackedIntDict, state=1)
    state > length(pd.data) && return nothing
    packed = pd.data[state]
    return (unpack_key(packed) => unpack_val(packed), state + 1)
end

"""
    TwoArrayDict

Simple parallel arrays without sorting - uses linear search.
Can be faster for very small collections (<20 elements).
"""
struct TwoArrayDict
    keys::Vector{Int32}
    vals::Vector{UInt8}
    
    TwoArrayDict() = new(Int32[], UInt8[])
end

function Base.getindex(td::TwoArrayDict, key::Int)
    idx = findfirst(==(Int32(key)), td.keys)
    idx === nothing && throw(KeyError(key))
    return Int(td.vals[idx])
end

function Base.setindex!(td::TwoArrayDict, value::Int, key::Int)
    0 <= value <= 255 || throw(ArgumentError("Value must be in range 0-255"))
    
    idx = findfirst(==(Int32(key)), td.keys)
    if idx !== nothing
        td.vals[idx] = UInt8(value)
    else
        push!(td.keys, Int32(key))
        push!(td.vals, UInt8(value))
    end
    return value
end

Base.haskey(td::TwoArrayDict, key::Int) = Int32(key) in td.keys
Base.length(td::TwoArrayDict) = length(td.keys)

function Base.iterate(td::TwoArrayDict, state=1)
    state > length(td.keys) && return nothing
    return (Int(td.keys[state]) => Int(td.vals[state]), state + 1)
end

# Memory comparison
function compare_memory()

    Random.seed!(42)
    
    # Test with different sizes
    for n in [10, 50, 100]
        println("\n=== Testing with $n elements ===")
        keys = sort(unique(rand(1:10000, n)))
        vals = rand(1:100, length(keys))
        
        # Standard Dict
        d1 = OrderedDict{Int32,Int8}()
        for (k, v) in zip(keys, vals)
            d1[k] = v
        end
        
        # CompactIntDict
        d2 = CompactIntDict()
        for (k, v) in zip(keys, vals)
            d2[k] = v
        end
        
        # PackedIntDict
        d3 = PackedIntDict()
        for (k, v) in zip(keys, vals)
            d3[k] = v
        end
        
        # TwoArrayDict
        d4 = TwoArrayDict()
        for (k, v) in zip(keys, vals)
            d4[k] = v
        end
        
        println("Dict{Int,Int}:  $(Base.summarysize(d1)) bytes")
        println("CompactIntDict: $(Base.summarysize(d2)) bytes")
        println("PackedIntDict:  $(Base.summarysize(d3)) bytes")
        println("TwoArrayDict:   $(Base.summarysize(d4)) bytes")
        
        # Calculate theoretical minimum
        theoretical = length(keys) * 5  # 4 bytes key + 1 byte value
        println("Theoretical min: $theoretical bytes (4+1 per entry)")
    end
end

# Performance comparison
function benchmark_implementations()

    
    Random.seed!(42)
    keys = sort(unique(rand(1:10000, 10)))
    vals = rand(1:100, length(keys))
    test_keys = rand(keys, 10)  # Random access pattern
    
    # Build dictionaries
    d1 = Dict{Int32,Int8}(k => v for (k,v) in zip(keys, vals))
    
    d2 = CompactIntDict()
    for (k, v) in zip(keys, vals)
        d2[k] = v
    end
    
    d3 = PackedIntDict()
    for (k, v) in zip(keys, vals)
        d3[k] = v
    end
    
    d4 = TwoArrayDict()
    for (k, v) in zip(keys, vals)
        d4[k] = v
    end
    
    println("\n=== Lookup Performance ($(length(test_keys)) random lookups) ===")
    
    @btime begin
        s = 0
        for k in $test_keys
            s += $d1[k]
        end
        s
    end
    
    @btime begin
        s = 0
        for k in $test_keys
            s += $d2[k]
        end
        s
    end
    
    @btime begin
        s = 0
        for k in $test_keys
            s += $d3[k]
        end
        s
    end
    
    @btime begin
        s = 0
        for k in $test_keys
            s += $d4[k]
        end
        s
    end

    println("\n=== Filling Dictionaries with $(length(keys)) keys ===")

    
    # Benchmark filling dics 
    @btime begin
        d1 = OrderedDict{Int32,Int32}()
        for (k,v) in zip($keys, $vals)
            d1[k] = v
        end
    end
    
    @btime begin
        d2 = CompactIntDict()
        for (k,v) in zip($keys, $vals)
            d2[k] = v
        end
    end
    
    @btime begin
        d3 = PackedIntDict()
        for (k,v) in zip($keys, $vals)
            d3[k] = v
        end
    end
    
    @btime begin
        d4 = TwoArrayDict()
        for (k,v) in zip($keys, $vals)
            d4[k] = v
        end
    end
    
    nothing
end

# Demo
function demo()
    println("=== CompactIntDict Demo ===")
    cd = CompactIntDict()
    cd[100] = 10
    cd[5000] = 50
    cd[300] = 30
    println(cd)
    println("Memory used: $(Base.summarysize(cd)) bytes")
    
    println("\n=== PackedIntDict Demo ===")
    pd = PackedIntDict()
    pd[100] = 10
    pd[5000] = 50
    pd[300] = 30
    println("Keys: ", [k for (k,v) in pd])
    println("Memory used: $(Base.summarysize(pd)) bytes")
    
    println("\n=== TwoArrayDict Demo ===")
    td = TwoArrayDict()
    td[100] = 10
    td[5000] = 50
    td[300] = 30
    println("Memory used: $(Base.summarysize(td)) bytes")
    
    compare_memory()
    
    # Uncomment to run benchmarks (requires BenchmarkTools)
    benchmark_implementations()
end

demo()
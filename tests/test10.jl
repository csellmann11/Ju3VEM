



using OrderedCollections

d = OrderedDict{Int,Int}(1=>2, 3=>4, 5=>6)



idx = OrderedCollections.ht_keyindex2(d, 2)

d[2] = 10

d.vals[4]

# d[1]







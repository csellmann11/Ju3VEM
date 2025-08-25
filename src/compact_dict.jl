using Dictionaries
using Chairmarks
using FixedSizeArrays

ids = distinct(rand(1:1000,100))

Base.summarysize(ids)

ids_dict = Dict{Int,Int}(i=>idsi for (i,idsi) in enumerate(ids))
Base.summarysize(ids_dict)



token_val = collect(ids)[30]


_,(_,token) = gettoken(ids,token_val)
test_token_val = gettokenvalue(ids,token)

@assert token_val == test_token_val

@b gettoken($ids,$token_val)
@b gettokenvalue($ids,$token)


function create_ids()
    ids = Indices{Int}()
    sizehint!(ids,100)
    for i in 1:100
        val = rand(1:1000)
        is_old,token = gettoken(ids,val)
        if !is_old
            insert!(ids,val)
        end
    end
    ids
end

b1 = @b create_ids()

create_ids()
function create_ids_dict()
    d = Dict{Int,Int}()
    sizehint!(d,10)
    for i in 1:10
        d[i] = rand(1:1000)
    end
    d
end

b2 = @b create_ids_dict()






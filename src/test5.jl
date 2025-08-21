





expr = Expr(:function)
push!(expr.args,:((x,)))
push!(expr.args,Expr(:if,:(x>0),:(x^2),:(x^3)))

x = 2 

eval(expr)(-2)










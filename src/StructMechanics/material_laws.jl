


module MaterialLaws

    using Symbolics
    using StaticArrays
    using LinearAlgebra

    struct Helmholtz{F<:Function,G<:Function,H<:Function,P}
        f::F                    # function to evaluate
        g::G                    # gradient
        h::H                    # hessian
        shape::Tuple{Int,Int}   # shape of argument
        pars::P                 # parameters
        
    end 


    """
        Helmholtz(f::F,dummy_arg::AbstractMatrix,pars...) where {F}

        Generates the gradient and the hessian function of the function f 
        w.r.t. the first argument of f.

        f must be a function of the form `f(∇u,pars...)` where ∇u is a matrix and `pars...` are the parameters.

    # Arguments
    - f: function to evaluate
    - D: dimension of the dummy argument
    - pars: parameters of the function f

    # Returns
    - Helmholtz object with the following fields:
    - f: function to evaluate
    - g: gradient function
    - h: hessian function
    - pars: parameters
    """
    function Helmholtz(f::F,U::I,D::I,pars...) where {F,I<:Integer}

        ∇u = Symbolics.variables(:x,1:U,1:D)
        sym_pars = Symbolics.variables.(pars)
        fsym = f(∇u,sym_pars...)[1]

        g = Symbolics.gradient(fsym,vec(∇u))
        h = Symbolics.hessian(fsym,vec(∇u))

        ffun = build_function(fsym,∇u,sym_pars...;expression = Val{false},cse = true,force_SA = true)
        gfun = build_function(g,∇u,sym_pars...;expression = Val{false},cse = true,force_SA = true)[1]
        hfun = build_function(h,∇u,sym_pars...;expression = Val{false},cse = true,force_SA = true)[1]
        
        return Helmholtz(ffun,gfun,hfun,(U,D),pars)
    end

    export Helmholtz

    @inline function eval_psi_fun(h::Helmholtz,∇u::SMatrix{U,D},
        pars = h.pars) where {U,D}
        @assert h.shape == (U,D) "argument shape mismatch"
        psi = h.f(∇u,pars...)
        return psi
    end

    @inline function eval_gradient(h::Helmholtz,∇u::SMatrix{U,D},
        pars = h.pars) where {U,D}
        @assert h.shape == (U,D) "argument shape mismatch"
        grad = SMatrix{Tuple{U,D}}(h.g(∇u,pars...)...)
        return grad
    end

    @inline function eval_hessian(h::Helmholtz,∇u::SMatrix{U,D},
        pars = h.pars) where {U,D}
        @assert h.shape == (U,D) "argument shape mismatch"
        hess = SArray{Tuple{U,D,U,D}}(h.h(∇u,pars...)...)
        return hess
    end

    @inline function eval_gradient_and_hessian(h::Helmholtz,
        ∇u::SMatrix,pars = h.pars)
        grad = eval_gradient(h,∇u,pars)
        hess = eval_hessian(h,∇u,pars)
        return grad,hess
    end

    export eval_psi_fun, eval_gradient, eval_hessian, eval_gradient_and_hessian


    @inline function Ψneo_hook(∇u::SMatrix{D,D},pars...) where D
        λ,μ = pars
        F = one(∇u) + ∇u 
        J = det(F)

        I1 = F ⋅ F
        W = μ/2*(I1-D) - μ*log(J) + λ/2*log(J)^2
        return W
    end


    @inline function Ψlin(∇u,pars...)
        λ,μ = pars
        ε = 1/2*(∇u + ∇u')
        W = λ/2 * tr(ε)^2 + μ*tr(ε*ε)
        return W 
    end

    @inline function Ψpoisson(∇u)
        energy = ∇u ⋅ ∇u
        return 1/2 * energy
    end


    export Ψneo_hook, Ψlin, Ψpoisson


    function E_ν_to_lame(E::Float64,ν::Float64)
        λ = E*ν/((1+ν)*(1-2ν))
        μ = E/(2(1+ν))
        return λ,μ
    end
    
    
    export E_ν_to_lame

    function lame_to_E_ν(λ::Float64,μ::Float64)
        ν = λ/(2(λ + μ))
        E = 2*(1+ν)*μ
        return E,ν
    end

    export lame_to_E_ν



end



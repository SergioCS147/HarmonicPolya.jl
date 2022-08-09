# HarmonicPolya

This is a package for homogeneous polynomial minimization on the sphere using harmonic hierarchies found in [Cristancho & Velasco](https://arxiv.org/abs/2202.12865) also implementing two required features: polynomial cubature/quadrature rules on the sphere using [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) and harmonic polynomial analysis on the sphere based on [Axler & Ramey](https://www.ams.org/journals/proc/1995-123-12/S0002-9939-1995-1277092-1/S0002-9939-1995-1277092-1.pdf). The implementation of polynomials uses [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) although [FixedPolynomials](https://github.com/JuliaAlgebra/FixedPolynomials.jl) is suggested for calculations using quadratures.

# Quadratures on the Sphere

The method `sphericalquadrature(n,deg)` defines the nodes and weights of a quadrature on the `n`-dimensional sphere $S^{n-1}$ for polynomials of degree $\leq$ `deg`. Example of usage:

```julia

julia> using HarmonicPolya

julia> z,w = sphericalquadrature(3,4) #this defines a quadrature on the 3-dimensional sphere for polynomials of degree =< 4

julia> poly(x) = x[1]^2*x[2]^2*x[3]^2 #defines a polynomial function on RR^3 

julia> I = sum(w*.(poly.(z))) #integrates f on the sphere using the quadrature rule

```

# Harmonic Analysis on the Sphere

Features three functions: `laplacian(poly,vars;power=1)` is the laplacian for a polynomial `poly` from [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) using variables `vars` to the power of `power` (default is 1);`generatebasissphere(deg,vars)` generates a basis for the space of homogeneous harmonic polynomials of degree `deg` in the variables `vars`; and `harmonicdecomposition(poly,vars)` produces a vector representing the decomposition in harmonic components (as a tuple (degree,component) in ascending degree) of the polynomial `poly` in the variables `vars`. 

```julia

julia> using HarmonicPolya, DynamicPolynomials

julia> @polyvar x[1:3]

julia> u = x[1]^2 + x[2]*x[3] # defines a polynomial in DynamicPolynomials

julia> ∇u = laplacian(u,x) # calculates the laplacian of u

julia> basis = generatebasissphere(3,x) # produces the basis of harmonic polynomials of degree 3 in the variables x

julia> desc = harmonicdecomposition(u,x) # produces a vector with the harmonic decomposition of u in ascending degree

julia> norm = sum(x .* x)

julia> u' = sum( map(v->norm^(floor(Int,v[1]/2))*v[2],desc)) #reconstructs u from its decomposition (working on a more elegant way)

```

# Polynomial Minimization on the Sphere

The primary function of the package. There are three different methods (as in [Cristancho & Velasco](https://arxiv.org/abs/2202.12865)): `upperbound(poly,vars,m)` produces an upper bound on the minimum of the homogeneous polynomial `poly` in variables `vars` using a quadrature of degree `m`; `lowerboundsquares(poly,vars,m)` produces a lower bound of `poly` in `vars` using the hierarchy of pure square powers of degree `m`; and `lowerboundfawzi(poly,vars,m)` produces a lower bound of `poly` in `vars` using the Fang-Fawzi hierarchy  of degree `m`.

```julia

julia> using HarmonicPolya, DynamicPolynomials

julia> @polyvar y[1:3]

julia> motzkin = y[1]^2*y[2]^4+ y[1]^4*y[2]^2+y[3]^6-3*y[1]^2*y[2]^2*y[3]^2 #defines the Motzkin polynomial (which is non-negative and homogeneous)

julia> u = upperbound(motzkin,y,20) #upper bound with quadrature of degree 20

julia> ls = lowerboundsquares(motzkin,y,20) #lower bound with squares hierarchy degree 20

julia> lf = lowerboundfawzi(motzkin,y,20) #lower bound with Fang-fawzi hierarchy degree 20

```

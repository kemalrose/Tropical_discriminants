
using Oscar

function Horn_param(data, v, λ)
    φλ = [prod(λ.^data.A[:,i]) for i in 1:n]
    φλ.*(data.B*v)
end

function get_Vdm(pts, mons; normalize = true)
    T = eltype(pts[1])
    V = zeros(T, (length(pts), length(mons)))
    for i in 1:length(pts)
        V[i, :] = [prod(pts[i].^mon) for mon in mons]
        if normalize
            V[i, :] = V[i, :]/norm(V[i, :])
        end
    end
    V
end


function newton_pol(data::Aux_data)
    n, d = data.n, data.d
    V_ambient = sampleRandom(data, 2(n - d))
    while Oscar.dim(convex_hull(hcat(V_ambient...)')) != n - d
        V_ambient = unique!([V_ambient; sampleRandom(data,(n - d))])
    end
    #V_ambient = [ [vi.num for vi in v] for v in V_ambient ]
    v_0 = V_ambient[1]
    V_0 = [data.Π*(v - v_0) for v in V_ambient]

    function getV(w)
        w_new = data.Π'*w
        w_new  = w_new * 10000 + rand(-100:100, size(w_new))
        vtx = getVertex(-w_new, data)
        monomial = data.Π * (vtx - v_0)
        monomial
    end
    vtcs, fcts, Pol = get_Polytope(V_0, getV);

    v_0, vtcs, fcts, Pol
end



function get_Polytope(V_0, getVtx)
    V_old, F_old = [], []
    V_new  = copy(V_0)
    F_new = []
    P = convex_hull(hcat(V_0...)')

    while size(V_old) != size(V_new)
        V_old = copy(V_new)
        F_old = copy(F_new)
        #println("in while loop")
        #println(V_old,F_old,getVtx)
        V_new, F_new, P = update(V_old, F_old, getVtx)
    end
    return copy(V_new), F_new, P
end

function update(V_old, F_old, getVtx)
    V_new = copy(V_old)
    P = convex_hull(hcat(V_old...)')
    F_new  = [Vector{Int64}(f.a)  for f in facets(P)]
    keep_facets = []

    for a in setdiff(F_new, F_old)
        vert = getVtx(a)
        if !in(vert, V_new)
        #if setdiff(V_old,vert) == V_old
            V_new = push!(V_new, vert)
        else

            index = findfirst(ℓ -> ℓ == a, F_new)
            keep_facets = push!(keep_facets, index)
        end
    end
    copy(V_new), [copy(F_old); F_new[keep_facets]], P
end


function verify(disc, mons, coeff, data; strategy = "substitution")
    A, B, d, n = data.A, data.B, data.d, data.n

    if strategy == "substitution"
        @var λ[1:d], v[1:n-d]
        φλ = [prod(λ.^A[:,i]) for i in 1:n]
        Φ = φλ.*(B*v)
        return expand(subs(disc, variables(disc) => Φ)) == 0
    elseif strategy == "interpolation"
        pts = [Horn_param(A,B,rand(-10:10, n-d)//10,rand(-10:10,d)//10) for j = 1:length(mons)*1.2]
        pts = Vector{Rational{BigInt}}.(pts)
        V = get_Vdm(pts, mons, normalize = false)
        S = MatrixSpace(QQ, size(V,1), size(V,2))
        if rank(S(V)) == size(V, 2) - 1
            return norm(V * coeff) == 0
        else
            print("Vdm is not of full rank!")
        end
    else
        print("Choose input strategy to either be substitution or interpolation")
    end
end




    data = data_from_matrix(A)
    d, n,  = data.d, data.n
    println("1. Compute the Newton polytope P of the discriminant")
    @time v_0, vtcs, fcts, Pol = newton_pol(data)
    println("-----------------------------------------------------------------")

    println("2. Find its lattice points")
    @time L_pts_proj = lattice_points(Pol)
    println("   P has $(length(L_pts_proj)) lattice points")
    println("   lifting the lattice points in $(n-d)-space to $(n)-space... ")
    @time mons = [ data.Π_rinv * v + v_0 for v in Vector{Int64}.(L_pts_proj)]
    mons = [ [m.num for m in mon] for mon in mons]
    println("-----------------------------------------------------------------")

    println("3. Construct the interpolation problem")
    if sample_method == "Horn"
        println("   points from Horn uniformization...")
        @time pts = [Horn_param(data,convert.(Complex{T},randn(ComplexF64, n-d)),convert.(Complex{T},exp.(2*pi*im*rand(d)))) for j = 1:convert(Int64,round(length(mons)*redundancy_factor))]
        pts = [pt/maximum(abs.(pt)) for pt ∈ pts]
    elseif sample_method == "monodromy"
        @time pts = sample_disc(A,convert(Int64,round(length(mons)*redundancy_factor)))
        pts = [pt/maximum(abs.(pt)) for pt ∈ pts]
    else println("choose sample_method = Horn or monodromy")
    end
    println("-----------------------------------------------------------------")

    println("   constructing Vandermonde matrix...")
    @time V = get_Vdm(pts, mons)
    println("   Constructed a Vandermonde matrix of size $(size(V))")
    println("-----------------------------------------------------------------")

    println("4. Find floating point coefficients...")
    if interpolation_method == "SVD"
        println("   ... using SVD")
        #println(typeof(V))
        @time U, sings, VV = svd(V)
        coeff = VV[:, end]
        #println(round.(log10.(sings)))
        err_tol = sings[end]/sings[end-1] * 100
        println(sings[end])
        println(sings[end-1])
        if sings[end] == 0.0
            err_tol = eps(T)/sings[end-1] * 100
        end
    elseif interpolation_method == "eigs"
        println("   ... using eigs")
        VtV = V' * V
        @time EVobj = eigs(VtV; nev = 2, which = :SM)
        coeff = EVobj[2][:, 1]
        println(EVobj[1])
        err_tol = abs(EVobj[1][1]/EVobj[1][2]) * 100
    else
        println("Choose interpolation_method to be either SVD or eigs.")
    end
    println("-----------------------------------------------------------------")

    println("5. Rationalizing with tolerance $(err_tol)")
    nonzinds = findall(ℓ->abs(ℓ)>err_tol, coeff)
    newcoeff = coeff/coeff[findfirst(ℓ->abs(ℓ) == minimum(abs.(coeff[nonzinds])),coeff)]
    ratcoeff = rationalize.(BigFloat.(real.(newcoeff)), tol = err_tol/minimum(abs.(coeff[nonzinds])))
    #intcoeff = coeff
    #coeff = coeff/coeff[findfirst(ℓ->abs(ℓ) == maximum(abs.(coeff)),coeff)]
    #ratcoeff = rationalize.(Float64.(real.(coeff)), tol = err_tol)

    @var a[1:n]
    monomials = [prod(a.^mon) for mon in mons]
    Δ = (ratcoeff'monomials)[1]
    #println(norm(System([Δ])(pts[1])))
    (ratcoeff'monomials)[1], ratcoeff, coeff, mons, err_tol, data
end



function sample_disc(A,nsamples)
    d,n = size(A)
    @var x[1:d] a[1:n] b[1:n-1,1:n] c[1:n-1]
    f = sum([a[i]*prod(x.^A[:,i]) for i = 1:n])
    L = b*a + c
    sys = System([subs(differentiate(f,x),x[1]=>1);L], variables = [a;x[2:end]], parameters = [b[:];c])
    # assume d < n, so that each point x can be a singularity
    xx = exp.(2*pi*im*rand(d-1))
    C = convert.(ComplexF64,jacobian(System(subs(differentiate(f,x),x=>[1;xx]),variables = a)))
    U,S,V = svd(C,full = true)
    #println(size(C))
    #println(U*diagm(S)*V-C)
    a_start = V[:,end]
    b_start = randn(n-1,n)
    c_start = - b_start*a_start
    monres = monodromy_solve(sys,[[a_start;xx]],[b_start[:];c_start])
    samples = [sol[1:length(a)] for sol ∈ solutions(monres)]
    startsols = solutions(monres)
    startparams = parameters(monres)
    while length(samples) < nsamples
        targetparams = randn(ComplexF64,length(startparams))
        R = HomotopyContinuation.solve(sys,startsols; target_parameters = targetparams, start_parameters = startparams)
        samples = vcat(samples, [sol[1:length(a)] for sol ∈ solutions(R)])
        println("found $(length(samples)) samples. ")
    end
    samples
end

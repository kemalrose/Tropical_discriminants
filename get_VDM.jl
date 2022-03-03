

struct Aux_data{T}
    A::Matrix{T}            # Input matrix. Full rank and all ones vector is in rowspan.
    d::Int              # number of rows of A
    n::Int                # number of columns of a
    B::Matrix{T}         # Gale dual to A
    U::Matrix{T}         # transform of the smith normal form S: U S V = B
    Π::Matrix{T}         # first n-d columns of U. Defines a projection onto the kernel of A.
    Π_rinv::Matrix{T}    # the first n-d rows of U^-1, a right-inverse to \Pi.
    matroid   # the matroid associated to B
    ΣB        # the Bergman fan associated to B
    Γ         # Inverses of the coefficient matrices for each cone σ and each standard basis vector e_i
    dets      # For each cone σ and each standard basis vector e_i, the determinant of the coefficient matrix of our linear system
    RR        # matrices with Z-bases for the lattices generated by all the cones
    Fσ        # matrices F such that σ = { F x <= 0 } for all cones σ
end

function data_from_matrix(A::Matrix{Int})
    d, n = size(A)
    S = MatrixSpace(ZZ, d, n)
    B = nullspace(S(A))[2]
    S,U,V = snf_with_transform(B)
    Π = Matrix{Int64}(U[1:n-d, :])
    B = Matrix{Int64}(B)
    Π_rinv  = Matrix{Int64}(inv(U))[:,1:n-d]
    U = Matrix{Int64}(U)
    matroid = Polymake.matroid.Matroid(VECTORS = B)
    ΣB = Polymake.tropical.matroid_fan{min}(matroid)
    Γ, dets, RR, Fσ = getConeData(A, ΣB)
    Aux_data(A, d, n, B, U, Π, Π_rinv, matroid, ΣB, Γ, dets, RR, Fσ)
end

function cone_containments(w, data::Aux_data)
    w_mat = MatrixSpace(QQ,data.n,1)(w)
    n = data.n
    is_contained = zeros(Int64,n, size(data.Γ, 1))
    is_in_face = zeros(Int64, n, size(data.Γ, 1))
    is_interiour = zeros(Int64, n, size(data.Γ, 1))
    for i = 1:n
        for j = 1:size(data.Γ,1)
            if data.dets[j, i] != 0
                x = Array(data.Γ[j,i] * w_mat)
                if x[end]≥0
                    vec = -data.Fσ[j]*(data.RR[j]*x[1:size(data.RR[j],2)])
                    if all(vec.≥0)
                        is_contained[i, j] = 1
                        if sum(vec.==0) > 0 || x[end] == 0
                            is_in_face[i, j] = 1
                        end
                    end
                end
            end
        end
    end
    is_contained, is_in_face, is_interiour
end


function getVertex(w, data::Aux_data)
    n = data.n

    monomial = zeros(fmpq,n)
    is_contained, is_in_face = cone_containments(w, data)
    is_generic = sum(is_in_face) == 0

    if !(is_generic)
        println("w is not generic: w = $w")
        w_new = 2000 * w + rand(-100:100, n)
        is_contained_new, is_in_face_new = cone_containments(w_new, data)
        while sum((is_contained_new - is_contained).>0) != 0
            w_new = 2000 * w + rand(-100:100, n)
            is_contained_new, is_in_face_new = cone_containments(w_new, data)
        end
        is_contained = is_contained_new
    end

    for i = 1:n
        for j = 1:size(data.Γ,1)
            if data.dets[j,i] != 0
                if is_contained[i, j] == 1
                    monomial[i] += data.dets[j,i]
                end
            end
        end
    end

    monomial
end




function getConeData(A, ΣB)
    n = length(ΣB.RAYS[1,:])-1
    maxpols = ΣB.MAXIMAL_POLYTOPES
    Ssq = MatrixSpace(QQ, n, n)
    E = Matrix{Int64}(I,n,n)
    Γ = fill(zero(Ssq),size(maxpols,1),n) # Inverses of the coefficient matrices for each cone σ and each standard basis vector e_i
    dets = fill(zero(fmpq),size(maxpols,1),n) # For each cone σ and each standard basis vector e_i, the determinant of the coefficient matrix of our linear system
    RR = []  # matrices with Z-bases for the lattices generated by all the cones
    Fσ = []  # matrices F such that σ = { F x <= 0 } for all cones σ
    for j = 1:size(maxpols,1)
        R = transpose(ΣB.RAYS[findall(ℓ->ℓ,maxpols[j,:])[1:end-1],2:end])
        S = MatrixSpace(QQ,size(R)...)
        R = S(convert(Array{Int64,2},R))
        if size(R,2) == ΣB.PROJECTIVE_DIM
            Rσ = R
        else
            SS, TT, U = snf_with_transform(R)
            Rσ = TT[:,1:ΣB.PROJECTIVE_DIM]
        end
        for i = 1:n
            mtx = Ssq(hcat(Array(Rσ),A',-E[:,i]))
            dmtx = det(mtx)
            if dmtx !=0
                Γ[j,i] = inv(mtx)
                dets[j,i] = abs(dmtx)
            end
        end
        fcts = collect(facets(positive_hull(Array(R'))))
        fctmtx = convert(Array{Int64,2},hcat([f.a for f ∈ fcts]...)')
        T = MatrixSpace(QQ,size(fctmtx)...)
        RR = push!(RR,Rσ)
        Fσ = push!(Fσ,T(fctmtx))
    end
    Γ, dets, RR, Fσ
end


function sampleRandom(data::Aux_data, nsamples)
    n = data.n
    monomials = []
    for i = 1:nsamples
        w = rand(-10000:10000,n)
        monomial = getVertex(w, data)
        monomials = unique!(push!(monomials,monomial))
    end
    monomials
end


function Horn_param(data, v, λ)
    φλ = [prod(λ.^data.A[:,i]) for i in 1:data.n]
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





function get_Vdm(A)
    T = BigFloat
    redundancy_factor = 1.2
    data = data_from_matrix(A)
    d, n  = data.d, data.n
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

    println("3. Constructing the VDM")
    println("   sampling points from Horn uniformization...")
    @time pts = [Horn_param(data,convert.(Complex{T},randn(ComplexF64, data.n-data.d)),convert.(Complex{T},exp.(2*pi*im*rand(data.d)))) for j = 1:convert(Int64,round(length(mons)*redundancy_factor))]
    pts = [pt/maximum(abs.(pt)) for pt ∈ pts]

    println("   constructing Vandermonde matrix...")
    @time V = get_Vdm(pts, mons)
    println("   Constructed a Vandermonde matrix of size $(size(V))")
    println("-----------------------------------------------------------------")
    V
end

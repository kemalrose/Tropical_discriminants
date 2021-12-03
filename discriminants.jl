

using Oscar, LinearAlgebra, HomotopyContinuation





function ishomogeneous(A)
    (n, m) = size(A)
    T = eltype(A)
    v = ones(Int64, m)
    rank(A) == rank([A; transpose(v)])
end

function tropicalized_discriminant(A)
    #if ! ishomogeneous(A)
    #    error("Not a homogeneous toric variety.")
    #end
    n, m = size(A)
    S = MatrixSpace(ZZ, n, m)
    B = nullspace(S(A))[2]
    A = Matrix{Int64}(A)
    B = Matrix{Int64}(B)

    matroid = Polymake.matroid.Matroid(VECTORS = B)
    #ΣB = Polymake.tropical.matroid_fan{min}(matroid)
    ΣB = Polymake.tropical.matroid_fan_from_flats{min}(matroid)

    max_cones = []
    raymatrix = Matrix{Int64}(ΣB.PROJECTIVE_VERTICES)[1:end-1,2:end]
    raymatrix = vcat(raymatrix, A, -A)
    n_rays = size(raymatrix, 1)

    incidence_max_cones = Matrix{Int64}(ΣB.MAXIMAL_POLYTOPES)[:, 1:end-1]
    n_cones = size(incidence_max_cones,1)
    incidence_max_cones = hcat(incidence_max_cones, ones(Int64, n_cones, 2*n))

    for ind in 1:size(incidence_max_cones, 1)
        ray_list = filter(λ -> Bool(incidence_max_cones[ind, λ]), 1:n_rays)
        cone = positive_hull(raymatrix[ray_list, :])
        push!(max_cones, Polyhedron(facets(cone),linear_span(cone)))
        #push!(max_cones,cone)
    end
    max_cones
end




A = [
1 1 1 1 1 1;
0 1 2 0 1 0;
0 0 0 1 1 2
]

A = [
1 1 1 1 0 0 0 0;
0 0 0 0 1 1 1 1;
2 3 5 7 11 13 17 19;
19 17 13 11 7 5 3 2;
]




n, m = size(A)
S = MatrixSpace(ZZ, n, m)
B = nullspace(S(A))[2]
B = Matrix{Int64}(B)

matroid = Polymake.matroid.Matroid(VECTORS = B)
#ΣB = Polymake.tropical.matroid_fan{min}(matroid)
ΣB = Polymake.tropical.matroid_fan_from_flats{min}(matroid)
maxpols = ΣB.MAXIMAL_POLYTOPES


S = MatrixSpace(QQ, m, m)
U = MatrixSpace(QQ,m,1)

E = Matrix{Int64}(I, m,m)



mlist = []
w = rand(-1000:1000,m)
w = U(w)

for k = 1:100
    w = rand(-1000:1000,m)

    w = [439, 464, 454, 360, 303, 279, 591, 583]


    w = U(w)
    monomial = []
    for i = 1:m
        mi = 0
        for j = 1:size(maxpols,1)
            #C = [ΣB.RAYS[findall(ℓ->ℓ,maxpols[j,:])[1:end-1],2:end];A;-E[i,:]']
            C = [ΣB.RAYS[findall(ℓ->ℓ,maxpols[j,:])[2:end],2:end];A;-E[i,:]']

            C = transpose(C)
            C = convert(Array{Int64,2},C)
            #dC = Int64(det(C))
            dC = det(S(C))
            number_rays = sum(ΣB.MAXIMAL_POLYTOPES[1,:]) - 1
            #println("det i = $i, j = $j: $(dC)")
            if dC != 0
                x = Oscar.solve(S(C),w)
                println(x)
                if prod([x[i] for i = 1:number_rays].>0) == 1 && x[m] > 0
                    mi += (abs(dC))
                end
            end
        end
        push!(monomial,mi)
    end
    println(monomial)
    mlist = push!(mlist,monomial)
    mlist = unique!(mlist)
end
mlist

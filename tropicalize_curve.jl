

using Oscar, IterTools


I = ideal(
x1 - (u1 - u0 * 5^(-2)),
x2 - (u1 - u0 * 5^(1)),
x3 - (u1 - u0 * 5^(3)),
x4 - u0
)

root_list = [QQ(5)^(-2), QQ(5)^(1), QQ(5)^(3)]



function trop_pols(root_list)
    A = hcat(ones(fmpq,length(root_list)), root_list)
    A = matrix(vcat(A, [0 1]))
    B = nullspace(A')[2]
    B = Matrix{Rational{Int64}}(Matrix(B))
    matroid = Polymake.matroid.Matroid(VECTORS = B)
    rk = matroid.RANK
    #bases = [ [b+=1 for b in base] for base in matroid.BASES]
    QQ5 = FlintPadicField(5,8)
    B_val = matrix(QQ5,fmpq.(B))
    #vals = [ valuation(det(B_val[base, :])) for base in bases]
    #vals = Vector{Rational{Int64}}(vals)
    #Polymake.Shell.vals = Polymake.Vector{Polymake.Rational}(vals)
    matroid = Polymake.matroid.Matroid(VECTORS = B)

#    circs = [[c+1 for c in circ] for circ in matroid.CIRCUITS]
#    matroid_val = Polymake.matroid.ValuatedMatroid{min}(BASES = matroid.BASES, VALUATION_ON_BASES = Polymake.Shell.vals)

 #   cycle = Polymake.tropical.MatroidRingCycle{min}(matroid_val)
 #   ΣB = Polymake.tropical.matroid_fan{min}(matroid)
 #   maxpols = ΣB.MAXIMAL_POLYTOPES
 #   rays = [transpose(ΣB.RAYS[findall(ℓ->ℓ,maxpols[j,:])[1:end-1],2:end-1]) for j in 1:length(root_list)]
 #   rays
    trop_circ_pol = []
    for j in 1:length(circs)
        circ = circs[j]
        pol = zeros(Int64, length(circ), size(B, 1)+1)
        for i in 1:length(circ)
            val = valuation(det(B_val[ circ[1:end.!=i],: ]))
            pol[i, circ[i]] += 1
            pol[i, end] += val
        end
        push!(trop_circ_pol, pol)
    end
    trop_circ_pol
end








A = hcat(ones(fmpq,length(root_list)), root_list)
A = matrix(vcat(A, [0 1]))
B = nullspace(A')[2]
B = Matrix{Rational{Int64}}(Matrix(B))

matroid = Polymake.matroid.Matroid(VECTORS = B)
bases = [ [b+=1 for b in base] for base in matroid.BASES]
QQ5 = FlintPadicField(5,8)
B_val = matrix(QQ5,fmpq.(B))
vals = [ valuation(det(B_val[base, :])) for base in bases]
Polymake.Shell.vals = Polymake.Vector{Polymake.Rational}(vals)
matroid = Polymake.matroid.Matroid(VECTORS = B)
matroid_val = Polymake.matroid.ValuatedMatroid{min}(BASES = matroid.BASES, VALUATION_ON_BASES = Polymake.Shell.vals, VECTORS = B)
cycle = Polymake.tropical.MatroidRingCycle{min}(matroid_val)


circs = [[c+1 for c in circ] for circ in matroid.CIRCUITS]
matroid_val = Polymake.matroid.ValuatedMatroid{min}(BASES = matroid.BASES, VALUATION_ON_BASES = Polymake.Shell.vals, VECTORS = B)
cycle = Polymake.tropical.MatroidRingCycle{min}(matroid_val)




function trop_lin_space(root_list)
    A = hcat(ones(fmpq,length(root_list)), root_list)
    A = matrix(vcat(A, [0 1]))
    B = nullspace(A')[2]
    BB = Matrix{Rational{Int64}}(
            [750//31   25//124;
            -781//31  -25//124;
                1//1     0//1;
                0//1     1//1]
            )
    BB = Matrix{Rational{Int64}}(Matrix(B))
    matroid = Polymake.matroid.Matroid(VECTORS = BB)
    matroid_val = Polymake.matroid.ValuatedMatroid{min}(BASES = matroid.BASES, VALUATION_ON_BASES = Polymake.Shell.vals, VECTORS=BB)
    cycle = Polymake.tropical.matroid_ring_cycle{min}(matroid_val)
    cycle
end

A = hcat(ones(fmpq,length(root_list)), root_list)
A = matrix(vcat(A, [0 1]))
B = nullspace(A')[2]
matroid = Polymake.matroid.Matroid(VECTORS = Matrix{Rational{Int64}}(Matrix(B)))
bases = [ [b+=1 for b in base] for base in matroid.BASES]
QQ5 = FlintPadicField(5,8)
B_val = matrix(QQ5,fmpq.(B))
vals = [ valuation(det(B_val[base, :])) for base in bases]
vals = rand(-10:10,6)
Polymake.Shell.vals = Polymake.Vector{Polymake.Rational}(vals)
BB = Matrix{Rational{Int64}}(Matrix(B))
matroid_val = Polymake.matroid.ValuatedMatroid{min}(BASES = matroid.BASES, VALUATION_ON_BASES = Polymake.Shell.vals, VECTORS=BB)
cycle = Polymake.tropical.matroid_ring_cycle{min}(matroid_val)
cycle



n_terms = [1:size(pol, 1) for pol in pols]
term_pairs = subsets.(n_terms, 2)
indices = IterTools.product(term_pairs...)

polys = []
for ind in indices
    A_eq = [(pols[i][ind[i][1],1:end-1] - pols[i][ind[i][2], 1:end-1])'  for i in 1:length(pols)]
    A_eq = vcat(A_eq...)
    b_eq = -[pols[i][ind[i][1],end] + pols[i][ind[i][2], end]  for i in 1:length(pols)]

    A_ineq = [pols[i] .- pols[i][ind[i][2], :]'  for i in 1:length(pols)]
    A_ineq = vcat(A_ineq...)[:, 1:end-1]
    b_ineq = -[pols[i] .- pols[i][ind[i][2], :]'  for i in 1:length(pols)]
    b_ineq = vcat(v_ineq...)[:, end]
    A = vcat([A_eq, -A_eq, A_ineq]...)
    b = vcat([b_eq, -b_eq, b_ineq]...)

    polyh = Oscar.Polyhedron(A, b)
    if isfeasible(polyh)
    push!(polys, polyh)
    end
end
A = [1 0; 0 1; -1 0; 0 -1]
b = [1, 1, 0, 0]
cube = Oscar.Polyhedron(A, b)

pols = trop_pols(root_list)
ind  = collect(indices)[10]
polyhedron_in_trop(pols, ind)

for index in indices
    @show index
end



function is_in_trop(pols, w)
    is_contained = true  
    for pol in pols
        a1, a2 = sort(pol * [w;1])[[1,2]]
        if a1 != a2
            return false
        end
    end
    is_contained
end



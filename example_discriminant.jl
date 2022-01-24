
using Arpack, Oscar, HomotopyContinuation, LinearAlgebra
import Pkg; Pkg.add("MultiFloats")
using MultiFloats
using GenericSVD
MultiFloats.use_bigfloat_transcendentals()
#Pkg.add("https://github.com/JuliaLinearAlgebra/GenericSVD.jl.git")


A =
[1 1 1 1 1 1;
0 1 2 0 1 0;
0 0 0 1 1 2]

A =
[
3 2 2 1 1 1 0 0 0 0;
0 1 0 2 1 0 3 2 1 0;
1 1 1 1 1 1 1 1 1 1]

A =
    [1 1 1 1 1 1  1 1 1 1 1 1;
    0 0 0 0 0 0  1 1 1 1 1 1;
    0 0 0 1 1 2 0 0 0 1 1 2;
    0 1 2 0 1 0 0 1 2 0 1 0]


#BBernd M
A =
    [ 1 1 1 1 0 0 0 0;
    0 0 0 0 1 1 1 1;
    2 3 5 7 11 13 17 19;
    19 17 13 11 7 5 3 2]


A =
    [1  1  1  1  0  0  0  0;
    0  0  0  0  1  1  1  1;
    2  3  1  0  1  0  2  0;
    2  0  0  1  0  2  2  1;]

A =
    [
    1   1   1   1   0   0   0   0 ;
    0   0   0   0   1   1   1   1;
    0   1   3   5  0 2 6 8;
    8  6 2 0   5   3   1   0;
    ]

A = A[:, 2:end]


A = [0 1 0 1 2 1; 0 0 1 1 1 2; 1 1 1 1 1 1]


AA =
[
 1 1 1 1 0 0 0;
 0 0 0 0 1 1 1;
 0 0 0 1 0 1 4;
 3 2 0 0 0 0 2;
]
A = AA

A =
    [ 1 1 1 1 0 0 0 0;
    0 0 0 0 1 1 1 1;
    2 3 5 7 11 13 17 19;
    19 17 13 11 7 5 3 2]

A = A[:, 1:end.!=4]



#A =[3 2 2 1 1 1 0 0 0 0;0 1 0 2 1 0 3 2 1 0;1 1 1 1 1 1 1 1 1 1]
Δ,ratcoeff, coeff, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x2, sample_method = "Horn", redundancy_factor = 1.5);
bool = verify(Δ, mons, coeff, data; strategy = "substitution")

coeff_fin, mons_fin, F = interpolate_discr_finitely(A, p = 373)
ratcoeff, mons, data = interpolate_discr_symbolically(A)

intcoeff = Matrix(integral_coeffs(ratcoeff))
intcoeff = [c.num for c in intcoeff]
primes = [373, 4289, 1697, 263, 941]
#How to turn finite-field element into integer?!
fin_coeffs = [ intcoeff.%p for p in primes]
fin_coeffs_pert = [ (rand(1:p-1).*(intcoeff.%p)).%p for p in primes]



inverses = zeros(Int64, (length(primes), length(primes)))
for i in 1:length(primes), j in 1:length(primes)
    if i != j
        inverses[i, j] = invmod(primes[i], primes[j])
    end
end

n = prod(primes)
idempots = ones(fmpz, length(primes))
for i in 1:length(primes)
    for j in 1:length(primes)
        if  j != i
            idempots[i] = (idempots[i] * primes[j] * inverses[j, i]) % n
        end
    end
end

function reconstruct_coeff(fin_coeffs, idempots)
    hcat(fin_coeffs ...) * idempots .%n
end

coeff_reconstr = reconstruct_coeff(fin_coeffs, idempots)
coeff_reconstr_pert = reconstruct_coeff(fin_coeffs_pert, idempots)

function integral_coeffs(ratcoeff)
    denom = lcm([BigInt(c.den) for c in ratcoeff])
    Matrix(ratcoeff.*denom)
end
function compare_results(ratcoeff, coeff_fin, F)
    intcoeff = integral_coeffs(ratcoeff)
    rats = [ F(coeff) ÷ F(intcoeff[1]) for coeff in intcoeff]
    rats2 = [ coeff ÷ coeff_fin[1] for coeff in coeff_fin]
    sum((rats-rats2).!=0) == 0
end

A =
    [ 1 1 1 1 0 0 0 0;
    0 0 0 0 1 1 1 1;
    2 3 5 7 11 13 17 19;
    19 17 13 11 7 5 3 2]

A = A[:, 1:end.!=4]


Δ,ratcoeff, coeff, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x8, sample_method = "Horn", redundancy_factor = 5)
bool = verify(Δ, mons, coeff, data; strategy = "substitution")

coeff, mons, F = interpolate_discr_finitely(A, p = 373)

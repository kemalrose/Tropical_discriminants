




using Oscar, Primes
prime_nrs = rand(primes(10000, 1000000), 30)
result_correct = fmpq.(rand(-10^30:10^30, 100000))


coeff_list = ...
fields = ...

function rational_rec(coeff_list, fields, prime_nrs)
    results_mod_p = [lift.(fields[i].(coeff_list[i]) .* fields[i](coeff_list[i][1])^(-1)) for i in 1:length(coeff_list)]
    env = Hecke.crt_env(fmpz.(prime_nrs))
    N = prod(fmpz.(prime_nrs))
    result_reconstr = zeros(fmpq, length(coeff_list[1]))
    for i in 1:length(coeff_list[1])
        ith_entries_mod_p = [res[i] for res in results_mod_p]
        rec = crt(ith_entries_mod_p , env) 
        bool, nom, denom = rational_reconstruction(rec, N)
        if !bool
            error("this is nonsense")
        end
        result_reconstr[i] = nom//denom
    end
end




indices_admissible = findall(F -> all([F(r.den) != 0 for r in result_correct]), fin_fields)
prime_nrs_ad = prime_nrs[indices_admissible]
env = Hecke.crt_env(fmpz.(prime_nrs_ad))
N = prod(fmpz.(prime_nrs_ad))

result_reconstr = zeros(fmpq, length(result_correct))
for i in 1:length(result_correct)
    entry = result_correct[i]
    entry_mod_primes = [lift(F(entry)) for F in fin_fields[indices_admissible]]
    rec = crt(entry_mod_primes, env) % N
    bool, nom, denom = rational_reconstruction(rec, N)
    if !bool
        error("this is nonsense")
    end
    result_reconstr[i] = nom//denom
end


F_applicables = [ all([ lift(F(nr.den)) !=0  for nr in result_correct  ]) for F in fin_fields]


env = Hecke.crt_env(fmpz.([2, 3, 5]))

rational_reconstruction(2, 5)

Hecke.solve_dixon



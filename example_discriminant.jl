


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

A =
    [ 1 1 1 1 0 0 0 0;
    0 0 0 0 1 1 1 1;
    2 3 5 7 11 13 17 19;
    19 17 13 11 7 5 3 2]

A = A[:, 1:end.!=4]


Δ,ratcoeff, coeffs, mons, err_tol, data = interpolate_discr(A; interpolation_method = "SVD", T = Float64x8, sample_method = "Horn", redundancy_factor = 5)
bool = verify(Δ, mons, coeffs, data; strategy = "substitution")

coeffs, mons, F = interpolate_discr_finitely(A, p = 373)
bool = verify(Δ, mons, coeffs, data; strategy = "substitution")

A = [
1 1 1 1 1 0 0 0 0 0
0 0 0 0 0 1 1 1 1 1
3 1 1 0 0 0 2 0 1 0
0 2 0 1 0 3 1 1 0 0
]
A = A[:, 2:end ]
P = newton_pol(A)

coeffs, mons, F = interpolate_discr_finitely(A, p = 373)


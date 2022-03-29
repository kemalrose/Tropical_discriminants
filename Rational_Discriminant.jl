
function interpolate_discr_symbolically(A; redundancy_factor = 1.2)
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

    println("3. Construct the interpolation problem")
        println("   points from Horn uniformization...")
        @time pts = [Horn_param(data, QQ.(rand(1:50,data.n-data.d).//rand(1:50,data.n-data.d)),QQ.(rand(1:50,data.d).//rand(1:50,data.d))) for j = 1:convert(Int64,round(length(mons)*redundancy_factor))]


    println("-----------------------------------------------------------------")

    println("   constructing Vandermonde matrix...")
    @time V = get_Vdm(pts, mons, normalize = false)
    println("   Constructed a Vandermonde matrix of size $(size(V))")
    println("-----------------------------------------------------------------")

    println("4. Find floating point coefficients...")
    println("Convert Vandermonde Matrix to OSCAR")
        @time V = MatrixSpace(QQ, size(V,1), size(V,2))(V)

        println("Compute rational coeffs")
        @time d, ratcoeff = nullspace(V)
        if d > 1
            print("V is not of corank 1!")
        end
    println("-----------------------------------------------------------------")

    @var a[1:data.n]
    #println(norm(System([Δ])(pts[1])))
    ratcoeff, mons, data
end




function construct_Vandermonde_rat(A; redundancy_factor = 1.2, mons = nothing, data = nothing)
    if data == nothing
        data = data_from_matrix(A)
    end

    

    if mons == nothing
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
    end

    println("3. Construct the interpolation problem")
    println("   points from Horn uniformization...")


    @time pts_int = [Horn_param(data, fmpz.(rand(-50:50,data.n-data.d)),fmpz.(rand(-50:50,data.d))) for j = 1:convert(Int64,round(length(mons)*redundancy_factor))]
    println("-----------------------------------------------------------------")
    println("   constructing Vandermonde matrix...")
    Mat_Space = Oscar.MatrixSpace(QQ, length(pts_rat), length(mons))
    V_rat = zero(Mat_Space)

    for i in 1:length(pts_rat)
        V_rat[i, :] = [prod(pts_rat[i].^mon) for mon in mons]
    end

    println("   Constructed a Vandermonde matrix of size $(size(V))")
    println("-----------------------------------------------------------------")

    V
end

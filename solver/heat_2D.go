package solver

// square plate
type Heat_2D_plate struct {
	N int // number of spatial points: N = M

	L     float64 // Dimmensions of plate (Width = Height)
	alpha float64 // thermal diffusivity (for the form du/dt = alpha * d2u/d2x)

	delx, delt float64 // time/space steps delx=dely

	// boundary conditions
	B  BoundarySet2D
	u0 plate_init_fcn

	U      Matrix
	u_last Matrix
	x, y   []float64

	CurrentTime float64
	K           int

	// internal params
	r float64

	// for BTCS
	F, F_inv, G Matrix
	Bp_inv, Cp  []Matrix
	Dp          []Row_vec

	/*
		a []float64
		b []float64
		c []float64
		d []float64

		cp []float64
		dp []float64
	*/
}

func (b *Heat_2D_plate) Create(num_points int, plate_width float64, B BoundarySet2D, alpha float64, u0 plate_init_fcn) {
	b.N = num_points
	b.L = plate_width
	b.delx = plate_width / float64(num_points)
	b.B = B
	b.u0 = u0

	// calc delt to enforce stability
	max_delt := 0.5 * b.delx * b.delx * b.delx * b.delx / (alpha * (b.delx*b.delx + b.delx*b.delx))
	b.delt = 100 * max_delt

	b.r = alpha * b.delt / (b.delx * b.delx)
	/*
		b.a = make([]float64, num_points+1)
		b.b = make([]float64, num_points+1)
		b.c = make([]float64, num_points+1)
		b.d = make([]float64, num_points+1)

		b.cp = make([]float64, b.N+1)
		b.dp = make([]float64, b.N+1)
	*/

	// allocate memory for data
	b.x = make([]float64, num_points+1)
	for i := 0; i <= num_points; i++ {
		b.x[i] = float64(i) * b.delx
	}
	b.y = make([]float64, num_points+1)
	for i := 0; i <= num_points; i++ {
		b.y[i] = float64(i) * b.delx
	}

	b.U = Make_matrix(num_points + 1)
	for n := 0; n <= num_points; n++ {
		b.U[n] = make([]float64, num_points+1)
	}

	b.u_last = Make_matrix(num_points + 1)

	// initialize using the user supplied function
	b.Reset()
}

func (b *Heat_2D_plate) Reset() {
	for r := 0; r <= b.N; r++ {
		y := b.y[r]
		for c := 0; c <= b.N; c++ {
			x := b.x[c]
			b.U[r][c] = b.u0(x, y, b.L)
		}
	}

	b.CurrentTime = 0
	b.K = 0

	//Print_matrix(b.U)
}

func (b *Heat_2D_plate) Update_FTCS() {
	b.CurrentTime += b.delt
	b.K++

	copy(b.u_last, b.U)

	switch b.B.Top.Type {
	case ConstantTemp:
		T := b.B.Top.Value
		for c := 1; c < b.N; c++ {
			b.U[0][c] = T
		}
	case ConstantFlux:
	}

	switch b.B.Right.Type {
	case ConstantTemp:
		T := b.B.Right.Value
		for r := 1; r < b.N; r++ {
			b.U[r][b.N] = T
		}
	case ConstantFlux:
	}

	switch b.B.Bot.Type {
	case ConstantTemp:
		T := b.B.Bot.Value
		for c := 1; c < b.N; c++ {
			b.U[b.N][c] = T
		}
	case ConstantFlux:
	}

	switch b.B.Left.Type {
	case ConstantTemp:
		T := b.B.Left.Value
		for r := 1; r < b.N; r++ {
			b.U[r][0] = T
		}
	case ConstantFlux:
	}

	for r := 1; r <= (b.N - 1); r++ {
		above := b.u_last[r-1]
		cent := b.u_last[r]
		below := b.u_last[r+1]
		for c := 1; c <= (b.N - 1); c++ {
			b.U[r][c] = cent[c] + b.r*(above[c]-2*cent[c]+below[c]) + b.r*(cent[c-1]-2*cent[c]+cent[c+1])
		}
	}
}

func calc_Bn(Bp, Cp_prev, F, G Matrix) {
	// G is diagonal (except 0 at the ends)
	N := len(Bp)
	for r := 0; r < N; r++ {
		k := G[r][r]
		for c := 0; c < N; c++ {
			Bp[r][c] = F[r][c] - Cp_prev[r][c]*k
		}
	}
}
func calc_Cn(Cp, Bp_inv, G Matrix) {
	// G is diagonal (except 0 at the ends)
	Matrix_mul(Cp, Bp_inv, G)
}

func make_matrix_arr(size, num int) []Matrix {
	arr := make([]Matrix, num)
	for n := 0; n < num; n++ {
		arr[n] = Make_matrix(size)
	}
	return arr
}

func make_vec_arr(size, num int) []Row_vec {
	arr := make([]Row_vec, num)
	for n := 0; n < num; n++ {
		arr[n] = make(Row_vec, size)
	}
	return arr
}

func copy_matrix(dst, src Matrix) {
	N := len(dst)
	for n := 0; n < N; n++ {
		copy(dst[n], src[n])
	}
}

func (b *Heat_2D_plate) Init_BTCS() {
	r := b.r
	N := b.N

	f := 1 + 4*r
	g := -r

	b.F = Make_matrix(N + 1)
	b.G = Make_matrix(N + 1)

	b.F[0][0] = 1
	b.F[N][N] = 1
	for n := 1; n < N; n++ {
		b.F[n][n-1] = g
		b.F[n][n] = f
		b.F[n][n+1] = g

		b.G[n][n] = g
	}

	b.F_inv = Make_matrix(N + 1)
	Invert_gauss(b.F, b.F_inv)

	b.Cp = make_matrix_arr(N+1, N-2)
	b.Bp_inv = make_matrix_arr(N+1, N-1)
	b.Dp = make_vec_arr(N+1, N-1)

	copy_matrix(b.Bp_inv[0], b.F_inv)
	calc_Cn(b.Cp[0], b.F_inv, b.G)
	for i := 1; i < (N - 2); i++ {
		calc_Bn(b.Bp_inv[i], b.Cp[i-1], b.F, b.G)
		Invert_gauss(b.Bp_inv[i], b.Bp_inv[i]) // TODO: invert in place with no alloc
		calc_Cn(b.Cp[i], b.Bp_inv[i], b.G)
	}
	calc_Bn(b.Bp_inv[N-2], b.Cp[N-3], b.F, b.G)
	Invert_gauss(b.Bp_inv[N-2], b.Bp_inv[N-2]) // TODO: invert in place with no alloc
}

func (b *Heat_2D_plate) Update_BTCS() {
	g := -b.r
	b.CurrentTime += b.delt
	b.K++

	copy_matrix(b.u_last, b.U)

	num_rows := b.N - 1
	tmp := make(Row_vec, b.N+1)

	for r := 0; r < num_rows; r++ {
		// calculate D for each row
		// first row: D_1 = R_1 - G*R_0
		// last row:  B_M = R_M - G*R_{M+1}

		//fmt.Println("row ", r+1)
		//fmt.Println(" R =", b.u_last[r+1])

		if r == 0 {
			b.U[r+1][0] = b.B.Left.Value
			b.U[r+1][b.N] = b.B.Right.Value
			b.Dp[r][0] = b.u_last[r+1][0]
			for c := 1; c < (b.N); c++ {
				b.Dp[r][c] = b.u_last[r+1][c] - b.B.Top.Value*g

				b.U[0][c] = b.B.Top.Value
			}
			b.Dp[r][b.N] = b.u_last[r+1][b.N]
		} else if r == (num_rows - 1) { // N-2
			b.U[r+1][0] = b.B.Left.Value
			b.U[r+1][b.N] = b.B.Right.Value
			b.Dp[r][0] = b.u_last[r+1][0]
			for c := 1; c < (b.N); c++ {
				b.Dp[r][c] = b.u_last[r+1][c] - b.B.Bot.Value*g

				b.U[b.N][c] = b.B.Bot.Value
			}
			b.Dp[r][b.N] = b.u_last[r+1][b.N]
		} else {
			b.Dp[r][0] = b.B.Left.Value
			for c := 1; c <= (b.N); c++ {
				b.Dp[r][c] = b.u_last[r+1][c]
			}
			b.Dp[r][b.N] = b.B.Right.Value

			b.U[r+1][0] = b.B.Left.Value
			b.U[r+1][b.N] = b.B.Right.Value
		}

		//fmt.Println(" U[] =", b.U[r+1])
		//fmt.Println(" D =", b.Dp[r])

		// calculate Dp for each row
		if r == 0 {
			//Dp_1 = Bp_inv_1*D_1
			Matrix_vec_mul(tmp, b.Bp_inv[r], b.Dp[r])
			copy(b.Dp[r], tmp)
		} else if r == (num_rows - 1) { // N-2
			//Dp_{M-1} = Bp_inv_{M-1}*D_{M-1}
			// G*Dp_{i-1} -> G is diagonal, so == g*Dp_{i-1}
			for c := 1; c < (b.N); c++ {
				b.Dp[r][c] = b.Dp[r][c] - g*b.Dp[r-1][c]
			}
			Matrix_vec_mul(tmp, b.Bp_inv[r], b.Dp[r])
			copy(b.Dp[r], tmp)
		} else {
			//Dp_i = Bp_inv_i*(D_i - G*Dp_{i-1})
			// G*Dp_{i-1} -> G is diagonal, so == g*Dp_{i-1}
			for c := 1; c < (b.N); c++ {
				b.Dp[r][c] = b.Dp[r][c] - g*b.Dp[r-1][c]
			}
			Matrix_vec_mul(tmp, b.Bp_inv[r], b.Dp[r])
			copy(b.Dp[r], tmp)
		}

		//fmt.Println(" Dp =", b.Dp[r])
	}

	// back-substitute to calculate X
	X := make(Row_vec, b.N+1)
	copy(X, b.Dp[num_rows-1]) // i=N-2
	// X_{M-1} = Dp_{M-1}
	for c := 1; c < b.N; c++ {
		b.U[b.N-1][c] = X[c]
	}

	for r := (num_rows - 2); r >= 0; r-- {
		// X_i = Dp_i - Cp_i*X_{i+1}
		Matrix_vec_mul(tmp, b.Cp[r], X) // X is X_{i+1}, tmp is X_i

		for c := 1; c < b.N; c++ {
			X[c] = b.Dp[r][c] - tmp[c]
			b.U[r+1][c] = X[c]
		}
	}

	//Print_matrix(b.U)
}

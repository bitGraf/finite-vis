package solver

import "math"

type Heat_2D_plate struct {
	N int // number of spatial points: x
	M int // number of spatial points: y

	Width, Height float64 // Dimmensions of plate
	alpha         float64 // thermal diffusivity (for the form du/dt = alpha * d2u/d2x)

	delx, dely, delt float64 // time/space steps

	// boundary conditions
	B  BoundarySet2D
	u0 plate_init_fcn

	U      [][]float64
	u_last [][]float64
	x, y   []float64

	CurrentTime float64
	K           int

	// internal params
	rx, ry float64

	/*
		a []float64
		b []float64
		c []float64
		d []float64

		cp []float64
		dp []float64
	*/
}

func (b *Heat_2D_plate) Create(num_points_x, num_points_y int, plate_width, plate_height float64, B BoundarySet2D, alpha float64, u0 plate_init_fcn) {
	b.N = num_points_x
	b.M = num_points_y
	b.Width = plate_width
	b.Height = plate_height
	b.delx = plate_width / float64(num_points_x)
	b.dely = plate_height / float64(num_points_y)
	b.B = B
	b.u0 = u0

	// calc delt to enforce stability
	max_delt_x := (0.5) * b.delx * b.delx / alpha
	max_delt_y := (0.5) * b.dely * b.dely / alpha
	b.delt = 0.7 * math.Min(max_delt_x, max_delt_y) // idk if this makes any sense...

	b.rx = alpha * b.delt / (b.delx * b.delx)
	b.ry = alpha * b.delt / (b.dely * b.dely)
	/*
		b.a = make([]float64, num_points+1)
		b.b = make([]float64, num_points+1)
		b.c = make([]float64, num_points+1)
		b.d = make([]float64, num_points+1)

		b.cp = make([]float64, b.N+1)
		b.dp = make([]float64, b.N+1)
	*/

	// allocate memory for data
	b.x = make([]float64, num_points_x+1)
	for i := 0; i <= num_points_x; i++ {
		b.x[i] = float64(i) * b.delx
	}
	b.y = make([]float64, num_points_y+1)
	for i := 0; i <= num_points_y; i++ {
		b.y[i] = float64(i) * b.dely
	}

	b.U = make([][]float64, num_points_x+1)
	for n := 0; n <= num_points_x; n++ {
		b.U[n] = make([]float64, num_points_y+1)
	}

	b.u_last = make([][]float64, num_points_x+1)
	for n := 0; n <= num_points_x; n++ {
		b.u_last[n] = make([]float64, num_points_y+1)
	}

	// initialize using the user supplied function
	b.Reset()
}

func (b *Heat_2D_plate) Reset() {
	for n := 0; n <= b.N; n++ {
		x := b.x[n]
		for m := 0; m <= b.M; m++ {
			y := b.y[m]
			b.U[n][m] = b.u0(x, y, b.Width, b.Height)
		}
	}

	b.CurrentTime = 0
	b.K = 0
}

func (b *Heat_2D_plate) Update_FTCS() {
	rx := b.rx
	ry := b.ry
	b.CurrentTime += b.delt
	b.K++

	copy(b.u_last, b.U)

	switch b.B.Top.Type {
	case ConstantTemp:
		T := b.B.Top.Value
		for n := 1; n <= b.N; n++ {
			b.U[n][0] = T
		}
	case ConstantFlux:
	}

	switch b.B.Right.Type {
	case ConstantTemp:
		T := b.B.Right.Value
		for m := 1; m <= b.M; m++ {
			b.U[b.N][m] = T
		}
	case ConstantFlux:
	}

	switch b.B.Bot.Type {
	case ConstantTemp:
		T := b.B.Bot.Value
		for n := b.N - 1; n >= 0; n-- {
			b.U[n][b.M] = T
		}
	case ConstantFlux:
	}

	switch b.B.Left.Type {
	case ConstantTemp:
		T := b.B.Left.Value
		for m := b.M - 1; m >= 0; m-- {
			b.U[0][m] = T
		}
	case ConstantFlux:
	}

	for n := 1; n <= (b.N - 1); n++ {
		left := b.u_last[n-1]
		cent := b.u_last[n]
		right := b.u_last[n+1]
		for m := 1; m <= (b.M - 1); m++ {
			b.U[n][m] = cent[m] + rx*(left[m]-2*cent[m]+right[m]) + ry*(cent[m-1]-2*cent[m]+cent[m+1])
		}
	}
}

/*
func (b *Heat_1D_bar) Update_BTCS() {
	r := b.r
	b.CurrentTime += b.delt

	copy(b.u_last, b.U)

	// construct tridagonal coefficients
	switch b.B1.Type {
	case ConstantTemp:
		//b.a[0] unused
		b.b[0] = 1
		b.c[0] = 0
		b.d[0] = b.B1.Value
	case ConstantFlux:
		k1 := b.B1.Value
		//b.a[0] unused
		b.b[0] = 1 + r
		b.c[0] = -r
		b.d[0] = b.u_last[0] - r*k1*b.delx
	}

	b.cp[0] = b.c[0] / b.b[0]
	b.dp[0] = b.d[0] / b.b[0]

	for i := 1; i < b.N; i++ {
		b.a[i] = -r
		b.b[i] = 1 + 2*r
		b.c[i] = -r
		b.d[i] = b.u_last[i]

		// reduced c coeffiecient
		b.cp[i] = b.c[i] / (b.b[i] - b.a[i]*b.cp[i-1])

		// reduced d coeffiecient
		b.dp[i] = (b.d[i] - b.a[i]*b.dp[i-1]) / (b.b[i] - b.a[i]*b.cp[i-1])
	}

	switch b.B2.Type {
	case ConstantTemp:
		b.a[b.N] = 0
		b.b[b.N] = 1
		// b.c[b.N] unused
		b.d[b.N] = b.B2.Value
	case ConstantFlux:
		k2 := b.B2.Value
		b.a[b.N] = -r
		b.b[b.N] = 1 + r
		// b.c[b.N] unused
		b.d[b.N] = b.u_last[b.N] + r*k2*b.delx
	}

	// b.cp[b.N] unused
	b.dp[b.N] = (b.d[b.N] - b.a[b.N]*b.dp[b.N-1]) / (b.b[b.N] - b.a[b.N]*b.cp[b.N-1])

	// substitue back using these new coefficients
	b.U[b.N] = b.dp[b.N]
	for i := (b.N - 1); i >= 0; i-- {
		b.U[i] = b.dp[i] - b.cp[i]*b.U[i+1]
	}
}

func (b *Heat_1D_bar) Update_CTCS() {
	r := b.r
	b.CurrentTime += b.delt

	copy(b.u_last, b.U)

	// construct tridagonal coefficients
	switch b.B1.Type {
	case ConstantTemp:
		//b.a[0] unused
		b.b[0] = 1
		b.c[0] = 0
		b.d[0] = b.B1.Value
	case ConstantFlux:
		k1 := b.B1.Value
		//b.a[0] unused
		b.b[0] = 2 + r
		b.c[0] = -r
		b.d[0] = (2-r)*b.u_last[0] + r*b.u_last[1] - 2*r*k1*b.delx
	}

	b.cp[0] = b.c[0] / b.b[0]
	b.dp[0] = b.d[0] / b.b[0]

	for i := 1; i < b.N; i++ {
		b.a[i] = -r
		b.b[i] = 2 + 2*r
		b.c[i] = -r
		b.d[i] = r*b.u_last[i-1] + (2-2*r)*b.u_last[i] + r*b.u_last[i+1]

		// reduced c coeffiecient
		b.cp[i] = b.c[i] / (b.b[i] - b.a[i]*b.cp[i-1])

		// reduced d coeffiecient
		b.dp[i] = (b.d[i] - b.a[i]*b.dp[i-1]) / (b.b[i] - b.a[i]*b.cp[i-1])
	}

	switch b.B2.Type {
	case ConstantTemp:
		b.a[b.N] = 0
		b.b[b.N] = 1
		// b.c[b.N] unused
		b.d[b.N] = b.B2.Value
	case ConstantFlux:
		k2 := b.B2.Value
		b.a[b.N] = -r
		b.b[b.N] = 2 + r
		// b.c[b.N] unused
		b.d[b.N] = r*b.u_last[b.N-1] + (2-r)*b.u_last[b.N] + 2*r*k2*b.delx
	}

	// b.cp[b.N] unused
	b.dp[b.N] = (b.d[b.N] - b.a[b.N]*b.dp[b.N-1]) / (b.b[b.N] - b.a[b.N]*b.cp[b.N-1])

	// substitue back using these new coefficients
	b.U[b.N] = b.dp[b.N]
	for i := (b.N - 1); i >= 0; i-- {
		b.U[i] = b.dp[i] - b.cp[i]*b.U[i+1]
	}
}
*/

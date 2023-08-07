package solver

type Heat_1D_bar struct {
	N int // number of spatial points

	L     float64 // length of bar
	alpha float64 // thermal diffusivity (for the form du/dt = alpha * d2u/d2x)

	delx, delt float64 // time/space steps

	// boundary conditions
	B1, B2 BoundaryCondition
	u0     bar_init_fcn

	U      []float64
	u_last []float64
	x      []float64

	CurrentTime float64

	// internal params
	r float64

	a []float64
	b []float64
	c []float64
	d []float64

	cp []float64
	dp []float64
}

func (b *Heat_1D_bar) Create(num_points int, bar_length float64, B1, B2 BoundaryCondition, alpha float64, u0 bar_init_fcn) {
	b.N = num_points
	b.L = bar_length
	b.delx = bar_length / float64(num_points)
	b.B1 = B1
	b.B2 = B2
	b.u0 = u0

	// calc delt to enforce stability
	max_delt := (0.5) * b.delx * b.delx / alpha
	b.delt = max_delt

	b.r = alpha * b.delt / (b.delx * b.delx)
	b.a = make([]float64, num_points+1)
	b.b = make([]float64, num_points+1)
	b.c = make([]float64, num_points+1)
	b.d = make([]float64, num_points+1)

	b.cp = make([]float64, b.N+1)
	b.dp = make([]float64, b.N+1)

	// allocate memory for data
	b.x = make([]float64, num_points+1)
	for i := 0; i <= num_points; i++ {
		b.x[i] = float64(i) * b.delx
	}

	b.U = make([]float64, num_points+1)
	b.u_last = make([]float64, num_points+1)

	// initialize using the user supplied function
	b.Reset()
}

func (b *Heat_1D_bar) Reset() {
	for n := 0; n <= b.N; n++ {
		x := b.x[n]
		b.U[n] = b.u0(x, b.L)
	}

	switch b.B1.Type {
	case ConstantTemp:
		//b.U[0] = b.B1.Value
	case ConstantFlux:
		//k1 := b.B1.Value
		//b.U[0] = b.U[1] + k1*b.delx
	}

	switch b.B2.Type {
	case ConstantTemp:
		//b.U[b.N] = b.B2.Value
	case ConstantFlux:
		//k2 := b.B2.Value
		//b.U[b.N] = b.U[b.N-1] + 3*k2*b.delx
	}

	b.CurrentTime = 0
}

func (b *Heat_1D_bar) Update_FTCS() {
	r := b.r
	b.CurrentTime += b.delt

	copy(b.u_last, b.U)

	switch b.B1.Type {
	case ConstantTemp:
		b.U[0] = b.B1.Value
	case ConstantFlux:
		k1 := b.B1.Value
		b.U[0] = (1-r)*b.u_last[0] + r*b.u_last[1] - r*k1*b.delx
	}

	for n := 1; n <= (b.N - 1); n++ {
		b.U[n] = r*b.u_last[n-1] + (1-2*r)*b.u_last[n] + r*b.u_last[n+1]
	}

	switch b.B2.Type {
	case ConstantTemp:
		b.U[b.N] = b.B2.Value
	case ConstantFlux:
		k2 := b.B2.Value
		b.U[b.N] = r*b.u_last[b.N-1] + (1-r)*b.u_last[b.N] + r*k2*b.delx
	}
}

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

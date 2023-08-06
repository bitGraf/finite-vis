package solver

type Method_type int

const (
	Forward Method_type = iota
	Backward
	CrankNicolson

	MaxTypes
)

type Heat_1D_bar struct {
	N int // number of spatial points

	L     float64 // length of bar
	alpha float64 // thermal diffusivity (for the form du/dt = alpha * d2u/d2x)

	delx, delt float64 // time/space steps

	// boundary conditions
	T1, T2 float64
	u0     bar_init_fcn

	U      []float64
	u_last []float64
	x      []float64

	CurrentTime float64

	// internal params
	internal struct {
		r float64

		a []float64
		b []float64
		c []float64
		d []float64
	}
}

// u0(x,L) = ...
type bar_init_fcn func(float64, float64) float64

func (b *Heat_1D_bar) Create(num_points int, bar_length, T1, T2, alpha float64, u0 bar_init_fcn) {
	b.N = num_points
	b.L = bar_length
	b.delx = bar_length / float64(num_points)
	b.T1 = T1
	b.T2 = T2
	b.u0 = u0

	// calc delt to enforce stability
	max_delt := (0.5) * b.delx * b.delx / alpha
	b.delt = max_delt

	b.internal.r = alpha * b.delt / (b.delx * b.delx)
	b.internal.a = make([]float64, num_points+1)
	b.internal.b = make([]float64, num_points+1)
	b.internal.c = make([]float64, num_points+1)
	b.internal.d = make([]float64, num_points+1)

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

	b.U[0] = b.T1
	b.U[b.N] = b.T2
	b.CurrentTime = 0
}

func (b *Heat_1D_bar) Update_FTCS() {
	r := b.internal.r
	b.CurrentTime += b.delt

	copy(b.u_last, b.U)

	b.U[0] = b.T1
	for n := 1; n <= (b.N - 1); n++ {
		b.U[n] = b.u_last[n] + r*(b.u_last[n-1]-2*b.u_last[n]+b.u_last[n+1])
	}
	b.U[b.N] = b.T2
}

func (b *Heat_1D_bar) Update_BTCS() {
	r := b.internal.r
	b.CurrentTime += b.delt

	copy(b.u_last, b.U)

	// construct tridagonal coefficients - these don't change
	b.internal.a[0] = 0 // unused
	for i := 1; i < b.N; i++ {
		b.internal.a[i] = -r
	}
	b.internal.a[b.N] = 0

	b.internal.b[0] = 1
	for i := 1; i < (b.N + 1); i++ {
		b.internal.b[i] = 1 + 2*r
	}
	b.internal.b[b.N] = 1

	b.internal.c[0] = 0
	for i := 1; i < b.N; i++ {
		b.internal.c[i] = -r
	}

	// reduced c coeffiecient
	cp := make([]float64, b.N+1)
	cp[0] = b.internal.c[0] / b.internal.b[0]
	for i := 1; i < b.N; i++ {
		cp[i] = b.internal.c[i] / (b.internal.b[i] - b.internal.a[i]*cp[i-1])
	}

	// d coefficient is constant in this method
	//b.internal.d[0] = b.T1
	for i := 0; i <= b.N; i++ {
		b.internal.d[i] = b.u_last[i]
	}
	//b.internal.d[b.N] = b.T2

	// reduced d coeffiecient (Thomas Method)
	dp := make([]float64, b.N+1)
	dp[0] = b.internal.d[0] / b.internal.b[0]
	for i := 1; i < b.N+1; i++ {
		dp[i] = (b.internal.d[i] - b.internal.a[i]*dp[i-1]) / (b.internal.b[i] - b.internal.a[i]*cp[i-1])
	}

	// substitue back using these new coefficients
	b.U[b.N] = dp[b.N]
	for i := (b.N - 1); i >= 0; i-- {
		b.U[i] = dp[i] - cp[i]*b.U[i+1]
	}
}

func (b *Heat_1D_bar) Update_CTCS() {
	r := b.internal.r
	b.CurrentTime += b.delt

	copy(b.u_last, b.U)

	// construct tridagonal coefficients - these don't change
	b.internal.a[0] = 0
	for i := 1; i < b.N; i++ {
		b.internal.a[i] = -1
	}
	b.internal.a[b.N] = 0

	b.internal.b[0] = 1
	for i := 1; i < (b.N + 1); i++ {
		b.internal.b[i] = 2 * (1 + r) / r
	}
	b.internal.b[b.N] = 1

	b.internal.c[0] = 0
	for i := 1; i < b.N; i++ {
		b.internal.c[i] = -1
	}

	// reduced c coeffiecient (Thomas Method)
	cp := make([]float64, b.N+1)
	cp[0] = b.internal.c[0] / b.internal.b[0]
	for i := 1; i < b.N; i++ {
		cp[i] = b.internal.c[i] / (b.internal.b[i] - b.internal.a[i]*cp[i-1])
	}

	// d coefficient
	b.internal.d[0] = b.u_last[0]
	for i := 1; i < b.N; i++ {
		b.internal.d[i] = b.u_last[i-1] + (2*(1-r)/r)*b.u_last[i] + b.u_last[i+1]
	}
	b.internal.d[b.N] = b.u_last[b.N]

	// reduced d coeffiecient (Thomas Method)
	dp := make([]float64, b.N+1)
	dp[0] = b.internal.d[0] / b.internal.b[0]
	for i := 1; i < b.N+1; i++ {
		dp[i] = (b.internal.d[i] - b.internal.a[i]*dp[i-1]) / (b.internal.b[i] - b.internal.a[i]*cp[i-1])
	}

	// substitue back using these new coefficients
	b.U[b.N] = dp[b.N]
	for i := (b.N - 1); i >= 0; i-- {
		b.U[i] = dp[i] - cp[i]*b.U[i+1]
	}
}

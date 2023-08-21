package solver

import "fmt"

type Row_vec []float64
type Matrix []Row_vec

func Make_matrix(size int) Matrix {
	m := make(Matrix, size)
	for n := 0; n < size; n++ {
		m[n] = make(Row_vec, size)
	}
	return m
}

// Inverts Matrix A, storing the result into A_inv without modifying A
func Invert_gauss(A, A_inv Matrix) {
	num_rows := len(A)
	num_cols := len(A[0])
	count := 0

	//fmt.Println("Inverting Matrix...")

	// create augmented half, as identity
	tmp := make(Matrix, num_rows)
	for r := 0; r < num_rows; r++ {
		tmp[r] = make(Row_vec, num_cols)
		copy(tmp[r], A[r])

		for c := 0; c < num_cols; c++ {
			A_inv[r][c] = 0 // Correct   :3
			//A_inv[r][r] = 0  // Incorrct >:(
		}
		A_inv[r][r] = 1
	}
	//copy(tmp, A)

	//Print_aug_matrix(tmp, A_inv)
	//fmt.Println("-------------------------")

	// perform gaussian elim. on (tmp|A_inv)
	for r := 0; r < num_rows; r++ {
		f := tmp[r][r]
		//fmt.Printf("Row %v, pivot=%v\n", r, f)
		for n := 0; n < num_cols; n++ {
			tmp[r][n] /= f
			A_inv[r][n] /= f

			count += 2
		}
		//Print_aug_matrix(A, A_inv)
		for n := 0; n < num_rows; n++ {
			if n == r {
				continue
			}
			g := tmp[n][r]
			for i := 0; i < num_cols; i++ {
				tmp[n][i] -= tmp[r][i] * g
				A_inv[n][i] -= A_inv[r][i] * g

				count += 4
			}
		}
		//Print_aug_matrix(tmp, A_inv)
		//fmt.Println("-------------------------")
	}

	//fmt.Println("Done!")
}

// C = AB <- no allocs done here
func Matrix_mul(C, A, B Matrix) {
	N := len(A)

	for r := 0; r < N; r++ {
		for c := 0; c < N; c++ {
			C[r][c] = 0
			for k := 0; k < N; k++ {
				C[r][c] += A[r][k] * B[k][c]
			}
		}
	}
}

// C = Av <- no allocs done here (v is a column_vec actually...)
func Matrix_vec_mul(C Row_vec, A Matrix, v Row_vec) {
	N := len(A)

	for r := 0; r < N; r++ {
		C[r] = 0
		for k := 0; k < N; k++ {
			C[r] += A[r][k] * v[k]
		}
	}
}

func Print_matrix(A Matrix) {
	num_rows := len(A)
	num_cols := len(A[0])
	for i := 0; i < num_rows; i++ {
		fmt.Printf("| ")
		for j := 0; j < num_cols; j++ {
			fmt.Printf("%6.3f, ", A[i][j])
		}
		fmt.Printf("|\n")
	}
	fmt.Printf("\n")
}

func Print_aug_matrix(A, aug Matrix) {
	num_rows := len(A)
	num_cols := len(A[0])
	for i := 0; i < num_rows; i++ {
		fmt.Printf("| ")
		for j := 0; j < num_cols; j++ {
			fmt.Printf("%6.3f, ", A[i][j])
		}
		fmt.Printf("| ")
		for j := 0; j < num_cols; j++ {
			fmt.Printf("%6.3f, ", aug[i][j])
		}
		fmt.Printf("|\n")
	}
	fmt.Printf("\n")
}

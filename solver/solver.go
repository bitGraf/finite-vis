package solver

type Method_type int

const (
	Forward Method_type = iota
	Backward
	CrankNicolson

	MaxMethodTypes
)

type Boundary_type int

const (
	ConstantTemp Boundary_type = iota
	ConstantFlux

	MaxBoundaryTypes
)

type BoundaryCondition struct {
	Type  Boundary_type
	Value float64
}

type BoundarySet1D struct {
	Left  BoundaryCondition
	Right BoundaryCondition
}

type BoundarySet2D struct {
	Top   BoundaryCondition
	Bot   BoundaryCondition
	Left  BoundaryCondition
	Right BoundaryCondition
}

// u0(x,L) = ...
type bar_init_fcn func(float64, float64) float64

// u0(x,y,W) = ...
type plate_init_fcn func(float64, float64, float64) float64

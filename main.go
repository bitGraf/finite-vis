package main

import (
	"flag"
	"fmt"
	"math"
	"runtime"
	"visualizer/solver"

	"github.com/fogleman/gg"
	"github.com/go-gl/gl/all-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
)

func init() {
	runtime.LockOSThread()
}

func write_to_image(filename string, u solver.Matrix, min_temp, max_temp float64, cell_width, cell_height int, map_type colormap_type, num_bands int) {
	N := len(u) - 1
	M := len(u[0]) - 1

	x_pixels := cell_width * (N + 1)
	y_pixels := cell_height * (M + 1)
	//bytes_per_pixel := 3 //rgb
	//num_pixels := x_pixels * y_pixels
	//num_bytes := num_pixels * bytes_per_pixel

	dc := gg.NewContext(x_pixels, y_pixels)

	//client_width := x_pixels
	//client_height := y_pixels

	dx := float64(cell_width)  //float64(1.0) / float64(N+1) * float64(client_width)
	dy := float64(cell_height) //float64(1.0) / float64(M+1) * float64(client_height)

	border_x := float64(0) // 1 pixels
	border_y := float64(0) // 1 pixels

	for n := 0; n <= N; n++ {
		xpos := dx * float64(n)

		col := u[n]
		for m := 0; m <= M; m++ {
			ypos := dy * float64(m)

			T := col[m]
			f := (T - min_temp) / (max_temp - min_temp)
			color := colormapN(f, map_type, num_bands)

			dc.SetRGB(color[0], color[1], color[2])
			dc.DrawRectangle(xpos, ypos, dx-border_x, dy-border_y)
			//renderer.DrawRect(xpos, start_y-ypos, dx-border_x, dy-border_y, float32(color[0]), float32(color[1]), float32(color[2]))
			dc.Fill() // fill the circle
		}
	}

	err := dc.SavePNG(filename) // save it
	if err != nil {
		panic(fmt.Errorf("Failed to write to image!"))
	}
}

func key_callback(key glfw.Key, scancode int, action glfw.Action, mods glfw.ModifierKey) {
	if (key == glfw.KeyR) && (action == glfw.Press) {
		bar.Reset()
		plate.Reset()
	}

	if (key == glfw.KeyF) && (action == glfw.Press) {
		method = solver.Forward
		fmt.Printf("Switching to forward-explicit method\n")
		bar.Reset()
		plate.Reset()
	}
	if (key == glfw.KeyB) && (action == glfw.Press) {
		method = solver.Backward
		fmt.Printf("Switching to backward-impicit method\n")
		bar.Reset()
		plate.Reset()
	}
	if (key == glfw.KeyC) && (action == glfw.Press) {
		method = solver.CrankNicolson
		fmt.Printf("Switching to Crank-Nicolson method\n")
		bar.Reset()
		plate.Reset()
	}

	if (key == glfw.KeySpace) && (action == glfw.Press) {
		auto_advance = !auto_advance
		if auto_advance {
			fmt.Printf("Playing automatically\n")
		} else {
			fmt.Printf("Use arrow kesys to advance\n")
		}
	}

	if !auto_advance {
		if (key == glfw.KeyRight) && (action != glfw.Release) {
			switch method {
			case solver.Forward:
				plate.Update_FTCS()
			case solver.Backward:
				plate.Update_BTCS()
			case solver.CrankNicolson:
				//plate.Update_CTCS()
			}
		}
	}

	if (key == glfw.KeyF4) && (action == glfw.Press) {
		filename := fmt.Sprintf("plate_%d.png", plate.K)
		fmt.Printf("Saving screenshot to %v...", filename)
		write_to_image(filename, plate.U, 0, 400, 8, 8, Plasma, 20)
		fmt.Printf("Done.\n")
	}
}

var (
	bar          solver.Heat_1D_bar
	plate        solver.Heat_2D_plate
	method       solver.Method_type = solver.Forward
	auto_advance bool               = true
	decimation   int                = 1
)

func main() {
	m := solver.Make_matrix(4)
	m[0][0] = 1
	m[0][1] = 2

	fmt.Println("m =", m)

	// set these by command-line args
	colormap_name := flag.String("cmap", "plasma", "Name of colormap to use. To get a list of supported maps, use -list-cmaps flag")
	colormap_list := flag.Bool("list-cmaps", false, "Print list of supported maps")
	tight := flag.Bool("tight", false, "Shrink window to min size")
	window_width := flag.Int("width", 800, "Window width")
	window_height := flag.Int("height", 840, "Window height")
	num_bands := flag.Int("bands", 10, "Number of discrete color bands. 0 means continuous")
	flag.IntVar(&decimation, "dec", 1, "Time decimation. How many timesteps occur for each visual update")
	flag.Parse()
	if *tight {
		*window_height = 104
	}

	others := flag.Args() // other args not parsed by flag
	if len(others) > 0 {
		fmt.Printf("Additional args:\n")
		for n := 0; n < len(others); n++ {
			fmt.Printf(" [%v]\n", others[n])
		}
	}

	if *colormap_list {
		print_map_types()
		return
	}

	map_type := map_type_from_string(colormap_name)

	fmt.Printf("Creating window [%vx%v], using %v colormap\n", *window_width, *window_height, map_type.ToString())
	//colormap_info()

	// Initialize GLFW
	var window Window
	if !window.InitGLFW(*window_width, *window_height, key_callback) {
		fmt.Println("Error initializing GLFW")
		return
	}
	// Create render data
	var renderer Renderer
	renderer.Init(window.width, window.height)

	defer func() {
		renderer.Shutdown()
		window.Close()
	}()

	DrawColormap := func(width, height int, min_temp, max_temp float64) {
		N := 255
		border_x := 2 / float32(width)  // 2 pixels
		border_y := 3 / float32(height) // 2 pixels
		dx := (float32(1.0) - 2*border_x) / float32(N+1)
		dy := float32(40) / float32(height)
		ypos := 1.0 - dy

		renderer.DrawRect(0, ypos, 1, dy, 0, 0, 0)
		for n := 0; n <= N; n++ {
			xpos := dx * float32(n)

			color := colormapN(float64(xpos), map_type, *num_bands)

			renderer.DrawRect(xpos+border_x, ypos+border_y, dx, dy-2*border_y, float32(color[0]), float32(color[1]), float32(color[2]))
		}

		//renderer.font.Printf(12, (1-ypos)*float32(height)-12, 0.5, "%.0fK", min_temp)
		//renderer.font.Printf(float32(width-80), (1-ypos)*float32(height)-12, 0.5, "%.0fK", max_temp)
		renderer.font.Printf(12, (1-ypos)*float32(height)-12, 0.5, "%.0fC", min_temp)
		renderer.font.Printf(float32(width-80), (1-ypos)*float32(height)-12, 0.5, "%.0fC", max_temp)
	}

	/*
		DrawHeatBar := func(u []float64, width, height int, min_temp, max_temp float64) {
			N := len(u) - 1
			dx := float32(1.0) / float32(N+1)
			dy := float32(64) / float32(height) // 100 pixels
			border_x := 1 / float32(width)      // 1 pixels
			var ypos float32
			if *tight {
				ypos = 0
			} else {
				ypos = float32(1.0/2.0) - dy/2.0
			}

			for n := 0; n <= N; n++ {
				xpos := dx * float32(n)

				T := u[n]
				f := (T - min_temp) / (max_temp - min_temp)

				color := colormapN(f, map_type, *num_bands)

				renderer.DrawRect(xpos, ypos, dx-border_x, dy, float32(color[0]), float32(color[1]), float32(color[2]))
			}
		}
	*/
	DrawHeatPlate := func(u solver.Matrix, screen_width, screen_height int, min_temp, max_temp float64) {
		N := len(u) - 1
		M := len(u[0]) - 1

		client_width := screen_width
		client_height := screen_height - 40

		dx := float32(1.0) / float32(N+1) * float32(client_width) / float32(screen_width)
		dy := float32(1.0) / float32(M+1) * float32(client_height) / float32(screen_height)

		border_x := 0 / float32(screen_width)  // 1 pixels
		border_y := 0 / float32(screen_height) // 1 pixels

		start_y := float32(client_height)/float32(screen_height) - dy

		for c := 0; c <= N; c++ {
			ypos := dy * float32(c)

			col := u[c]
			for r := 0; r <= M; r++ {
				xpos := dx * float32(r)

				T := col[r]
				f := (T - min_temp) / (max_temp - min_temp)
				color := colormapN(f, map_type, *num_bands)

				renderer.DrawRect(xpos, start_y-ypos, dx-border_x, dy-border_y, float32(color[0]), float32(color[1]), float32(color[2]))
			}
		}
	}

	// Create heat bar
	B1 := solver.BoundaryCondition{Type: solver.ConstantTemp, Value: 300}
	B2 := solver.BoundaryCondition{Type: solver.ConstantFlux, Value: 0}
	//dist := func(x, L float64) float64 { return 400 * (1 - math.Abs(x-(L/2))) }
	dist := func(x, L float64) float64 { return 0 }
	bar.Create(20, 1.0, B1, B2, 111.0, dist)

	// Create heat plate
	B := solver.BoundarySet2D{
		Top:   solver.BoundaryCondition{Type: solver.ConstantTemp, Value: 300},
		Right: solver.BoundaryCondition{Type: solver.ConstantTemp, Value: 200},
		Bot:   solver.BoundaryCondition{Type: solver.ConstantTemp, Value: 100},
		Left:  solver.BoundaryCondition{Type: solver.ConstantTemp, Value: 400},
	}
	peaks := func(x, y, W float64) float64 {
		// x and y come in on the ranges [0,W] and [0,H]
		// transform them into the range [-3,3] and [-3,3]
		x = (x - (W / 2)) * (6 / W)
		y = (y - (W / 2)) * (6 / W)

		// peaks function
		p1 := 3 * (1 - x) * (1 - x) * math.Exp(-x*x-(y+1)*(y+1))
		p2 := -10 * (x/5 - x*x*x - y*y*y*y*y) * math.Exp(-x*x-y*y)
		p3 := math.Exp(-(x+1)*(x+1)-y*y) / 3
		peaks := p1 + p2 + p3 // on the range [-5, 5]

		// transform to the range [0 300]
		return (peaks + 5) * 300 / 10
	}
	sinc_2d := func(x, y, W float64) float64 {
		// x and y come in on the ranges [0,W] and [0,H]
		// transform them into the range [-5,5]
		x = (x - (W / 2)) * (2 / W) * 15
		y = (y - (W / 2)) * (2 / W) * 15

		// sinc function = sin(r)/r where r = sqrt(x^2 + y^2)
		// at r=0, sinc(0) = 1
		r := math.Sqrt(x*x + y*y)
		var sinc float64
		if math.Abs(r) < 1e-6 {
			sinc = 1
		} else {
			sinc = math.Sin(r) / r // on the range [~, 1]
		}

		// transform to the range [~, 300]
		return 1000 * ((sinc)*5/6 + 1/6)
	}
	zero_2d := func(x, y, W float64) float64 { return 0 }
	_ = peaks(0, 0, 0)
	_ = sinc_2d(0, 0, 0)
	_ = zero_2d(0, 0, 0)
	plate.Create(101, 1, B, 1.0, zero_2d)
	plate.Init_BTCS()

	//gl.ClearColor(0.4, 0.2, 0.5, 1.0)
	gl.ClearColor(0.3, 0.3, 0.3, 1.0)
	for !window.ShouldClose() {
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

		DrawColormap(window.width, window.height, 0, 500)
		//DrawHeatBar(bar.U, window.width, window.height, 0, 400)
		DrawHeatPlate(plate.U, window.width, window.height, 0, 500)

		renderer.font.SetColor(1.0, 1.0, 1.0, 1.0)
		renderer.font.Printf(12, -12+float32(*window_height), 0.5, "t = %5.3f ms [%d]", plate.CurrentTime*1000, plate.K)

		if auto_advance {
			for d := 0; d < decimation; d++ {
				/*
					switch method {
					case solver.Forward:
						bar.Update_FTCS()
					case solver.Backward:
						bar.Update_BTCS()
					case solver.CrankNicolson:
						bar.Update_CTCS()
					}
				*/
				switch method {
				case solver.Forward:
					plate.Update_FTCS()
				case solver.Backward:
					plate.Update_BTCS()
				case solver.CrankNicolson:
					//plate.Update_CTCS()
				}

			}

			/*
				filename := fmt.Sprintf("out/plate_%d.png", plate.K)
				fmt.Printf("Saving screenshot to %v...", filename)
				write_to_image(filename, plate.U, 0, 400, 8, 8, Plasma, 20)
				fmt.Printf("Done.\n")
			*/

			//if plate.K >= 1000 {
			//	window.window.SetShouldClose(true)
			//}
		}

		//window.window.SetShouldClose(true)

		glfw.PollEvents()
		window.SwapBuffers()
	}
}

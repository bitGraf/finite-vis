package main

import (
	"flag"
	"fmt"
	"runtime"
	"visualizer/solver"

	"github.com/go-gl/gl/all-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
)

func init() {
	runtime.LockOSThread()
}

func key_callback(key glfw.Key, scancode int, action glfw.Action, mods glfw.ModifierKey) {
	if (key == glfw.KeyR) && (action == glfw.Press) {
		bar.Reset()
	}

	if (key == glfw.KeyF) && (action == glfw.Press) {
		method = solver.Forward
		fmt.Printf("Switching to forward-explicit method\n")
		bar.Reset()
	}
	if (key == glfw.KeyB) && (action == glfw.Press) {
		method = solver.Backward
		fmt.Printf("Switching to backward-impicit method\n")
		bar.Reset()
	}
	if (key == glfw.KeyC) && (action == glfw.Press) {
		method = solver.CrankNicolson
		fmt.Printf("Switching to Crank-Nicolson method\n")
		bar.Reset()
	}
}

var (
	bar    solver.Heat_1D_bar
	method solver.Method_type = solver.Forward
)

func main() {
	// set these by command-line args
	colormap_name := flag.String("cmap", "plasma", "Name of colormap to use. To get a list of supported maps, use -list-cmaps flag")
	colormap_list := flag.Bool("list-cmaps", false, "Print list of supported maps")
	tight := flag.Bool("tight", false, "Shrink window to min size")
	window_width := flag.Int("width", 1024, "Window width")
	window_height := flag.Int("height", 240, "Window height")
	num_bands := flag.Int("bands", 10, "Number of discrete color bands. 0 means continuous")
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

		renderer.font.Printf(12, (1-ypos)*float32(height)-12, 0.5, "%.0fK", min_temp)
		renderer.font.Printf(float32(width-80), (1-ypos)*float32(height)-12, 0.5, "%.0fK", max_temp)
	}
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

	// Create heat bar
	bar.Create(100, 1.0, 300, 400, 111.0, func(x, L float64) float64 { return 0 })

	//gl.ClearColor(0.4, 0.2, 0.5, 1.0)
	gl.ClearColor(0.3, 0.3, 0.3, 1.0)
	for !window.ShouldClose() {
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

		DrawColormap(window.width, window.height, 0, 400)
		DrawHeatBar(bar.U, window.width, window.height, 0, 400)
		renderer.font.SetColor(1.0, 1.0, 1.0, 1.0)
		renderer.font.Printf(12, -12+float32(*window_height), 0.5, "t = %5.3f ms", bar.CurrentTime*1000)

		switch method {
		case solver.Forward:
			bar.Update_FTCS()
		case solver.Backward:
			bar.Update_BTCS()
		case solver.CrankNicolson:
			bar.Update_CTCS()
		}

		glfw.PollEvents()
		window.SwapBuffers()
	}
}

package main

import (
	"fmt"

	"github.com/go-gl/gl/all-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
)

type Window struct {
	window *glfw.Window
	width  int
	height int
}

type key_callback_t func(key glfw.Key, scancode int, action glfw.Action, mods glfw.ModifierKey)

func (this *Window) InitGLFW(window_width, window_height int, key_callback key_callback_t) bool {
	// Initialize GLFW
	if err := glfw.Init(); err != nil {
		return false
	}

	this.width = window_width
	this.height = window_height

	glfw.WindowHint(glfw.Resizable, glfw.False)
	glfw.WindowHint(glfw.ContextVersionMajor, 4)
	glfw.WindowHint(glfw.ContextVersionMinor, 6)
	glfw.WindowHint(glfw.OpenGLProfile, glfw.OpenGLCoreProfile)
	glfw.WindowHint(glfw.OpenGLForwardCompatible, glfw.True)

	var err error
	this.window, err = glfw.CreateWindow(this.width, this.height, "Visualizer", nil, nil)
	if err != nil {
		return false
	}
	this.window.MakeContextCurrent()

	// Create OpenGL context
	if err := gl.Init(); err != nil {
		return false
	}
	version := gl.GoStr(gl.GetString(gl.VERSION))
	fmt.Println("OpenGL version", version)

	// Set callbacks
	this.window.SetKeyCallback(func(_ *glfw.Window, key glfw.Key, scancode int, action glfw.Action, mods glfw.ModifierKey) {
		if (key == glfw.KeyEscape) && (action == glfw.Press) {
			this.window.SetShouldClose(true)
		}

		key_callback(key, scancode, action, mods)
	})

	return true
}

func (this *Window) Close() {
	this.window = nil
	glfw.Terminate()
}

func (this *Window) ShouldClose() bool {
	return this.window.ShouldClose()
}

func (this *Window) SwapBuffers() {
	this.window.SwapBuffers()
}

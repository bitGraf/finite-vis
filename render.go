package main

import (
	"fmt"
	"log"
	"strings"

	"github.com/go-gl/gl/all-core/gl"
	"github.com/nullboundary/glfont"
)

type Renderer struct {
	vbo, vao uint32
	shader   uint32

	// shader locations
	xpos_loc   int32
	ypos_loc   int32
	xscale_loc int32
	yscale_loc int32

	color_loc int32

	font *glfont.Font
}

func (this *Renderer) Init(width, height int) {
	// create rect buffer
	rect := []float32{
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		1.0, 1.0, 0.0,

		0.0, 0.0, 0.0,
		1.0, 1.0, 0.0,
		0.0, 1.0, 0.0,
	}

	gl.GenBuffers(1, &this.vbo)
	gl.BindBuffer(gl.ARRAY_BUFFER, this.vbo)
	gl.BufferData(gl.ARRAY_BUFFER, 4*len(rect), gl.Ptr(rect), gl.STATIC_DRAW)

	gl.GenVertexArrays(1, &this.vao)
	gl.BindVertexArray(this.vao)
	gl.EnableVertexAttribArray(0)
	gl.BindBuffer(gl.ARRAY_BUFFER, this.vbo)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 0, nil)

	// create simple shader
	vertSrc := `
		#version 150
		in vec3 vp;

		// coords in range [0,1] from bottom left corner
		uniform float xpos;
		uniform float ypos;
		uniform float xscale;
		uniform float yscale;

		void main() {
			gl_Position = vec4(2*(vp.x*xscale+xpos) - 1, 2*(vp.y*yscale+ypos) - 1, vp.z, 1.0);
		}
	` + "\x00" // null-terminate so C-bindings can handle it
	fragSrc := `
		#version 150
		out vec4 frag_colour;

		uniform vec3 color;

		void main() {
			frag_colour = vec4(color, 1);
		}
	` + "\x00" // null-terminate so C-bindings can handle it
	this.shader = CreateShader(vertSrc, fragSrc)

	this.xpos_loc = gl.GetUniformLocation(this.shader, gl.Str("xpos"+"\x00"))
	this.ypos_loc = gl.GetUniformLocation(this.shader, gl.Str("ypos"+"\x00"))
	this.xscale_loc = gl.GetUniformLocation(this.shader, gl.Str("xscale"+"\x00"))
	this.yscale_loc = gl.GetUniformLocation(this.shader, gl.Str("yscale"+"\x00"))
	this.color_loc = gl.GetUniformLocation(this.shader, gl.Str("color"+"\x00"))

	gl.UseProgram(this.shader)
	gl.Uniform1f(this.xpos_loc, 0)
	gl.Uniform1f(this.ypos_loc, 0)
	gl.Uniform1f(this.xscale_loc, 1)
	gl.Uniform1f(this.yscale_loc, 1)
	gl.Uniform3f(this.color_loc, 1, 1, 1)

	//load font (fontfile, font scale, window width, window height
	var err error
	this.font, err = glfont.LoadFont("FiraCode-Regular.ttf", int32(52), width, height)
	if err != nil {
		log.Panicf("LoadFont: %v", err)
	}
}

func (this *Renderer) Shutdown() {
	//if this.font != nil {
	//	this.font.Release()
	//}
}

func (this *Renderer) DrawRect(posx, posy float32, scalex, scaley float32, r, g, b float32) {
	gl.UseProgram(this.shader)

	gl.Uniform1f(this.xpos_loc, posx)
	gl.Uniform1f(this.ypos_loc, posy)
	gl.Uniform1f(this.xscale_loc, scalex)
	gl.Uniform1f(this.yscale_loc, scaley)

	gl.Uniform3f(this.color_loc, r, g, b)

	gl.BindVertexArray(this.vao)
	gl.DrawArrays(gl.TRIANGLES, 0, 6)
}

func CreateShader(vertSrc, fragSrc string) uint32 {
	vert := CompileShaderComponent(vertSrc, gl.VERTEX_SHADER)
	if vert == 0 {
		return 0
	}
	frag := CompileShaderComponent(fragSrc, gl.FRAGMENT_SHADER)
	if frag == 0 {
		return 0
	}

	prog := gl.CreateProgram()
	gl.AttachShader(prog, vert)
	gl.AttachShader(prog, frag)
	gl.LinkProgram(prog)

	// check for errors
	var status int32
	gl.GetProgramiv(prog, gl.LINK_STATUS, &status)
	if status == gl.FALSE {
		var logLength int32
		gl.GetProgramiv(prog, gl.INFO_LOG_LENGTH, &logLength)

		log := strings.Repeat("\x00", int(logLength+1))
		gl.GetProgramInfoLog(prog, logLength, nil, gl.Str(log))

		gl.DeleteProgram(prog)
		gl.DeleteShader(vert)
		gl.DeleteShader(frag)

		fmt.Printf("failed to link: %v", log)

		return 0
	}

	return prog
}

func CompileShaderComponent(src string, shader_type uint32) uint32 {
	shader := gl.CreateShader(shader_type)

	csources, free := gl.Strs(src)
	gl.ShaderSource(shader, 1, csources, nil)
	free()
	gl.CompileShader(shader)

	// check for errors
	var status int32
	gl.GetShaderiv(shader, gl.COMPILE_STATUS, &status)
	if status == gl.FALSE {
		var logLength int32
		gl.GetShaderiv(shader, gl.INFO_LOG_LENGTH, &logLength)

		log := strings.Repeat("\x00", int(logLength+1))
		gl.GetShaderInfoLog(shader, logLength, nil, gl.Str(log))

		fmt.Printf("failed to compile %v: %v", src, log)

		return 0
	}

	return shader
}

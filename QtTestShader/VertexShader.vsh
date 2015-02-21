#version 120

uniform mat4 mat;
in vec4 vertex;

void main(void)
{
	gl_Position = vertex * mat;
	gl_TexCoord[0] = gl_MultiTexCoord0;
}
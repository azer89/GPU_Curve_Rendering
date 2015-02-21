#version 120


uniform vec4 insideColor;

uniform vec4 outsideColor;


void main(void)
{
	vec3 uv = gl_TexCoord[0].xyz;
	float val = pow(uv.x, 3) - uv.y * uv.z;

	vec4 mycol = vec4(outsideColor.xyz, 1.0);

	if(val > 0.0 )  mycol = vec4(insideColor.xyz, 1.0);
	
	gl_FragColor = mycol;
}


#version 120

// convex or concave (not yet implemented)
uniform int flag;

void main(void)
{
	vec2 p = gl_TexCoord[0].st;
	
	if( (p.s * p.s - p.t) > 0.0 )
		discard;

	gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);
}

 
#version 330

uniform bool bLighting;
uniform mat4 projection, modelview;
uniform mat3 normalMatrix;

in vec3 position;
in vec3 normal;
in vec4 color;

out vec4 fragColor;

void main()
{
	// Transform position from pixel coordinates to clipping coordinates
	gl_Position = projection * modelview * vec4(position, 1.0);
	
	fragColor = color;
	if(bLighting)
		fragColor = fragColor * normalize(normalMatrix * normal).z;
	
	// Modulate color with lighting and apply gamma correction
	fragColor = pow(fragColor, vec4(1.0 / 2.1));
}


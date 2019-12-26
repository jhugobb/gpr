#version 330

uniform mat4 projection, modelview;
uniform vec4 color;

in vec3 position;

out vec4 fragColor;

void main()
{
	// Transform position from pixel coordinates to clipping coordinates
	gl_Position = projection * modelview * vec4(position, 1.0);
	
	fragColor = pow(color, vec4(1.0 / 2.1));
}


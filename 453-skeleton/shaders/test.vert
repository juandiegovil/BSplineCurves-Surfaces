#version 330 core
layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 color;

uniform mat4 M;
uniform mat4 V;
uniform mat4 P;
out vec3 fragColor;

void main() {
	fragColor = color;
	gl_Position = P * V * M * vec4(pos, 1.0);
	//gl_Position = vec4(pos, 1.0);
}

#version 460 core

in vec4 color;
in vec2 uv;

out vec4 fragColor;
uniform sampler2D flagSampler;

void main(void)
{
	// Set the color of this fragment to the interpolated color
	// value computed by the rasterizer.
	
	vec4 kdTex= texture(flagSampler, uv);
	fragColor = kdTex;
}

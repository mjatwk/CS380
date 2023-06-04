// Input vertex data, different for all executions of this shader.
attribute vec3 vertexPosition_modelspace;

// Values that stay constant for the whole mesh.
uniform mat4 depthMVP;

varying float fragmentdepth;

void main(){
	gl_Position =  depthMVP * vec4(vertexPosition_modelspace, 1.0);
	fragmentdepth = gl_Position.z; // Store the depth value in the varying variable.
}
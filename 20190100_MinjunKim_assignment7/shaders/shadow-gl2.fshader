varying float fragmentdepth;

void main(){
	gl_FragColor = vec4(fragmentdepth, fragmentdepth, fragmentdepth, 1.0); // Use fragmentdepth to output the color.
}

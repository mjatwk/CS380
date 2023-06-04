// Varying values from the vertex shaders
varying vec2 UV;
varying vec4 ShadowCoord;

// Values that stay constant for the whole mesh.
uniform sampler2D myTextureSampler;
uniform sampler2D shadowMap;

void main(){

    // Light emission properties
    vec3 LightColor = vec3(1,1,1);

    // Material properties
    vec3 MaterialDiffuseColor = texture2D( myTextureSampler, UV ).rgb;

    vec3 ShadowXYZ = (ShadowCoord.xyz / ShadowCoord.w);

    float visibility = 1.0;
    float offset = 0.005;
    if (texture2D( shadowMap, ShadowXYZ.xy ).r < ShadowXYZ.z - offset) {
        visibility = 0.0;
    }

    gl_FragColor = vec4(visibility * MaterialDiffuseColor * LightColor, 1.0);
}

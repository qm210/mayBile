#version 300 es
precision highp float;
out vec4 fragColor;

uniform float iTexSize;
uniform float iBlockOffset;
uniform float iSampleRate;

#define PI radians(180.)
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }

float Tsample;

void main()
{
    Tsample = 1./iSampleRate;
    float t = (iBlockOffset + gl_FragCoord.x + gl_FragCoord.y*iTexSize) * Tsample;
    vec2 s = vec2(_sin(210. * t), _sin(210. * .985 * t));
    vec2 v  = floor((0.5+0.5*s)*65535.0);
    vec2 vl = mod(v,256.0)/255.0;
    vec2 vh = floor(v/256.0)/255.0;
    fragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}

"use strict";

var gl;

var sampleRate = 44100;
var texSize = 512;
var duration = 100; // enforces nBlocks = 1

var vertexShaderSource = `#version 300 es
in vec2 position;
//uniform vec2 u_resolution;
void main() {
    gl_Position = vec4(position, 0, 1);
}
`;

function main() {
    var canvas = document.getElementById("canvas");
    gl = canvas.getContext("webgl2");
    if (!gl) {
        return;
    }

    var fragmentShaderSource = null;
    var request = new XMLHttpRequest();
    request.open('GET', './feelit.glsl', false);
    request.onload = function() {
        fragmentShaderSource = request.responseText
            .replace('#version 130', '#version 300 es\nprecision highp float;\nout vec4 fragColor;\n')
            .replace('gl_FragColor', 'fragColor');
    }
    request.send();
    if (!fragmentShaderSource) {
        document.body.innerHTML = "<b>No Fragment Shader Found</b><br><br>" + document.body.innerHTML;
        return;
    }

    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);

    var program = webglUtils.createProgramFromSources(gl, [vertexShaderSource, fragmentShaderSource]);

    var lower = -1;
    var upper = 1;
    var vertices = [
        lower, lower,
        upper, lower,
        lower, upper,
        upper, lower,
        upper, upper,
        lower, upper,
    ];

    var vertexBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);

    var vertexAttributeLocation = gl.getAttribLocation(program, 'position');
    var resolutionLocation = gl.getUniformLocation(program, 'resolution');
    var iTexSizeLocation = gl.getUniformLocation(program, 'iTexSize');
    var iBlockOffsetLocation = gl.getUniformLocation(program, 'iBlockOffset');
    var iSampleRateLocation = gl.getUniformLocation(program, 'iSampleRate');

    webglUtils.resizeCanvasToDisplaySize(gl.canvas);
    gl.viewport(0, 0, gl.canvas.width, gl.canvas.height);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.useProgram(program);

//    gl.uniform2f(resolutionUniformLocation, texSize, texSize);
    gl.uniform1f(iTexSizeLocation, texSize);
    gl.uniform1f(iSampleRateLocation, sampleRate);

    gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);

    gl.vertexAttribPointer(vertexAttributeLocation, 2, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(vertexAttributeLocation);
    gl.drawArrays(gl.TRIANGLES, 0, 6);
    //gl.disableVertexAttribArray(vertexAttributeLocation);

    var frameBuffer = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, frameBuffer);
    gl.pixelStorei(gl.PACK_ALIGNMENT, 4);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);

    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, texSize, texSize, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST)
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);

    var blockSize = texSize ** 2;
    var nBlocks = Math.ceil(duration * sampleRate / blockSize);

    console.log(duration, blockSize, nBlocks);

    //gl.viewport(0, 0, texSize, texSize);

    var musicBlock = new Uint8Array(4 * blockSize);
    var music = new Uint8Array(4 * blockSize * nBlocks);
    for (let i = 0; i < nBlocks; i++) {
        console.log("block" + i);
        gl.uniform1f(iBlockOffsetLocation, i)

        gl.drawArrays(gl.TRIANGLES, 0, 6);
        gl.flush();
        gl.readPixels(0, 0, texSize, texSize, gl.RGBA, gl.UNSIGNED_BYTE, music, 4 * blockSize * i);
        //for (let k = 0; k < 4 * blockSize; k++) {
        //    music[4*i + k] = musicBlock[k];
        //}
    }
    // k, Probleme:
    // entweder es is zu kurz oder es repeatet zu schnell
    // mit texSize größer überlagert sich alle scheiße..??
    // vllt an viewport und resolution arbeiten, ich war froh genug dass irgendwas schonmal tat =o

    console.log(duration, sampleRate);
    var buffer = audioContext.createBuffer(2, duration * sampleRate, sampleRate);

    var channelL = buffer.getChannelData(0);
    var channelR = buffer.getChannelData(1);
    for(let i = 0; i < buffer.length / 4; i++) {
        var LL = music[4*i];
        var LH = music[4*i+1];
        var RL = music[4*i+2];
        var RH = music[4*i+3];
        var scale = 1 << 16;
        channelL[i] = ((LH << 8) + LL - scale) / scale;
        channelR[i] = ((RH << 8) + RL - scale) / scale;
    }
    console.log(scale, channelR);
    playBuffer(buffer);
}

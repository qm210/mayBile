"use strict";

var audioContext = null;

initAudio();
function initAudio() {

    window.addEventListener('load', init, false);
    function init() {
        try {
            window.AudioContext = window.AudioContext || window.webkitAudioContext;
            audioContext = new AudioContext();
            console.log("WebAudio is supported :)");
            console.log(audioContext.destination);
        }
        catch (e) {
            alert('Web Audio API is not supported');
        }
    }
}

function testAudio() {
    var length = 20;
    var buffer = audioContext.createBuffer(2, length * sampleRate, sampleRate);

    var freq = 210;

    var channelL = buffer.getChannelData(0);
    var channelR = buffer.getChannelData(1);
    for(let i = 0; i < buffer.length; i++) {
        channelL[i] = Math.sin(2. * Math.PI * freq * i / sampleRate);
        channelR[i] = Math.sin(2. * Math.PI * freq * 0.99 * i / sampleRate);
    }
    playBuffer(buffer);
}

function playBuffer(buffer) {
    var source = audioContext.createBufferSource();
    console.log(buffer);
    source.buffer = buffer;
    source.connect(audioContext.destination);
    source.start(0);
}

function mainAudioOld() {
    // XML request for loading sound
    var request = new XMLHttpRequest();
    request.open('GET', './somedistortedstrings.ogg', true);
    request.responseType = "arraybuffer";

    request.onload = function() {
        audioContext.decodeAudioData(request.response, function(buffer) {
            console.log("XML request completed.")
            if (!buffer) {
                alert('error decoding file: ' + url);
                return;
            }
            playBuffer(buffer);
        }, function() {
            console.log("Error in XML request");
        });
    }
    request.send();
}
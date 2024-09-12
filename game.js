const GLOBAL_Z = 5;
const GRAVITY = 0.1;
let CELLSIZE = 32;
let GameDimR = 10;
let GameDimC = 10;
const ENABLECAMERA = false;
var songBgm = {songData: [{ i: [0, 0, 140, 0, 0, 0, 140, 0, 0, 255, 158, 158, 158, 0, 0, 0, 0, 51, 2, 1, 2, 58, 239, 0, 32, 88, 1, 157, 2 ],p: [1,1,1,1],c: [{n: [161,,,,,,,,,,,,,,,,163,,,,,,,,159],f: []}]},{ i: [0, 91, 128, 0, 0, 95, 128, 12, 0, 0, 12, 0, 72, 0, 0, 0, 0, 0, 0, 0, 2, 255, 0, 0, 32, 83, 3, 130, 4 ],p: [1,1,2,1],c: [{n: [144,,151,,149,,147,,146,,147,,146,,144,,144,,151,,149,,147,,146,,147,,146,,144],f: []},{n: [156,,163,,161,,159,,158,,159,,158,,156,,156,,163,,161,,159,,158,,159,,158,,168],f: []}]},{ i: [0, 16, 133, 0, 0, 28, 126, 12, 0, 0, 2, 0, 60, 0, 0, 0, 0, 0, 0, 0, 2, 91, 0, 0, 32, 47, 3, 157, 2 ],p: [1,2,1,2],c: [{n: [144,,151,,149,,147,,146,,147,,146,,144,,144,,151,,149,,147,,146,,147,,146,,144],f: []},{n: [168,,175,,173,,171,,170,,171,,170,,168,,168,,175,,173,,171,,170,,171,,170,,168],f: []}]},{ i: [0, 255, 116, 79, 0, 255, 116, 0, 83, 0, 4, 6, 69, 52, 0, 0, 0, 0, 0, 0, 2, 14, 0, 0, 32, 0, 0, 0, 0 ],p: [1,1,1,1],c: [{n: [144,,151,,149,,147,,146,,147,,146,,144,,144,,151,,149,,147,,146,,147,,146,,144,,,159,,,,159,,,,159,,,,,,,,,,,,159,,159],f: []}]},],rowLen: 8269,   patternLen: 32,  endPattern: 3,  numChannels: 4  };
// var songBgm2 = {songData: [{ i: [2, 100, 128, 0, 3, 201, 128, 0, 0, 0, 5, 6, 58, 0, 0, 0, 0, 195, 6, 1, 2, 135, 0, 0, 32, 147, 6, 121, 6 ],p: [2,2],c: [{n: [155,,160,,155,,160,,162,,160,,162,,163,,162,,160,,158,,155,,160,,155,,160,,155,,160,,158,,157,,155,,,,160,,,,160,,,,,,,,,,,,,,,,,,160,,,,160],f: []},{n: [157,,162,,157,,162,,164,,162,,164,,165,,164,,162,,160,,157,,162,,157,,162,,157,,162,,160,,159,,157,,,,162,,,,162,,,,,,,,,,,,,,,,,,162,,,,162],f: []}]},{ i: [0, 214, 104, 64, 0, 204, 104, 0, 64, 229, 4, 40, 43, 51, 0, 0, 0, 231, 6, 1, 3, 183, 15, 0, 32, 232, 4, 74, 6 ],p: [2,1],c: [{n: [144,,,,,,,,,,,,146,146,,,,,,,135,,,,135,135,,,,,,,135],f: []},{n: [157,,162,,157,,162,,,,,,,,,,,,,,,,,,162,,157,,162,,,,,,,,,,,,,,,,,,162],f: []}]},{ i: [0, 255, 106, 64, 0, 255, 106, 0, 64, 0, 5, 7, 164, 0, 0, 0, 0, 0, 0, 0, 2, 255, 0, 2, 32, 83, 5, 25, 1 ],p: [2,1],c: [{n: [143,,,,,,,,,,,,,,,,,,,,,,143],f: []},{n: [135,,,,,,,,135,,,,,,,,135,,,,,,,,135,,,,,,,,135],f: []}]},{ i: [0, 214, 104, 64, 0, 204, 104, 0, 64, 229, 4, 40, 43, 51, 0, 0, 0, 231, 6, 1, 3, 183, 15, 0, 32, 232, 4, 74, 6 ],p: [1,1],c: [{n: [,,,,,,,,176,,174,,176,,177,,176,,174,,172,,169,,,,,,,,169,,174,,172,,171,,169],f: []}]},{ i: [1, 255, 128, 0, 1, 154, 128, 9, 0, 0, 7, 5, 52, 0, 0, 0, 0, 0, 0, 0, 2, 255, 0, 0, 32, 47, 3, 146, 2 ],p: [2,2],c: [{n: [,,,,,,,,170,,,,,,,,166,,164,,,,,,,,,,170,,,,166,,164,,163],f: []},{n: [,,162,,,,162,,,,,,,,165,,,,,,,,,,162,,,,162,,,,162,,,,159,,,,,,162,,,,162,,,,,,,,,,,,,,,,,,162,,,,162],f: []}]},],rowLen: 5513,   patternLen: 40,  endPattern: 1,  numChannels: 5  };
class CPlayer {
    constructor() {
        this.mOscillators = [
            this.osc_sin,
            this.osc_square,
            this.osc_saw,
            this.osc_tri
        ];
        this.mSong = null;
        this.mLastRow = 0;
        this.mCurrentCol = 0;
        this.mNumWords = 0;
        this.mMixBuf = null;
    }
    osc_sin(value) {
        return Math.sin(value * 6.283184);
    }
    osc_saw(value) {
        return 2 * (value % 1) - 1;
    }
    osc_square(value) {
        return (value % 1) < 0.5 ? 1 : -1;
    }
    osc_tri(value) {
        const v2 = (value % 1) * 4;
        if (v2 < 2) return v2 - 1;
        return 3 - v2;
    }
    getnotefreq(n) {
        return 0.003959503758 * (2 ** ((n - 128) / 12));
    }
    createNote(instr, n, rowLen) {
        const osc1 = this.mOscillators[instr.i[0]];
        const o1vol = instr.i[1];
        const o1xenv = instr.i[3] / 32;
        const osc2 = this.mOscillators[instr.i[4]];
        const o2vol = instr.i[5];
        const o2xenv = instr.i[8] / 32;
        const noiseVol = instr.i[9];
        const attack = instr.i[10] * instr.i[10] * 4;
        const sustain = instr.i[11] * instr.i[11] * 4;
        const release = instr.i[12] * instr.i[12] * 4;
        const releaseInv = 1 / release;
        const expDecay = -instr.i[13] / 16;
        let arp = instr.i[14];
        const arpInterval = rowLen * (2 ** (2 - instr.i[15]));
        const noteBuf = new Int32Array(attack + sustain + release);
        let c1 = 0, c2 = 0;
        let o1t = 0;
        let o2t = 0;
        for (let j = 0, j2 = 0; j < attack + sustain + release; j++, j2++) {
            if (j2 >= 0) {
                arp = (arp >> 8) | ((arp & 255) << 4);
                j2 -= arpInterval;
                o1t = this.getnotefreq(n + (arp & 15) + instr.i[2] - 128);
                o2t = this.getnotefreq(n + (arp & 15) + instr.i[6] - 128) * (1 + 0.0008 * instr.i[7]);
            }
            let e = 1;
            if (j < attack) {
                e = j / attack;
            } else if (j >= attack + sustain) {
                e = (j - attack - sustain) * releaseInv;
                e = (1 - e) * (3 ** (expDecay * e));
            }
            c1 += o1t * e ** o1xenv;
            let rsample = osc1(c1) * o1vol;
            c2 += o2t * e ** o2xenv;
            rsample += osc2(c2) * o2vol;
            if (noiseVol) {
                rsample += (2 * Math.random() - 1) * noiseVol;
            }
            noteBuf[j] = (80 * rsample * e) | 0;
        }
        return noteBuf;
    }
    initGenBuffer(song,context,callback){
        this.init(song);
        var loop = ()=>{
            var done = this.generate();
            if(done == 1){
                var buffer = this.createAudioBuffer(context);
                return callback(buffer);
            }
            else{
                requestAnimationFrame(loop);
            }
        }
        requestAnimationFrame(loop);
    }
    init(song) {
        this.mSong = song;
        this.mLastRow = song.endPattern;
        this.mCurrentCol = 0;
        this.mNumWords = song.rowLen * song.patternLen * (this.mLastRow + 1) * 2;
        this.mMixBuf = new Int32Array(this.mNumWords);
    }
    generate() {
        let i, j, b, p, row, col, n, cp, k, t, lfor, e, x, rsample, rowStartSample, f, da;
        const chnBuf = new Int32Array(this.mNumWords);
        const instr = this.mSong.songData[this.mCurrentCol];
        const rowLen = this.mSong.rowLen;
        const patternLen = this.mSong.patternLen;
        let low = 0, band = 0, high;
        let lsample, filterActive = false;
        const noteCache = [];
        for (p = 0; p <= this.mLastRow; ++p) {
            cp = instr.p[p];
            for (row = 0; row < patternLen; ++row) {
                const cmdNo = cp ? instr.c[cp - 1].f[row] : 0;
                if (cmdNo) {
                    instr.i[cmdNo - 1] = instr.c[cp - 1].f[row + patternLen] || 0;
                    if (cmdNo < 17) {
                        noteCache.length = 0;
                    }
                }
                const oscLFO = this.mOscillators[instr.i[16]];
                const lfoAmt = instr.i[17] / 512;
                const lfoFreq = (2 ** (instr.i[18] - 9)) / rowLen;
                const fxLFO = instr.i[19];
                const fxFilter = instr.i[20];
                const fxFreq = instr.i[21] * 43.23529 * 3.141592 / 44100;
                const q = 1 - instr.i[22] / 255;
                const dist = instr.i[23] * 1e-5;
                const drive = instr.i[24] / 32;
                const panAmt = instr.i[25] / 512;
                const panFreq = 6.283184 * (2 ** (instr.i[26] - 9)) / rowLen;
                const dlyAmt = instr.i[27] / 255;
                const dly = instr.i[28] * rowLen & ~1;  
                rowStartSample = (p * patternLen + row) * rowLen;
                for (col = 0; col < 4; ++col) {
                    n = cp ? instr.c[cp - 1].n[row + col * patternLen] : 0;
                    if (n) {
                        if (!noteCache[n]) {
                            noteCache[n] = this.createNote(instr, n, rowLen);
                        }
                        const noteBuf = noteCache[n];
                        for (j = 0, i = rowStartSample * 2; j < noteBuf.length; j++, i += 2) {
                          chnBuf[i] += noteBuf[j];
                        }
                    }
                }
                for (j = 0; j < rowLen; j++) {
                    k = (rowStartSample + j) * 2;
                    rsample = chnBuf[k];
                    if (rsample || filterActive) {
                        f = fxFreq;
                        if (fxLFO) {
                            f *= oscLFO(lfoFreq * k) * lfoAmt + 0.5;
                        }
                        f = 1.5 * Math.sin(f);
                        low += f * band;
                        high = q * (rsample - band) - low;
                        band += f * high;
                        rsample = fxFilter == 3 ? band : fxFilter == 1 ? high : low;
                        if (dist) {
                            rsample *= dist;
                            rsample = rsample < 1 ? rsample > -1 ? this.osc_sin(rsample * .25) : -1 : 1;
                            rsample /= dist;
                        }
                        rsample *= drive;
                        filterActive = rsample * rsample > 1e-5;
                        t = Math.sin(panFreq * k) * panAmt + 0.5;
                        lsample = rsample * (1 - t);
                        rsample *= t;
                    } else {
                        lsample = 0;
                    }
                    if (k >= dly) {
                        lsample += chnBuf[k - dly + 1] * dlyAmt;
                        rsample += chnBuf[k - dly] * dlyAmt;
                    }
                    chnBuf[k] = lsample | 0;
                    chnBuf[k + 1] = rsample | 0;
                    this.mMixBuf[k] += lsample | 0;
                    this.mMixBuf[k + 1] += rsample | 0;
                }
            }
        }
        this.mCurrentCol++;
        return this.mCurrentCol / this.mSong.numChannels;
    }
    createAudioBuffer(context) {
        const buffer = context.createBuffer(2, this.mNumWords / 2, 44100);
        for (let i = 0; i < 2; i++) {
            const data = buffer.getChannelData(i);
            for (let j = i; j < this.mNumWords; j += 2) {
                data[j >> 1] = this.mMixBuf[j] / 65536;
            }
        }
        return buffer;
    }
    createWave() {
        const headerLen = 44;
        const l1 = headerLen + this.mNumWords * 2 - 8;
        const l2 = l1 - 36;
        const wave = new Uint8Array(headerLen + this.mNumWords * 2);
        wave.set([
            82, 73, 70, 70, 
            l1 & 255, (l1 >> 8) & 255, (l1 >> 16) & 255, (l1 >> 24) & 255,
            87, 65, 86, 69, 
            102, 109, 116, 32, 
            16, 0, 0, 0, 
            1, 0, 
            2, 0, 
            68, 172, 0, 0, 
            16, 177, 2, 0, 
            4, 0, 
            16, 0, 
            100, 97, 116, 97, 
            l2 & 255, (l2 >> 8) & 255, (l2 >> 16) & 255, (l2 >> 24) & 255
        ]);
        for (let i = 0, idx = headerLen; i < this.mNumWords; ++i) {
            let y = this.mMixBuf[i];
            y = y < -32767 ? -32767 : (y > 32767 ? 32767 : y);
            wave[idx++] = y & 255;
            wave[idx++] = (y >> 8) & 255;
        }
        return wave;
    }
    getData(t, n) {
        const i = 2 * Math.floor(t * 44100);
        const d = new Array(n);
        for (let j = 0; j < 2 * n; j += 1) {
            const k = i + j;
            d[j] = t > 0 && k < this.mMixBuf.length ? this.mMixBuf[k] / 32768 : 0;
        }
        return d;
    }
}
class SoundSystem{
    constructor(autostart = true){
        this.audioContext = new (window.AudioContext || window.webkitAudioContext)();
        this.audioContextSingleFire = new (window.AudioContext || window.webkitAudioContext)();
        this.buffer1 = this.generateShootingSound();
        this.buffer2 = this.generateExplosion();
        var cplayer = new CPlayer();
        var cplayer2 = new CPlayer();
        this.bgmTime = 0;
        this.pausedTime = 0;
        this.startTime = 0;
        cplayer.initGenBuffer(songBgm, this.audioContext,(buffer)=>{
            this.bgmBuffer = buffer;
            if(autostart) this.startBgm();
        });
        // cplayer2.initGenBuffer(songBgm2, this.audioContext,(buffer)=>{
        //     this.bgm2Buffer = buffer;
        // });
    }
    generateShootingSound() {
        const sampleRate = this.audioContext.sampleRate;
        const duration = 0.3; 
        const buffer = this.audioContext.createBuffer(1, sampleRate * duration, sampleRate);
        const data = buffer.getChannelData(0);
        for (let i = 0; i < data.length; i++) {
            
            data[i] = (Math.random() - 0.5) * 2;
        }
        const attackTime = 0.01; 
        const decayTime = 0.1;  
        const sustainLevel = 0.2; 
        const releaseTime = duration - attackTime - decayTime; 
        for (let i = 0; i < data.length; i++) {
            let time = i / sampleRate;
            if (time < attackTime) {
                data[i] *= time / attackTime; 
            } else if (time < attackTime + decayTime) {
                data[i] *= 1 - (time - attackTime) / decayTime * (1 - sustainLevel); 
            } else if (time > duration - releaseTime) {
                data[i] *= (duration - time) / releaseTime; 
            }
        }
        for (let i = 0; i < data.length; i++) {
            let time = i / sampleRate;
            
            data[i] *= Math.sin(2 * Math.PI * time * (440 + Math.random() * 100)); 
        }
        return buffer;
    }
    generateSound() {
        const sampleRate = this.audioContext.sampleRate;
        const duration = 0.01; 
        const frequency = 10; 
        const buffer = this.audioContext.createBuffer(1, sampleRate * duration, sampleRate);
        const data = buffer.getChannelData(0);
        for (let i = 0; i < data.length; i++) {
          data[i] = Math.sin(2 * Math.PI * frequency * i / sampleRate);
        }
        return buffer;
    }
    generateExplosion() {
        const sampleRate = this.audioContext.sampleRate;
        const duration = 0.5; 
        const buffer = this.audioContext.createBuffer(1, sampleRate * duration, sampleRate);
        const data = buffer.getChannelData(0);
        for (let i = 0; i < data.length; i++) {
            data[i] = Math.random() * 2 - 1; 
        }
        const attackTime = 0.05; 
        const decayTime = 0.2; 
        const sustainLevel = 0.0; 
        const releaseTime = duration - attackTime - decayTime; 
        for (let i = 0; i < data.length; i++) {
            let time = i / sampleRate;
            if (time < attackTime) {
                data[i] *= time / attackTime; 
            } else if (time < attackTime + decayTime) {
                data[i] *= 1 - (time - attackTime) / decayTime * (1 - sustainLevel); 
            } else if (time > duration - releaseTime) {
                data[i] *= (duration - time) / releaseTime; 
            }
        }
        return buffer;
    }
    playS1(){
        const source = this.audioContextSingleFire.createBufferSource();
        source.buffer = this.buffer1;
        source.connect(this.audioContextSingleFire.destination);
        source.start();
    }
    playS2(){
        const source = this.audioContextSingleFire.createBufferSource();
        source.buffer = this.buffer2;
        source.connect(this.audioContextSingleFire.destination);
        source.start();
    }
    startBgm(id = 1){
        if(this.bgmsource){
            this.bgmsource.stop();
            this.bgmsource = null;
        }
        if(this.bgmBuffer){
            this.bgmsource = this.audioContext.createBufferSource();
            this.bgmsource.buffer = id==1 ? this.bgmBuffer : this.bgm2Buffer;
            this.bgmsource.connect(this.audioContext.destination);
            this.bgmsource.loop = true;
            this.bgmsource.start(0, this.pausedTime);
            this.startTime = this.audioContext.currentTime - this.pausedTime;
        }
    }
    stopBgm(id){
        if(this.bgmsource){
            this.pausedTime = this.audioContext.currentTime - this.startTime;
            this.bgmsource.stop();
            this.bgmsource = null;
        }
    }
}
class Webgl2 {
    constructor() {
        this.models = {};
        this.terminate = false;
        this.addModels();
    }
    reset(canvas) {
        let t;
        this.canvas = canvas;
        this.objs = 0;
        this.current = {};
        this.next = {};
        this.textures = {};
        if(canvas == null ) return;
        this.gl = canvas.getContext('webgl2');
        this.gl.blendFunc(770, 771);
        this.gl.activeTexture(33984);
        this.program = this.gl.createProgram();
        this.gl.enable(2884);
        this.gl.shaderSource(t = this.gl.createShader(35633), `#version 300 es
        precision highp float;                        
        in vec4 pos, col, uv, normal;                 
        uniform mat4 pv, eye, m, im;                  
        uniform vec4 bb;                              
        out vec4 v_pos, v_col, v_uv, v_normal;        
        void main() {                                 
        gl_Position = pv * (                        
        v_pos = bb.z > 0.                         
        ? m[3] + eye * (pos * bb)                 
        : m * pos                                 
        );                                          
        v_col = col;                                
        v_uv = uv;
        v_normal = transpose(inverse(m)) * normal;  
        }`);
        this.gl.compileShader(t);
        this.gl.attachShader(this.program, t);
        this.gl.shaderSource(t = this.gl.createShader(35632), `#version 300 es
        precision highp float;                  
        in vec4 v_pos, v_col, v_uv, v_normal;   
        uniform vec3 light;                     
        uniform vec4 o;                         
        uniform sampler2D sampler;              
        out vec4 c;                             
        void main() {
        c = mix(texture(sampler, v_uv.xy), v_col, o[3]);  
        if(o[1] > 0.){                                    
        c = vec4(                                       
            c.rgb * (max(0., dot(light, -normalize(       
            o[0] > 0.                                   
            ? vec3(v_normal.xyz)                        
            : cross(dFdx(v_pos.xyz), dFdy(v_pos.xyz))   
            )))
            + o[2]),                                      
            c.a                                           
        );
        }
        }`);
        this.gl.compileShader(t);
        this.gl.attachShader(this.program, t);
        this.gl.linkProgram(this.program);
        this.gl.useProgram(this.program);
        this.gl.clearColor(1, 1, 1, 1);
        this.clearColor("fff");
        this.gl.enable(2929);
        this.light({
            'y': -1
        });
        this.camera({
            'fov': 30
        });
        setTimeout( () => {
            this.draw();
        }
        , 16);
    }
    clearColor(color) {
        var c = this.col(color);
        this.gl.clearColor(...c);
    }
    getSettingsEmpty(_) {
        var s = {};
        s['w'] = 1;
        s['h'] = 1;
        s['d'] = 1;
        s['x'] = 0;
        s['y'] = 0;
        s['z'] = 0;
        s['rx'] = 0;
        s['ry'] = 0;
        s['rz'] = 0;
        s['b'] = '888';
        s['mode'] = 4;
        s['mix'] = 0;
        return s;
    }
    setState(state, type, texture, i, normal=[], A, B, C, Ai, Bi, Ci, AB, BC) {
        state[`n`] ||= 'o' + this.objs++;
        if (state[`size`])
            state[`w`] = state[`h`] = state[`d`] = state[`size`];
        if (state[`t`] && state[`t`][`width`] && !this.textures[state[`t`][`id`]]) {
            texture = this.gl.createTexture();
            this.gl.pixelStorei(37441, true);
            this.gl.bindTexture(3553, texture);
            this.gl.pixelStorei(37440, 1);
            this.gl.texImage2D(3553, 0, 6408, 6408, 5121, state[`t`]);
            this.gl.generateMipmap(3553);
            this.textures[state[`t`][`id`]] = texture;
        }
        if (state[`fov`]) {
            this.projection = new DOMMatrix([(1 / Math.tan(state[`fov`] * Math.PI / 180)) / (this.canvas.width / this.canvas.height), 0, 0, 0, 0, (1 / Math.tan(state[`fov`] * Math.PI / 180)), 0, 0, 0, 0, -1001 / 999, -1, 0, 0, -2002 / 999, 0]);
        }
        var s = this.getSettingsEmpty();
        state = {
            type,
            ...(this.current[state[`n`]] = this.next[state[`n`]] || s),
            ...state,
            'f': 0
        };
        if (this.models[state[`type`]] && this.models[state[`type`]]['vertices'] && 
                !this.models?.[state[`type`]]['verticesBuffer']) {
            this.gl.bindBuffer(34962, this.models[state[`type`]]['verticesBuffer'] = this.gl.createBuffer());
            this.gl.bufferData(34962, new Float32Array(this.models[state[`type`]]['vertices']), 35044);
            if (!this.models[state[`type`]][`normals`] && this.smooth)
                this.smooth(state);
            if (this.models[state[`type`]][`normals`]) {
                this.gl.bindBuffer(34962, this.models[state[`type`]][`normalsBuffer`] = this.gl.createBuffer());
                this.gl.bufferData(34962, new Float32Array(this.models[state[`type`]][`normals`].flat()), 35044);
            }
        }
        if ((this.models[state[`type`]] && this.models[state[`type`]][`uv`]) && !this.models[state[`type`]]['uvBuffer']) {
            this.gl.bindBuffer(34962, this.models[state[`type`]][`uvBuffer`] = this.gl.createBuffer());
            this.gl.bufferData(34962, new Float32Array(this.models[state[`type`]][`uv`]), 35044);
        }
        if ((this.models[state[`type`]] && this.models[state[`type`]][`indices`]) && !this.models[state[`type`]][`indicesBuffer`]) {
            this.models[state[`type`]][`indicesBuffer`] = this.gl.createBuffer();
            this.gl.bindBuffer(34963, this.models[state[`type`]][`indicesBuffer`]);
            this.gl.bufferData(34963, new Uint16Array(this.models[state[`type`]][`indices`]), 35044);
        }
        if (!state[`t`]) {
            state[`mix`] = 1;
        } else if (state[`t`] && !state[`mix`]) {
            state[`mix`] = 0;
        }
        this.next[state[`n`]] = state;
    }
    draw(now, dt, v, i, transparent=[]) {
        if(this.terminate) return;
        if(this.canvas == null) return;
        dt = now - this.lastFrame;
        this.lastFrame = now;
        requestAnimationFrame( (...t) => {
            this.draw(...t)
        }
        );
        if (this.next['camera']['g']) {
            this.render(this.next[this.next['camera']['g']], dt, 1);
        }
        v = this.animation('camera');
        if (this.next['camera']?.['g']) {
            v.preMultiplySelf(this.next[this.next['camera']['g']][`M`] || this.next[this.next['camera']['g']][`m`]);
        }
        this.gl.uniformMatrix4fv(this.gl.getUniformLocation(this.program, 'eye'), false, v.toFloat32Array());
        v['invertSelf']();
        v['preMultiplySelf'](this.projection);
        this.gl.uniformMatrix4fv(this.gl.getUniformLocation(this.program, 'pv'), false, v.toFloat32Array());
        this.gl.clear(16640);
        for (i in this.next) {
            if (this.next[i]['type'] == undefined)
                continue;
            if (this.next[i][`b`] == undefined) {
                debugger ;
            }
            if (!this.next[i][`t`] && this.col(this.next[i][`b`])[3] == 1) {
                if (this.next[i]['type'] == undefined) {
                    debugger ;
                }
                this.render(this.next[i], dt);
            } else {
                transparent.push(this.next[i]);
            }
        }
        transparent.sort( (a, b) => {
            return this.dist(b) - this.dist(a);
        }
        );
        this.gl.enable(3042);
        for (i of transparent) {
            if (["plane", "billboard"].includes(i[`type`]))
                this.gl.depthMask(0);
            this.render(i, dt);
            this.gl.depthMask(1);
        }
        this.gl.disable(3042);
        this.gl.uniform3f(this.gl.getUniformLocation(this.program, 'light'), this.lerp('light', 'x'), this.lerp('light', 'y'), this.lerp('light', 'z'));
    }
    render(object, dt, just_compute=['camera', 'light', 'group'].includes(object[`type`]), buffer) {

        if (object[`t`]) {
            this.gl.bindTexture(3553, this.textures[object[`t`][`id`]]);
            this.gl.uniform1i(this.gl.getUniformLocation(this.program, 'sampler'), 0);
        }
        if (object[`f`] < object[`a`])
            object[`f`] += dt;
        if (object[`f`] > object[`a`])
            object[`f`] = object[`a`];
        this.next[object[`n`]][`m`] = this.animation(object['n']);
        if (this.next[object[`g`]]) {
            this.next[object[`n`]][`m`].preMultiplySelf(this.next[object[`g`]][`M`] || this.next[object[`g`]][`m`]);
        }
        this.gl.uniformMatrix4fv(this.gl.getUniformLocation(this.program, 'm'), false, (this.next[object[`n`]][`M`] || this.next[object[`n`]][`m`]).toFloat32Array());
        try {
            this.gl.uniformMatrix4fv(this.gl.getUniformLocation(this.program, 'im'), false, (new DOMMatrix(this.next[object[`n`]][`M`] || this.next[object[`n`]][`m`]))['invertSelf']().toFloat32Array());
        } catch (e) {
            console.log(e,this.next[object[`n`]]);
        }
        if (!just_compute) {
            if (this.models[object[`type`]] == undefined) {
                debugger ;
            }
            this.gl.bindBuffer(34962, this.models[object[`type`]][`verticesBuffer`]);
            this.gl.vertexAttribPointer(buffer = this.gl.getAttribLocation(this.program, 'pos'), 3, 5126, false, 0, 0)
            this.gl.enableVertexAttribArray(buffer);
            if (this.models[object[`type`]][`uvBuffer`]) {
                this.gl.bindBuffer(34962, this.models[object[`type`]][`uvBuffer`]);
                this.gl.vertexAttribPointer(buffer = this.gl.getAttribLocation(this.program, 'uv'), 2, 5126, false, 0, 0);
                this.gl.enableVertexAttribArray(buffer);
            }
            if ((object[`s`] || this.models[object[`type`]][`customNormals`]) && this.models[object[`type`]][`normalsBuffer`]) {
                this.gl.bindBuffer(34962, this.models[object[`type`]][`normalsBuffer`]);
                this.gl.vertexAttribPointer(buffer = this.gl.getAttribLocation(this.program, 'normal'), 3, 5126, false, 0, 0);
                this.gl.enableVertexAttribArray(buffer);
            }
            this.gl.uniform4f(this.gl.getUniformLocation(this.program, 'o'), object[`s`], ((object[`mode`] > 3) || (this.gl[object[`mode`]] > 3)) && !object[`ns`] ? 1 : 0, this.ambientLight || 0.2, object[`mix`]);
            this.gl.uniform4f(this.gl.getUniformLocation(this.program, 'bb'), object[`w`], object[`h`], object[`type`] == 'billboard', 0);
            if (this.models[object[`type`]][`indicesBuffer`]) {
                this.gl.bindBuffer(34963, this.models[object['type']][`indicesBuffer`]);
            }
            this.gl.vertexAttrib4fv(this.gl.getAttribLocation(this.program, 'col'), this.col(object['b']));
            if (this.models[object[`type`]][`indicesBuffer`]) {
                this.gl.drawElements(+object[`mode`] || this.gl[object[`mode`]], this.models[object[`type`]][`indices`].length, 5123, 0);
            } else {
                this.gl.drawArrays(+object[`mode`] || this.gl[object[`mode`]], 0, this.models[object[`type`]][`vertices`].length / 3);
            }
        }
    }
    lerp(item, property) {
        try {
            if (this.next[item] && this.next[item][`a`]) {
                return this.current[item][property] + (this.next[item][property] - this.current[item][property]) * (this.next[item][`f`] / this.next[item][`a`]);
            } else {
                return this.next[item][property];
            }
        } catch (e) {
            return 0;
        }
    }
    animation(item, m=new DOMMatrix) {
        var xitem = this.next[item];
        if (!xitem)
            return m;
        var x = this.lerp(item, 'x');
        var y = this.lerp(item, 'y');
        var z = this.lerp(item, 'z');
        var rx = this.lerp(item, 'rx');
        var ry = this.lerp(item, 'ry');
        var rz = this.lerp(item, 'rz');
        var w = this.lerp(item, 'w');
        var h = this.lerp(item, 'h');
        var d = this.lerp(item, 'd');
        var m2 = m.translateSelf(x, y, z).rotateSelf(rx, ry, rz).scaleSelf(w, h, d);
        return m2;
    }
    dist(a, b=null) {
        if (b == null)
            b = this.next.camera;
        return (a && a[`m`]) && (b && b[`m`]) ? (b[`m`][`m41`] - a[`m`][`m41`]) ** 2 + (b[`m`][`m42`] - a[`m`][`m42`]) ** 2 + (b[`m`][`m43`] - a[`m`][`m43`]) ** 2 : 0
    }
    ambient(a) {
        this.ambientLight = a
    }
    col(c) {
        return [...c.replace("#", "").match(c.length < 5 ? /./g : /../g).map(a => ('0x' + a) / (c.length < 5 ? 15 : 255)), 1]
    }
    add(name, objects) {
        this.models[name] = objects;
        if (objects['normals']) {
            this.models[name]['customNormals'] = 1;
        }
        this[name] = (settings) => {
            this.setState(settings, name);
        }
    }
    group(t) {
        this.setState(t, 'group')
    }
    move(t, delay) {
        return setTimeout( () => {
            this.setState(t)
        }
        , delay || 1);
    }
    delete(t, delay) {
        return setTimeout( () => {
            delete this.next[t]
        }
        , delay || 1);
    }
    camera(t, delay) {
        return setTimeout( () => {
            this.setState(t, t[`n`] = 'camera')
        }
        , delay || 1);
    }
    light(t, delay) {
        return delay ? setTimeout( () => {
            this.setState(t, t[`n`] = 'light')
        }
        , delay) : this.setState(t, t[`n`] = 'light');
    }
    smooth(state, dict={}, vertices=[], iterate, iterateSwitch, i, j, A, B, C, Ai, Bi, Ci, normal) {
        this.models[state[`type`]][`normals`] = [];
        for (i = 0; i < this.models[state[`type`]][`vertices`].length; i += 3) {
            vertices.push(this.models[state[`type`]][`vertices`].slice(i, i + 3));
        }
        if (iterate = this.models[state[`type`]][`indices`])
            iterateSwitch = 1;
        else
            iterate = vertices,
            iterateSwitch = 0;
        for (i = 0; i < iterate.length * 2; i += 3) {
            j = i % iterate.length;
            A = vertices[Ai = iterateSwitch ? this.models[state[`type`]][`indices`][j] : j];
            B = vertices[Bi = iterateSwitch ? this.models[state[`type`]][`indices`][j + 1] : j + 1];
            C = vertices[Ci = iterateSwitch ? this.models[state[`type`]][`indices`][j + 2] : j + 2];
            let AB = [B[0] - A[0], B[1] - A[1], B[2] - A[2]];
            let BC = [C[0] - B[0], C[1] - B[1], C[2] - B[2]];
            normal = i > j ? [0, 0, 0] : [AB[1] * BC[2] - AB[2] * BC[1], AB[2] * BC[0] - AB[0] * BC[2], AB[0] * BC[1] - AB[1] * BC[0]];
            dict[A[0] + "_" + A[1] + "_" + A[2]] ||= [0, 0, 0];
            dict[B[0] + "_" + B[1] + "_" + B[2]] ||= [0, 0, 0];
            dict[C[0] + "_" + C[1] + "_" + C[2]] ||= [0, 0, 0];
            this.models[state[`type`]][`normals`][Ai] = dict[A[0] + "_" + A[1] + "_" + A[2]] = dict[A[0] + "_" + A[1] + "_" + A[2]].map( (a, i) => a + normal[i]);
            this.models[state[`type`]][`normals`][Bi] = dict[B[0] + "_" + B[1] + "_" + B[2]] = dict[B[0] + "_" + B[1] + "_" + B[2]].map( (a, i) => a + normal[i]);
            this.models[state[`type`]][`normals`][Ci] = dict[C[0] + "_" + C[1] + "_" + C[2]] = dict[C[0] + "_" + C[1] + "_" + C[2]].map( (a, i) => a + normal[i]);
        }
    }
    addModels() {
        var plane_settings = {};
        var pyramid_settings = {};
        var cube_settings = {};
        plane_settings['vertices'] = [.5, .5, 0, -.5, .5, 0, -.5, -.5, 0, .5, .5, 0, -.5, -.5, 0, .5, -.5, 0];
        plane_settings['uv'] = [1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0];
        cube_settings['vertices'] = [.5, .5, .5, -.5, .5, .5, -.5, -.5, .5, .5, .5, .5, -.5, -.5, .5, .5, -.5, .5, .5, .5, -.5, .5, .5, .5, .5, -.5, .5, .5, .5, -.5, .5, -.5, .5, .5, -.5, -.5, .5, .5, -.5, -.5, .5, -.5, -.5, .5, .5, .5, .5, -.5, -.5, .5, .5, .5, .5, .5, -.5, .5, .5, -.5, .5, -.5, -.5, -.5, -.5, -.5, .5, .5, -.5, -.5, -.5, -.5, -.5, .5, -.5, .5, -.5, .5, .5, -.5, .5, -.5, -.5, -.5, .5, -.5, .5, -.5, -.5, -.5, -.5, -.5, .5, -.5, .5, -.5, -.5, .5, -.5, -.5, -.5, .5, -.5, .5, -.5, -.5, -.5, .5, -.5, -.5];
        cube_settings['uv'] = [1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0];
        pyramid_settings['vertices'] = [-.5, -.5, .5, .5, -.5, .5, 0, .5, 0, .5, -.5, .5, .5, -.5, -.5, 0, .5, 0, .5, -.5, -.5, -.5, -.5, -.5, 0, .5, 0, -.5, -.5, -.5, -.5, -.5, .5, 0, .5, 0, .5, -.5, .5, -.5, -.5, .5, -.5, -.5, -.5, .5, -.5, .5, -.5, -.5, -.5, .5, -.5, -.5];
        pyramid_settings['uv'] = [0, 0, 1, 0, .5, 1, 0, 0, 1, 0, .5, 1, 0, 0, 1, 0, .5, 1, 0, 0, 1, 0, .5, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0];
        this.add("plane", plane_settings);
        this.add("billboard", this.models['plane']);
        this.add("cube", cube_settings);
        this.add("pyramid", pyramid_settings);
    }
    deleteAll() {
        this.terminate = true;
        var keys = Object.keys(this.current).filter(x => x!= 'camera' && x != 'light');
        for(let i = 0 ; i < keys.length;i++){
            var n = keys[i];
            this.delete(n,0);
        }
        var numTextureUnits = this.gl.getParameter(this.gl.MAX_TEXTURE_IMAGE_UNITS);
        for (var unit = 0; unit < numTextureUnits; ++unit) {
            this.gl.activeTexture(this.gl.TEXTURE0 + unit);
            this.gl.bindTexture(this.gl.TEXTURE_2D, null);
            this.gl.bindTexture(this.gl.TEXTURE_CUBE_MAP, null);
        }
        this.gl.bindBuffer(this.gl.ARRAY_BUFFER, null);
        this.gl.bindBuffer(this.gl.ELEMENT_ARRAY_BUFFER, null);
        this.gl.bindRenderbuffer(this.gl.RENDERBUFFER, null);
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
        for(let i in this.textures){
            this.gl.deleteTexture(this.textures[i]);
        }
    }
}
class MazeGenerator {
    constructor(rows, cols) {
        this.rows = rows;
        this.cols = cols;
        this.grid = new Array(rows).fill(null).map( () => new Array(cols).fill(true));
        this.generateMaze();
        for (var i = 0; i < this.rows; i++) {
            this.grid[i][cols - 1] = false;
        }
    }
    invertGrid() {
        var grid2 = new Array(this.rows).fill(null).map( () => new Array(this.cols).fill(false));
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.cols; j++) {
                grid2[i][j] = !this.grid[i][j];
            }
        }
        return grid2;
    }
    generateMaze() {
        this.clearMaze();
        this.carvePassage(0, 0);
        return this.grid;
    }
    clearMaze() {
        for (let row = 0; row < this.rows; row++) {
            for (let col = 0; col < this.cols; col++) {
                this.grid[row][col] = true;
            }
        }
    }
    carvePassage(row, col) {
        this.grid[row][col] = false;
        const directions = this.shuffleDirections();
        for (const direction of directions) {
            const newRow = row + direction[0];
            const newCol = col + direction[1];
            if (this.isValidCell(newRow, newCol) && this.grid[newRow][newCol]) {
                const betweenRow = row + direction[0] / 2;
                const betweenCol = col + direction[1] / 2;
                this.grid[betweenRow][betweenCol] = false;
                this.carvePassage(newRow, newCol);
            }
        }
    }
    isValidCell(row, col) {
        return row >= 0 && row < this.rows && col >= 0 && col < this.cols;
    }
    shuffleDirections() {
        const directions = [[-2, 0], [2, 0], [0, -2], [0, 2]];
        for (let i = directions.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [directions[i],directions[j]] = [directions[j], directions[i]];
        }
        return directions;
    }
}
class G {
    static makeCanvas(w=0, h=0, gl=false) {
        let c = document.createElement('canvas');
        c.width = w;
        c.height = h;
        c.w = w;
        c.h = h;
        if (gl) {// c.ctx = c.getContext('webgl');
        } else {
            c.ctx = c.getContext('2d');
        }
        c.center = {
            x: w / 2,
            y: h / 2
        }
        c.clear = () => {
            c.ctx.clearRect(0, 0, w, h);
        }
        c.fill = (color) => {
            c.ctx.fillStyle = color;
            c.ctx.fillRect(0, 0, w, h);
        }
        c.fillPatern = (img) => {
            const pattern = c.ctx.createPattern(img, "repeat");
            c.ctx.fillStyle = pattern;
            c.ctx.fillRect(0, 0, w, h);
        }
        c.id = Math.random();
        return c;
    }
    static GenTable(rows, cols) {
        var html = ``;
        for (let i = 0; i < rows; i++) {
            html += `<tr>`;
            for (let j = 0; j < cols; j++) {
                html += `<td></td>`;
            }
            html += `</tr>`;
        }
        var table = document.createElement('table');
        table.innerHTML = html;
        var entities = [];
        var trs = table.querySelectorAll('tr');
        for (let i = 0; i < trs.length; i++) {
            var tds = trs[i].querySelectorAll('td');
            entities[i] = [...tds];
        }
        table.entities = entities;
        return table;
    }
    static Point(pos) {
        return new Point(pos);
    }
    static getEmojiSprite(emoji, size, factor=1.3, color='#000') {
        let canvas = G.makeCanvas(size, size);
        var ctx = canvas.ctx;
        ctx.font = `${size / factor}px sans-serif`;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillStyle = color;
        ctx.fillText(emoji, size / 2, size * 1.1 / 2);
        canvas.id = Math.random();
        return canvas;
    }
    static getTextSprite(text, size, color, factor=0.8) {
        text = text.toUpperCase();
        let canvas = G.makeCanvas(size * text.length, size);
        for (let i = 0; i < text.length; i++) {
            var ls = G.getEmojiSprite(text[i], size, factor, color);
            canvas.ctx.drawImage(ls, i * size, 0);
        }
        return canvas;

    }
    static fuseImage(canvas, canvas2, composite='source-atop') {
        let buffer = G.makeCanvas(canvas.width, canvas.height);
        let ctx = buffer.ctx;
        ctx.drawImage(canvas, 0, 0);
        ctx.globalCompositeOperation = composite;
        for (let i = 0; i < canvas.width / canvas2.width; i++) {
            for (let j = 0; j < canvas.height / canvas2.height; j++) {
                ctx.drawImage(canvas2, i * canvas2.width, j * canvas2.height);
            }
        }
        return buffer;
    }
    static rotateCanvas(_image, deg) {
        var image = (deg % 90 != 0) ? G.prepForRotate(_image) : _image;
        var canvas = G.makeCanvas(image.width, image.height);
        var ctx = canvas.ctx;
        ctx.save();
        ctx.translate(canvas.width / 2, canvas.height / 2);
        ctx.rotate(deg * Math.PI / 180);
        ctx.drawImage(image, -image.width / 2, -image.height / 2);
        ctx.restore();
        return canvas;
    }
    static prepForRotate(image) {
        let d = Math.sqrt(Math.pow(image.width, 2) + Math.pow(image.height, 2));
        let buffer = G.makeCanvas(d, d);
        buffer.ctx.drawImage(image, (d - image.width) / 2, (d - image.height) / 2);
        return buffer;
    }
    static mirror(canvas, hor=true) {
        let buffer = G.makeCanvas(canvas.width, canvas.height);
        let context = buffer.ctx;
        context.save();
        if (hor) {
            context.scale(-1, 1);
            context.drawImage(canvas, 0, 0, canvas.width * -1, canvas.height);
        } else {
            context.scale(1, -1);
            context.drawImage(canvas, 0, 0, canvas.width, canvas.height * -1);
        }
        context.restore();
        return buffer;
    }
    static gridBG(color1="lightgrey", color2=null, scale=8, width=1) {
        var canvas = G.makeCanvas(scale, scale);
        var ctx = canvas.ctx;
        ctx.fillStyle = color1;
        ctx.fillRect(0, 0, scale, scale);
        if (color2 == null) {
            ctx.clearRect(0, 0, scale - width, scale - width);
        } else {
            ctx.fillStyle = color2;
            ctx.fillRect(0, 0, scale - width, scale - width);
        }
        return canvas;
    }
    static Lightify(canvas, opacity) {
        let buffer = G.makeCanvas(canvas.width, canvas.height);
        buffer.ctx.globalAlpha = opacity;
        buffer.ctx.drawImage(canvas, 0, 0);
        buffer.ctx.globalAlpha = 1;
        return buffer;
    }
    static makeDom(html) {
        var h = document.createElement('div');
        h.innerHTML = html;
        return h.firstChild;
    }
    static shuffleArray(array) {
        for (let i = array.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [array[i],array[j]] = [array[j], array[i]];
        }
        return array;
    }
    static repeatCanvas(canvas, r, c=0) {
        if (c == 0)
            c = r;
        var buffer = G.makeCanvas(canvas.width * c, canvas.height * r);
        for (let i = 0; i < r; i++) {
            for (let j = 0; j < c; j++) {
                buffer.ctx.drawImage(canvas, j * canvas.width, i * canvas.height);
            }
        }
        return buffer;
    }
    static merge(list, w, h) {
        var c = G.makeCanvas(w, h);
        for (let i in list) {
            c.ctx.drawImage(list[i], 0, 0);
        }
        return c;
    }
    static brickPattern(color1="#fff", color2="#000", r=1) {
        var canvas = G.makeCanvas(8, 8);
        var ctx = canvas.ctx;
        ctx.fillStyle = color1;
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = color2;
        ctx.fillRect(7, 0, 1, 4);
        ctx.fillRect(0, 3, 8, 1);
        ctx.fillRect(4, 4, 1, 4);
        ctx.fillRect(0, 7, 8, 1);
        if (r > 1) {
            return G.repeatCanvas(canvas, r, r);
        }
        return canvas;
    }
    static tilePattern(color1, color2, r) {
        var canvas = G.makeCanvas(8, 8);
        canvas.fill(color1);
        canvas.ctx.fillStyle = color2;
        canvas.ctx.fillRect(1, 1, 6, 6);
        canvas.ctx.fillStyle = color1;
        canvas.ctx.fillRect(2, 2, 4, 4);
        canvas.ctx.fillStyle = color2;
        canvas.ctx.fillRect(3, 3, 2, 2);

        return G.repeatCanvas(canvas, r, r);
    }
    static randomPattern(color1, color2, bias=0.3, w=8, h=8) {
        var canvas = G.makeCanvas(w, h);
        var ctx = canvas.ctx;
        ctx.fillStyle = color1;
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = color2;
        for (let i = 0; i < h; i++) {
            for (let j = 0; j < w; j++) {
                if (Math.random() < bias)
                    ctx.fillRect(i, j, 1, 1);
            }
        }
        return canvas;
    }
    static MakeCircle(r, stroke=null, fill=null) {
        var s = G.makeCanvas(r * 2 + 2, r * 2 + 2);
        var ctx = s.ctx;
        ctx.beginPath();
        ctx.arc(s.width / 2, s.height / 2, r, 0, Math.PI * 2, false);
        if (stroke != null) {
            ctx.strokeStyle = stroke;
            ctx.stroke();
        }
        if (fill != null) {
            ctx.fillStyle = fill;
            ctx.fill();
        }
        return s;
    }
    static movePointToward(pos, rotation, distance) {
        const rRad = rotation * (Math.PI / 180);
        const vx = distance * Math.cos(rRad);
        const vy = distance * Math.sin(rRad);
        return {
            x: pos.x + vx,
            y: pos.y + vy
        }
    }
    static loadImage(url, callback) {
        var img = new Image();
        img.src = url;
        img.addEventListener('load', () => {
            callback(img);
        }
        );
    }
    static getColor(r, g, b, a) {
        if (r + g + b + a == 0) {
            return null;
        } else if (r + g + b == 0) {
            return '#000000';
        } else if (r > 255 || g > 255 || b > 255) {
            return '#000000';
        }
        return '#' + ((r << 16) + (g << 8) + b).toString(16).padStart(6, '0');
    }
    static getColorMatrix(canvas, changefct) {
        var context = canvas.getContext('2d');
        var width = canvas.width;
        var height = canvas.height;
        var imageData = context.getImageData(0, 0, width, height);
        var data = imageData.data;
        var colorMatrix = [];
        for (var i = 0; i < data.length; i += 4) {
            colorMatrix.push(G.getColor(data[i], data[i + 1], data[i + 2], data[i + 3]));
        }
        var matrix = [];
        for (let i = 0; i < canvas.height; i++) {
            matrix[i] = [];
        }
        let c = 0
          , r = 0;
        for (let i = 0; i < colorMatrix.length; i++) {
            if (c >= canvas.width) {
                r++;
                c = 0
            }
            matrix[r][c] = colorMatrix[i];
            if (changefct)
                matrix[r][c] = changefct(matrix[r][c]);
            c++;
        }
        return matrix;
    }
    static imgToCanvas(img) {
        var c = G.makeCanvas(img.width, img.height);
        c.ctx.drawImage(img, 0, 0);
        return c;
    }
    static colorsMatrixToSprite(matrix, scale=1, deform=null) {
        let height = matrix.length;
        let width = Math.max(...matrix.map( (row) => row.length));
        var buffer = G.makeCanvas(width * scale, height * scale);
        var ctx = buffer.ctx;
        for (let i = 0; i < height; i++) {
            for (let j = 0; j < width; j++) {
                var color = matrix[i][j];
                if (deform)
                    color = deform(color);
                if (!color || color == '')
                    continue;
                ctx.fillStyle = color;
                ctx.fillRect(j * scale, i * scale, scale, scale);
            }
        }
        return buffer;
    }
    static crop(canvas, x, y, width, height) {
        let buffer = G.makeCanvas(width, height);
        buffer.ctx.drawImage(canvas, x, y, width, height, 0, 0, width, height);
        return buffer;
    }
    static randomColor() {
        var letters = "0123456789ABCDEF";
        var color = "#";
        for (var i = 0; i < 6; i++) {
            color += letters[Math.floor(Math.random() * 16)];
        }
        return color;
    }
    static rand(a=1, b=0) {
        return b + (a - b) * Math.random();
    }
    static randInt(a=1, b=0) {
        return G.rand(a, b) | 0;
    }
}
class Point {
    constructor(pos) {
        this.x = pos.x;
        this.y = pos.y;
        this.z = pos.z;
    }
    moveToward(p2, dist=1) {
        var vx = this.x == p2.x ? 0 : this.x < p2.x ? dist : -dist;
        var vy = this.y == p2.y ? 0 : this.y < p2.y ? dist : -dist;
        var vz = this.z == p2.z ? 0 : this.z < p2.z ? dist : -dist;
        this.x += vx;
        this.y += vy;
        this.z += vz;
    }
    distance(p2) {
        let distance = 0;
        distance += Math.pow((this.x - p2.x), 2);
        distance += Math.pow((this.y - p2.y), 2);
        distance += Math.pow((this.z - p2.z), 2);
        distance = Math.sqrt(distance);
        return distance;
    }
    //only xy
    getAngleTo(target) {
        let dx = target.x - this.x;
        let dy = target.y - this.y;
        let angleRadians = Math.atan2(dy, dx);
        return angleRadians * 180 / Math.PI;
    }
    moveByAngle(rotation, distance) {
        const rRad = rotation * (Math.PI / 180);
        const vx = distance * Math.cos(rRad);
        const vy = distance * Math.sin(rRad);
        this.x = this.x + vx;
        this.y = this.y + vy;
    }
}
class SpriteEngine {
    constructor(img) {
        var imgCanvas = G.imgToCanvas(img);
        var mat = G.getColorMatrix(imgCanvas, (r) => {
            if (r == '')
                return null;
            if (r == '#fff')
                return null;
            if (r == '#ffffff')
                return null;
            return r;
        }
        );
        var MULT = CELLSIZE / 16;
        var cvs = G.colorsMatrixToSprite(mat, MULT);

        this.sprite0 = G.crop(cvs, (MULT * 16 * 0), 0, MULT * 16, MULT * 16);
        this.sprite1 = G.crop(cvs, (MULT * 16 * 1), 0, MULT * 16, MULT * 16);
        this.sprite2 = G.crop(cvs, (MULT * 16 * 2), 0, MULT * 16, MULT * 16);
        this.sprite3 = G.crop(cvs, (MULT * 16 * 3), 0, MULT * 16, MULT * 16);
    }
}
class Camera {
    constructor(game) {
        // this.game = game;
        this.n = "camera";
        this.x = 0;
        this.y = 2;
        this.z = 0;

        this.rx = -20;
        this.ry = 0;
        this.rz = 0;

        this.fov = 25;

        this.vz = 0;
        this.vx = 0;
        this.vy = 0;
    }
    update() {

        //use rotation to update camera pos

        const rad = this.ry * Math.PI / 180
        // Adjust the movement based on the rotation
        const cos = Math.cos(rad);
        const sin = Math.sin(rad);
        // Update z and x based on the rotated direction
        this.z += this.vz * cos - this.vx * sin;
        this.x += this.vz * sin + this.vx * cos;

        //forward backward
        // var rotatedX = this.rx;
        // this.z += this.vz;
        // this.x += this.vx;

        this.y += this.vy;

        // this.z = (this.z > 0 ? -1 : 1) * Math.ceil(Math.abs(this.z));
        // this.x = (this.x > 0 ? -1 : 1) * Math.ceil(Math.abs(this.x));
        // this.y = (this.y > 0 ? -1 : 1) * Math.ceil(Math.abs(this.y));

        // this.y -= GRAVITY;
        if (this.y < 0)
            this.y = 0;
        if (this.y > 100)
            this.y = 100;

    }
    getTable() {
        var attrib = "n,fov,x,y,z,rx,ry,rz,vx,vy,vz".split(',');
        var t = G.GenTable(2, attrib.length);
        t.classList.add('table-view');
        var e = t.entities;
        for (let i = 0; i < attrib.length; i++) {
            e[0][i].innerHTML = `${attrib[i]}`;
            e[1][i].innerHTML = `${isNaN(this[attrib[i]]) ? this[attrib[i]] : this[attrib[i]].toFixed(2)}`;
        }
        return t;
    }
    getAttrib(n) {
        switch (n) {
        case 'x':
            return this.x;
        case 'y':
            return this.y;
        case 'z':
            return this.z;
        case 'rx':
            return this.rx;
        case 'ry':
            return this.ry;
        case 'rz':
            return this.rz;
        case 'size':
            return this.size;
        case 'fov':
            return this.fov;
        }
        return 0;
    }
    genS() {
        var s2 = {};
        s2[`x`] = this.getAttrib('x');
        s2[`y`] = this.getAttrib('y');
        s2[`z`] = this.getAttrib('z');
        s2[`rx`] = this.getAttrib('rx');
        s2[`ry`] = this.getAttrib('ry');
        s2[`rz`] = this.getAttrib('rz');
        s2[`size`] = this.getAttrib('size');
        s2[`fov`] = this.getAttrib('fov');
        return s2;
    }
}
class Player {
    constructor(game,pos) {
        this.game = game;

        this.r = 0;
        this.n = 'player';
        this.b = '#0ff';

        this.x = pos.x;
        this.y = pos.y;
        this.z = pos.z;

        this.w = 0.1;
        this.h = 0.1;
        this.d = 0.2;

        this.rx = 0;
        this.ry = 0;
        this.rz = 0;

        this.vz = 0;
        this.vx = 0;
        this.vy = 0;

        this.keys = {
            w: false,
            a: false,
            s: false,
            d: false,
            e: false,
            q: false,
        }
        this.pos = G.Point(pos);
        this.destination = G.Point(pos);
        this.speed = 0.4;
        this.or = 'U';
        this.enableEventListners();
    }
    getAttrib(n) {
        switch (n) {
        case 'n':
            return this.n;
        case 'x':
            return this.x;
        case 'y':
            return this.y;
        case 'z':
            return this.z;
        case 'rx':
            return this.rx;
        case 'ry':
            return this.ry;
        case 'rz':
            return this.rz;
        case 'w':
            return this.w;
        case 'h':
            return this.h;
        case 'd':
            return this.d;
        case 'b':
            return this.b;
        }
        return 0;
    }
    genS() {
        var s2 = {};
        s2[`n`] = this.getAttrib('n');
        s2[`x`] = this.getAttrib('x');
        s2[`y`] = this.getAttrib('y');
        s2[`z`] = this.getAttrib('z');
        s2[`rx`] = this.getAttrib('rx');
        s2[`ry`] = this.getAttrib('ry');
        s2[`rz`] = this.getAttrib('rz');
        s2[`w`] = this.getAttrib('w');
        s2[`h`] = this.getAttrib('h');
        s2[`d`] = this.getAttrib('d');
        s2[`b`] = this.getAttrib('b');
        return s2;
    }
    keyup(code) {
        switch (code.toLowerCase()) {
        case 'a':
            this.keys.a = false;
            break;
        case 'w':
            this.keys.w = false;
            break;
        case 's':
            this.keys.s = false;
            break;
        case 'd':
            this.keys.d = false;
            break;
        case 'e':
            this.keys.e = false;
            break;
        case 'q':
            this.keys.q = false;
            break;
        }
    }
    keydown(code) {
        console.log(`keydown`,code);
        switch (code.toLowerCase()) {
        case 'a':
            this.keys.a = true;
            this.keys.d = false;
            break;
        case 'w':
            this.keys.w = true;
            this.keys.s = false;
            break;
        case 's':
            this.keys.s = true;
            this.keys.w = false;
            break;
        case 'd':
            this.keys.d = true;
            this.keys.a = false;
            break;
        case 'q':
            this.keys.q = true;
            this.keys.e = false;
            break;
        case 'e':
            this.keys.e = true;
            this.keys.q = false;
            break;
        }
    }
    resetKeys() {
        this.keys = {
            w: false,
            a: false,
            s: false,
            d: false,
            e: false,
            q: false,
        }
    }
    getOrientation() {
        switch (this.ry) {
        case 0:
            return 'U';
        case 90:
        case -270:
            return 'R';
        case 180:
        case -180:
            return 'B';
        case 270:
        case -90:
            return 'L';
        }
    }
    update(t) {
        this.or = this.getOrientation();
        var d = this.destination.distance(this.pos);
        this.vx = 0;
        this.vz = 0;
        if (this.keys.d) this.vx = -1;
        if (this.keys.a) this.vx = 1;

        if (this.keys.e) this.vx = -1;
        if (this.keys.q) this.vx = 1;


        if (this.keys.w)this.vz = 1;
        if (this.keys.s)this.vz = -1;
        

        // var moveside = this.keys.q || this.keys.e;
        var moveside = this.keys.a || this.keys.d;

        if (this.vx != 0 || this.vz != 0) {}
        
        if (this.vx != 0 && ! moveside) {
            this.ry += this.vx < 0 ? -90 : 90;
            this.vx = 0;
            if (this.ry >= 360)
                this.ry = 0;
            if (this.ry <= -360)
                this.ry = 0;
            this.updateCamera();
            this.keys.d = this.keys.a = false;
            this.keys.e = this.keys.q = false;
        }

        if (d > this.speed) {
            this.pos.moveToward(this.destination, this.speed);
            this.x = this.pos.x;
            this.y = this.pos.y;
            this.z = this.pos.z;
            this.updateCamera();
        } else if (this.vz != 0 || this.vx != 0) {
            this.pos = G.Point(this.destination);
            var dest = G.Point(this.destination);
            
            if (this.ry == 0) {dest.z -= this.vz;} 
            else if (this.ry == 90 * 1 || this.ry == -90 * 3) {dest.x -= this.vz;} 
            else if (this.ry == 90 * 2 || this.ry == -90 * 2) {dest.z += this.vz;} 
            else if (this.ry == 90 * 3 || this.ry == -90 * 1) {dest.x += this.vz;}
            
            if(moveside){
                if (this.ry == 0) {dest.x -= this.vx;} 
                else if (this.ry == 90 * 2 || this.ry == -90 * 2) {dest.x += this.vx;} 
                else if (this.ry == 90 * 1 || this.ry == -90 * 3) {dest.z += this.vx;} 
                else if (this.ry == 90 * 3 || this.ry == -90 * 1) {dest.z -= this.vx;}
            }

            if (this.game.validLocation(dest)) {
                this.destination = G.Point(dest);
            }
        }

    }
    updateCamera() {

        var povPly = {
            y: 0.7,
            z: 0.8,
            x: 0.8
        }
        
        var or = this.getOrientation();
        this.game.camera.y = this.y + povPly.y;
        this.game.camera.rx = 0;
        var ox = 0;
        var oz = 0;
        if (or == 'U')
            oz = +povPly.z;
        if (or == 'B')
            oz = -povPly.z;
        if (or == 'R')
            ox = +povPly.x;
        if (or == 'L')
            ox = -povPly.x;
        if(this.game.config.birdview){
            this.game.camera.x = this.x;
            this.game.camera.z = this.z;
            this.game.camera.ry = this.ry;
            this.game.camera.ry = this.ry;
            this.game.camera.y = this.y + 12;
            this.game.camera.rx = -90;
        }
        else{
            this.game.camera.x = this.x + ox;
            this.game.camera.z = this.z + oz;
            this.game.camera.ry = this.ry;
            this.game.camera.rx = -20;
        }
        


    }
    getTable() {
        var attrib = "or,x,z,rx,ry,rz".split(',');
        var t = G.GenTable(2, attrib.length);
        t.classList.add('table-view');
        var e = t.entities;
        for (let i = 0; i < attrib.length; i++) {
            e[0][i].innerHTML = `${attrib[i]}`;
            e[1][i].innerHTML = `${isNaN(this[attrib[i]]) ? this[attrib[i]] : this[attrib[i]].toFixed(2)}`;
        }
        return t;
    }
    enableEventListners() {
        window.addEventListener('keydown', (e) => {
            if(!ENABLECAMERA) this.keydown(e.key);
            switch (e.key.toLowerCase()) {
            case 'arrowup':
                this.keydown('w');
                break;
            case 'arrowdown':
                this.keydown('s');
                break;
            case 'arrowright':
                this.keydown('d');
                break;
            case 'arrowleft':
                this.keydown('a');
                break;
            }
        }
        );
        window.addEventListener('keyup', (e) => {
            if(!ENABLECAMERA) this.keyup(e.key);
            switch (e.key.toLowerCase()) {
            case 'arrowup':
                this.keyup('w');
                break;
            case 'arrowdown':
                this.keyup('s');
                break;
            case 'arrowright':
                this.keyup('d');
                break;
            case 'arrowleft':
                this.keyup('a');
                break;
            }
        }
        );
    }
    initHandlers() {

        document.addEventListener('keydown', (e) => {
            switch (e.key.toLowerCase()) {
            case 'arrowup':
                this.player.vz = +1 * 1.1;
                break;
            case 'arrowdown':
                this.player.vz = -1 * 1.1;
                break;
            case 'arrowright':
                this.player.vx = +1 * 1.1;
                break;
            case 'arrowleft':
                this.player.vx = -1 * 1.1;
                break;
            }
        }
        )
        document.addEventListener('keyup', (e) => {
            switch (e.key.toLowerCase()) {
            case 'arrowup':
                this.player.vz = 0;
                break;
            case 'arrowdown':
                this.player.vz = 0;
                break;
            case 'arrowright':
                this.player.vx = 0;
                break;
            case 'arrowleft':
                this.player.vx = 0;
                break;
            }
        }
        )
    }
}
class Cookie{
    constructor(game,pos){
        this.game = game;
        this.pos = G.Point(pos);
        this.id = Math.random();
        this.rx = 0;
        this.ry = 0;
        this.rz = 0;
        this.size = 0.5;
        var cookieSprite = G.getEmojiSprite('🍪',128,1.1);
        this.sprite = cookieSprite;
        var sides = this.getSides();
        sides.forEach(s=>{
            this.game.W['plane'](s);
        })
    }
    getSides(){
        var f = {
            'n' : `C${this.id}_a`,
            'x' : this.pos.x,
            'y' : this.pos.y,
            'z' : this.pos.z + 0.01,
            'rx' :  this.rx,
            'ry':   this.ry,
            'rz':   this.rz,
            'size': this.size,
            't':    this.sprite
        }
        var b = {
            'n' : `C${this.id}_b`,
            'x' : this.pos.x,
            'y' : this.pos.y,
            'z' : this.pos.z - 0.01,
            'rx' :this.rx,
            'ry': this.ry-180,
            'rz': this.rz,
            'size':this.size,
            't': this.sprite
        }
        return [f,b];
    }
    update(t){
        this.ry += 15;
        if(this.ry == 360) this.ry = 0;
        
        try{
            if(this.pos.distance(this.game.player.pos) <= 1){
                this.game.cookies += 1;
                this.game.objects = this.game.objects.filter(x=>x!=this);
                this.game.W.delete(`C${this.id}_a`,0);
                this.game.W.delete(`C${this.id}_b`,0);
            }
        }
        catch(e){

        }
        var sides = this.getSides();
        sides.forEach(s=>{
            this.game.W.move(s);
        })

    }
}
class Game {
    constructor() {
        this.config = {
            music : false,
            sound : false,
            song : 1,
            controls:false,
            birdview:false,
        };
        GameDimR = Math.floor(window.innerHeight/CELLSIZE) - 2;
        GameDimC = Math.floor(window.innerWidth/CELLSIZE) - 1;
        document.body.innerHTML = ``;
        this.resetBody();
        this.preLoading();
        this.objects = [];
        this.mainScene();
    }
    mainScene(){
        this.gameover = true;
        this.gamePased = true;
        this.resetBody();
        var canvas = G.makeCanvas(GameDimC*CELLSIZE,GameDimR*CELLSIZE);
        canvas.fill('#000');
        this.getMainMenuBg(canvas);
        this.body.append(canvas);
        this.showMenu();
    }
    showMenu(){
        this.gamePased = true;
        if(this.dialog != null){this.dialog.remove();}
        this.dialog = Object.assign(document.createElement('div'), { className: 'menuDialog'});
        
        var navItems = [];
        if(this.gameover){
            navItems.push({html : '<button >New Game</button>', f:'newgame'});
        }
        else{
            navItems.push({html : '<button >Resume</button>', f:'resume'});
        }
        navItems.push(...[
            {html : '<button >Help</button>',   f:'help'},
            {html : `<button >Music ${this.config.music ? 'ON': 'OFF'}</button>`,   f:'music'},
            {html : `<button >Controls ${this.config.controls ? 'ON': 'OFF'}</button>`,   f:'controls'},
            {html : `<button >Bird View ${this.config.birdview ? 'ON': 'OFF'}</button>`,   f:'birdview'},
        ]);
        if(!this.gameover){
            navItems.push({html : '<button >Quit</button>',   f:'quit'},);
        }
        var nav = G.GenTable(navItems.length,1);
        for(let i in navItems){
            var dom = G.makeDom(navItems[i].html)
            dom.style.width = `${GameDimC*CELLSIZE * 0.9}px`;
            dom.style.fontSize = `24pt`;

            nav.entities[i][0].append(dom);
            nav.entities[i][0].onclick = ()=>{
                this.ApplyMenuItem(navItems[i].f);
            }
        }
        this.dialog.append(nav);
        this.body.append(this.dialog);
    }
    ApplyMenuItem(item){
        if(item == 'newgame'){
            this.gamePased = false;
            this.gameover = false;
            this.newGame();
        }
        else if(item == 'resume'){
            this.gamePased = false;
            this.dialog.remove();
            this.update(this.time);
        }
        else if(item == 'controls'){
            this.config.controls = !this.config.controls;
            this.prepFootercontrols();
            this.showMenu();
        }
        else if(item == 'birdview'){
            this.config.birdview = !this.config.controls;
            this.showMenu();
        }
        else if(item == 'music'){
            if(!this.SoundSystem){
                this.SoundSystem = new SoundSystem();
            }
            var currentval = this.config.music;
            if(currentval){
                this.SoundSystem.stopBgm();
            }
            else{
                this.SoundSystem.startBgm();
            }
            this.config.music = !this.config.music;
            this.dialog.remove();
            this.showMenu();
        }
        else if(item == 'quit'){
            if(document.webkitIsFullScreen) document.exitFullscreen();
            this.gamePased = true;
            this.gameover = true;
            this.dialog.remove();
            this.mainScene();
        }
        else if(item == `help`){
            this.gamePased = true;
            if(this.dialog != null){this.dialog.remove();}
            this.dialog = Object.assign(document.createElement('div'), { className: 'menuDialog'});
            this.dialog.style.width = `${GameDimC*CELLSIZE}px`;
            var h2 = `
                <div class="helpDiv">
                    <h2>Help</h2>
                    <p>You are trapped in a maze </p>
                    <p>Gather cookies (🍪) before time run out</p>
                    <p>You have 13 minutes</p>
                    <p><b>Controls</b></p>
                    <table style="width:200px">
                        <tr><td><b>W</b></td><td>move forward</td></tr> 
                        <tr><td><b>S</b></td><td>move backward</td></tr> 
                        <tr><td><b>A</b></td><td>lean left</td></tr> 
                        <tr><td><b>D</b></td><td>lean right</td></tr> 
                        <tr><td><b>Q</b></td><td>rotate left</td></tr> 
                        <tr><td><b>R</b></td><td>rotate right</td></tr> 
                    </table>
                    <p><b>Game By Mhmd Jawad ZD</b></p>
                </div>
            `;
            var mdom = G.makeDom('<button>Menu</button>');
            mdom.onclick = ()=>{
                this.gamePased = false;
                this.dialog.remove();
                this.showMenu();
                // this.update(this.time);
            }
            this.dialog.innerHTML += h2;
            var helpDiv = this.dialog.querySelector('.helpDiv');
            helpDiv.style['overflow-y'] = `auto`;
            // helpDiv.style.height = GameDimR*CELLSIZE * 0.8  + `px`;
            this.dialog.append(mdom);
            this.body.append(this.dialog);
        }
    }
    prepheader(){
        var headerTable = G.GenTable(2,6);
        headerTable.style.width = GameDimC * CELLSIZE + "px";
        var entities = headerTable.entities;
        this.leveldom = document.createElement('div');
        this.cookiedom = document.createElement('div');
        this.timedom = document.createElement('div');
        this.playerhealthdom = document.createElement('div');
        this.menuDom = document.createElement('div');
        entities[0][0].append(G.getEmojiSprite('🍪',32,1.4));
        entities[1][0].append(this.cookiedom);
        // entities[0][1].append(G.getEmojiSprite('💓',32,1.4));
        // entities[1][1].append(this.playerhealthdom);
        entities[0][2].append(G.getEmojiSprite('⌛',32,1.4));
        entities[1][2].append(this.timedom);
        entities[0][4].append(`Level`);
        entities[1][4].append(this.leveldom);
        entities[0][5].rowSpan = 2;
        entities[0][5].append(G.getEmojiSprite('📋',40,1.4));
        entities[1][5].remove();
        entities[0][5].onclick = ()=>{this.showMenu();}
        this.header.append(headerTable);
    }
    newGame(){
        this.resetBody();
        this.prepheader();
        this.prepFootercontrols();
        this.newLevel(1);
    }
    newLevel(level){
        this.gamePased = false;
        this.gameover = false;
        this.input = {};
        this.level = level;
        this.timeup = false;
        this.leveldom.innerHTML = this.level;
        this.cookies = 0;
        this.levelCookies = 13;
        this.timeremaining = 13*60;
        this.initMainCanvas();
        this.body.innerHTML = ``;
        this.body.append(this.canvas);
        this.initHandlers();
        this.objects = [];
        this.WItems = [];
        this.start();
    }
    prepFootercontrols(){
        if(this.config.controls == false){
            this.footer.innerHTML = '';
            return;
        }
        this.footer.innerHTML = '';
        var table = G.GenTable(2,3);
        table.classList.add('gamecontrolstable');
        table.style.width = GameDimC * CELLSIZE + "px";
        var entities = table.entities;
        var keys = [
            {html : '<span> <h1>w</h1> </span>', f : 'w' , r : 0 , c : 1},
            {html : '<span> <h1>s</h1> </span>', f : 's' , r : 1 , c : 1},
            {html : '<span> <h1>a</h1> </span>', f : 'a' , r : 1 , c : 0},
            {html : '<span> <h1>d</h1> </span>', f : 'd' , r : 1 , c : 2},
            {html : '<span> <h1>d</h1> </span>', f : 'q' , r : 0 , c : 0},
            {html : '<span> <h1>d</h1> </span>', f : 'e' , r : 0 , c : 2},
        ]
        keys.forEach(k=>{
            var dom = G.makeDom(k.html);
            entities[k.r][k.c].addEventListener('touchstart',(e)=> this.player.keydown(k.f));
            entities[k.r][k.c].addEventListener('touchend',(e)=> this.player.keyup(k.f));
            entities[k.r][k.c].addEventListener('mousedown',(e)=> this.player.keydown(k.f));
            entities[k.r][k.c].addEventListener('mouseup',(e)=> this.player.keyup(k.f));
            entities[k.r][k.c].append(dom) ;
            entities[k.r][k.c].style.border = '2px solid black';
            entities[k.r][k.c].style.background = 'blue';
            entities[k.r][k.c].style.color = '#fff';
        })
        // entities[0][2].rowSpan = 2;
        // entities[1][2].remove();
        // entities[0][0].rowSpan = 2;
        // entities[1][0].remove();
        this.footer.appendChild(table);
    }
    resetBody(){
        var div_w_class = `<div class='_class_'></div>`;
        this.layout = G.makeDom(div_w_class.replace('_class_','layout'));
        this.header = G.makeDom(div_w_class.replace('_class_','header'));
        this.body = G.makeDom(div_w_class.replace('_class_','body'));
        this.footer = G.makeDom(div_w_class.replace('_class_','footer'));
        this.layout.appendChild(this.header);
        this.layout.appendChild(this.body);
        this.layout.appendChild(this.footer);
        document.body.innerHTML = ``;
        document.body.appendChild(this.layout);
    }
    preLoading(){
        var about = G.makeDom(`<div>Loading....</div>`);
        this.body.append(about);
    }
    initMainCanvas() {
        if(this.W){
            // this.W.reset(null);
            this.W.deleteAll();
            this.W = null;
        }
        var s = Math.min(GameDimC*CELLSIZE,GameDimR*CELLSIZE);
        this.canvas = G.makeCanvas(s,s, true);
        this.W = new Webgl2();
        this.W.reset(this.canvas);
        this.W.light({x:0,y:0.5,z:0});
        this.W.ambient(0.5);
    }
    colorMatrixTo3D(matrix, scale=1, deform=null) {
        let height = matrix.length;
        let width = Math.max(...matrix.map( (row) => row.length));
        var buffer = G.makeCanvas(width * scale, height * scale);
        var ctx = buffer.ctx;
        this.W.group({
            n: 'mob',
            x: 0,
            y: 0,
            z: 0
        });
        for (let i = 0; i < height; i++) {
            for (let j = 0; j < width; j++) {
                var color = matrix[i][j];
                if (deform)
                    color = deform(color);
                if (!color || color == '')
                    continue;
                ctx.fillStyle = color;
                ctx.fillRect(j * scale, i * scale, scale, scale);
                var s2 = 0.05;
                var ox = j * s2;
                var oy = -i * s2;
                this.W.cube({
                    g: 'mob',
                    x: ox,
                    y: oy,
                    z: 0,
                    size: s2,
                    b: color
                })
            }
        }
        return buffer;
    }
    start() {
        this.objects = [];
        this.genMaze(this.W);
        var pos = G.Point({
            x: -this.mazeDim.r/2,
            y: 0,
            z: -this.mazeDim.c/2
        });
        this.player = new Player(this,pos);
        this.camera = new Camera(this);
        this.W[`cube`](this.player.genS());
        this.player.updateCamera();
        this.timestart = 0;
        this.update(0);
    }
    update(t) {
        if(this.gameover == true) return this.gameOverScene();
        if(this.gamePased == true){return;}
        var tins = parseInt(t/1000);
        var timeremaining = this.timeremaining - tins;
        this.timedom.innerHTML = `${this.parseTime(timeremaining)}`;
        if(timeremaining <=0){
            this.gameover == true;
            this.timeup == true;
            this.gameOverScene()
        }
        this.cookiedom.innerHTML = `${this.cookies}/${this.levelCookies}`;
        if(this.cookies >= this.levelCookies){
            this.LevelEndScene();
        }
        this.player.update();
        this.camera.update();
        this.W.camera(this.camera.genS());
        this.W.move(this.player.genS());
        this.objects.forEach(x=> x.update(t));

        // this.footer.innerHTML = `time ${t}`;
        // this.logDom.append(this.camera.getTable())
        // this.footer.append(this.player.getTable())
        // W.light(this.camera);
        setTimeout(() => {
            requestAnimationFrame( (time) => this.update(time));
        }, 16);
    }
    initHandlers() {
        if (!ENABLECAMERA)
            return;
        document.addEventListener('keydown', (e) => {
            switch (e.key.toLowerCase()) {
            case 'w':
                this.camera.vz = -1 * 0.1;
                break;
            case 's':
                this.camera.vz = +1 * 0.1;
                break;
            case 'a':
                this.camera.vx = -1 * 0.1;
                break;
            case 'd':
                this.camera.vx = +1 * 0.1;
                break;
            case 'c':
                this.camera.vy = -1 * 1.1;
                break;
            case ' ':
                this.camera.vy = +1 * 1.1;
                break;

                // case 'arrowup':     this.player.vz = +1 * 1.1 ;break;
                // case 'arrowdown':   this.player.vz = -1 * 1.1 ;break;
                // case 'arrowright':  this.player.vx = +1 * 1.1 ;break;
                // case 'arrowleft':   this.player.vx = -1 * 1.1 ;break;

            case 't':
                document.exitPointerLock();
                break;
            case 'r':
                this.camera = new Camera();
                break;

            }
        }
        )
        document.addEventListener('keyup', (e) => {
            switch (e.key.toLowerCase()) {
            case 'w':
                this.camera.vz = 0;
                break;
            case 's':
                this.camera.vz = 0;
                break;
            case 'a':
                this.camera.vx = 0;
                break;
            case 'd':
                this.camera.vx = 0;
                break;
            case 'c':
                this.camera.vy = 0;
                break;
            case ' ':
                this.camera.vy = 0;
                break;

                // case 'arrowup':     this.player.vz = 0 ;break;
                // case 'arrowdown':   this.player.vz = 0 ;break;
                // case 'arrowright':  this.player.vx = 0 ;break;
                // case 'arrowleft':   this.player.vx = 0 ;break;

            }
        }
        )
        this.canvas.onclick = (e) => {
            this.canvas.requestPointerLock();
        }
        document.addEventListener('wheel', (e) => {
            this.camera.fov += e.deltaY > 0 ? 1 : -1;
        }
        )
        document.addEventListener('mousemove', (e) => {
            if (document.pointerLockElement === this.canvas || document.mozPointerLockElement === this.canvas) {
                this.camera.ry -= (e.movementX / 20 > 0 ? 1 : -1) * 1.5;
                this.camera.rx -= (e.movementY / 20 > 0 ? 1 : -1) * 1.5;

                // if (this.camera.rx > 45)
                //     this.camera.rx = 45;
                // if (this.camera.rx < -45)
                //     this.camera.rx = -45;

                if (this.camera.ry >= 360)
                    this.camera.ry = 0;
                // if(this.camera.ry > 45) this.camera.ry = 45;
                // if(this.camera.ry < -45) this.camera.ry = -45;
            }
        }
        )
    }
    validLocation(pos) {
        try {
            if (pos.x < -this.mazeDim.r / 2)
                return false;
            if (pos.z < -this.mazeDim.c / 2)
                return false;
            if (pos.x >= this.mazeDim.r / 2)
                return false;
            if (pos.z >= this.mazeDim.c / 2)
                return false;
            var r = pos.z + this.mazeDim.r / 2;
            var c = pos.x + this.mazeDim.c / 2;
            var gridval = this.maze.grid[r][c] || false;
            return !gridval;
        } catch (e) {
            console.log(e);
            return false;
        }
    }
    getMiniMap(grid, rows, cols) {
        var cellSize = 2;
        var grass = G.randomPattern('#2d7d00', '#509e26', 0.1, cellSize, cellSize);
        var dirt = G.randomPattern('#924200', '#3d1c00', 0.01, cellSize, cellSize);
        var canvas = G.makeCanvas(cols * cellSize, rows * cellSize);
        var ctx = canvas.ctx;
        for (let i = 0; i < rows; i++) {
            var row = grid[i];
            for (let j = 0; j < cols; j++) {
                var s = row[j] == path ? grass : dirt;
                ctx.drawImage(s, j * s.width, i * s.height);
            }
        }
        return canvas;
    }
    genMaze() {
        this.mazeDim = {
            r: 20 + this.level*2,
            c: 20 + this.level*2,
        }
        var rows = this.mazeDim.r;
        var cols = this.mazeDim.c;
        var cellSize = 64;
        this.maze = new MazeGenerator(rows,cols);
        var color = G.randomColor();
        var grass = G.randomPattern('#00ab00', '#008000', 0.2, cellSize, cellSize);
        var dirt = G.randomPattern('#b85900', '#874100', 0.2, cellSize, cellSize);
        var tile1 = G.tilePattern(color, '#fff', cellSize / 8);
        for (let i = 0; i < this.mazeDim.r; i++) {
            for (let j = 0; j < this.mazeDim.c; j++) {
                if (this.maze.grid[i][j] == true) {
                    if (Math.random() < 0.3) {
                        this.maze.grid[i][j] = false;
                    }
                }
            }
        }
        var openareas = [];
        for (let i = 0; i < this.mazeDim.r; i++) {
            for (let j = 0; j < this.mazeDim.c; j++) {
                if (this.maze.grid[i][j] == true) {
                    // this.W['plane'](this.genS(j - cols/2, 0 , i-rows/2, -90,null,null, 1,null,null,null,null,dirt));
                    this.W['cube'](this.genS(j - cols / 2, 0.5, i - rows / 2, -90, null, null, 1, null, null, null, null, tile1));
                } else {
                    var ss = Math.random() > 0.5 ? grass : dirt;
                    this.W['plane'](this.genS(j - cols / 2, 0, i - rows / 2, -90, null, null, 1, null, null, null, null, ss));
                    // this.W['plane'](this.genS(j - cols / 2, 1, i - rows / 2, 90, null, null, 1, null, null, null, null, tile1));
                    openareas.push({
                        x : j - cols / 2,
                        y : 0.5,
                        z : i - rows / 2
                    });
                }
            }
        }
        for (let i = 0; i < rows + 2; i++) {
            this.W['cube'](this.genS(i - cols / 2 - 1, 0.5, -rows / 2 - 1, null, null, null, 1, null, null, null, null, tile1));
            this.W['cube'](this.genS(-rows / 2 - 1, 0.5, i - cols / 2 - 1, null, null, null, 1, null, null, null, null, tile1));
            this.W['cube'](this.genS(i - cols / 2 - 1, 0.5, rows / 2, null, null, null, 1, null, null, null, null, tile1));
            this.W['cube'](this.genS(rows / 2, 0.5, i - cols / 2 - 1, null, null, null, 1, null, null, null, null, tile1));
        }
        var locations = G.shuffleArray(openareas).splice(0,this.levelCookies);
        locations.forEach(loc=>{
            this.objects.push(
                new Cookie(this,loc)
            )
        })
    }
    genS(x=null, y=null, z=null, rx=null, ry=null, rz=null, size=null, w=null, h=null, d=null, b=null, t=null, fov=null) {
        var s2 = {};
        if (x != null && x != undefined)
            s2[`x`] = x;
        if (y != null && y != undefined)
            s2[`y`] = y;
        if (z != null && z != undefined)
            s2[`z`] = z;
        if (rx != null && rx != undefined)
            s2[`rx`] = rx;
        if (ry != null && ry != undefined)
            s2[`ry`] = ry;
        if (rz != null && rz != undefined)
            s2[`rz`] = rz;
        if (size != null && size != undefined)
            s2[`size`] = size;
        if (w != null && w != undefined)
            s2[`w`] = w;
        if (h != null && h != undefined)
            s2[`h`] = h;
        if (d != null && d != undefined)
            s2[`d`] = d;
        if (b != null && b != undefined)
            s2[`b`] = b;
        if (t != null && t != undefined)
            s2[`t`] = t;
        if (fov != null && fov != undefined)
            s2[`fov`] = fov;
        return s2;
    }
    parseTime(s){
        let m = Math.floor(s / 60);
        let h = Math.floor(m / 60);
        h = h == 0 ? '' : h < 10 ? `0${h}:` : `${h}:`;
        m = Math.floor(m % 60);
        m = m == 0 ? '' : m < 10 ? `0${m}:` : `${m}:`;
        s = Math.floor(s % 60);
        return `${h}${m}${s}`;
    }
    LevelEndScene(){
        this.gamePased = true;
        if(this.dialog != null){this.dialog.remove();}
        this.dialog = Object.assign(document.createElement('div'), { className: 'menuDialog'});
        this.dialog.innerHTML = `<h1>Well Done<h1><h2> Level ${this.level} finished</h2>`;
        var nextLevelButton = G.makeDom(`<button id="nextLevel"><h3>Next Level</h3></button>`);
        nextLevelButton.onclick = ()=>{
            this.newLevel(this.level+1);
        }
        this.dialog.append(nextLevelButton);
        this.body.append(this.dialog);
    }
    gameOverScene(){
        this.gamePased = true;
        this.gameover = true;
        if(this.dialog != null){this.dialog.remove();}
        this.dialog = Object.assign(document.createElement('div'), { className: 'menuDialog'});
        this.dialog.style.width = `${GameDimC*CELLSIZE}px`;
        this.dialog.style.height = `${GameDimR*CELLSIZE * 0.8}px`;
        this.dialog.innerHTML = `<h1>Game Over</h1>`;
        if(this.timeup){
            this.dialog.innerHTML += `<h2>Time Up</h2>`;
        }
        var button = G.makeDom(`<button id="nextLevel"><h2>New Game</h2></button>`);
        button.style.width = `${GameDimC*CELLSIZE * 0.90}px`;

        button.onclick = ()=>{
            this.newGame();
        }
        this.dialog.append(button);
        this.body.append(this.dialog);
    }
    getMainMenuBg(canvas){
        var ctx = canvas.ctx;
        var grass = G.randomPattern('#2d7d00','#509e26',0.1,CELLSIZE,CELLSIZE);
        var cookie = G.getEmojiSprite(`🍪`,CELLSIZE*2,1.1)
        var cookiedim = G.Lightify(cookie,0.3);
        var cookielocs = [];
        for(let i = 0 ; i < 50;i++){
            cookielocs.push(
                {x : G.randInt(cookiedim.w/2,canvas.w) , y : -G.randInt(CELLSIZE,CELLSIZE*4) ,}
            )
        }
        var slogan = G.getTextSprite(`The   MAZE `,   50, `#fff`, 1.5, '');
        var t3d = G.getTextSprite(`3d`,   50, `#f00`, 1.5, '');
        var credit = G.getTextSprite(`BY MHMDJAWADZD`,   50, `#fff`, 1.5, '');

        function upd(t){
            canvas.fill('#000');
            canvas.fillPatern(grass);
            // ctx.clearRect(0,0,canvas.w,canvas.h);
            ctx.fillStyle = '#fff';
            // ctx.fillText(t,20,20);
            ctx.drawImage(t3d, canvas.w/2 + slogan.w/2 - t3d.w,canvas.h/2-slogan.h);
            ctx.drawImage(cookie,canvas.w/2 - cookie.w,canvas.h/2-slogan.h/2);
            ctx.drawImage(slogan,canvas.w/2-slogan.w/2,canvas.h/2-slogan.h/2);
            ctx.drawImage(credit,canvas.w/2-credit.w/2,canvas.h-credit.h);
            
            cookielocs.forEach(c=>{
                ctx.drawImage(cookiedim,
                    c.x,c.y
                );
                c.y += G.randInt(1,CELLSIZE/4)
                if(c.y > canvas.h) c.y = -CELLSIZE;
                
            })

            if(t > 10000) return;
            requestAnimationFrame((t)=>upd(t));
            setTimeout(() => {
            }, 500);
        }
        upd(0);
    }
}
document.addEventListener('DOMContentLoaded', function() {
    window.game = new Game("body");
}, false);

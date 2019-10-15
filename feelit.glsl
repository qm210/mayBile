#version 300 es
precision highp float;
out vec4 fragColor;

uniform float iTexSize;
uniform float iBlockOffset;
uniform float iSampleRate;

#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0.,1e-3,clamp(x,0.,1e-3)); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sin_(float a, float p) { return sin(2. * PI * mod(a,1.) + p); }
float freqC1(float note){ return 32.7 * pow(2., note/12.); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;
float linmix(float x, float a, float b, float y0, float y1) { return mix(y0,y1,clamp(a*x+b,0.,1.)); }
float lofi(float time, float quant) { return floor(quant*time+.5)/quant; }

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

#define NTIME 2
const float pos_B[2] = float[2](0.,8.);
const float pos_t[2] = float[2](0.,13.7143);
const float pos_BPS[1] = float[1](.5833);
const float pos_SPB[1] = float[1](1.7144);
float BPS, SPB, BT;

float Tsample;

const float filterthreshold = 1e-3;

const float sequence_texture[816] = float[816](0.,2.,9.,10.,6.,35.,12.,.89990234375,.7998046875,.300048828125,.1500244140625,1.,.7001953125,0.,0.,0.,.25,0.,.114990234375,0.,4.,1.,2.,3.,4.,5.,6.,7.,4.,4.,8.,2.,3.,4.,5.,6.,7.,8.,8.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-12.,0.,26.,52.,56.,60.,64.,68.,72.,76.,80.,105.,0.,.125,.25,.375,.625,.75,1.,1.125,1.25,1.375,1.625,1.75,1.875,2.,2.125,2.25,2.375,2.625,2.75,2.875,3.,3.25,3.375,3.5,3.625,3.75,0.,.125,.25,.375,.625,.75,1.,1.125,1.25,1.375,1.625,1.75,1.875,2.,2.125,2.25,2.375,2.625,2.75,2.875,3.,3.25,3.375,3.5,3.625,3.75,0.,.25,.5,.75,0.,.25,.5,.75,0.,.25,.5,.75,0.,.25,.5,.75,0.,.25,.5,.75,0.,.25,.5,.75,0.,.25,.5,.75,0.,.125,.25,.375,.5,.75,1.,1.125,1.25,1.375,1.5,1.75,2.125,2.25,2.375,2.5,2.625,2.75,3.,3.125,3.25,3.375,3.5,3.625,3.75,.125,.25,.375,.5,.75,1.,1.125,1.25,1.375,1.5,1.75,1.875,2.,2.125,2.25,2.375,2.5,2.75,2.875,3.,3.125,3.375,3.5,3.625,3.75,4.,.125,.25,.375,.5,.75,1.,1.125,1.25,1.375,1.5,1.75,1.875,2.,2.125,2.25,2.375,2.5,2.75,2.875,3.,3.125,3.375,3.5,3.625,3.75,4.,.0625,.375,.5625,.875,.0625,.375,.5625,.875,.0625,.375,.5625,.875,.0625,.375,.5625,.875,.0625,.375,.5625,.875,.0625,.375,.5625,.875,.0625,.375,.5625,.875,.125,.25,.375,.5,.75,1.,1.125,1.25,1.375,1.5,1.75,2.125,2.25,2.375,2.5,2.625,2.75,3.,3.125,3.25,3.375,3.5,3.625,3.75,4.,14.,26.,29.,14.,14.,17.,14.,26.,29.,14.,14.,17.,29.,14.,26.,29.,14.,26.,29.,14.,17.,17.,19.,16.,14.,12.,14.,26.,29.,14.,14.,17.,14.,26.,29.,14.,14.,17.,29.,14.,26.,29.,14.,26.,29.,14.,17.,17.,19.,16.,14.,12.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,45.,50.,53.,45.,48.,57.,45.,48.,52.,53.,52.,43.,50.,53.,45.,52.,55.,52.,45.,48.,41.,45.,43.,41.,40.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,.89990234375,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,-3.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.,0.,0.,0.,0.,0.,-3.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,0.,0.,-5.,0.,-1.,-2.,0.,1.,4.,0.,-1.,0.,-2.,-1.,1.,1.,-1.,3.,-2.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.);


float s_atan(float a) { return 2./PI * atan(a); }

float drop_phase(float time, float t1, float f0, float f1)
{
    float t = min(time, t1);
    float phi = f0*t + .5*(f1-f0)/t1*t*t;

    if(time > t1)
    {
        phi += f1 * (time - t1);
    }
    return phi;
}

float drop2_phase(float time, float t1, float t2, float f0, float f1, float f2)
{
    float t = min(time, t1);
    float phi = f0*t + .5*(f1-f0)/t1*t*t;

    if(time > t1)
    {
        t = min(time, t2);
        phi += f1*t - t1*f2 + .5*(f2-f1)/(t2-t1)*(t*t - t1*t1);
    }
    if(time > t2)
    {
        phi += f2 * (time-t2);
    }
    return phi;
}

float lpnoise(float t, float fq)
{
    t *= fq;
    float tt = fract(t);
    float tn = t - tt;
    return mix(pseudorandom(floor(tn) / fq), pseudorandom(floor(tn + 1.0) / fq), smoothstep(0.0, 1.0, tt));
}

float env_AHDSR(float x, float L, float A, float H, float D, float S, float R)
{
    return (x<A ? x/A : x<A+H ? 1. : x<A+H+D ? (1. - (1.-S)*(x-H-A)/D) : x<=L-R ? S : x<=L ? S*(L-x)/R : 0.);
}

float waveshape(float s, float amt, float A, float B, float C, float D, float E)
{
    float w;
    float m = sign(s);
    s = abs(s);

    if(s<A) w = B * smoothstep(0.,A,s);
    else if(s<C) w = C + (B-C) * smoothstep(C,A,s);
    else if(s<=D) w = s;
    else if(s<=1.)
    {
        float _s = (s-D)/(1.-D);
        w = D + (E-D) * (1.5*_s*(1.-.33*_s*_s));
    }
    else return 1.;

    return m*mix(s,w,amt);
}

float sinshape(float x, float amt, float parts)
{
    return (1.-amt) * x + amt * sign(x) * 0.5 * (1. - cos(parts*PI*x));
}

float QFM_FB(float PH, float FB) // my guessing of feedback coefficients, FB>0 'saw', FB<0 'sq'
{
    if(FB > 0.) return abs(FB) * .8*_sin(PH + .35*_sin(PH));
    else return abs(FB) * _sin(PH + .5*PI);
}

float QFM(float t, float f, float phase, float LV1, float LV2, float LV3, float LV4, float FR1, float FR2, float FR3, float FR4, float FB1, float FB2, float FB3, float FB4, float ALGO)
{
    int iALGO = int(ALGO);
    float PH1 = FR1 * f * t + phase;
    float PH2 = FR2 * f * t + phase;
    float PH3 = FR3 * f * t + phase;
    float PH4 = FR4 * f * t + phase;

    float LINK41 = 0., LINK42 = 0., LINK43 = 0., LINK32 = 0., LINK31 = 0., LINK21 = 0.;
    if(iALGO == 1)       {LINK43 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 2)  {LINK42 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 3)  {LINK41 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 4)  {LINK42 = 1.; LINK43 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 5)  {LINK41 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 6)  {LINK43 = 1.; LINK32 = 1.;}
    else if(iALGO == 7)  {LINK43 = 1.; LINK32 = 1.; LINK31 = 1.;}
    else if(iALGO == 8)  {LINK21 = 1.; LINK43 = 1.;}
    else if(iALGO == 9)  {LINK43 = 1.; LINK42 = 1.; LINK41 = 1.;}
    else if(iALGO == 10) {LINK43 = 1.; LINK42 = 1.;}
    else if(iALGO == 11) {LINK43 = 1.;}

    float OP4 = LV4 * _sin(PH4 + QFM_FB(PH4, FB4));
    float OP3 = LV3 * _sin(PH3 + QFM_FB(PH3, FB3) + LINK43*OP4);
    float OP2 = LV2 * _sin(PH2 + QFM_FB(PH2, FB2) + LINK42*OP4 + LINK32*OP3);
    float OP1 = LV1 * _sin(PH1 + QFM_FB(PH1, FB1) + LINK41*OP4 + LINK31*OP3 + LINK32*OP2);

    float sum = OP1;
    if(LINK21 > 0.) sum += OP2;
    if(LINK31 + LINK32 > 0.) sum += OP3;
    if(LINK41 + LINK42 + LINK43 > 0.) sum += OP4;

    return s_atan(sum);
}

float bandpassBPsaw1(float time, float f, float tL, float vel, float fcenter, float bw, float M)
{
    float y = 0.;

    float facM = 2.*PI/M;
    float facL = 2.*PI*Tsample * (fcenter - bw);
    float facH = 2.*PI*Tsample * (fcenter + bw);

    if(facL < 0.) facL = 0.;
    if(facH > PI) facH = PI;

    float _TIME, mm, w, h;

    M--;
    for(float m=1.; m<=M; m++)
    {
        mm = m - .5*M;
        w = .42 - .5 * cos(mm*facM) - .08 * cos(2.*mm*facM);
        h = 1./(PI*mm) * (sin(mm*facH) - sin(mm*facL));

        _TIME = time - m*Tsample;
        y += w*h*(2.*fract(f*_TIME)-1.);
    }

    return s_atan(M*M*y); // I DO NOT CARE ANYMORE
}

float _SHKENV0(float t){return t <=.009? linmix(t,111.1111,0.,0.,1.):t <=.062? linmix(t,18.8679,-.1698,1.,.492):t <=.271? linmix(t,4.7847,-.2967,.492,.945):t <=.467? linmix(t,5.102,-1.3827,.945,.293):t <=.746? linmix(t,3.5842,-1.6738,.293,.508):t <=1.2? linmix(t,2.2026,-1.6432,.508,0.):0.;}
float _SHKENV1(float t){return t <=.029? linmix(t,34.4828,0.,0.,1.):t <=.046? linmix(t,58.8235,-1.7059,1.,.008):t <=1.146? linmix(t,.9091,-.0418,.008,0.):0.;}
float _SHKENV2(float t){return t <=.022? linmix(t,45.4545,0.,0.,1.):t <=.632? linmix(t,1.6393,-.0361,1.,.239):t <=1.2? linmix(t,1.7606,-1.1127,.239,0.):0.;}
float _SHKENV3(float t){return t <=.05? linmix(t,20.,0.,0.,1.):t <=.5? linmix(t,2.2222,-.1111,1.,.5):t <=1.2? linmix(t,1.4286,-.7143,.5,0.):0.;}
float _SHKENV5(float t){return t <=.05? linmix(t,20.,0.,25.658,3742.273):3742.273;}
float _SHKENV8(float t){return t <=.016? linmix(t,62.5,0.,.92,11.682):11.682;}
float _SHKENV9(float t){return t <=.016? linmix(t,62.5,0.,.92,11.682):11.682;}
float _SHKENV10(float t){return t <=.08? linmix(t,12.5,0.,1.667,3.337):t <=.8? linmix(t,1.3889,-.1111,3.337,.5):t <=1.2? linmix(t,2.5,-2.,.5,.2):.2;}

float rfloat(int off){return sequence_texture[off];}

#define NTRK 3
#define NMOD 10
#define NPTN 10
#define NNOT 105
#define NDRM 11

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}
float trk_pre(int index)    {return     rfloat(index+1+4*NTRK);}
float trk_slide(int index)  {return     rfloat(index+1+5*NTRK);} // idea for future: change to individual note_slide_time
float mod_on(int index)     {return     rfloat(index+1+6*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+6*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+6*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+6*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+6*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_pan(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+3*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+4*NNOT);}
float note_slide(int index) {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+5*NNOT);}
float note_aux(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+6*NNOT);}
float drum_rel(int index)   {return     rfloat(index+2+6*NTRK+4*NMOD+NPTN+7*NNOT);}

vec2 mainSynth(float time)
{
    float sL = 0.;
    float sR = 0.;
    float dL = 0.;
    float dR = 0.;

    time = mod(time, 14.9143);

    int _it;
    for(_it = 0; _it < NTIME - 2 && pos_t[_it + 1] < time; _it++);
    BPS = pos_BPS[_it];
    SPB = pos_SPB[_it];
    BT = pos_B[_it] + (time - pos_t[_it]) * BPS;

    float time2 = time - .0002;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, pre, f, amtL, amtR, env, slide, aux;
    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, syn, drum;

    for(int trk = 0; trk < NTRK; trk++)
    {
        tsep0 = trk_sep(trk);
        tsep1 = trk_sep(trk + 1);

        syn = trk_syn(trk);
        rel = trk_rel(trk);
        pre = trk_pre(trk);

        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1) - pre); _modU++);
        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + rel); _modL++);

        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            B = BT - mod_on(_mod);

            ptn   = mod_ptn(_mod);
            psep0 = ptn_sep(ptn);
            psep1 = ptn_sep(ptn + 1);

            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1) - pre); _noteU++);
            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + rel); _noteL++);

            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                if(syn == 35)
                {
                    drum = int(note_pitch(_note));
                    rel = drum_rel(drum);
                }

                amaysynL  = 0.;
                amaysynR  = 0.;
                amaydrumL = 0.;
                amaydrumR = 0.;

                Bon   = note_on(_note) - pre;
                Boff  = note_off(_note) + rel;
                L     = Boff - Bon;
                tL    = L * SPB;
                Bprog = B - Bon;
                Bproc = Bprog / L;
                _t    = Bprog * SPB;
                _t2   = _t - .0002;
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.);
                slide = note_slide(_note);
                aux   = note_aux(_note);

                if(syn == 35)
                {
                    env = trk_norm(trk) * theta(Bprog) * theta(L - Bprog);
                    if(drum == 0) { sidechain = min(sidechain, 1. - vel * (clamp(1e4 * Bprog,0.,1.) - pow(Bprog/(L-rel),8.))); }
                    else if(drum == 3){
                        amaydrumL = vel*(vel*_SHKENV0(_t)*(1.66*s_atan(2.1*_SHKENV8(_t)*(.5*_sin_(drop_phase(_t,.053,250.,57.878),.6*_SHKENV0(_t)*(.5*lpnoise(_t,325.282)+.5*lpnoise(_t,.745*325.282)))+.5*_sin_(.494*drop_phase(_t,.053,250.,57.878),.6*_SHKENV0(_t)*(.5*lpnoise(_t,325.282)+.5*lpnoise(_t,.745*325.282)))))+1.*_sin(drop_phase(_t,.053,250.,57.878)))+vel*_SHKENV1(_t)*(.13*sinshape((.5*pseudorandom(_SHKENV5(_t)*_t)+.5*pseudorandom(.45*_SHKENV5(_t)*_t)),_SHKENV9(_t),3.))+vel*_SHKENV2(_t)*(1.4*(.5*lpnoise(_t,325.282)+.5*lpnoise(_t,.745*325.282)))+vel*_SHKENV3(_t)*(1.03*_sin(drop2_phase(lofi(_t,(892.*_SHKENV10(_t))),.15,.4,901.313,1000.,200.))));
                        amaydrumR = vel*(vel*_SHKENV0(_t)*(1.66*s_atan(2.1*_SHKENV8(_t)*(.5*_sin_(drop_phase((_t-.00039),.053,250.,57.878),.6*_SHKENV0(_t)*(.5*lpnoise((_t-.00039),325.282)+.5*lpnoise((_t-.00039),.745*325.282)))+.5*_sin_(.494*drop_phase((_t-.00039),.053,250.,57.878),.6*_SHKENV0(_t)*(.5*lpnoise((_t-.00039),325.282)+.5*lpnoise((_t-.00039),.745*325.282)))))+1.*_sin(drop_phase(_t,.053,250.,57.878)))+vel*_SHKENV1(_t)*(.13*sinshape((.5*pseudorandom(_SHKENV5(_t)*(_t-.00168))+.5*pseudorandom(.45*_SHKENV5(_t)*(_t-.00168))),_SHKENV9(_t),3.))+vel*_SHKENV2(_t)*(1.4*(.5*lpnoise((_t-.00149),325.282)+.5*lpnoise((_t-.00149),.745*325.282)))+vel*_SHKENV3(_t)*(1.03*_sin(drop2_phase(lofi(_t,(892.*_SHKENV10(_t))),.15,.4,901.313,1000.,200.))));
                    }

                    dL += amtL * s_atan(env * amaydrumL);
                    dR += amtR * s_atan(env * amaydrumR);
                }
                else
                {
                    f = freqC1(note_pitch(_note) + mod_transp(_mod));

                    if(abs(slide) > 1e-3) // THIS IS SLIDEY BIZ
                    {
                        float Bslide = trk_slide(trk);
                        float fac = slide * log(2.)/12.;
                        if (Bprog <= Bslide)
                        {
                            float help = 1. - Bprog/Bslide;
                            f *= Bslide * (fhelp(fac) - help * fhelp(fac*help*help)) / Bprog;
                        }
                        else
                        {
                            f *= 1. + (Bslide * (fhelp(fac)-1.)) / Bprog;
                        }
                    }

                    env = theta(Bprog) * (1. - smoothstep(Boff-rel, Boff, B));
                    if(syn == 0){amaysynL = _sin(f*_t); amaysynR = _sin(f*_t2);}
                    else if(syn == 6){

                        amaysynL = .8*env_AHDSR(Bprog,L,.001,0.,.1,1.,.3)*(waveshape(clip(1.6*QFM((_t-0.0*(1.+2.*_sin(.15*_t))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8)
      +waveshape(clip(1.6*QFM((_t-2.0e-03*(1.+2.*_sin(.15*_t))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8)
      +waveshape(clip(1.6*QFM((_t-4.0e-03*(1.+2.*_sin(.15*_t))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8));
                        amaysynR = .8*env_AHDSR(Bprog,L,.001,0.,.1,1.,.3)*(waveshape(clip(1.6*QFM((_t2-0.0*(1.+2.*_sin(.15*_t2))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8)
      +waveshape(clip(1.6*QFM((_t2-2.0e-03*(1.+2.*_sin(.15*_t2))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8)
      +waveshape(clip(1.6*QFM((_t2-4.0e-03*(1.+2.*_sin(.15*_t2))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8));
                    }
                    else if(syn == 12){

                        amaysynL = env_AHDSR(Bprog,L,.002,0.,.15,.25,.13)*bandpassBPsaw1(_t,f,tL,vel,(2000.+(1500.*_sin(.25*BT))),10.,100.);
                        amaysynR = env_AHDSR(Bprog,L,.002,0.,.15,.25,.13)*bandpassBPsaw1(_t2,f,tL,vel,(2000.+(1500.*_sin(.25*BT))),10.,100.);
                    }

                    sL += amtL * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynL);
                    sR += amtR * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynR);
                }
            }
        }
    }
    return vec2(s_atan(.8 * sidechain * sL + .8 * dL), s_atan(.8 * sidechain * sR + .8 * dR));
}

void main()
{
    Tsample = 1./iSampleRate;
    float t = (iBlockOffset + gl_FragCoord.x + gl_FragCoord.y*iTexSize) * Tsample;
    vec2 s = mainSynth(t);
    vec2 v  = floor((0.5+0.5*s)*65535.0);
    vec2 vl = mod(v,256.0)/255.0;
    vec2 vh = floor(v/256.0)/255.0;
    fragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}

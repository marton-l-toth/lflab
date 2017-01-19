#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nz.h"
#include "glob.h"

#define NZGB 7
static double divpts[129] = {
0.0,0.00177486136336095,0.00354981218576197,0.00532503143037788,0.00710069822279239,0.00887699197759647,
0.0106540925254837,0.0124321802410763,0.0142114361717174,0.0159920421674704,0.0177741810125668,
0.0195580365585536,0.0213437938593926,0.0231316393087741,0.0249217607799136,0.0267143477681104,
0.0285095915363538,0.0303076852642789,0.0321088242007799,0.0339132058206074,0.0357210299852891,
0.0375324991087279,0.0393478183278528,0.0411671956787144,0.0429908422784388,0.0448189725134781,
0.0466518042346192,0.0484895589592411,0.05033246208134,0.052180743089875,0.0540346357960209,
0.0558943785699553,0.0577602145878462,0.0596323920897549,0.0615111646492153,0.0633967914553076,
0.0652895376081007,0.0671896744284037,0.0690974797828351,0.0710132384252957,0.0729372423560146,
0.0748697911994269,0.0768111926022446,0.0787617626531887,0.0807218263259707,0.0826917179472444,
0.0846717816913903,0.0866623721041566,0.0886638546573518,0.0906766063369793,0.0927010162674127,
0.0947374863744475,0.0967864320903217,0.0988482831040836,0.100923484161004,0.103012495915082,
0.105115795839077,0.107233879196952,0.10936726008408,0.111516472541113,0.113682071748023,
0.115864635305504,0.118064764611654,0.120283086342767,0.122520254047956,0.12477694986847,
0.12705388639375,0.129351808667664,0.131671496359938,0.134013766119565,0.136379474128998,0.138769518880265,
0.141184844196752,0.143626442527465,0.146095358544019,0.148592693074634,0.151119607414033,0.153677328053459,
0.156267151881278,0.158890451911817,0.16154868360855,0.164243391877595,0.166976218819101,0.16974891233777,
0.172563335729935,0.175421478383799,0.178325467752331,0.181277582785674,0.18428026904279,0.187336155741733,
0.190448075056031,0.193619084023164,0.196852489502846,0.200151876710912,0.203521141963796,0.20696453040428,
0.21048667964931,0.214092670514977,0.217788086245888,0.221579082024098,0.225472466981219,0.229475801520236,
0.233597513517888,0.237847037990487,0.242234986159776,0.24677335168629,0.251475764343689,0.256357804882,
0.261437399712891,0.266735321024276,0.272275828055423,0.278087500231326,0.28420433543561,0.290667221538737,
0.297525944406742,0.304841985271995,0.312692510915287,0.321176222295754,0.330422203262448,0.340603818163979,
0.351961538175002,0.364842534915935,0.379774190860012,0.397613016389257,0.419883453025696,0.449684875926625,
0.494798582564708,0.584003173991606,0.99999999999997};

#define NZG (1<<NZGB)

//static double divpts[33] = {   0.0, 0.00752119373899836,0.0150491979445374, 0.0225976953249547,
//         0.0301805926584789, 0.0378122020923033, 0.0455074352696795, 0.0532820172644446,
//         0.0611527285765055, 0.0691376853377064, 0.0772566706091478, 0.0855315335468761,
//         0.0939866787762359, 0.102649676324104,  0.111552034133378,  0.120730192477159,
//         0.130226825691108,  0.140092576898522,  0.150388415032167,  0.161188906866031,
//         0.172586870221355,  0.184700176146081,  0.197682014725453,  0.211736979953816,
//         0.227147426800969,  0.2443190828995,    0.263865537029438,  0.286779013127202,
//         0.314818117067482,  0.351544376555854,  0.405873479341178,  0.510820679386546, 1.0 };

double divptsx[33] = { 0.0, 0.00217093382031652, 0.00441859982834941, 0.00674856892484449,
        0.00916703747595322, 0.0116809240880007, 0.0142979857682489, 0.017026958307811,
        0.019877727186392, 0.0228615372893915, 0.0259912524692444, 0.0292816798076564,
        0.0327499788486178, 0.0364161838537321, 0.0403038785099428, 0.0444410794721183,
        0.0488614109067238, 0.0536056923279233, 0.0587241260785742, 0.0642793760652808,
        0.0703510080153575, 0.0770420762329545, 0.0844892204077339, 0.0928787540455465,
        0.102473517456856, 0.11366030822135, 0.127039790208231, 0.143613108166267,
        0.16521902355074, 0.19574755744358, 0.245502923888483, 0.355802952884758, 1.0};

static const double nz_1mul = 1.0 / (double)(1u<<31u);
static const double nz_2mul = 1.0 / (double)(1u<<30u);
static const double nz_3mul = 1.0 / (double)(0x2ffffe);
static const double nz_xmul = 1.0 / (65536.0*65536.0*32768.0*32768.0);
static const double nz_ymul = 1.0 / (32768.0*32768.0*32768.0*32768.0);
static const double nz_qmul[2] = {.999999999/(32767.0*32767.0*32767.0*32767.0),
      	                         -.999999999/(32767.0*32767.0*32767.0*32767.0) };
static const double nz_hmul[2] = {.999999999/(32767.0*32767.0*32767.0*32767.0*32767.0*32767.0),
      	                         -.999999999/(32767.0*32767.0*32767.0*32767.0*32767.0*32767.0) };
typedef struct { int iacc,rsrv; double x0, xu, y0, yu; } zgent_t;
static zgent_t zgent[2<<NZGB];
static zgent_t zgentx[64];

static inline double f(double x) { return exp(-16.0*x*x); }
static inline double fx(double x) { return exp(-16.0*fabs(x)); }

static struct random_data r_dat; static char r_buf[256];
 
void nz_init() {
	initstate_r(zero_sec, r_buf, 245, &r_dat);
	double xmul = 1.0/(double)((1u<<31u)), ymul = 1.0/(double)((1<<(30-NZGB)));
	for (int i=0; i<NZG; i++) {
		zgent_t * p = zgent + i;
		p->x0 = divpts[i]; p->xu = xmul * (divpts[i+1]-p->x0); p->x0 += .5*p->xu;
		p[NZG].x0 = - p->x0; p[NZG].xu = - p->xu;
		p[NZG].yu = p->yu = ymul * (f(p->x0));
		p[NZG].iacc = p->iacc = (int)lround(f(divpts[i+1]) / p->yu);
	}
	for (int i=0; i<32; i++) {
		zgent_t * p = zgentx + i;
		p->x0 = divptsx[i]; p->xu = xmul * (divptsx[i+1]-p->x0); p->x0 += .5*p->xu;
		p[32].x0 = - p->x0; p[32].xu = - p->xu;
		p[32].yu = p->yu = ymul * (fx(p->x0));
		p[32].iacc = p->iacc = (int)lround(fx(divptsx[i+1]) / p->yu);
	}
}

static int iacc_ny[2];
void mk_gwn(double * to, int n) {
        int i=0; while (i<n) {
		int xi, yi; random_r(&r_dat, &xi); random_r(&r_dat, &yi);
                int k = yi>>(30-NZGB);
                yi &= ( (1<<(30-NZGB)) - 1);
                zgent_t * zp = zgent + k;
                double x = zp->x0 + (double)(xi) * zp->xu;
		if(
#ifdef NZ_DEBUG
		(k=(yi<zp->iacc), ++iacc_ny[k], k)
#else
                (yi<zp->iacc)
#endif
				|| (double)yi * zp->yu < f(x)) to[i++] = x;
        }}

void gwn_debug(int step) {
	int n = iacc_ny[0], y = iacc_ny[1], tot = n+y; double prc = 100.0 / (double)tot;
	log("perfstat: iacc: tot:%d y:%d(%g%%) n:%d(%g%%)", tot, y, prc*(double)y,
								 n, prc*(double)n); }

static void mk_exp(double * to, int n) {
        int i=0; while (i<n) {
                int xi = random(), yi = random();
                int k = yi>>25; yi &= 0x1ffffff;
                zgent_t * zp = zgentx + k;
                double x = zp->x0 + (double)(xi) * zp->xu;
                if ((yi<zp->iacc) || (double)yi * zp->yu < fx(x)) to[i++] = x;
        }}

static void mk_uexp(double * to, int n) {
        int i=0; while (i<n) {
                int xi = random(), yi = random();
                int k = yi>>26; yi &= 0x1ffffff;
                zgent_t * zp = zgentx + k;
                double x = zp->x0 + (double)(xi) * zp->xu;
                if ((yi<zp->iacc) || (double)yi * zp->yu < fx(x)) to[i++] = x;
        }}

static void mk_1x31(double *to, int n) {
        for (int k,i=0; i<n; i++) k=random(), k += k-0x7fffffff, to[i] = nz_1mul * (double)k; }

static void mk_2x31(double *to, int n) {
        for (int k,i=0; i<n; i++) k=random()-0x7fffffff, k += random(), to[i] = nz_1mul * (double)k; }

static void mk_3x20(double *to, int n) {
        for (int k,x,y,i=0; i<n; i++) {
                x=random(); y=random(); k = (unsigned int)(x+y) >> 11u;
                k += (x&1023) + 1024*(y&1023); k += k - 0x2ffffd;
                to[i] = nz_3mul * (double)k;
        }}

static void mk_mul2(double *to, int n) {
        for (int x,y,i=0; i<n; i++) {
                x=random(); y=random();
                x += (x-0x7fffffff); y += (y-0x7fffffff);
                to[i] = nz_xmul * (double)x * (double)y;
        }}

static void mk_mul3(double *to, int n) {
        for (int x,y,z,i=0; i<n; i++) {
                x=random(); y=random(); z = 1024*(x&1023) + (y&1023);
                x>>=10; x&=~1; x-=0xfffff;
                y>>=10; y&=~1; y-=0xfffff;
                z += z - 0xfffff;
                to[i] = nz_ymul * (double)x * (double)y * (double)z;
        }}

static void mk_mul4(double *to, int n) {
        for (int x,y,z,v,i=0; i<n; i++) {
                x=random(); y=random();
                z = (x&32767) * (y&32767);
                v = (x>>16) * (y>>16);
                to[i] = nz_qmul[((x^y)>>15)&1] * (double)z * (double)v;
        }}

static void mk_mul6(double *to, int n) {
        for (int x,y,z,i=0; i<n; i++) {
                x=random(); y=random(); z=random();
		int sg = ((x^y^z)>>15) & 1;
		x = (x&32767) * (x>>16);
		y = (y&32767) * (y>>16);
		z = (z&32767) * (z>>16);
                to[i] = nz_hmul[sg] * (double)x * (double)y * (double)z;
        }}

int rndbits(int * st2, int n) {
	int r = 0;
	if (n > st2[1]) n -= st2[1], r = *st2<<n, *st2 = random(), st2[1] = 31;
	r += *st2 & ((1<<n)-1), *st2 >>= n, st2[1] -= n; return r;
}

int nz_bv_c(unsigned int * to, double v, int n) {
	int r = 0, hl = (int)(lround(v * 1048576.0)), hl5 = hl >> 15;
	int n2 = (n+31) >> 5; memset(to, 0, 4*n2);
	int rs[2]; rs[0] = rs[1] = 0;
	for (int i=0; i<n; i++) {
		int k0 = rndbits(rs,5), k = k0 - hl5;
		unsigned int m = (unsigned int) (k ? k>>6 : ((k0<<15) + rndbits(rs,15) - hl) >> 31);
		r -= (int) m; to[i>>5] |= (m & (1u << (i&31u))); }
	return r;
}

nz_fun0_t nz_fun0[NZ_N_FUN0] = {mk_gwn, mk_1x31, mk_2x31, mk_3x20, mk_mul2, mk_mul3, mk_mul4, mk_mul6, mk_exp, mk_uexp } ;


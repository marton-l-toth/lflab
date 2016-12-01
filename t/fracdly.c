#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define FP0(S) 		fprintf(stderr, S "\n")
#define FP1(S,X) 	fprintf(stderr, S "\n", (X))
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y))
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z))
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W))

double * store64 = NULL;
int nnn;
#define BLK(I,J) (store64 + 64*(nnn*((I)+1)+J))

void ini64(double *to, double dly) {
	if (fabs(dly)<1e-99) 	 return (void)(to[0]=1.0, 	     memset(to+1, 0, 504));
	if (fabs(dly-1.0)<1e-14) return (void)(to[0]=0.0, to[1]=1.0, memset(to+2, 0, 496));
	double z = (1.0-dly)/(1.0+dly), y = 1.0-z*z;
	to[0]=z; int i; for(i=1; i<64; i++) to[i] = y, y *= -z;
}

void filt64(double *to, double *in, double dly) {
	int i; double x, z = (1.0-dly)/(1.0+dly), x1 = 0.0, y1 = 0.0;
	for (i=0; i<64; i++) x=in[i], to[i]=y1=z*(x-y1)+x1, x1=x;
}

#define DIFF(P,Q,J) (dtmp=(P)[J]-(Q)[J], dtmp*dtmp)
int diffsq_lim(double *plim, double *p, double *q) {
	double dtmp, acc = DIFF(p,q,0); if ((acc += DIFF(p,q,1)) >=* plim) return 0;
	if (acc+=DIFF(p,q,2), acc+=DIFF(p,q,3), acc+=DIFF(p,q,4), acc+=DIFF(p,q,5), acc>=*plim) return 0;
	int i; for (i=6;  i<16; i++) acc += DIFF(p,q,i);   if (acc>=*plim) return 0;
	for 	   (i=16; i<64; i++) acc += DIFF(p,q,i);   if (acc>=*plim) return 0;
	*plim = acc; return 1;
}

#define T(J) (un*(double)(J))
int main(int ac, char** av) {
	if (ac!=5) FP1("usage: %s op t0 t1 div", *av), exit(1);
	int div = atoi(av[4]);
	double t0d = atof(av[2]), t1d = atof(av[3]), ddiv = (double)div, un = 1.0/div;
	int k0 = (int)lround(t0d*ddiv), j, k, m, n, kl, kc, jm=-1, km=-1,
	    k1 = (int)lround(t1d*ddiv);
	nnn = k1-k0+1; store64 = malloc((nnn*nnn+nnn+1)*512);
	for (j=0; j<nnn; j++) ini64(BLK(-1,j), T(k0+j));
	for (j=0; j<nnn; j++) for (k=0; k<nnn; k++) filt64(BLK(j,k),BLK(-1,j),T(k0+k));

	double diffsq = 1e99, *q = BLK(nnn,0);
	for (j=k0; j<=k1; j++) { fprintf(stderr, "(%.15g) ", T(j)); fflush(stderr); for (k=k0; k<=k1; k++) {
		double *q2 = BLK(j-k0,k-k0);
		int t3 = j+k, mc, nc;
		for (m=k0; (mc=t3-m)>=2*k0; m++) for (n=k0; n<=k1&&(nc=mc-n)>=k0; n++)
			if (filt64(q, BLK(n-k0, nc-k0), T(m)), diffsq_lim(&diffsq, q, q2))
				fprintf(stderr, "%.15g + %.15g = %.15g + %.15g + %.15g (%.15g)\n",
						T(j),T(k),T(m),T(n),T(nc),sqrt(diffsq));
	}}
	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define FP0(S) 		fprintf(stderr, S "\n")
#define FP1(S,X) 	fprintf(stderr, S "\n", (X))
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y))
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z))
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W))

#define CLKDIF(P,Q) ((int)(1000000000*((P)->tv_sec-(Q)->tv_sec)+(P)->tv_nsec-(Q)->tv_nsec))

clockid_t clkt[3] = { CLOCK_MONOTONIC, CLOCK_MONOTONIC_COARSE, CLOCK_MONOTONIC_RAW };

inline int clk_upd(struct timespec *q, clockid_t ty) {
	struct timespec t; memcpy(&t, q, sizeof(t));
	clock_gettime(ty, q); return CLKDIF(q, &t); }

int main(int ac, char** av) {
	if (ac != 3) return FP1("usage: %s flg cnt", *av);
	int j,flg = atoi(av[1]), n = atoi(av[2]);
	clockid_t ty = (flg&1) ? CLOCK_MONOTONIC_RAW : CLOCK_MONOTONIC;
	if (flg&2) {
		struct timespec ts1; clk_upd(&ts1, ty);
		int td[n]; for (j=0; j<n; j++) td[j] = clk_upd(&ts1, ty);
		for (j=0; j< n; j++) printf(" %d", td[j]);
	} else {
		struct timespec ts[n+1];
		for (j=0; j<=n; j++) clock_gettime(ty, ts+j);
		for (j=0; j< n; j++) printf(" %d", CLKDIF(ts+j+1, ts+j));
	}
	putchar(10); return 0;
}

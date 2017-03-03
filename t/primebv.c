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

int r() { char buf[64]; return fgets(buf, 63, stdin) ? atoi(buf) : -1; }

static unsigned int bv[2048];

int main(int ac, char** av) {
	int i,k;
	while((k=r())>=0) k>>=1, bv[(k>>5)&2047] |= 1u<<(k&31);
	for (i=0; i<2048; i++) printf("\n% 10u,"+!!(i%10),bv[i]);
	puts("");
	return 0;
}

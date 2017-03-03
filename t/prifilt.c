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


static char buf[256];
static int bufl = 0, cnt = 0;

int r() { char buf[64]; return fgets(buf, 63, stdin) ? atoi(buf) : -1; }
void fl() { int c=buf[bufl]; buf[bufl]=10; write(1, buf, bufl+1); buf[bufl]=c; }
void w(unsigned int x) {
	int j=bufl,n=sprintf(buf+j,"%u,",x);++cnt; bufl=j+n<112?j+n:(fl(),memcpy(buf,buf+j,n),n); }

static unsigned int bv[2048];

int main(int ac, char** av) {
	int i,k;
	double z = 0.0, mm = ac>1 ? atof(av[1]) : 1.0;
	while ((k=r())>=0) { if ((double)k>=z) w(k>>0), z=(double)k*mm; }
	fl(); printf("// %d\n", cnt);
}

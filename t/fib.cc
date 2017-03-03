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
void w(double x) {
	int j=bufl,n=sprintf(buf+j,"%.15g,",x);++cnt; bufl=j+n<112?j+n:(fl(),memcpy(buf,buf+j,n),n); }

static unsigned int bv[2048];

int main(int ac, char** av) {
	int i,k; double x,y,z;
	int n = ac>1 ? atoi(av[1]) : 20, 
	    m = ac>2 ? atoi(av[2]) : 0;
	double pre[m+!!m];
	for (i=0,x=0.0,y=1.0; i<m; i++) z=x, pre[i]=x=y-x, y=z;
	for (i=0,x=0.0,y=1.0; i<n; i++) w(x), z=y, y+=x, x=z;
	for (i=m-1; i>=0; i--) w(pre[i]);
	fl();
}

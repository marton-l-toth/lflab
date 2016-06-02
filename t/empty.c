#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define FP0(S) 		fprintf(stderr, S "\n");
#define FP1(S,X) 	fprintf(stderr, S "\n", (X));
#define FP2(S,X,Y) 	fprintf(stderr, S "\n", (X),(Y));
#define FP3(S,X,Y,Z) 	fprintf(stderr, S "\n", (X),(Y),(Z));
#define FP4(S,X,Y,Z,W) 	fprintf(stderr, S "\n", (X),(Y),(Z),(W));

int main(int ac, char** av) {
	FP2("hello %s %d", *av, ac);
	return 0;
}

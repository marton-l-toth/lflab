#include <stdio.h>
#include <stdlib.h>

static int osb_hash(const char *s) { // od|sed|bc
	int r=0; while (*s) r *= (*s>99)?1000:100, r += *s, s++;   return r&0x3fffffff; }

#include "cfg2.inc"

int main(int ac, char**av) {
	int i; for (i=1; i<ac; i++) printf(" %d", osb_find(av[i]));
	putchar(10); return 0; }


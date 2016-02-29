#ifndef __qwe_wrap_h__
#define __qwe_wrap_h__

class BoxGen; class AReader;
static inline int wr_ixtr(int ix) { return (ix&95) + ((0x106>>((ix>>3)&12))&15); }
static inline int wr_ixtr_r(int ix) { int k = ix>>6; return ix + (ix<6+60*k ? 32 : 5*k-6); }
// these fns *assume* that arg is a wrbox
const char * wrap_rgb(BoxGen * bx);
int wrap_2mx_txdkl(BoxGen * bx, int trg, int xflf, int dly, int mxky, int lim);
int wrap_qdiff(BoxGen * b1, BoxGen * b2);
int wrap_dump_keycfg(char * to, unsigned int * bv, short ** pk);
int wrap_key_op(BoxGen * bx, int ky, int op, const char *s, int nof);
void wrap_set_trec(BoxGen * bx, int j);
AReader * wrap_avreader(BoxGen * bx);

#endif // __qwe_wrap_h__

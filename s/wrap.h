#ifndef __qwe_wrap_h__
#define __qwe_wrap_h__

class BoxGen; class AReader; class ABoxNode;
static inline int wr_ixtr(int ix) { return (ix&95) + ((0x106>>((ix>>3)&12))&15); }
static inline int wr_ixtr_r(int ix) { int k = ix>>6; return ix + (ix<6+60*k ? 32 : 5*k-6); }
// these fns *assume* that arg is a wrbox
int wrap_2mx_txdlv(BoxGen * bx, int trg, int xflg, int dly, int lim, double *v);
int wrap_nd_2mx(ABoxNode * bn, int trg, double bpm, int dly);
int wrap_qdiff(BoxGen * b1, BoxGen * b2);
int wrap_dump_keycfg(char * to, unsigned int * bv, short ** pk);
int wrap_key_op(BoxGen * bx, int ky, int op, const char *s, int nof, int dly_tg=-1);
void wrap_set_trec(BoxGen * bx, int j);
AReader * wrap_avreader(BoxGen * bx, int cflg);
int swrap_grab_c(BoxGen *bx, int f);
int wrap_midi_ev(unsigned int j5i20o7, int ky, int val, const unsigned int * blk);
void wrap_set_trg(BoxGen * bx, int trg);

#endif // __qwe_wrap_h__

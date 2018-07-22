#ifndef __qwe_trk_h__
#define __qwe_trk_h__

class BoxGen; class ABoxNode;
extern BoxGen * trk_rec_trg;
BoxGen * trk_mk(ABoxNode * nd, ANode * g0, ANode * g1);
ANode * trk_bkm_find(BoxGen * abx, int j);
void trk_bkm_add(BoxGen * abx, ANode * nd);
void trk_bkm_rm(BoxGen * abx, ANode * nd);
int trk_cond_pm(BoxGen * abx, ANode * nd, int pm);
int trk_rec(ANode * wnd, int mxbi);
int trk_glob_paste(BoxGen * bx, int nof);
void trk_w_gna(int id);
void trk_sane(BoxGen * abx, int j);
int trk_mx99(int k);
int trk_keyop_2q(int k, int op, const char *xys, int id, int trg);


#ifdef SANCHK_TRK
#define TRK_SANE1(J) sane(J)
#define TRK_SANE2(P,J) trk_sane(P,J)
#define TRK_CO2(P,Q,S) chk_ord(P,Q,S)
#define TRK_CO2B(B,P,Q,S) (static_cast<TBoxNode*>(B)->chk_ord(P,Q,S))
#define TRK_CO1(P,S) (chk_ord((P)->cth()->pv,P,S "1"),chk_ord(P,(P)->next(),S "2"))
#else
#define TRK_SANE1(J)	  ((void)0)
#define TRK_SANE2(P,J)	  ((void)0)
#define TRK_CO2(P,Q,S)	  ((void)0)
#define TRK_CO2B(B,P,Q,S) ((void)0)
#define TRK_CO1(P,S)	  ((void)0)
#endif

#endif // __qwe_trk_h__

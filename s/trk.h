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

#endif // __qwe_trk_h__

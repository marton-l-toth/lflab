#ifndef __qwe_midi_h__
#define __qwe_midi_h__

extern unsigned int midi_bv;
extern int midi_fd[32];
void midi_input(int i);
int midi_cmd(const char *s);
int midi_grab(int id, int ix, int dev, int ch, int kc, const unsigned char *kv, int flg);
int midi_w_ln(char *to, int j),
    midi_w_slg(char *to, int j),
    midi_w_flg();

inline const unsigned int * mi_getblk(int ldev, int ch) {
	extern unsigned int**mi_root[32]; extern char mi_tr_l2p[32]; return mi_root[(int)mi_tr_l2p[ldev]][ch];}

#endif // __qwe_midi_h__

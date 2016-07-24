#include <math.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "uc0.h"
#include "uc1.h"
#include "glob.h"
#include "util.h"
#include "errtab.inc"
#include "cfgtab.inc"
#include "wrap.h"

#define DBGC (debug_flags&DFLG_MIDIEV)

// .0078740157480315
int rawmidi_desc(char *to, int id, int maxlen); //asnd.cc
static unsigned int mi_dfltc[256], *mi_keyspd;
static unsigned int*mi_dflti[16] = {mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,
				     mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc};
unsigned int ** mi_root[32] = { mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti};
static char *mi_dsc[32], mi_devname[24], keytrans[128];
char mi_tr_l2p[32], mi_tr_p2l[32];
int midi_fd[32];
unsigned int midi_bv;
static const char ktrr_02_27[26] = {15,2,16,3,17,4,18,19,6,20,7,21,8,22,23,10,24,11,25,26,13,27,14,5,9,12},
	     	  ktrr_31_54[24] = {44,31,45,32,46,47,34,48,35,49,36,50,51,38,52,39,53,54,40,33,37,41,42,43};
static int ktr_flg;

static int devnm_cmp(const void *p, const void *q) { 
	return strcmp(mi_dsc[(int)mi_tr_l2p[(int)*(const char*)p]],
		      mi_dsc[(int)mi_tr_l2p[(int)*(const char*)q]]); }

static int set_kflg(int f) {
	int d = (ktr_flg ^ f);
	if(d&1) { if(f&1) keytrans[69]=87, keytrans[70]=88, keytrans[87]=69, keytrans[88]=70;
		  else	  keytrans[69]=69, keytrans[70]=70, keytrans[87]=87, keytrans[88]=88; }
	if(d&2) { if(f&2) for(int i=0; i<26; i++) keytrans[(int)ktrr_02_27[i]]  = i+2;
		  else    for(int i=0; i<26; i++) keytrans[i+2]			= i+2; }
	if(d&4) { if(f&4) for(int i=0; i<24; i++) keytrans[(int)ktrr_31_54[i]]  = i+31;
		  else    for(int i=0; i<24; i++) keytrans[i+31]		= i+31; }
	return ktr_flg = f, 0;
}

static void mi_i_ini(int dev) { memcpy( mi_root[dev] = (unsigned int**)nf_alloc(16*sizeof(int*)),
					mi_dflti, 16*sizeof(int*)); }
inline static unsigned int * mi_ro(int dev, int ch, int j) { return mi_root[dev][ch]+j; }
inline static unsigned int * mi_rw(int dev, int ch, int j) {  unsigned int **pp=mi_root[dev];
	return (pp[ch]==mi_dfltc ? (pp[ch]=(unsigned int*)nf_alloc(1024)) : pp[ch])+j; }

static int midi_wr(int ix, const char * s) {
	unsigned char buf[64];  int ec, n=0;
	for (int i=0; i<64; i++) {
		while (*s==32) ++s;
		if (!*s) break;
		if (!s[1]) return MDE_PARSE;
		buf[n++] = hex2(s); s += 2;
	}
	if ((ec=write(midi_fd[ix], buf, n))<=0) return ec ? EEE_ERRNO : EEE_ZEROLEN;
	return 0;
}
static void mi_xchg(int i, int j){ int I=mi_tr_l2p[i], J=mi_tr_l2p[j];  mi_tr_l2p[i]=J; mi_tr_l2p[j]=I; 
									mi_tr_p2l[I]=j; mi_tr_p2l[J]=i; }
static void mi_log(int i, const char *s, const unsigned char *p, int n) {
	char buf[1024]; int bl = sprintf(buf, "midi: %02d (%02d,%s) -- %s %db:",
					      mi_tr_p2l[i], i, mi_dsc[i], s, n);
	for (int x,j=0; j<n; j++) x=p[j], buf[bl++]=32,buf[bl++]=hexc1(x>>4),buf[bl++]=hexc1(x&15);
	buf[bl++]=10; log_sn("", buf, bl);
}

static void mi_err(int i, int f, int ec) {
	log("midi(%d): %s, %s", i, err_str(ec), f?"removed":"msg ignored"); }

static int midi_f0(int i, const unsigned char *p, int r) {
	int j; for (j=1; j<r; j++) if (p[j]==0xf7) goto found;
	mi_log(i, "unterm.sys.", p, r); return MDE_FRAGF0;
found:	mi_log(i, "system", p, ++j); return j; }

static void midi_kc(int i, int ch, int k, int v) {
 	unsigned int *p0 = mi_rw(i, ch, 0), *p = p0+k, x = *p, x25 = x & ~127u;
	if (DBGC) log("midi_kc: i=(p%d l%d) ch=%d kc=%d v=%d oldv=%d", i, mi_tr_p2l[i], ch, k, v, *p);
	if (!x25) return (void) (*p = v);  else *p = x25|v;
	int ec = wrap_midi_ev(x, k, v, p0); if (ec>=0) return;
	if (ec==MDE_KEEPV) return (void) (*p = x);
	if (ec!=EEE_NOEFF) gui_errq_add(ec,"midi/swr");
	else if (DBGC) log("midi_ev: 0x%x -> 0x%x - no effect",x25,v);
}

int midi_grab(int id, int ix, int dev, int ch, int kc, const unsigned char *kv, int flg) {
	if (DBGC) log_n("midi_grab: id=0x%x, ix=%d, d=%d ch=%d kc=%d [", id, ix, dev, ch, kc);
	int dev2 = mi_tr_l2p[dev]; if (!((midi_bv|0x80000000)&(1u<<dev2))) return MDE_GRABWHAT;
	unsigned int m = ((id<<7)|(ix<<27)) & ((flg&1)-1), *p = mi_rw(dev2, ch, 0);
	if (DBGC) for (int j, i=0; i<kc; i++) j = kv[i], p[j] = (p[j]&127)|m, log_n(" %02x",kv[i]);
	else      for (int j, i=0; i<kc; i++) j = kv[i], p[j] = (p[j]&127)|m;
	if (DBGC) log(" ]");
	return 0;
}

#define CHK1 if (r<2 || p[1]>127) return mi_err(i,0,MDE_FRAG12)
#define CHK2 if (r<3 || ((p[1]|p[2])&128)) return mi_err(i,0,MDE_FRAG23-(r==1||p[2]>127))
void midi_input(int i) {
	unsigned char buf[512], *p = buf;
	int v, c, r = read(midi_fd[i], buf, 512), x=*p;
	if (r<0) return close(midi_fd[i]), midi_bv &= ~(1<<i), mi_err(i,1,EEE_ERRNO);
	while (r) {
		switch(x>>4) {
			case 8: v=c=0; goto kc;
			case 9: *mi_keyspd = v = p[2]; c=0; goto kc;
			case 10: case 14: goto l2;
			case 11: v=p[2]; c=128; *mi_keyspd = 127; goto kc;
			case 12: case 13: goto l1;
			case 15: switch(x) { case 0xf0: goto f0;
					 case 0xf1: case 0xf3: goto l1;
					 case 0xf2: goto l2;
					 default: goto l0; }
			default: return mi_err(i,0,MDE_SEVEN); 
		}
l0:		mi_log(      i, "unused", p, 1); ++p;  --r;  continue;
l1:		CHK1; mi_log(i, "unused", p, 2); p+=2; r-=2; continue;
l2:		CHK2; mi_log(i, "unused", p, 3); p+=3; r-=3; continue;
kc:		CHK2; midi_kc(i, x&15, p[1]+c, v); p+=3, r-=3; continue;
f0:		if ((v = midi_f0(i,p,r))<0) return mi_err(i,0,v);  p+=v; r-=v; continue;
	}}

int midi_open(int ix, int id) {
	char nmb[256]; int l = 1; nmb[0] = 63; nmb[1] = 0;
	mi_devname[14] = 48+(id>>4); mi_devname[16] = 48+(id&15);
	int fd = open(mi_devname, ((debug_flags&DFLG_RWMIDI) ? O_RDWR : O_RDONLY)|O_NONBLOCK);
	if (fd<0) goto err; else close(fd);
	l = rawmidi_desc(nmb, id, 255)+1;  memcpy(mi_dsc[ix]=nf_alloc(l), nmb, l);
	mi_devname[14] = 48+(id>>4); mi_devname[16] = 48+(id&15);
	if ((fd = open(mi_devname, (debug_flags&DFLG_RWMIDI) ? O_RDWR : O_RDONLY)) < 0) goto err;
	log("midi%x: %s(%s): opened(%d)", ix, mi_devname , nmb, fd); return fd;
err:	gui_errq_add(EEE_ERRNO); gui_errq_add(MDE_OFAIL);
	log("midi%x: %s(%s): %s", ix, mi_devname, nmb, strerror(errno)); return -1;
}


int midi_cmd(const char *s) {
	if (!s || !*s) return MDE_PARSE;
	int j, c = *(s++), ix = b32_to_i(c);
	if (ix<0) { switch(c) {
		case '@': return set_kflg(atoi_h(s));
		default: return MDE_PARSE; }}
	switch(*s) {
		case 'K': return *mi_keyspd=127, j=hex2(s+1), midi_kc(31,0, keytrans[j&127], 127&(-(j>>7))), 0;
		case 'W': return midi_wr(ix, s+1);
		case 'X': if (ix==31) return MDE_VXCHG;
			  j = b32_to_i(s[1]);
			  if (j<0) { if ((s[1]-43)&~2) return MDE_PARSE; else j=44-s[1]; }
			  return ((j&=31) == 31) ? MDE_VXCHG : (mi_xchg(ix, j), 0);
		default:  return MDE_PARSE;
	}}

void midi_init() {
	for (int i=0; i<32; i++) mi_tr_l2p[i] = i;   memcpy(mi_tr_p2l, mi_tr_l2p, 32);
	mi_i_ini(31); mi_keyspd = mi_rw(31, 0, 255);
	memcpy(mi_devname, "/dev/snd/midiCxDy", 18);
	unsigned char buf[32]; 
	int n = 0, n0 = find_dev(buf, 1, 31);
	for (int i=0; i<n0; i++) if ((midi_fd[n]=midi_open(n,buf[i]))>=0) midi_bv|=1u<<n, mi_i_ini(n), n++;
	if (n>1) qsort(mi_tr_l2p, n, 1, &devnm_cmp);
	for (int i=0; i<n; i++) mi_tr_p2l[(int)mi_tr_l2p[i]] = i;
	for (int i=0; i<128; i++) keytrans[i] = i;   set_kflg(CFG_MIDI_KTR.i);
}


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

#define DBGC (debug_flags&DFLG_MIDIEV)

// .0078740157480315
int rawmidi_desc(char *to, int id, int maxlen); //asnd.cc
static unsigned int mi_dfltc[256];
static unsigned int*mi_dflti[16] = {mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,
				     mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc};
unsigned int ** mi_root[32] = { mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti};
static char *mi_dsc[32], mi_devname[24];
int midi_fd[32];
unsigned int midi_bv;

char mi_tr_l2p[32]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},
     mi_tr_p2l[32]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};

static void mi_i_ini(int dev) { memcpy( mi_root[dev] = (unsigned int**)nf_alloc(16*sizeof(int*)),
					mi_dflti, 16*sizeof(int*)); }
inline static unsigned int * mi_ro(int dev, int ch, int j) { return mi_root[dev][ch]+j; }
inline static unsigned int * mi_rw(int dev, int ch, int j) {  unsigned int **pp=mi_root[dev];
	return (pp[ch]==mi_dfltc ? (pp[ch]=(unsigned int*)nf_alloc(1024)) : pp[ch])+j; }

static void mi_xchg(int i, int j){ int I=mi_tr_l2p[i], J=mi_tr_l2p[j];  mi_tr_l2p[i]=J; mi_tr_l2p[j]=I; 
									mi_tr_p2l[I]=j; mi_tr_p2l[J]=i; }
int midi_open(int ix, int id) {
	char nmb[256]; int l = rawmidi_desc(nmb, id, 255)+1;
	memcpy(mi_dsc[ix]=nf_alloc(l), nmb, l);
	mi_devname[14] = 48+(id>>4); mi_devname[16] = 48+(id&15);
	int fd = open(mi_devname, (debug_flags&DFLG_RWMIDI) ? O_RDWR : O_RDONLY); 
	if (fd<0) log("midi%x: %s(%s): %s", ix, mi_devname, nmb, strerror(errno));
	else 	  log("midi%x: %s(%s): opened(%d)", ix, mi_devname , nmb, fd);
	return fd;
}

int midi_cmd(const char *s) { // TODO
	if (!s || !*s) return MDE_PARSE;
	unsigned char buf[64];
	int ec, ix = b32_to_i(*(s++)), n = 0;
	for (int i=0; i<64; i++) {
		while (*s==32) ++s;
		if (!*s) break;
		if (!s[1]) return MDE_PARSE;
		buf[n++] = hex2(s); s += 2;
	}
	if ((ec=write(midi_fd[ix], buf, n))<=0) return ec ? EEE_ERRNO : EEE_ZEROLEN;
	return 0;
}

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

static void midi_kc(int i, int ch, int kc, int v) {
 	unsigned int *p = mi_rw(i, ch, kc);
	if (DBGC) log("midi_kc: i=%d ch=%d kc=%d v=%d oldv=%d", i, ch, kc, v, *p);
	*p = v; // TODO 
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
			case 9: v=p[2]; c=0; goto kc;
			case 10: case 14: goto l2;
			case 11: v=p[2]; c=128;  goto kc;
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

void midi_init() {
	if (!CFG_DEVEL.i) return;
	memcpy(mi_devname, "/dev/snd/midiCxDy", 18);
	unsigned char buf[32]; 
	int n = find_dev(buf, 1, 32);
	for (int i=0; i<n; i++) if ((midi_fd[i]=midi_open(i,buf[i]))>=0) midi_bv|=1u<<i, mi_i_ini(i);
}


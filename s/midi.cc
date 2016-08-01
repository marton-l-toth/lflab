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
#include "util2.h"
#include "errtab.inc"
#include "cfgtab.inc"
#include "wrap.h"

#define DBGC (debug_flags&DFLG_MIDIEV)

void gui_midi(int flg); 			//guistub.cc
int rawmidi_desc(char *to, int id, int maxlen); //asnd.cc
static unsigned int mi_dfltc[256], *mi_keyspd;
static unsigned int*mi_dflti[16] = {mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,
				     mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc,mi_dfltc};
unsigned int ** mi_root[32] = { mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,
				mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti,mi_dflti};
static char *mi_dsc[32], keytrans[128], mi_devname[24] = "/dev/snd/midiCxDy";

char mi_tr_l2p[32], mi_tr_p2l[32], mi_dscl[32], midi_dvid[32];
int midi_fd[32];
unsigned int midi_bv, midi_err_bv;
static const char ktrr_02_27[26] = {15,2,16,3,17,4,18,19,6,20,7,21,8,22,23,10,24,11,25,26,13,27,14,5,9,12},
	     	  ktrr_31_54[24] = {44,31,45,32,46,47,34,48,35,49,36,50,51,38,52,39,53,54,40,33,37,41,42,43};
static int midi_dispflg = 0, midi_vsl_ix[2] = {128, 136};

static int devnm_cmp(const void *p, const void *q) {
	int i = *(const char*)p, j = *(const char*)q, sc = strcmp(mi_dsc[i], mi_dsc[j]);
	return sc ? sc : midi_dvid[i] - midi_dvid[j]; }

static int set_kflg(int f) {
	int d = f<0 ? (f=CFG_MIDI_KTR.i, -1) : (CFG_MIDI_KTR.i^f); if (!d) return 0;
	if(d&1) { if(f&1) keytrans[69]=87, keytrans[70]=88, keytrans[87]=69, keytrans[88]=70;
		  else	  keytrans[69]=69, keytrans[70]=70, keytrans[87]=87, keytrans[88]=88; }
	if(d&2) { if(f&2) for(int i=0; i<26; i++) keytrans[(int)ktrr_02_27[i]]  = i+2;
		  else    for(int i=0; i<26; i++) keytrans[i+2]			= i+2; }
	if(d&4) { if(f&4) for(int i=0; i<24; i++) keytrans[(int)ktrr_31_54[i]]  = i+31;
		  else    for(int i=0; i<24; i++) keytrans[i+31]		= i+31; }
	return CFG_MIDI_KTR.i = f, 4096;
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
static void mi_xchg(int i, int j){ int t, I=mi_tr_l2p[i], J=mi_tr_l2p[j]; mi_tr_l2p[i]=J; mi_tr_l2p[j]=I; 
	char *q = CFG_MIDI_PRM.s; t=q[i],q[i]=q[j],q[j]=t; 		  mi_tr_p2l[I]=j; mi_tr_p2l[J]=i; }
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
	if (DBGC|midi_dispflg) log("midi: dev%02d(%02d) ch%02d ky%03d %03d=>%03d",mi_tr_p2l[i],i,ch,k,*p,v);
	if (!x25) return (void) (*p = v);  else *p = x25|v;
	int ec = wrap_midi_ev(x, k, v, p0); if (ec>=0) return;
	if (ec==MDE_KEEPV) return (void) (*p = x);
	if (ec!=EEE_NOEFF) gui_errq_add(ec,"midi/swr");
	else if (DBGC) log("midi_ev: 0x%x -> 0x%x - no effect",x25,v);
}

static void midi_rls_all(int ix) {
	unsigned int *q, **qq = mi_root[ix]; if (qq==mi_dflti) return;
	for (int i=0; i<16; i++) { if ((q=qq[i])!=mi_dfltc) {
		for (int j=0; j<256; j++) q[j] &= 127u;
		if (DBGC) log("rls_all: ixP=%d, ch=%d",ix,i); }}}

static int midi_c_slg(int j, const char *s) {
	if (intv_cmd(midi_vsl_ix+(j&=1),s,128,246,0x400801)) gui_midi(1024<<j);   return 0; }

static int midi_open(int ix, int id) {
	char nmb[256]; nmb[0] = 63; nmb[1] = 0;
	mi_devname[14] = 48+(id>>4); mi_devname[16] = 48+(id&15);
	int l, fd = open(mi_devname, ((debug_flags&DFLG_RWMIDI) ? O_RDWR : O_RDONLY)|O_NONBLOCK|O_CLOEXEC);
	if (fd<0) goto err; else close(fd);
	l = rawmidi_desc(nmb, id, 127) + 1; 
	memcpy(mi_dscl[ix]<l ? (mi_dsc[ix]=nf_alloc(l)) : mi_dsc[ix], nmb, l); mi_dscl[ix] = l;
	mi_devname[14] = 48+(id>>4); mi_devname[16] = 48+(id&15);
	if ((fd = open(mi_devname, ((debug_flags&DFLG_RWMIDI) ? O_RDWR : O_RDONLY)|O_CLOEXEC)) < 0) goto err;
	log("midi%x: %s(%s): opened(%d)", ix, mi_devname , nmb, fd); return fd;
err:	gui_errq_add(EEE_ERRNO); gui_errq_add(MDE_OFAIL); midi_err_bv |= 1<<ix;
	log("midi%x: %s(%s): %s", ix, mi_devname, nmb, strerror(errno)); return -1;
}

static int chk_perm(char *to, const char *s) {
	unsigned int m, bv = 0u;
	for (int j,i=0; i<32; i++) 
		if ((j=b32_to_i(s[i]))<0 || (bv&(m=1u<<j))) return 0; else to[i] = j, bv |= m;
	return !~bv;
}

static void midi_i_devls() {
	static const char *idprm = "0123456789abcdefghijklmnopqrstuv";
	unsigned char buf[32]; int n = 0, n0 = find_dev(buf, 1, 31);
	for (int i=0; i<n0; i++) if ((midi_fd[n]=midi_open(n,buf[i]))>=0) 
		midi_bv|=1u<<n, mi_i_ini(n), midi_dvid[n]=buf[i], n++;
	if (n>1) qsort(mi_tr_l2p, n, 1, &devnm_cmp);
	char po[32], pb[32], *q = CFG_MIDI_PRM.s;
	int k = CFG_MIDI_PRM.i; if (k!=32) { if(k) goto err; else goto mk; }
	if (!memcmp(q, idprm, 32)) goto inv;
	if (!chk_perm(pb, q) || pb[31]!=31) goto err;
	memcpy(po, mi_tr_l2p, 32); for (int i=0; i<32; i++) mi_tr_l2p[i] = po[(int)pb[i]];   goto inv;
err:    log("midi: invalid MIDI_PERM[%d]=\"%s\"", CFG_MIDI_PRM.i, CFG_MIDI_PRM.s);
mk:     memcpy(CFG_MIDI_PRM.s, idprm, 33); CFG_MIDI_PRM.i = 32;
inv:	for (int i=0; i<31; i++) mi_tr_p2l[(int)mi_tr_l2p[i]] = i;
}

//////////////////////////////////////////////////////////////////////////////

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
	if (r<0) return close(midi_fd[i]), midi_fd[i]=-1, midi_bv &= ~(1<<i), midi_err_bv |= (1<<i),
			mi_err(i,1,EEE_ERRNO), gui_midi(32*mi_tr_p2l[i] + 1);
	while (r) {
		switch(x>>4) {
			case 8: v=c=0; goto kc;
			case 9: goto kv;
			case 10: case 14: goto l2;
			case 11: v=p[2]; c=128; goto kc;
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
kv:		CHK2; c=*mi_keyspd; midi_kc(i, x&15, p[1], *mi_keyspd=p[2]); p+=3,r-=3; *mi_keyspd=c; continue;
f0:		if ((v = midi_f0(i,p,r))<0) return mi_err(i,0,v);  p+=v; r-=v; continue;
	}}

int midi_cmd(const char *s) {
	if (!s || !*s) return MDE_PARSE;
	int j, m, c = *(s++), ix = b32_to_i(c);
	if (ix<0) { switch(c) {
		case 'W': return gui_midi(0x4000fc1e), 0;
		case '<': m = 1<<(*s&7); j = m &- (s[1]&1);
			  return gui_midi(set_kflg((CFG_MIDI_KTR.i&~m)|j)), 0;
		case '>': midi_dispflg = *s&1; gui_midi(4096); return 0;
 		case 'x': return midi_c_slg(0, s);
 		case 'y': return midi_c_slg(1, s);
		case 'z': BVFOR_JM(midi_bv|(1u<<31)) midi_rls_all(j);  return 0;
		default: return MDE_PARSE; }}
	switch(*s) {
		case 'o': j = mi_tr_l2p[ix]; m = 1<<j; if (!((midi_bv|midi_err_bv)&m)) return MDE_NODEV;
			  if (midi_bv&m) close(midi_fd[j]), midi_bv &= ~m;
			  if ((midi_fd[j] = midi_open(j, midi_dvid[j]))>=0) midi_bv|=m, midi_err_bv &= ~m;
			  return gui_midi(32*ix+1), 0;
		case 'z': return midi_rls_all(mi_tr_l2p[ix]), 0;
		case 'S': return (s[1]&&s[2]) ? (midi_kc(ix,0,hex2(s+1),atoi_h(s+3)&127),0) : MDE_PARSE;
		case 'K': return j=hex2(s+1), midi_kc(ix,0, keytrans[j&127], 127&(-(j>>7))), 0;
		case 'W': return midi_wr(ix, s+1);
		case 'X': if (ix==31) return MDE_VXCHG;
			  j = b32_to_i(s[1]);
			  if (j<0) { if ((s[1]-43)&~2) return MDE_PARSE; else j=ix+44-s[1]; }
			  if (!(~j&31)) return MDE_VXCHG;   if (ix==j) return EEE_NOEFF; 
			  mi_xchg(ix, j); gui_midi(ix<j?31*ix+j+1:31*j+ix+1); return 0;
		default:  return MDE_PARSE;
	}}

int midi_w_flg() { return 8*midi_dispflg + CFG_MIDI_KTR.i; }

int midi_w_ln(char *to, int j0) {
	int j = mi_tr_l2p[j0], id = midi_dvid[j], n = mi_dscl[j]-1, m = (1<<j);
	if (!(midi_bv&m)) return (midi_err_bv&m) ? (*to='!', to[1]=48+(id>>4), to[2]=48+(id&15), 3) 
						 : (*to='-', 1);
	to[0] = 48+(id>>4); to[1] = 48+(id&15); memcpy(to+2, mi_dsc[j], n); return n+2; }

int midi_w_slg(char *to, int j) {
	int buf[6], k0 = midi_vsl_ix[j];
	unsigned int *p = mi_root[31][0];
	buf[0] = qh4(256*k0 + (p[255-j]&127)); p += k0;
	for (int i=0; i<4; i++) buf[i+1] = qh4( 256*(p[2*i]&127) + (p[2*i+1]&127) );
	memcpy(to, buf, 20); return 20;
}

void midi_init() {
	for (int i=0; i<32; i++) mi_tr_l2p[i] = i;   memcpy(mi_tr_p2l, mi_tr_l2p, 32);
	mi_i_ini(31); *(mi_keyspd = mi_rw(31, 0, 255)) = 127;  
	midi_i_devls();
	for (int i=0; i<128; i++) keytrans[i] = i;   set_kflg(-1);
}

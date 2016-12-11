#include <alsa/asoundlib.h>

#define CHK1(X,S) if ((e=(X))<0) { es=(S); goto err2; }
inline int bitfill(int x) { return x|=(x>>16), x|=(x>>8), x|=(x>>4), x|=(x>>2), x|(x>>1); }
static const char* au_op_cfg(snd_pcm_t ** pp_h, const char *nm, int *bufsiz, int *nchan, int *srate, int flg) {
        snd_pcm_hw_params_t *hwpar = 0;
        snd_pcm_uframes_t bs = *bufsiz;
        unsigned int rate = *srate, nc = *nchan; 
        const char *es = "BUG!!!";
        int x, e = snd_pcm_open(pp_h, nm, SND_PCM_STREAM_PLAYBACK, 0); if (e<0) { es="0open"; goto err1; }
        CHK1(snd_pcm_hw_params_malloc(&hwpar), "1allocpar");
        CHK1(snd_pcm_hw_params_any(*pp_h,hwpar), "2initpar");
        CHK1(snd_pcm_hw_params_set_access(*pp_h,hwpar,SND_PCM_ACCESS_RW_INTERLEAVED), "3access");
        CHK1(snd_pcm_hw_params_set_format (*pp_h,hwpar,SND_PCM_FORMAT_S16_LE), "4format");
        CHK1(snd_pcm_hw_params_set_rate_near (*pp_h,hwpar,&rate,0), "5rate");
        CHK1(snd_pcm_hw_params_set_channels_min (*pp_h, hwpar, &nc), "6#chan");
        CHK1(snd_pcm_hw_params_set_buffer_size_min(*pp_h, hwpar, &bs), "7bufsiz/m");
        if (snd_pcm_hw_params_set_buffer_size_near(*pp_h, hwpar, &bs) < 0) {
                int bs2 = bitfill(bs-1)+1;
                while ((bs=bs2)<4444 && snd_pcm_hw_params_set_buffer_size_near(*pp_h, hwpar, &bs)<0) bs2<<=1; }
        CHK1(snd_pcm_hw_params(*pp_h, hwpar),"setpar");
        CHK1(snd_pcm_hw_params_get_channels(hwpar, &nc), "get_nchan");
        CHK1(snd_pcm_hw_params_get_rate(hwpar, &rate, &x), "get_rate");
        CHK1(snd_pcm_hw_params_get_buffer_size(hwpar, &bs), "get_bufsiz");
        snd_pcm_hw_params_free(hwpar); hwpar = 0;
        CHK1(snd_pcm_nonblock(*pp_h,flg&1),"nonblock");
        CHK1(snd_pcm_prepare(*pp_h),"prepare");
        return *bufsiz = bs, *nchan = nc, *srate = rate, (char*)0;
err2:   snd_pcm_close(*pp_h); *pp_h = 0; if (hwpar) snd_pcm_hw_params_free(hwpar);
err1:   *bufsiz = e; return es;
}

#ifndef __cplusplus
#define FP2(S,X,Y)      fprintf(stderr, S "\n", (X),(Y))
#define FP4(S,X,Y,Z,W)  fprintf(stderr, S "\n", (X),(Y),(Z),(W))
static inline int hxd2i(int c) { return 15&(c+9*(c>>6)); }
static int atoi_h(const char *s) { int r=0; while(*s&80) r=16*r+hxd2i(*(s++)); return r; }
static int report(int x) { int v=(x<<16)+0xa307052; return write(1, &v, 4), x; }

int main(int ac, char** av) {
        if (ac!=5) return FP2("%s: ac=%d (exp.5 (dev flg/nch nf bufsz), BUG)", *av, ac), 1;
        snd_pcm_t * hnd;
        int nc0=atoi_h(av[2]), nf=atoi_h(av[3]), hbsiz = atoi_h(av[4]),
            nc = nc0&255, flg = nc0>>8, rate = 44100;
        const char *nm = av[1], *es = au_op_cfg(&hnd, nm, &hbsiz, &nc, &rate, flg);
        if (es) return FP2("error: %s: %s", es, snd_strerror(hbsiz)), report(9), 1;
        FP4("%s started, #ch=%d, hwbs=%d, rate=%d", *av, nc, hbsiz, rate);
        int k, fsiz = 2*nc, asiz = fsiz*nf, errtemp = 0;
        char buf[asiz]; memset(buf, 0, asiz);
        report(16+nc);
        for(;;){k=0;do{ int r2,r=snd_pcm_writei(hnd,buf+k*fsiz,nf-k);
			if (r>0) { k+=r; errtemp -= !!errtemp; continue; }
                        if (!r) return report(7);
                        if (errtemp+=16>256) return report(6);
                        if ((r2=snd_pcm_recover(hnd, r, 1))<0) return report(8); else report(3); } while(k<nf);
                k=0;do{ int r=read(0,buf+k,asiz-k); if (r<=0) return report(4+!!r); k+=r; } while(k<asiz);
                report(1);
        }}
#else // cplusplus

#include "cfgtab.inc"
#include "util.h"
#include "mx.h"
#include "guistub.h"
#include "util2.h"
#include "asnd.h"
#include "pt.h"
#include "box0.h"

#define IFDBG IFDBGX(AUDIO)
#define CDEBUG(S) IFDBG debug(S)
#define SAMP2NS(X) (((X)*m_ns_16f+8)>>4)
#define ADDSAMP(X) (m_clk_nbuf += SAMP2NS(X))
static fa_writer fa_wr;

int ASnd::close() { clk0.cls(); p_close(&m_pump_ofd); if (m_hnd) snd_pcm_close(m_hnd);
		    return m_hnd=0, m_cur_cpf = &cpf_mute, m_n_chan=0, m_pump_st&=~3, w(1024), 0; }
int ASnd::err(int k, const char *s, int ec) { log("alsa/%s: %s(%d)", s, snd_strerror(k), k); close();
					      if (ec) gui2.errq_add(ec); return k; }

#define FUNCHK(N) if (m_cur_cpf==&cpf_##N) fnm = #N; else
void ASnd::debug(const char *s) {
	const char *fnm; FUNCHK(mute) FUNCHK(true) FUNCHK(liar) FUNCHK(pump) fnm = "???";
	log("au/%s: fun=%s, pump_cmd=%d, ahnd=%p flg=0x%x pump_st=0x%x", 
			s?s:"debug", fnm, *m_pump_pcp, m_hnd, m_flg, m_pump_st); }

int ASnd::start(int flg, int mxid, int *ppcp) {
	int n = CFG_AU_TRY_N.i;  if (ppcp) m_pump_pcp = ppcp, m_pump_et = 0;
	if (flg&1) n = (n+1)>>1; else if (mxid>=0) m_mxid = mxid;
	cfg_pre(CFG_AU_SPD.i, CFG_AU_RSRV.i, CFG_AU_CMODE.i);
	if (m_flg&2) return pump_launch(n);
	double clf = (n<2) ? 1.0 : exp(M_LN10 / (double)(n-1)), sclim = 1e3*(double)CFG_AU_CLKLIM.i;
	int e=-1, us = 1000 * CFG_AU_TRY_MS.i;   m_flg &= ~2;
	
	for (int i=0; i<n; u_sleep(us), i++, sclim*=clf)
		if ((e=start1(min_i((1<<27)-1, (int)lround(sclim))))>=0) return w(1024), e;
	return err(e, "start", AUE_START), w(-1), 0;
}

void ASnd::cfg_pre(int spd, int rsrv, int mode) {
	int bs  = m_bs  = (int)lround(2048.0 * ipow(0.965936328924846, spd)),
	    bs2 = m_bs2 = (bs*rsrv+50)/100;
	clk0.bcfg(sample_rate, bs, bs2, (3*bs)>>1);
	m_flg &= ~3; m_flg |= mode;
	if (mode) m_hwbs_trg = bs+bs2, m_cur_cpf = (mode&2) ? &cpf_mute : &cpf_liar;
	else 	  m_hwbs_trg = bitfill(1024|bs*(500+3*rsrv)/200) + 1, m_cur_cpf = &cpf_true;
}

int ASnd::pump_launch(int nt) {
	int aa[5]; aa[0]=qh4(CFG_AU_CHCFG.i)&65535; aa[1]=qh4(m_bs); aa[3]=qh4(m_bs+2*m_bs2); aa[2]=aa[4]=0;
	int pid = launch(QENV('q'), "!><k", &m_pump_ofd, m_pump_pcp,
			 CFG_AU_NAME.s, (char*)aa, (char*)(aa+1), (char*)(aa+3), (char*)0);
	IFDBG log("pump/cp=%d, args: %s %s %s %s", *m_pump_pcp, CFG_AU_NAME.s, (char*)aa, (char*)aa+1, (char*)aa+3);
	return w(1024), (pid<0) ? EEE_ERRNO : (m_pump_st = 4 + (nt ? 256*nt : m_pump_st&0xff00), 0); }

int ASnd::pump_op(int x) {
	if (debug_flags&DFLG_AUDIO&-(x!=1)) log_n("pump_op %d ", x), debug("p/o");
	if (x==1) return (m_pump_st&1) ? AUE_P_STATE : (++m_pump_st,clk0.pump_1(),0);
	int e,k; switch(x) {
		case 0: e = AUE_P_CRASH; goto down;
		case 1: m_pump_st |= 1; clk0.pump_1(); return 0;
		case 4: e = 0;		 goto down0;
		case 5: case 6: case 7: case 8: case 9: e = AUE_P_ERREX; goto down;
		case 3: return AUE_P_RECOV;
		default: break; }
	if ((unsigned int)(x-17)>15u) return AUE_P_WHAT;
	if ((m_pump_st&7)!=4) return AUE_P_STATE; else m_pump_st = 3, clk0.pump_1();
	mx_au16_cfg(&m_cfg, m_n_chan=x-16, CFG_AU_CHCFG.s); m_cur_cpf = &cpf_pump;
	log("pump process started, #chan=%d", m_n_chan); clk0.pump_cfg(1);
	return w(1024), 0;
down:	IFDBG log("pump/dn: et=%d", m_pump_et); 
	if ((m_pump_et+=32) < 999) gui_errq_add(AUE_P_RE), m_pump_st |= 8; else m_pump_et = 0;
down0:	p_close(m_pump_pcp), p_close(&m_pump_ofd); k = m_pump_st&8; m_pump_st &= 0xff00;
	clk0.pump_cfg(0); m_cur_cpf = &cpf_mute;
	e = k ? start() : (m_pump_st ? (m_pump_st-=256, pump_launch()) : e); return w(1024), e;
}

int ASnd::start_buf(int tlim) {
	int wr0=64, lr=m_flg&1; // TODO: config
	if (m_bufsiz<m_hwbs_trg) return log("alsa/bsiz: ret(%d) < rq(%d)", m_bufsiz, m_hwbs_trg), close(), -1;
	if(lr){ m_bs2 = 2*(m_bufsiz-m_bs);
		clk0.bcfg(sample_rate, m_bs, m_bs2, (3*m_bs)>>1); m_cur_cpf = &cpf_liar;
		int sb1 = 2 * m_bufsiz * m_cfg.nch;
		if (sb1>m_swb_siz) { if(m_swb_siz)free(m_swb);  m_swb = (short*)malloc(m_swb_siz=sb1); }
		memset(m_swb, 0, sb1); clk0.set_f(m_bufsiz+wr0); snd_pcm_writei(m_hnd,m_swb,wr0);
		return 0; }
	else  { clk0.bcfg(sample_rate, m_bs, m_bs2, (3*m_bs)>>1);  m_cur_cpf = &cpf_true;
		return play_and_adj((short*)zeroblkD, wr0, tlim|(1<<30)); }}

int ASnd::start1(int sc_lim) {
	if (m_hnd) close(); clk0.ev('o');
	if (CFG_AU_KILL_PA.i) gui2.errq_add(pt_kill_pa(CFG_AU_KILL_PA.i));
	int rate = 44100, chan = CFG_AU_CHCFG.i, bsiz = m_hwbs_trg;
	const char *es = au_op_cfg(&m_hnd, CFG_AU_NAME.s, &bsiz, &chan, &rate, m_flg&1);
	if (es) return log("snd/start1: %s failed, %s", es+1, snd_strerror(bsiz)), bsiz;
	m_n_chan = chan; sample_rate = rate; *clk0.pa() = m_bufsiz = bsiz;
	mx_au16_cfg(&m_cfg, m_n_chan, CFG_AU_CHCFG.s);
	return start_buf(sc_lim);
}

int ASnd::hcp_start(int t) { return m_hcp ? JQE_DUP : ((fa_start(&fa_wr, 2)<0) ? EEE_A20 : 
		(m_hcp = t, m_hcp_s0 = 0x7fffffff, memcpy(m_hcp_lbl, "R --:--", 8), 0)); }
int ASnd::hcp_end(int f) { return !(f|m_hcp) ? EEE_NOEFF : (m_hcp = 0, gui_acv_op(fa_wr.id),
	gui2.setwin(7,'.'), gui2.wupd_s('W',"rec"), fa_end(&fa_wr)<0 ? EEE_ERRNO : 0); }

int ASnd::e_msg_re(int e1, const char *s, int re) {
	if (e1>=0) return e1; if (!re) return log("audio/%s: %s (giving up)", s, snd_strerror(e1)), e1;
	int e2 = snd_pcm_recover(m_hnd, e1, 1);
	return log("audio/%s: %s (R: %s)", s, snd_strerror(e1), e2<0 ? snd_strerror(e2) : "OK"), e2;
}

int ASnd::adj2(int dif, unsigned int *to) { int adj = ivlim(dif/16,-960000,960000);
					    clk0.add(adj); *to = (unsigned int)adj; clk0.gcond(); return 0; }

void ASnd::err3(int flg) { if (flg&  1) gui_errq_add(AUE_LTNOSEE);
			   if (flg&  2) gui_errq_add(AUE_MCFLY);
			   if (flg&508) gui_errq_add(AUE_RECOVER); }

#define TWICE(X,S) ( (ec1=(X))<0 && (++rcnt, e_msg_re(ec1,S,1)<0 || e_msg_re((X),S,0)<0) )
int ASnd::play_and_adj(short *buf, int nf, int opt) {
	int wr,ec1, av0=999999,av1=999999, t2,t1,t0=clk0.ev('a'), rcnt=0, cs=opt&(1<<30), oldce=clk0.err();
	if ((!cs && TWICE(av0=snd_pcm_avail(m_hnd),"avail1")) 
		 || TWICE((t1=clk0.ev('w'), wr=snd_pcm_writei(m_hnd,buf,nf)),"writei")
		 || TWICE((t2=clk0.ev('v'), av1=snd_pcm_avail(m_hnd)),"avail2")) return AUE_UNRECOV;
	int ce=clk0.err(), clke = ce|oldce, err2 = clke+4*rcnt;
	if (err2|cs) {  if (err2) err3(err2);
			return (cs && ((ce|rcnt) || (t0-t1)>(opt&0xfffffff))) ? AUE_CLOCK
				: (clk0.set_f(m_bufsiz - av1), 0); }
	if (wr<nf) log("alsa/wr: %d/%d written", wr, nf), gui2.errq_add(wr?AUE_LESSWR:AUE_ZEROWR); // TODO
	int t3 = clk0.ev('p'), /*tea = t0-t1,*/ teb = t2-t3, tb1 = clk0.f2ns(m_bufsiz-av1), tb0 = tb1-teb,
	    taz = clk0.f2ns(m_bufsiz-av0+wr)+t3, ta0 = taz-t0, ta1 = taz-t1;
	unsigned int *paa = clk0.pa();
	*clk0.pa(-1) = (unsigned int) (m_bufsiz-av1); 
	*clk0.pa(-2) = (unsigned int) (m_bufsiz-av0); 
	int t4 = clk0.add_f(wr), avg = (ta0+ta1+tb0+tb1+2)>>2; return adj2(avg-t4, paa);
}

void ASnd::play2(short *buf, int nf) {
	int ec = mx_calc_int(m_mxid, buf, &m_cfg, m_hcp?&fa_wr : 0, nf);
	m_total_played += nf; if (ec<0) gui_errq_add(ec);
	if (m_hcp) { if (ec==MXE_HCPFAIL || (m_hcp-=nf)<=0) { if ((ec=hcp_end(1))<0) gui_errq_add(ec); }
		     else if (m_hcp<m_hcp_s0) { int t = m_hcp/44100; m_hcp_s0 = 44100*t;
			     			d99(m_hcp_lbl+2, t/60), d59(m_hcp_lbl+5, t%60),
						gui2.setwin(7,'.'), gui2.wupd_s('W', m_hcp_lbl); }}}

void ASnd::cpf_mute(ASnd *p) {
	int nf = clk0.f2play(1); if (!nf) return; else p->play2((short*)junkbufC, nf);
	glob_flg|=GLF_SILENCE, clk0.add_f(nf), clk0.ev('q'), clk0.gcond(); }

void ASnd::cpf_pump(ASnd *p) {
	if (!(p->m_pump_st&1)) return clk0.pump_n(); else p->m_pump_st &= ~1, p->m_pump_et -= !!p->m_pump_et;
	int nf = p->m_bs, bsc = nf*p->m_cfg.nch;  clk0.ev2('P', nf); clk0.ev2('w');
	short buf[bsc]; p->play2(buf, nf); write(p->m_pump_ofd, buf, 2*bsc); clk0.pump_y(); }

void ASnd::cpf_true(ASnd *p) {
	int nf = clk0.f2play(2); if (!nf) return;
	short buf[nf*p->m_cfg.nch]; p->play2(buf, nf);
	int ec = p->play_and_adj(buf, nf, 0); if (ec<0)  gui_errq_add(ec), p->start(1); }

void ASnd::cpf_liar(ASnd *p) {
	int nf0 = clk0.f2play(0); if (!nf0) return;
	int nch = p->m_cfg.nch, bs3 = p->m_bufsiz, t0 = clk0.t();
	short *buf = p->m_swb;
	int wr = snd_pcm_writei(p->m_hnd,buf,bs3), re = (wr<0);
	if (re && (p->e_msg_re(wr,"writei",1)<0 || 
		   (wr=p->e_msg_re(snd_pcm_writei(p->m_hnd,buf,bs3),"writei",0))<0 ))
		return wr==-EAGAIN ? (void)(clk0.ev2('f'), clk0.ev('p'), p->adj2(clk0.f2ns(nf0>>1), clk0.pa()))
				   : (void)(log("liar: nf0=%d",nf0), gui_errq_add(AUE_UNRECOV), p->start(1));
	int t1 = clk0.ev2('f', wr), rs = bs3-wr, rs2 = nch*rs, tx = clk0.f2ns(bs3+rs);
	if ((re*=4, re|=clk0.err())) p->err3(re);
	memmove(buf, buf+nch*wr, 2*rs2); p->play2(buf+rs2, wr); 
	clk0.add_f(wr); clk0.ev('p'); p->adj2(tx-((t0+t1)>>1), clk0.pa());
}

int ASnd::w(int flg) {
	if (flg&   1) gui2.cre(0x67, 'S'); else gui2.setwin(0x67, 'S');
	if (flg&   2) gui2.wupd_i2('s', CFG_AU_SPD.i);
	if (flg&   4) gui2.wupd_i2('r', CFG_AU_RSRV.i);
	if (flg&   8) gui2.wupd_s('n', CFG_AU_NAME.s);
	if (flg&  16) gui2.wupd_s('o', CFG_AU_CHCFG.s);
	if (flg&  32) gui2.wupd_c48('c', CFG_AU_CMODE.i);
	if (flg&  64) gui2.wupd_i2('t', CFG_AU_TRY_N.i);
	if (flg& 128) gui2.wupd_i2('w', CFG_AU_TRY_MS.i);
	if (flg& 256) gui2.wupd_i1('0', CFG_AU_KILL_PA.i&1);
	if (flg& 512) gui2.wupd_i1('1', CFG_AU_KILL_PA.i>>1);
	if (flg&1024) { 
		char buf[8]; memcpy(buf, "#out: 0", 8); 
		if (m_n_chan>9) buf[5]='1', buf[6] = 38+m_n_chan; else buf[6] += m_n_chan;
		gui2.wupd_s('#', buf);
		if (m_flg&2) {  gui2.wupd_s('C', "%%%XXX(off)\0%%%zz%pump \0%%%zz%pu...\0zz%z%%BUG!\0 "
						 + 6*(m_pump_st&6)); }
		else if (!m_hnd) { gui2.wupd_s('C', "%%%XXX(no audio output)"); }
		else {  snd_pcm_info_t * info;
			if (snd_pcm_info_malloc(&info)<0 || snd_pcm_info(m_hnd, info)<0) {
				gui2.wupd_s('C', "zz%z%%(getinfo failed)"); }
			else {  char buf[24]; const char *s = snd_pcm_info_get_name(info);
				int l = strlen(s); 
				if (l<21) memcpy(buf, s, l+1);
				else 	  memcpy(buf, s, 20), memcpy(buf+20, "...", 4);
				gui2.wupd_s('C', "%%%zz%"); gui2.sz(buf);
			}}}
	return 0;
}

int ASnd::cmd(const char *s) { CDEBUG(s); switch(*s) {
	case 'Y': if (is_up()) close(); else start(); CDEBUG("Y1");  return 0;
	case 'R': if (m_flg&2) return (m_pump_st&2) ? m_pump_st|=8, close() : start();
		  if (m_hnd) close(); return start();
	case 'W': return w(-1);
	case 's': cfg_setint(&CFG_AU_SPD,   atoi_h(s+1)); return 0; 
	case 'r': cfg_setint(&CFG_AU_RSRV,  atoi_h(s+1)); return 0; 
	case 'c': CFG_AU_CMODE.i = s[1]&3; return w(32);
	case 't': cfg_setint(&CFG_AU_TRY_N, atoi_h(s+1)); return 0; 
	case 'w': cfg_setint(&CFG_AU_TRY_MS,atoi_h(s+1)); return 0;
	case '0': CFG_AU_KILL_PA.i = 3*(s[1]&1); return w(768);
	case '1': return (CFG_AU_KILL_PA.i) ? (CFG_AU_KILL_PA.i=1+2*((s[1]&1)), w(512)) : EEE_NOEFF; 
	case 'n': cfg_setstr(&CFG_AU_NAME,  s+1); return w(8);
	case 'o': cfg_setstr(&CFG_AU_CHCFG, s+1); return w(16);
	case 'N': cfg_setstr(&CFG_AU_NAME,  s+1); return 0;
	case 'O': cfg_setstr(&CFG_AU_CHCFG, s+1); return 0;
	case '?': debug(); return 0;
	default: return GCE_PARSE; }}

int rawmidi_desc(char *to, int id, int maxlen) {
	static char nm[8]; if (!nm[0]) memcpy(nm, "hw:x,y\0", 8);
	nm[3] = 48 + (id>>4); nm[5] = 48 + (id&15);
	snd_rawmidi_t * dev = 0;
	snd_rawmidi_info_t * info = 0;
	int len; const char * s;
	if (snd_rawmidi_open(&dev, 0, nm, 0)    < 0) len= 9, s="open fail";
	else if (snd_rawmidi_info_malloc(&info) < 0) len=11, s="malloc fail";
	else if (snd_rawmidi_info (dev, info)   < 0) len=12, s="getinfo fail";
	else s = snd_rawmidi_info_get_name(info),  len = min_i(strlen(s), maxlen-1);
	memcpy(to, s, len); to[len] = 0;
	if (info) snd_rawmidi_info_free(info);
	if (dev)  snd_rawmidi_close(dev);
	return len;
}
#endif // cplusplus

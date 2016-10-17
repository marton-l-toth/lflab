#include <alsa/asoundlib.h>

#include "cfgtab.inc"
#include "util.h"
#include "mx.h"
#include "guistub.h"
#include "util2.h"
#include "asnd.h"
#include "pt.h"

#ifdef DUMMY
#include "t/dummybox.h"
#else
#include "box0.h"
#endif

#define IFDBG if (debug_flags&DFLG_AUDIO)
#define SAMP2NS(X) (((X)*m_ns_16f+8)>>4)
#define ADDSAMP(X) (m_clk_nbuf += SAMP2NS(X))
static fa_writer fa_wr;

int ASnd::close() { clk0.cls(); if (m_hnd) snd_pcm_close(m_hnd); 
		    return w(1024), m_hnd=0, m_cur_cpf = &cpf_mute, m_n_chan=0; }
int ASnd::err(int k, const char *s, int ec) { log("alsa/%s: %s(%d)", s, snd_strerror(k), k); close();
					      if (ec) gui2.errq_add(ec); return k; }
int ASnd::start(int flg, int mxid) {
	int n = CFG_AU_TRY_N.i;
	if (flg&1) n = (n+1)>>1; else if (mxid>=0) m_mxid = mxid;
	double clf = (n<2) ? 1.0 : exp(M_LN10 / (double)(n-1)), sclim = 1e3*(double)CFG_AU_CLKLIM.i;
	int e=-1, us = 1000 * CFG_AU_TRY_MS.i;   m_flg &= ~2;
	for (int i=0; i<n; u_sleep(us), i++, sclim*=clf)
		if ((e=start1(min_i((1<<27)-1, (int)lround(sclim))))>=0) return w(1024), e;
	return err(e, "start", AUE_START), w(-1), AUE_START;
}

void ASnd::cfg_pre(int spd, int rsrv, int liar) {
	int bs  = m_bs  = (int)lround(2048.0 * ipow(0.965936328924846, spd)),
	    bs2 = m_bs2 = (bs*rsrv+50)/100;
	m_flg &= ~1; m_flg |= liar;
	if (liar) { m_cur_cpf = &cpf_liar, m_hwbs_trg = bs+bs2; }
	else {	int n = 4096, l = bs*(500+3*rsrv)/200 + 1;
		while (n<l) n+=n; m_cur_cpf = &cpf_true; m_hwbs_trg = n; }}

int ASnd::start_buf(int tlim) {
	int wr0=64, lr=m_flg&1; // TODO: config
	if (m_bufsiz<m_hwbs_trg) return log("alsa/bsiz: ret(%d) < rq(%d)", m_bufsiz, m_hwbs_trg), close(), -1;
	if(lr){ m_bs2 = 2*(m_bufsiz-m_bs);
		clk0.bcfg(sample_rate, m_bs, m_bs2, (3*m_bs)>>1); m_cur_cpf = &cpf_liar;
		int sb1 = 2 * m_bufsiz * m_cfg.nch;
		if (sb1>m_swb_siz) { if(m_swb_siz)free(m_swb);  m_swb = (short*)malloc(m_swb_siz=sb1); }
		memset(m_swb, 0, sb1); clk0.set_f(2*m_bufsiz); u_sleep(10000); return 0; }
	else  { clk0.bcfg(sample_rate, m_bs, m_bs2, (3*m_bs)>>1);  m_cur_cpf = &cpf_true;
		return play_and_adj((short*)zeroblkD, wr0, tlim|(1<<30)); }}

int ASnd::start1(int sc_lim) {
	if (m_hnd) close(); clk0.ev('o');
	if (!(m_flg&2)) cfg_pre(CFG_AU_SPD.i, CFG_AU_RSRV.i, CFG_AU_LIAR.i);
	if (CFG_AU_KILL_PA.i) gui2.errq_add(pt_kill_pa(CFG_AU_KILL_PA.i));
        int qw, e = snd_pcm_open(&m_hnd, CFG_AU_NAME.s, SND_PCM_STREAM_PLAYBACK, 0);
        if (e<0) return err(e, CFG_AU_NAME.s); 
        snd_pcm_hw_params_t *hwpar;
        if ((e=snd_pcm_hw_params_malloc(&hwpar))<0) return err(e, "allocpar");
        if ((e=snd_pcm_hw_params_any(m_hnd,hwpar))<0) return err(e, "initpar");
        if ((e=snd_pcm_hw_params_set_access(m_hnd,hwpar,SND_PCM_ACCESS_RW_INTERLEAVED))<0) return err(e, "access");
        if ((e=snd_pcm_hw_params_set_format (m_hnd,hwpar,SND_PCM_FORMAT_S16_LE))<0) return err(e, "format");
        unsigned int rate = 44100, chan = CFG_AU_CHCFG.i; snd_pcm_uframes_t bsiz = m_hwbs_trg;
        if ((e=snd_pcm_hw_params_set_rate_near (m_hnd,hwpar,&rate,0))<0) return err(e, "rate");
	if ((e=snd_pcm_hw_params_set_channels_min (m_hnd, hwpar, &chan))<0) return err(e, "#chan");
        if ((e=snd_pcm_hw_params_set_buffer_size_min(m_hnd, hwpar, &bsiz))<0) return err(e, "bufsiz/m");
        if ((e=snd_pcm_hw_params_set_buffer_size_near(m_hnd, hwpar, &bsiz))<0) {
		log("bufsiz/near(%d) failed", (int)bsiz);
		if (m_flg&2) { if (bsiz<8192) bsiz *= 2; }
		else         { m_flg |= 2; for (bsiz=1024; bsiz<m_hwbs_trg; bsiz+=bsiz); }
		m_hwbs_trg = bsiz; return err(e, "bufsiz/n");
	}
        if ((e=snd_pcm_hw_params(m_hnd, hwpar))<0) return err(e,"setpar");
	chan=rate=999888777; bsiz = 0;
	if ((e=snd_pcm_hw_params_get_channels(hwpar, &chan))<0) err(e, "get_nchan");
	if ((e=snd_pcm_hw_params_get_rate(hwpar, &rate, &qw))<0) err(e, "get_rate");
	if (chan>16) return log("snd/init2: #chan=%d, WTF???"), -1;
	if ((e=snd_pcm_hw_params_get_buffer_size(hwpar, &bsiz))<0) err(e, "get_bufsiz");
	if ((int)bsiz < m_hwbs_trg) return log("alsa/bsiz: ret(%d) < rq(%d)", bsiz, m_hwbs_trg), -1;
	m_n_chan=(int)chan, sample_rate=(int)rate; mx_au16_cfg(&m_cfg, m_n_chan, CFG_AU_CHCFG.s);
        snd_pcm_hw_params_free(hwpar); *clk0.pa() = m_bufsiz = bsiz;
	log("snd: open(%s) OK, #chan=%d, rate=%d, bufsiz=%d (trg:%d)", 
		CFG_AU_NAME.s, m_n_chan, sample_rate, m_bufsiz,m_hwbs_trg);
        if ((e=snd_pcm_nonblock(m_hnd,1))<0) return err(e,"nonblock");
        if ((e=snd_pcm_prepare(m_hnd)) < 0) return err(e,"prepare");
	return start_buf(sc_lim);
	return m_cur_cpf = &cpf_true, play_and_adj((short*)zeroblkD, 64, sc_lim|(1<<30));
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

int ASnd::adj2(int dif) { int adj = ivlim(dif/16,-960000,960000);
			  clk0.add(adj); *clk0.pa() = (unsigned int)adj; clk0.gcond(); return 0; }

#define TWICE(X,S) ( (ec1=(X))<0 && (++rcnt, e_msg_re(ec1,S,1)<0 || e_msg_re((X),S,0)<0) )
int ASnd::play_and_adj(short *buf, int nf, int opt) {
	int wr,ec1, av0=999999,av1=999999, t2,t1,t0=clk0.ev('a'), rcnt=0, cs=opt&(1<<30), oldce=clk0.err();
	if ((!cs && TWICE(av0=snd_pcm_avail(m_hnd),"avail1")) 
		 || TWICE((t1=clk0.ev('w'), wr=snd_pcm_writei(m_hnd,buf,nf)),"writei")
		 || TWICE((t2=clk0.ev('v'), av1=snd_pcm_avail(m_hnd)),"avail2")) return AUE_UNRECOV;
	int ce=clk0.err(), clke = ce|oldce, err2 = clke|rcnt;
	if (err2) { if (clke&1) gui_errq_add(AUE_LTNOSEE); if (clke&2) gui_errq_add(AUE_MCFLY);
		    if (rcnt)   gui_errq_add(AUE_RECOVER); }
	if (err2|cs) return (cs && ((ce|rcnt) || (t0-t1)>(opt&0xfffff))) ? AUE_CLOCK
		: (clk0.set_f(m_bufsiz - av1), 0);
	if (wr<nf) log("alsa/wr: %d/%d written", wr, nf), gui2.errq_add(wr?AUE_LESSWR:AUE_ZEROWR); // TODO
	int t3 = clk0.ev('p'), tea = t0-t1, teb = t2-t3, tb1 = clk0.f2ns(m_bufsiz-av1), tb0 = tb1-teb,
	    taz = clk0.f2ns(m_bufsiz-av0+wr)+t3, ta0 = taz-t0, ta1 = taz-t1;
	*clk0.pa(-1) = (unsigned int) (m_bufsiz-av1); 
	*clk0.pa(-2) = (unsigned int) (m_bufsiz-av0); 
	if (clke|rcnt) { if (rcnt) gui_errq_add(AUE_RECOVER);
			 if (clke&1) gui_errq_add(AUE_MCFLY);
			 if (clke&2) gui_errq_add(AUE_LTNOSEE); }
	int t4 = clk0.add_f(wr), avg = (ta0+ta1+tb0+tb1+2)>>2; return adj2(avg-t4);
}

void ASnd::play2(short *buf, int nf) {
	int ec = mx_calc_int(m_mxid, buf, &m_cfg, m_hcp?&fa_wr : 0, nf);
	m_total_played += nf; if (ec<0) gui_errq_add(ec);
	if (m_hcp) { if (ec==MXE_HCPFAIL || (m_hcp-=nf)<=0) { if ((ec=hcp_end(1))<0) gui_errq_add(ec); }
		     else if (m_hcp<m_hcp_s0) { int t = m_hcp/44100; m_hcp_s0 = 44100*t;
			     			d99(m_hcp_lbl+2, ec/60), d59(m_hcp_lbl+5, ec%60),
						gui2.setwin(7,'.'), gui2.wupd_s('W', m_hcp_lbl); }}}

void ASnd::cpf_mute(ASnd *p) {
	int nf = clk0.f2play(1); if (!nf) return; else p->play2((short*)junkbufC, nf);
	glob_flg|=GLF_SILENCE, clk0.add_f(nf), clk0.ev('q'), clk0.gcond(); }

void ASnd::cpf_true(ASnd *p) {
	int nf = clk0.f2play(2); if (!nf) return;
	short buf[nf*p->m_cfg.nch]; p->play2(buf, nf);
	int ec = p->play_and_adj(buf, nf, 0); if (ec<0)  gui_errq_add(ec), p->start(1); }

void ASnd::cpf_liar(ASnd *p) {
	int nf0 = clk0.f2play(0); if (!nf0) return;
	int nch = p->m_cfg.nch, bs3 = p->m_bufsiz, t0 = clk0.t();
	short *buf = p->m_swb;
	int wr = snd_pcm_writei(p->m_hnd,buf,bs3);
	if (wr<0 && (p->e_msg_re(wr,"writei",1)<0 || 
		     p->e_msg_re(wr=snd_pcm_writei(p->m_hnd,buf,bs3),"writei",0)<0 ))
		return (void)(log("liar: nf0=%d", nf0), gui_errq_add(AUE_UNRECOV), p->start(1));
	int t1 = clk0.ev2('f', wr), rs = bs3-wr, rs2 = nch*rs, tx = clk0.f2ns(bs3+rs);
	memmove(buf, buf+nch*wr, 2*rs2); p->play2(buf+rs2, wr); 
	clk0.add_f(wr); clk0.ev('p'); p->adj2(tx-((t0+t1)>>1));
}

int ASnd::w(int flg) {
	if (flg&   1) gui2.cre(0x67, 'S'); else gui2.setwin(0x67, 'S');
	if (flg&   2) gui2.wupd_i2('s', CFG_AU_SPD.i);
	if (flg&   4) gui2.wupd_i2('r', CFG_AU_RSRV.i);
	if (flg&   8) gui2.wupd_s('n', CFG_AU_NAME.s);
	if (flg&  16) gui2.wupd_s('o', CFG_AU_CHCFG.s);
	if (flg&  32) gui2.wupd_i1('c', CFG_AU_LIAR.i);
	if (flg&  64) gui2.wupd_i2('t', CFG_AU_TRY_N.i);
	if (flg& 128) gui2.wupd_i2('w', CFG_AU_TRY_MS.i);
	if (flg& 256) gui2.wupd_i1('0', CFG_AU_KILL_PA.i&1);
	if (flg& 512) gui2.wupd_i1('1', CFG_AU_KILL_PA.i>>1);
	if (flg&1024) { 
		char buf[8]; memcpy(buf, "#out: 0", 8); 
		if (m_n_chan>9) buf[5]='1', buf[6] = 38+m_n_chan; else buf[6] += m_n_chan;
		gui2.wupd_s('#', buf);
		if (!m_hnd) { gui2.wupd_s('C', "%%%XXX(no audio output)"); }
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

int ASnd::cmd(const char *s) { switch(*s) {
	case 'Y': if (m_hnd) close(); else start();   return 0;
	case 'R': if (m_hnd) close();      start();   return 0;
	case 'W': return w(-1);
	case 's': cfg_setint(&CFG_AU_SPD,   atoi_h(s+1)); return 0; 
	case 'r': cfg_setint(&CFG_AU_RSRV,  atoi_h(s+1)); return 0; 
	case 'c': CFG_AU_LIAR.i = s[1]&1; return w(32);
	case 't': cfg_setint(&CFG_AU_TRY_N, atoi_h(s+1)); return 0; 
	case 'w': cfg_setint(&CFG_AU_TRY_MS,atoi_h(s+1)); return 0;
	case '0': CFG_AU_KILL_PA.i = 3*(s[1]&1); return w(768);
	case '1': return (CFG_AU_KILL_PA.i) ? (CFG_AU_KILL_PA.i=1+2*((s[1]&1)), w(512)) : EEE_NOEFF; 
	case 'n': cfg_setstr(&CFG_AU_NAME,  s+1); return w(8);
	case 'o': cfg_setstr(&CFG_AU_CHCFG, s+1); return w(16);
	case 'N': cfg_setstr(&CFG_AU_NAME,  s+1); return 0;
	case 'O': cfg_setstr(&CFG_AU_CHCFG, s+1); return 0;
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

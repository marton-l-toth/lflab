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

int ASnd::cond_clk(int * p, int min_ms) { 
	int cs = m_clk_ts.tv_sec, cu = m_clk_ts.tv_nsec, sd = (cs - p[0]) & 0x7fffffff;
	return (sd>(63<<15) || (1000*sd+(cu-p[1])/1000000 >= min_ms)) && (p[0]=cs, p[1]=cu, 1); }

int ASnd::err(int k, const char *s, int ec) { log("alsa/%s: %s(%d)", s, snd_strerror(k), k); close();
					      if (ec) gui2.errq_add(ec); return k; }

void ASnd::cfg(int tpipe_fd, int mxid) {
	if (tpipe_fd>-2) m_tlogpipe = tpipe_fd;
	if (mxid>0) m_mxid = mxid;
	int bs = (int)lround(2048.0 * ipow(0.965936328924846, CFG_AU_SPD.i)),
	    bs2 = (CFG_AU_RSRV.i * bs + 50) / 100;
	double nspf = 1e9 * sample_length;
	m_nspf   = (int)lround(	      nspf);   m_t_empty = (int)lround(nspf * (double)(  bs+bs2));
	m_ns_16f = (int)lround(16.0 * nspf);   m_t_full  = (int)lround(nspf * (double)(2*bs+bs2)); 
	int f_e = m_t_full - m_t_empty;
	m_t_half = m_t_empty + (f_e>>1);  m_sel_min = m_t_empty + (f_e>>3);
}

int ASnd::start() {
	cfg(-2, -2); m_clk_type = (CFG_AU_CLK.i&1) ? CLOCK_MONOTONIC_RAW : CLOCK_MONOTONIC;
	if (clock_gettime(m_clk_type, &m_clk_ts)<0) return gui2.errq_add(EEE_ERRNO), AUE_CLOCK;
	if (m_ev_arg) tlog(m_ev_arg, 256000000); m_ev_arg = 'o';
	return try_start(CFG_AU_TRY_N.i);
}

int ASnd::try_start(int n) {
	double clf = (n<2) ? 1.0 : exp(M_LN10 / (double)(n-1)), sclim = 1e3*(double)CFG_AU_CLK.i;
	int e=-1, us = 1000 * CFG_AU_TRY_MS.i;
	for (int i=0; i<n; u_sleep(us), i++, sclim*=clf)
		if ((e=start1(min_i((1<<27)-1, (int)lround(sclim))))>=0) return m_state = 1, w(1024), e;
	return err(e, "start", AUE_START), w(-1), AUE_START;
}

int ASnd::start1(int sc_lim) {
	if (m_state||m_hnd) close();
	if (CFG_AU_KILL_PA.i) gui2.errq_add(pt_kill_pa(CFG_AU_KILL_PA.i));
        int qw, e = snd_pcm_open(&m_hnd, CFG_AU_NAME.s, SND_PCM_STREAM_PLAYBACK, 0);
        if (e<0) return err(e, CFG_AU_NAME.s); 
        snd_pcm_hw_params_t *hwpar;
        if ((e=snd_pcm_hw_params_malloc(&hwpar))<0) return err(e, "allocpar");
        if ((e=snd_pcm_hw_params_any(m_hnd,hwpar))<0) return err(e, "initpar");
        if ((e=snd_pcm_hw_params_set_access(m_hnd,hwpar,SND_PCM_ACCESS_RW_INTERLEAVED))<0) return err(e, "access");
        if ((e=snd_pcm_hw_params_set_format (m_hnd,hwpar,SND_PCM_FORMAT_S16_LE))<0) return err(e, "format");
        unsigned int rate = 44100, chan = CFG_AU_CHCFG.i; snd_pcm_uframes_t bsiz = 8192;
        if ((e=snd_pcm_hw_params_set_rate_near (m_hnd,hwpar,&rate,0))<0) return err(e, "rate");
	if ((e=snd_pcm_hw_params_set_channels_min (m_hnd, hwpar, &chan))<0) return err(e, "#chan");
        if ((e=snd_pcm_hw_params_set_buffer_size_min(m_hnd, hwpar, &bsiz))<0) return err(e, "bufsiz/m");
        if ((e=snd_pcm_hw_params(m_hnd, hwpar))<0) return err(e,"setpar");
	chan=rate=999888777;
	if ((e=snd_pcm_hw_params_get_channels(hwpar, &chan))<0) err(e, "get_nchan");
	if ((e=snd_pcm_hw_params_get_rate(hwpar, &rate, &qw))<0) err(e, "get_rate");
	if (chan>16) return log("snd/init2: #chan=%d, WTF???"), -1;
	log("snd: open(%s) OK, #chan=%d, rate=%d", CFG_AU_NAME.s, m_n_chan=(int)chan, sample_rate=(int)rate);
	mx_au16_cfg(&m_cfg, m_n_chan, CFG_AU_CHCFG.s);
	if ((e=snd_pcm_hw_params_get_buffer_size(hwpar, &bsiz))<0) err(e, "get_bufsiz");
	if (bsiz < 8192) return log("alsa/bsiz: ret(%d) < rq(8192)", bsiz), -1;
        snd_pcm_hw_params_free(hwpar); m_bsiz_ref = bsiz;
        if ((e=snd_pcm_nonblock(m_hnd,1))<0) return err(e,"nonblock");
        if ((e=snd_pcm_prepare(m_hnd)) < 0) return err(e,"prepare");
	return play_and_adj((short*)zeroblkD, 64, sc_lim|(1<<30));
}

int ASnd::hcp_start(int t) { return m_hcp ? JQE_DUP : ((fa_start(&fa_wr, 2)<0) ? EEE_A20 : 
		(m_hcp = t, m_hcp_s0 = 0x7fffffff, memcpy(m_hcp_lbl, "R --:--", 8), 0)); }
int ASnd::hcp_end(int f) { return !(f|m_hcp) ? EEE_NOEFF : (m_hcp = 0, gui_acv_op(fa_wr.id),
	gui2.setwin(7,'.'), gui2.wupd_s('W',"rec"), fa_end(&fa_wr)<0 ? EEE_ERRNO : 0); }

int ASnd::close() { if (m_hnd) snd_pcm_close(m_hnd); m_hnd = 0; m_n_chan = 0;
	mark('x'); m_clk_nbuf = m_t_half-1; w(1024); return m_state = 0; }

void ASnd::tlog(unsigned int arg, int t) {
	static int pl_cnt = 0;
	unsigned int t2 = (unsigned int)((t+512)>>10); if (t2>262141u) t2 = 262142u + (t<0); 
	m_l_p[m_l_i] = (arg<<18u) + t2; m_l_i++; if (m_l_i>1023) log("ERROR: dropping tlog"), m_l_i = 0;
	if ((arg|1)=='q' && ++pl_cnt>=4) pl_cnt = 0, write(m_tlogpipe, m_l_p, 4*m_l_i), m_l_i = 0;
}

int ASnd::nclk() {
	struct timespec ts; memcpy(&ts, &m_clk_ts, sizeof(ts));
	clock_gettime(m_clk_type, &m_clk_ts);
	int td = m_clk_ts.tv_nsec - ts.tv_nsec, ds = m_clk_ts.tv_sec - ts.tv_sec;
	int e = ds ? ( (ds==1) ? (td+=1000000000)>0xff00000 : 1+(ds<1) )
		   : ( ((unsigned int)(td)>0xff00000u) ? 1+(td<0) : 0 );
	if (e) m_clk_err |= e, td = m_clk_nbuf&(e-2);
	if (m_ev_arg=='j') {
		if (m_jt_n >= 2000) jt_flush();
		if (!m_jt_n) m_jt_sum = m_jt_sum1 = m_jt_Mi = m_jt_Mv = 0;
		if (td>m_jt_Mv) m_jt_Mi = m_jt_n, m_jt_sum1 = m_jt_sum, m_jt_Mv = td;
		++m_jt_n; m_jt_sum += td; m_ev_arg = '?';
	} else if (m_ev_arg) {
		if (m_jt_n) jt_flush();
		tlog(m_ev_arg, td); m_ev_arg = '?';
	}
 	return m_clk_nbuf -= td;
}

int ASnd::e_msg_re(int e1, const char *s, int re) {
	if (e1>=0) return e1; if (!re) return log("audio/%s: %s (giving up)", s, snd_strerror(e1)), e1;
	int e2 = snd_pcm_recover(m_hnd, e1, 1);
	return log("audio/%s: %s (R: %s)", s, snd_strerror(e1), e2<0 ? snd_strerror(e2) : "OK"), e2;
}

#define TWICE(X,S) ( (ec1=(X))<0 && (++rcnt, e_msg_re(ec1,S,1)<0 || e_msg_re((X),S,0)<0) )
int ASnd::play_and_adj(short *buf, int nf, int opt) {
	int wr,ec1, av0,av1, t2,t1,t0=mark('p'), rcnt=0, cs=opt&(1<<30), oldce=m_clk_err;    m_clk_err = 0;
	if ((!cs && TWICE(av0=snd_pcm_avail(m_hnd),"avail1")) 
		 || TWICE((t1=mark('w'), wr=snd_pcm_writei(m_hnd,buf,nf)),"writei")
		 || TWICE((t2=mark('v'), av1=snd_pcm_avail(m_hnd)),"avail2")) return AUE_UNRECOV;
	if (rcnt|m_clk_err|oldce|cs) return (cs && ((m_clk_err|rcnt) || (t0-t1)>(opt&0xfffff))) ? AUE_CLOCK
		: (m_clk_nbuf = SAMP2NS(m_buftot - av1), m_clk_err = 0);
	int t3 = nclk(), tea = t0-t1, teb = t2-t3, tb1 = SAMP2NS(m_buftot-av1), tb0 = tb1-teb,
	    taz = SAMP2NS(m_buftot-av0+nf)+t3, ta0 = taz-t0, ta1 = taz-t1;
	ADDSAMP(nf);
	if (++m_ca_cnt>99) { IFDBG log("snd/adj100: cur=%d ta0=%d ta1=%d tb0=%d tb1=%d", 
						m_clk_nbuf,ta0,   ta1,	 tb0,	tb1); m_ca_cnt = 0; }
	int avg = (ta0+ta1+tb0+tb1+2)>>2, adj = ivlim((avg-m_clk_nbuf)/16,-960000,960000);
	m_clk_nbuf += adj; m_ev_arg = (2181*512 + adj)>>10;
	return 0;
}

void ASnd::c_play() {
	int t = nclk(); if (!m_hnd && (t<0 || t>2*m_t_full)) return (void) (mark('y'), m_clk_nbuf = m_t_half);
	if (t>m_t_half) return (void) (m_ev_arg = 'P');
	int nf = (m_t_full-t)/m_nspf; if ((unsigned int)nf>4095u) nf = 4095;
	short buf[nf*m_cfg.nch];
	int ec = 0, mxr = mx_calc_int(m_mxid, buf, &m_cfg, m_hcp?&fa_wr : 0, nf);
	if (mxr<0) gui_errq_add(mxr);
	m_ev_arg = 4096 + nf; m_total_played += nf;
	if (m_hcp) {
		if (mxr==MXE_HCPFAIL || (m_hcp -= nf)<=0) ec = hcp_end(1);
		else if (m_hcp<m_hcp_s0) ec = m_hcp / 44100, m_hcp_s0 = ec * 44100,
					 d99(m_hcp_lbl+2, ec/60), d59(m_hcp_lbl+5, ec%60),
					 gui2.setwin(7,'.'), gui2.wupd_s('W', m_hcp_lbl), ec = 0;
	}
	if (!m_state) return glob_flg|=GLF_SILENCE, ADDSAMP(nf), (void)mark('q');
	if ((ec=play_and_adj(buf, nf, 0))<0) gui_errq_add(ec), try_start((CFG_AU_TRY_N.i+1)>>1);
}

void ASnd::jt_flush() {
	if (m_jt_Mi) tlog(2048+m_jt_Mi, m_jt_sum1);
	tlog(2049, m_jt_Mv);
	int k = m_jt_n - m_jt_Mi - 1;
	if (k>0) tlog(2048+k, m_jt_sum - m_jt_Mv - m_jt_sum1);
	m_jt_n = 0;
}

inline int d2s16(short * p, double v) { int x = (int)lround(v); 
	return (unsigned int)(x+32767)>65534u ? (*p = x<0 ? -32767 : 32767, 1) : (*p=x, 0); }

int ASnd::w(int flg) {
	if (flg&   1) gui2.cre(0x67, 'S'); else gui2.setwin(0x67, 'S');
	if (flg&   2) gui2.wupd_i2('s', CFG_AU_SPD.i);
	if (flg&   4) gui2.wupd_i2('r', CFG_AU_RSRV.i);
	if (flg&   8) gui2.wupd_s('n', CFG_AU_NAME.s);
	if (flg&  16) gui2.wupd_s('o', CFG_AU_CHCFG.s);
	if (flg&  32) gui2.wupd_s('c', "avail\0delay\0retBS"+6*CFG_AU_CLK.i);
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
	case 'Y': if (m_state) close(); else start();   return 0;
	case 'R': if (m_state) close();      start();   return 0;
	case 'W': return w(-1);
	case 's': cfg_setint(&CFG_AU_SPD,   atoi_h(s+1)); return 0; 
	case 'r': cfg_setint(&CFG_AU_RSRV,  atoi_h(s+1)); return 0; 
	case 'c': cfg_setint(&CFG_AU_CLK,   atoi_h(s+1)); return w(32);
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

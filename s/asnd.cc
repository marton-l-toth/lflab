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
static fa_writer fa_wr;

int ASnd::cond_clk(int * p, int min_us) { int sd = m_clk_sec - p[0];
	return (sd>1800 || (1000000*sd+m_clk_usec-p[1]) >= min_us) && (p[0]=m_clk_sec, p[1]=m_clk_usec, 1); }

int ASnd::err(int k, const char *s, int ec) { log("alsa/%s: %s(%d)", s, snd_strerror(k), k); close();
					      if (ec) gui2.errq_add(ec); return k; }

void ASnd::cfg(int tpipe_fd, int mxid) {
	if (tpipe_fd>-2) m_tlogpipe = tpipe_fd;
	if (mxid>0) m_mxid = mxid;
	int bs = (int)lround(2048.0 * ipow(0.965936328924846, CFG_AU_SPD.i)),
	    bs2 = (CFG_AU_RSRV.i * bs + 50) / 100;
	m_uspf = 1e6 * sample_length;
	m_buf64 = (int)lround(64.0 * m_uspf);
	m_fp32m = ((sample_rate<<13) + 125000) / 250000;
	m_t_empty = (int)lround(m_uspf * (double)(  bs+bs2));
	m_t_full  = (int)lround(m_uspf * (double)(2*bs+bs2));
	m_t_half = (m_t_empty + m_t_full) >> 1;
	m_sel_min = (7*m_t_empty + m_t_full) >> 3;
}

void ASnd::clk_samp(int n) {
	double x = m_clk_frac + m_uspf * (double)n;
	int k = (int)lround(x);
	m_clk_buf += k; m_clk_frac = x - (double)k;
}

int ASnd::start() {
	cfg(-2, -2);
	if (!m_clk_sec) { struct timeval tv; gettimeofday(&tv, 0); 
			  m_clk_sec=tv.tv_sec; m_clk_usec=tv.tv_usec; m_ev_arg='0'; }
	int e=-1, us = 1000 * CFG_AU_TRY_MS.i, sclim = CFG_AU_CLK.i ? 999999 : CFG_AU_CLKLIM.i;
	for (int i=CFG_AU_TRY_N.i; i>0; u_sleep(us), i--, sclim*=5, sclim>>=2)
		if ((e=start1(sclim))>=0) return m_state = 1, w(1024), e;
	return err(e, "start", AUE_START), w(-1), close();
}

int ASnd::start1(int sc_lim) {
	if (CFG_AU_KILL_PA.i) gui2.errq_add(pt_kill_pa(CFG_AU_KILL_PA.i));
        int e = snd_pcm_open(&m_hnd, CFG_AU_NAME.s, SND_PCM_STREAM_PLAYBACK, 0);
        if (e<0) return err(e, CFG_AU_NAME.s); else log("snd: \"%s\" opened OK", CFG_AU_NAME.s);
        snd_pcm_hw_params_t *hwpar;
        if ((e=snd_pcm_hw_params_malloc(&hwpar))<0) return err(e, "allocpar");
        if ((e=snd_pcm_hw_params_any(m_hnd,hwpar))<0) return err(e, "initpar");
        if ((e=snd_pcm_hw_params_set_access(m_hnd,hwpar,SND_PCM_ACCESS_RW_INTERLEAVED))<0) return err(e, "access");
        if ((e=snd_pcm_hw_params_set_format (m_hnd,hwpar,SND_PCM_FORMAT_S16_LE))<0) return err(e, "format");
        unsigned int rate = 44100; snd_pcm_uframes_t bsiz = 8192;
        if ((e=snd_pcm_hw_params_set_rate_near (m_hnd,hwpar,&rate,0))<0) return err(e, "rate");
        log("rate=%d", sample_rate = (int)rate);
	unsigned int chan = strlen(CFG_AU_CHCFG.s);
	if ((e=snd_pcm_hw_params_set_channels_min (m_hnd, hwpar, &chan))<0) return err(e, "#chan");
        if ((e=snd_pcm_hw_params_set_buffer_size_min(m_hnd, hwpar, &bsiz))<0) return err(e, "bufsiz/m");
//        if ((e=snd_pcm_hw_params_set_buffer_size_near(m_hnd, hwpar, &bsiz))<0) return err(e, "bufsiz/n");
        if ((e=snd_pcm_hw_params(m_hnd, hwpar))<0) return err(e,"setpar");
	if ((e=snd_pcm_hw_params_get_channels(hwpar, &chan))<0) err(e, "get_nchan");
	log("#chan=%d",m_n_chan = (int)chan); if (chan>16) return log("snd/init2: #chan=%d, WTF???"), -1;
	mx_au16_cfg(&m_cfg, m_n_chan, CFG_AU_CHCFG.s);
	if ((e=snd_pcm_hw_params_get_buffer_size(hwpar, &bsiz))<0) err(e, "get_bufsiz");
	if (bsiz < 8192) return log("alsa/bsiz: ret(%d) < rq(8192)", bsiz), -1;
        snd_pcm_hw_params_free(hwpar);
        if ((e=snd_pcm_nonblock(m_hnd,1))<0) return err(e,"nonblock");
        if ((e=snd_pcm_prepare(m_hnd)) < 0) return err(e,"prepare");
        return m_bsiz_ref = CFG_AU_CLK.i==2 ? bsiz : -CFG_AU_CLK.i, set_clk(sc_lim);
}

int ASnd::hcp_start(int t) { return m_hcp ? JQE_DUP : ((fa_start(&fa_wr, 2)<0) ? EEE_ERRNO : 
		(m_hcp = t, m_hcp_s0 = 0x7fffffff, memcpy(m_hcp_lbl, "R --:--", 8), 0)); }
int ASnd::hcp_end(int f) { return !(f|m_hcp) ? EEE_NOEFF : (m_hcp = 0, gui2.acv_open(fa_wr.id),
	gui2.setwin(7,'.'), gui2.wupd_s('W',"rec"), fa_end(&fa_wr)<0 ? EEE_ERRNO : 0); }

int ASnd::set_clk(int lim) {
	struct timeval tv0, tv1;
	if (gettimeofday(&tv0, 0)<0) return log("gettimeofday failed"), -1;
	int e1 = snd_pcm_writei(m_hnd, zeroblkC, 64); if (e1<0) return err(e1, "write 64*zero"), -1;
	int av = avail();
	if (gettimeofday(&tv1, 0)<0) return log("gettimeofday failed"), -1;
	m_buftot = m_bsiz_ref ? m_bsiz_ref : av+64;  m_clk_buf = m_buf64;
	m_clk_sec = tv0.tv_sec; m_clk_usec = tv0.tv_usec; m_clk_frac = 0.0;
	m_ca_min = m_ca_max = m_ca_acc = m_ca_cnt = 0; m_ev_arg = 'K';
	int t = 1000000 * (tv1.tv_sec - tv0.tv_sec) + tv1.tv_usec - tv0.tv_usec;
	log("bufsiz = %d (%d), t(wr64+avail) = %d us, clk_buf=%d, %s",
			av + 64, m_bsiz_ref, t, m_clk_buf, t<lim ? "OK":"fail"); 
	return t<lim ? 0 : (close(), -1);
}

void ASnd::adj_clk() {
	int t = clk(), ibuf = m_buftot<0 ? delay(AUE_DLY) : m_buftot - avail(AUE_AVAIL); if (!m_hnd) return;
	int dif = (int)lround(m_uspf * (double)ibuf) - t, dif2 = ivlim((dif+8)>>4, -950, 950);
	if (++m_ca_cnt>99) {
		IFDBG log("snd/adj100: min=%d max=%d acc=%d   ibuf=%d",
			m_ca_min, m_ca_max, m_ca_acc, ibuf);
		m_ca_min = m_ca_max = m_ca_acc = m_ca_cnt = 0; }
	if (dif<m_ca_min) m_ca_min = dif; else if (dif>m_ca_max) m_ca_max = dif;
	m_ca_acc += dif2; m_clk_buf += dif2; m_ev_arg = 1090 + dif2;
}

int ASnd::close() { if (m_hnd) snd_pcm_close(m_hnd); m_hnd = 0; m_n_chan = 0;
	mark('x'); m_clk_buf = m_t_half-1; w(1024); return m_state = 0; }

void ASnd::tlog(unsigned int arg, int t) {
	static int pl_cnt = 0;
	unsigned int t2 = (unsigned int)t; if (t2>262141u) t2 = 262142u + (t<0); 
	m_l_p[m_l_i] = (arg<<18u) + t2; m_l_i++; if (m_l_i>1023) log("ERROR: dropping tlog"), m_l_i = 0;
	if ((arg|1)=='q' && ++pl_cnt>=4) pl_cnt = 0, write(m_tlogpipe, m_l_p, 4*m_l_i), m_l_i = 0;
}

int ASnd::clk() {
	struct timeval tv;
	if (gettimeofday(&tv, 0)<0) log("gettimeofday failed");
	int t = 1000000*(tv.tv_sec - m_clk_sec) + tv.tv_usec - m_clk_usec;
	if (t<0) log("WTF: time just went backwards %d usec", -t), gui_errq_add(AUE_MCFLY), t = 0;
	m_clk_sec = tv.tv_sec; m_clk_usec = tv.tv_usec;
	if (m_ev_arg=='j') {
		if (m_jt_n >= 2000) jt_flush();
		if (!m_jt_n) m_jt_sum = m_jt_sum1 = m_jt_Mi = m_jt_Mv = 0;
		if (t>m_jt_Mv) m_jt_Mi = m_jt_n, m_jt_sum1 = m_jt_sum, m_jt_Mv = t;
		++m_jt_n; m_jt_sum += t; m_ev_arg = '?';
	} else if (m_ev_arg) {
		if (m_jt_n) jt_flush();
		tlog(m_ev_arg, t); m_ev_arg = '?';
	}
	return m_clk_buf -= t;
}

void ASnd::jt_flush() {
	if (m_jt_Mi) tlog(2048+m_jt_Mi, m_jt_sum1);
	tlog(2049, m_jt_Mv);
	int k = m_jt_n - m_jt_Mi - 1;
	if (k>0) tlog(2048+k, m_jt_sum - m_jt_Mv - m_jt_sum1);
	m_jt_n = 0;
}

void ASnd::c_play() {
	int t = clk(); if (t>m_t_half) return (void) (m_ev_arg = 'P');
	int nf = ((m_t_full - t) * m_fp32m) >> 15; if (nf>4095) nf = 4095;
	short buf[nf*m_cfg.nch];
	int ec = 0, mxr = mx_calc_int(m_mxid, buf, &m_cfg, m_hcp?&fa_wr : 0, nf);
	if (mxr<0) gui_errq_add(mxr);
	m_ev_arg = 4096 + nf; mark('q' - m_state);
	if (m_hcp) {
		if (mxr==MXE_HCPFAIL || (m_hcp -= nf)<=0) ec = hcp_end(1);
		else if (m_hcp<m_hcp_s0) ec = m_hcp / 44100, m_hcp_s0 = ec * 44100,
					 d59(m_hcp_lbl+2, ec/60), d59(m_hcp_lbl+5, ec%60),
					 gui2.setwin(7,'.'), gui2.wupd_s('W', m_hcp_lbl), ec = 0;
	}
	int ae = m_state ? snd_pcm_writei(m_hnd, buf, nf) : nf;
	if (ae<0) { if (recover(ae, "play1", AUE_PLAYRCV)) return;
		    ae = snd_pcm_writei(m_hnd, buf, nf);
		    if (ae<0) return err(ae, "play1/2", AUE_WRITEI), (void)close();
	}
	if (ae<nf) log("alsa/wr: %d/%d written", ae, nf), gui2.errq_add(ae?AUE_LESSWR:AUE_ZEROWR); // TODO
	m_total_played += nf; clk_samp(nf); if (m_state) adj_clk();
}

int ASnd::recover(int ec, const char * from, int ec2) {
	const char * s0 = snd_strerror(ec);
	if ((ec=snd_pcm_recover(m_hnd, ec, 1))<0)
		return log("snd: %s: %s (cannot recover from: %s)", from, snd_strerror(ec), s0), close(),
		       gui2.errq_add(ec2), -1;
	log("snd: %s: recovered from %s", from, s0); gui2.errq_add(AUE_RECOVER);
	int e = set_clk(10000); if (e<0) { return err(e, "recover/clk", AUE_RE_CLK), close(), start(), 0; }
	return 0;
}

inline int d2s16(short * p, double v) { int x = (int)lround(v); 
	return (unsigned int)(x+32767)>65534u ? (*p = x<0 ? -32767 : 32767, 1) : (*p=x, 0); }

#define TWICE(NM,FUN) int ASnd::NM(int ec) { \
	if (!m_hnd) return 0xdeadbeef; \
	int ret = (int) FUN (m_hnd); if (ret<0) { \
		if (recover(ret, #NM, ec)) return 0xdeadbeef;  \
		ret = (int) FUN (m_hnd);  \
		if (ret<0) return err(ret, #NM "/2", ec), close(), 0xdeadbeef; \
	} return ret; }

int ASnd::delay2(snd_pcm_t * q) {
	snd_pcm_sframes_t av, dly; int ec;
	return (ec=snd_pcm_avail_delay(q, &av, &dly))<0 ? ec : dly; }
TWICE(avail, snd_pcm_avail)
TWICE(delay, delay2)

int ASnd::stat(unsigned int * to, int n) { return GCE_SORRY; }

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

#ifndef __qwe_asnd_h__
#define __qwe_asnd_h__

struct _snd_pcm; // i know, i know...
typedef struct _snd_pcm snd_pcm_t;

class ASnd {
        public: 
                ASnd() : m_hnd(0), m_ev_arg(0), m_hcp(0) {}
                void cfg(int tpipe_fd, int mxid);
                int start(), close();
		void c_play();
		int time4job() { int r = (nclk() > m_t_half); return m_ev_arg = 'J'+32*r, r; }
		int time4sel(int zero) { int t = nclk(), zf = zero || (t<m_sel_min);
			return (zf) ? (m_ev_arg='S',0) : (m_ev_arg='s',(t-m_t_empty)>>10); }
		inline int mark(int arg14) { int t = nclk(); m_ev_arg = arg14 & 16383u; return t; }
		long long total_played() const { return m_total_played; }
		void set_vol(int x) { m_cfg.vol = x; }
		int vol() const { return m_cfg.vol; }
		int hcp_start(int t), hcp_end(int f = 0);
		int hcp() const { return m_hcp; }
		int cmd(const char *s), w(int flg);
		int cond_clk(int * trg /*int[2]*/, int min_ms = 0x80000000); 
        protected:
		int start1(int sc_lim), try_start(int n);
		void tlog(unsigned int arg, int t);
                int err(int k, const char *s, int ec = 0);
		void jt_flush();
		int e_msg_re(int e1, const char *s, int re);
		int nclk();
		int timer_test(int ty, int n);
		int play_and_adj(short *buf, int nf, int opt);

                snd_pcm_t * m_hnd;
		int m_state, m_n_chan, m_mxid;

		int m_avu_trg, m_sel_min, m_bsiz_ref;
		int m_job_dl;
		int m_nspf, m_ns_16f, m_buftot, m_t_full, m_t_half, m_t_empty;
		int m_ca_cnt, m_ca_min, m_ca_max, m_ca_acc;

		clockid_t m_clk_type;
		struct timespec m_clk_ts;
		int m_ev_arg, m_clk_nbuf, m_clk_err;
		int m_jt_n, m_jt_sum, m_jt_sum1, m_jt_Mi, m_jt_Mv;

		unsigned int m_l_p[1024];
		long long m_total_played;
		int m_l_i, m_f_min, m_f_max;

		int m_tlogpipe, m_hcp, m_hcp_s0;
		char m_hcp_lbl[8];
		au16w_t m_cfg;
};

extern ASnd snd0;

#endif // __qwe_asnd_h__

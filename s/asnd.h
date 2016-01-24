#ifndef __qwe_asnd_h__
#define __qwe_asnd_h__

struct _snd_pcm; // i know, i know...
typedef struct _snd_pcm snd_pcm_t;

class ASnd {
        public: 
                ASnd() : m_hnd(0), m_clk_sec(0), m_hcp(0) {}
                void cfg(int tpipe_fd, int mxid);
                int start();
		void c_play();
                int avail(int ec = 0), delay(int ec = 0), delay2(snd_pcm_t * q);
		int time4job() { int r = (clk() > m_t_half); return m_ev_arg = 'J'+32*r, r; }
		int time4sel(int zero) { int t = clk(), zf = zero || (t<m_sel_min);
			return (zf) ? (m_ev_arg='S',0) : (m_ev_arg='s',t-m_t_empty); }
		void mark(int arg14) { clk(); m_ev_arg = arg14 & 16383u; }
		int t0();
		long long total_played() const { return m_total_played; }
		int stat(unsigned int * to, int n);
		int close();
		void set_vol(int x) { m_cfg.vol = x; }
		int vol() const { return m_cfg.vol; }
		int hcp_start(int t);
		int hcp_end(int f = 0);
		int hcp() const { return m_hcp; }
		int cmd(const char *s);
		int w(int flg);
		int cond_clk(int * trg /*int[2]*/, int min_us); 
        protected:
		int start1(int sc_lim);
		void tlog(unsigned int arg, int t);
                int err(int k, const char *s, int ec = 0);
		void jt_flush();
		int recover(int ec, const char * from, int ec2 = 0);
		int set_clk(int lim);
		void adj_clk();
		int clk();
		void clk_samp(int n);
		int timer_test(int ty, int n);

                snd_pcm_t * m_hnd;
		int m_state, m_n_chan, m_mxid;

		int m_avu_trg, m_sel_min, m_bsiz_ref;
		int m_job_dl;
		int m_buf64, m_nspf, m_fp32m, m_buftot, m_t_full, m_t_half, m_t_empty;
		int m_ca_cnt, m_ca_min, m_ca_max, m_ca_acc;
		double m_uspf, m_clk_frac;
		int m_ev_arg, m_clk_sec, m_clk_usec, m_clk_buf;
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

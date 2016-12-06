#ifndef __qwe_asnd_h__
#define __qwe_asnd_h__

struct _snd_pcm; // i know, i know...
typedef struct _snd_pcm snd_pcm_t;

class ASnd {
	public:
                ASnd() : m_hnd(0), m_hcp(0), m_swb_siz(0) {}
                void cfg(int mxid);
                int start(int qre = 0, int mxid = -1, int *ppcp = 0), close();
		inline void c_play() { (*m_cur_cpf)(this); }
		long long total_played() const { return m_total_played; }
		void set_vol(int x) { m_cfg.vol = x; }
		int vol() const { return m_cfg.vol; }
		int hcp_start(int t), hcp_end(int f = 0);
		int hcp() const { return m_hcp; }
		int cmd(const char *s), w(int flg);
		int pump_op(int x);
        protected:
		static void cpf_mute(ASnd*), cpf_true(ASnd*), cpf_liar(ASnd*), cpf_pump(ASnd*);
		int start1(int sc_lim), try_start(int n);
		void cfg_pre(int spd, int rsrv, int liar);
		int start_buf(int tlim);
		int adj2(int dif, unsigned int *to);
		void err3(int flg);
                int err(int k, const char *s, int ec = 0);
		void play2(short *buf, int nf);
		int e_msg_re(int e1, const char *s, int re);
		int play_and_adj(short *buf, int nf, int opt);

                snd_pcm_t * m_hnd;
		void (*m_cur_cpf)(ASnd*);
		long long m_total_played;
		au16w_t m_cfg;
		char m_hcp_lbl[8];
		int m_bufsiz, m_n_chan, m_mxid, m_hcp, m_hcp_s0;
		int m_bs, m_bs2, m_hwbs_trg, m_swb_siz, m_flg; // flg: 1:liar 2:skipbufnear 4:pump
		int m_pump_st, m_pump_ofd, *m_pump_pcp; // st: 1-rdy2read 2-running 4-wt4up 8:re
		short * m_swb;
};
extern ASnd snd0;

#endif // __qwe_asnd_h__

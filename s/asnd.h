#ifndef __qwe_asnd_h__
#define __qwe_asnd_h__

struct _snd_pcm; // i know, i know...
typedef struct _snd_pcm snd_pcm_t;

class ASnd {
	public:
                ASnd() : m_hnd(0), m_hcp(0) {}
                void cfg(int mxid);
                int start(int qre = 0, int mxid = -1), close();
		inline void c_play() { (*m_cur_cpf)(this); }
		long long total_played() const { return m_total_played; }
		void set_vol(int x) { m_cfg.vol = x; }
		int vol() const { return m_cfg.vol; }
		int hcp_start(int t), hcp_end(int f = 0);
		int hcp() const { return m_hcp; }
		int cmd(const char *s), w(int flg);
        protected:
		static void cpf_mute(ASnd*), cpf_true(ASnd*);
		int start1(int sc_lim), try_start(int n);
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
};
extern ASnd snd0;

#endif // __qwe_asnd_h__

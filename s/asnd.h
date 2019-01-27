#ifndef __qwe_asnd_h__
#define __qwe_asnd_h__

struct _snd_pcm; // i know, i know...
typedef struct _snd_pcm snd_pcm_t;
struct cfg_ent;

#define ASND_KF(NM) inline cfg_ent * k##NM() const { return &CFG_AU_##NM + m_cfg_offs; }
class ASnd {
	public:
                ASnd() : m_hnd(0), m_hcp(0), m_swb_siz(0), m_flg(0), m_pump_st(0), m_pump_ofd(-1) {}
                void cfg(int mxid);
                int start(int flg = 0, int mxid = -1, int *ppcp = 0), close(int x=0); // 32:aux 64:prep
		inline void c_play() { (*m_cur_cpf)(this); }
		inline int flg(int m=-1) const { return m_flg&m; }
		long long total_played() const { return m_total_played; }
		void set_vol(int x) { m_cfg.vol = x; }
		int vol() const { return m_cfg.vol; }
		int hcp_start(int t), hcp_end(int f = 0);
		int hcp() const { return m_hcp; }
		int cmd(const char *s), w(int flg);
		int pump_op(int x), pump_launch(int nt = 0, int flg = 0);
		int apmp_op(int x), aux_pump();
		int is_up() const { return (m_flg&2) ? m_pump_st&2 : !!m_hnd; }
		void debug(const char *s = 0);
        protected:
		static void cpf_mute(ASnd*), cpf_true(ASnd*), cpf_liar(ASnd*), cpf_pump(ASnd*), cpf_bug(ASnd*);
		int start1(int sc_lim), try_start(int n);
		void cfg_pre(int spd, int rsrv, int adjlim, int liar);
		int start_buf(int tlim);
		int adj2(int dif, unsigned int *to);
		void err3(int flg);
                int err(int k, const char *s, int ec = 0);
		void play2(short *buf, int nf);
		int e_msg_re(int e1, const char *s, int re);
		int play_and_adj(short *buf, int nf, int opt);
		int retry(int re);
		ASND_KF(NAME)  ASND_KF(CHCFG)  ASND_KF(SPD)     ASND_KF(RSRV)  ASND_KF(CLKLIM)
		ASND_KF(TRY_N) ASND_KF(TRY_MS) ASND_KF(SUSP_PA) ASND_KF(CMODE) ASND_KF(BFILL) ASND_KF(ADJLIM)

                snd_pcm_t * m_hnd;
		void (*m_cur_cpf)(ASnd*);
		long long m_total_played;
		au16w_t m_cfg;
		char m_hcp_lbl[8];
		int m_bufsiz, m_n_chan, m_mxid, m_hcp, m_hcp_s0, m_adj_lim;
		int m_bs, m_bs2, m_hwbs_trg, m_swb_siz, m_cfg_offs, m_flg; // flg: 1:liar 2:pump 8:fill 32:aux
		int m_pump_st, m_pump_ofd, m_pump_et, *m_pump_pcp; // st: 1-r2r 2-runn 4-wt4up 8:re 256*#try
		int m_re_usec, m_sclim;
		struct timespec m_re_stamp;
		int m_aux_jid;
		short * m_swb;
};

#define snd0 snd01[0]
#define snd1 snd01[1]
extern ASnd snd01[2];	struct jobq_ent_t;
int asnd_mkjob(jobq_ent_t * ent, char * arg), asnd_gcmd(const char *s);
#endif // __qwe_asnd_h__

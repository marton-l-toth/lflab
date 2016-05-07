#ifndef __qwe_job_h__
#define __qwe_job_h__

class Job; class CmdBuf; class ANode;

class JobQ {
        public: 
		typedef struct { Job * p; int nid; short st, xst, jid, i4f; int plttwwii; } ent_t;
		void init();
		int nj() const { return m_nj; }
		Job *  job(ANode * nd, int ix);
		int    jst(ANode * nd, int ix);
		int launch(ANode * nd, int ix, char * arg);
		bool run();
		void upd_gui(bool force = false);
		void upd_gui_1(ANode * nd, int ix) { if ((ix=lu(nd,ix))>=0) upd_gui_1p(m_ent+ix); }
		int kill(ANode * nd, int ix, int ec = JQE_KILL) { return kill5(lu(nd,ix), ec); }
		int purge();
		int cmd(ANode * nd, int ix, char * arg);
		void set_upd_t(int nsamp) { m_upd_t = nsamp; }
		void debug();
        protected:
		int lu(ANode * nd, int ix);
		int kill5(int j, int ec = JQE_KILL);
		void upd_gui_1p(ent_t * p);
		bool need_upd();
		unsigned char m_px[32];
		ent_t m_ent[32];
		int m_nj;
		unsigned int m_ent_bv;
		long long m_last_upd;
		int m_upd_t;
};

class Job {
        public: 
		Job(JobQ::ent_t * ent) : m_ent(ent) {}
		virtual ~Job() {}
		virtual int run1() = 0;
		virtual void abort() {}
		virtual int cmd(char * arg) { return JQE_NOCMD; }
	protected:
		JobQ::ent_t * m_ent;
	private:
		Job() {}
};

extern JobQ jobq;

#endif // __qwe_job_h__

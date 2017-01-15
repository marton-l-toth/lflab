#ifndef __qwe_job_h__
#define __qwe_job_h__

class Job; class CmdBuf; class ANode;

struct jobq_ent_t { Job * p; int nid; short st, xst, jid, i4f; int plttwwii; };
class JobQ {
        public: 
		typedef jobq_ent_t ent_t;
		typedef int (*jlfun_t)(ent_t *, char*);
		static int jl_dummy(ent_t * e, char* a) { return JQE_UNDEF; }		
		void init();
		int nj() const { return m_nj; }
		Job *  job(ANode * nd, int ix);
		int    jst(ANode * nd, int ix);
		int launch(ANode * nd, int ix, char * arg);
		bool run();
		void upd_gui(int force = 0);
		void upd_gui_1(ANode * nd, int ix) { if ((ix=lu(nd,ix))>=0) upd_gui_1p(m_ent+ix); }
		int wake(ANode * nd, int ix, int k);
		int kill(ANode * nd, int ix, int ec = JQE_KILL);
		int unblock(ANode * nd, int ix);
		int purge();
		int cmd(ANode * nd, int ix, char * arg);
		void set_upd_t(int nsamp) { m_upd_t = nsamp; }
		void debug();
        protected:
		int lu(ANode * nd, int ix);
		int kill_j(int j, int ec);
		int kill5(int j);
		void upd_gui_1p(ent_t * p);
		bool need_upd();
		unsigned char m_px[32];
		ent_t m_ent[32];
		int m_nj;
		unsigned int m_ent_bv;
		long long m_last_upd;
		int m_upd_t, m_upd_flg;
		jlfun_t m_jlfun[16];
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

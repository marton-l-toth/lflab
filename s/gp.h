#ifndef __qwe_gp_wrap_h__
#define __qwe_gp_wrap_h__

class PlotPar { public: virtual void dummy() {} };

class PlotPar_x : public PlotPar { public: double x; };
class PlotPar_xy : public PlotPar { public: double x,y ; };
class PlotPar_ix : public PlotPar { public: int i; double x; };
class PlotPar_ixy : public PlotPar { public: int i; double x,y ; };
class PlotPar_arr : public PlotPar {
	public:
		double *p;
		double x0, mul;
		PlotPar_arr(double *pp, int n, double x_0, double x_1) :
			p(pp), x0(x_0), mul((double)(n-1)/(x_1-x_0)) {}
		double y(double x) {
			return p[ (int)round( mul * (x - x0) ) ]; }
};

typedef double (*plotfun1_t) (double, PlotPar*);
typedef double (*plotfun2_t) (double, double, PlotPar*);

double arrfun1(double x, PlotPar* par);

class Gnuplot {
        public:
                Gnuplot() : m_pid(0), m_wid(1), m_y2bits(0) {}
                ~Gnuplot() { if (m_pid) stop(); }

                int start();
                void stop();
		int restart() { stop(); return start(); }
                int cmd(const char * fmt, ...);
                int cmd_n(const char * fmt, ...);
                int read_file(const char * path);
                FILE * get_inpipe() const { return m_inpipe; }
		void setfun1(int ix, plotfun1_t fun, PlotPar * arg = 0, bool y2 = false,
			     const char * nm = 0, const char * c = 0) {
			m_fun1[ix] = fun, m_funarg[ix] = arg;
			m_funname[ix] = nm, m_funcolor[ix] = c;
			ix = 1<<ix; y2 ? (m_y2bits |= ix) : (m_y2bits &= ~ix);
		}
		void setfun2(int ix, plotfun2_t fun, PlotPar * arg = 0, 
			     const char * nm = 0, const char * c = 0) {
			m_fun2[ix] = fun, m_funarg[ix] = arg;
			m_funname[ix] = nm, m_funcolor[ix] = c; }
		void plot1(int sel, double x0, double x1, int n, bool nwin = true);
		void plot2(int sel, double x0, double x1, double y0, double y1, int n, int m, bool nwin = true);
		static Gnuplot * sg();
		static int gp_dead(int pid, int stat, int td);
        protected:
		static Gnuplot * m0_sg;
                static const char * m0_path;
                static const char m0_gpini[];
                static const int m0_gpini_len;

                int m_pid;
                FILE * m_inpipe;
		int m_wid;
		plotfun1_t m_fun1[8];
		plotfun2_t m_fun2[8];
		PlotPar * m_funarg[8];
		const char * m_funname[8];
		const char * m_funcolor[8];
		int m_y2bits;
};

#endif //  __qwe_gp_wrap_h__


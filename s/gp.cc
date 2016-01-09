#include <signal.h>
#include <fcntl.h>
#include <sys/types.h>

#include "util.h"
#include "gp.h"
#include "glob.h"
#include "pt.h"

double arrfun1(double x, PlotPar* par) {
	PlotPar_arr * par1 = dynamic_cast<PlotPar_arr*> (par);
	return par1 -> y(x);
}

Gnuplot * Gnuplot::m0_sg = 0;
const char * Gnuplot::m0_path = "gnuplot";
const char Gnuplot::m0_gpini[] = "set style data line\n"
                                 "set autoscale fix\n";

const int Gnuplot::m0_gpini_len = sizeof(m0_gpini)-1;

Gnuplot * Gnuplot::sg() {
	return m0_sg ? m0_sg : (m0_sg = new Gnuplot(), m0_sg->start(), m0_sg); }

int Gnuplot::gp_dead(int pid, int stat, int td) {
	if (!m0_sg) log("BUG: unexp. gp_dead() (!sg) ");
	else if (m0_sg->m_pid!=pid) log("BUG: unexp. gp_dead() (%d!=%d)", m0_sg->m_pid, pid);
	else fclose(m0_sg->m_inpipe), m0_sg->m_pid = 0, delete(m0_sg), m0_sg = 0; 
	return 0; 
}

int Gnuplot::start() {
        if (m_pid) return fprintf(stderr,"gnuplot already started: pid %d", m_pid), -2;
	int rv, inpipe = -1;
	m_pid = launch(m0_path, "!>p1", &inpipe, (char*)0);
        if ((rv = write(inpipe, m0_gpini, m0_gpini_len))!=m0_gpini_len) perror("gnuplot/inpipe");
        m_inpipe = fdopen(inpipe, "w");
	pt_reg(PT_GPLT, m_pid, &gp_dead);
	return 0;
}

void Gnuplot::stop() {
	if (!m_pid) return;
        fclose(m_inpipe);
        kill(m_pid, 9);
        m_pid = 0;
}

int Gnuplot::cmd(const char * fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        int rv = vfprintf(m_inpipe, fmt, ap);
        if (rv<0) perror("gnuplot/inpipe");
        va_end(ap);
        fprintf(m_inpipe,"\n");
        fflush(m_inpipe);
        return rv;
}

int Gnuplot::cmd_n(const char * fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        int rv = vfprintf(m_inpipe, fmt, ap);
        if (rv<0) perror("gnuplot/inpipe");
        va_end(ap);
        return rv;
}


void Gnuplot::plot1(int sel, double x0, double x1, int n, bool nwin)
{
	if (nwin) cmd("set term x11 %d", m_wid++);
	cmd("set y2tics nomirror format \"%s\"",
			(sel&m_y2bits) ? "%g" : "");

	cmd_n("plot ");
	int ch = ' ';
	for (int i=0; i<8; i++) {
		if ( !(sel&(1<<i)) ) continue;
		cmd_n("%c'-'", ch);
		if (m_y2bits & (1<<i)) cmd_n("axes x1y2");
		if (m_funname[i]) cmd_n(" title \"%s\"", m_funname[i]);
		if (m_funcolor[i]) cmd_n(" lc rgbcolor \"%s\"", m_funcolor[i]);
		ch = ',';
	}
	cmd("");
	double xstep = (x1-x0) / (double)(n-1);
	for (int i=0; i<8; i++) {
		if ( !(sel&(1<<i)) ) continue;
		double x = x0;
		for (int j=0; j<n; j++) {
			double v = (*m_fun1[i]) (x , m_funarg[i]);
			if (v==v)
				cmd("%g %g", x, v );
			x += xstep;
		}
		cmd("e");
	}
}

void Gnuplot::plot2(int sel, double x0, double x1, double y0, double y1, int n, int m, bool nwin)
{
	if (debug_flags & DFLG_PLOT2) log("plot2 sel=%d",sel);
	if (nwin) cmd("set term x11 %d", m_wid++);
	cmd_n("splot ");
	int ch = ' ';
	for (int i=0; i<8; i++) {
		if ( !(sel&(1<<i)) ) continue;
		if (debug_flags & DFLG_PLOT2) 
			log("i=%d title='%s' color='%s'",i,str0(m_funname[i]),str0(m_funcolor[i]));
		cmd_n("%c'-'",ch); ch = ',';
		if (m_funname[i]) cmd_n(" title \"%s\"", m_funname[i]);
		if (m_funcolor[i]) cmd_n(" lc rgbcolor \"%s\"", m_funcolor[i]);
	}
	cmd("");
	double xstep = (x1-x0) / (double)(n-1);
	double ystep = (y1-y0) / (double)(m-1);
	for (int i=0; i<8; i++) {
		if ( !(sel&(1<<i)) ) continue;
		double x = x0;
		for (int j=0; j<n; j++) {
			double y = y0;
			for (int k=0; k<m; k++) {
				double v = (*m_fun2[i]) (x, y, m_funarg[i]);
				if (v==v)
					cmd("%g %g %g", x, y, v );
				y += ystep;
			}
			cmd("");
			x += xstep;
		}
		cmd("e");
	}
}


#include <fcntl.h>
#include <sys/stat.h>
#include <signal.h>

#define QWE_DEFINE_ERRTAB
#include "glob.h"
#include "cfgtab.inc"
#include "cmd.h"
#include "util.h"
#include "util2.h"
#include "guistub.h"
#include "asnd.h"
#include "midi.h"
#include "pt.h"

int glob_flg = GLF_EMPTY | GLF_INI0, debug_flags = 0, sample_rate = 44100, killer_fd = -1;
double sample_length = 1.0/44100.0, natural_bpm = 65.625, natural_bpm_r = 1.0/65.625, junkbuf[4096];
char zeroblkC[32768], save_file_name[1024], debug_utxt_buf[1024];

GuiStub gui2;
JobQ jobq;
BufClock clk0;
ASnd snd0;
QuickStat qstat;

#define N_SLCMD 4
#define PFD(J) sl_cmd[J].pfd()
static CmdBuf sl_cmd[N_SLCMD];
static struct timespec asv_ts;

void hi() { log("lflab: linear filter based audio lab %d.%02d (%s)\n%s\n%s\n%s", v_major, v_minor, pt_hello,
	    "Copyright (C) 2014-2016 Marton Laszlo Toth","This is free software with ABSOLUTELY NO WARRANTY.",
	    "see the file \"COPYING\" for details"); }

static void rf(const char ** ppf) {
	for (int ec, i=0, ro=1; ppf[i]; i++) {
		if (!ppf[i][0]) continue;
		if (ppf[i][0]=='/' && !ppf[i][1]) { ro = 0; continue; }
		log("reading %s... (r%c)", ppf[i], 119-8*ro);
		if ((ec = CmdBuf::read_batch(ppf[i], NOF_BATCHF|(ro&&ppf[i+1]))) < 0) gui_errq_add(ec);
		if (ec<0) log("%s: %s(%d)", ppf[i], err_str(ec), ec);
	}}

#define INI_LIST hi(), nd0_init(), nz_init(), calc_init(), graph_init(), nd_init(), mx_init(), wrap_init(), \
		 track_init(), itb_init(), midi_init()
static void ini(const char ** ppf) {
	extern void INI_LIST; INI_LIST;
	CmdBuf::st_init();  jobq.init();  glob_flg ^= (GLF_INI0|GLF_INI1); snd0.set_vol(92); // :)
	if (CFG_INI_ORDER.i) gui2.start(PFD(1)), rf(ppf); else rf(ppf), gui2.start(PFD(1));
	sl_cmd[0].init( -1,  NOF_ERRMSG|NOF_FGUI);
	sl_cmd[1].init( -1,  NOF_GUI, 126);
	sl_cmd[2].init( -1 , NOF_ERRMSG);
	sl_cmd[3].init(-'A', NOF_ERRMSG);
	ADirNode *btin = static_cast<ADirNode*>(ANode::lookup_n_q(1)),
		 *hlp  = static_cast<ADirNode*>(ANode::lookup_n_q(2));
	gui2.root_expand(); gui2.tree_expand(1, btin); gui2.tree_expand(1, hlp);
	log("### i2m=%d, gcp=%d", *PFD(0), *PFD(1));
	if (glob_flg&GLF_HITHERE) {const char *s="getting started"; snd0.w(-1); hlp->sn(&s)->draw_window(11);}
	pt_wrk_start(0); if (CFG_DEVEL.i) pt_con_op("-1");
	snd0.start(0, 0); clk0.tcond(&asv_ts); glob_flg &= ~GLF_INI1; glob_flg |= GLF_SAVED;
}

#define FOR_SLC for (int k,i=0; i<N_SLCMD; i++) if ((k=sl_cmd[i].fd()) >= 0)
static void sel_loop() {
	fd_set rset; struct timeval tv; int nj = 1;
	while (1) {
		int nfd = 0; FD_ZERO(&rset); FOR_SLC { FD_SET(k, &rset); if (k>=nfd) nfd = k+1; }
		BVFOR_JM(midi_bv) { int k=midi_fd[j];  FD_SET(k, &rset); if (k>=nfd) nfd = k+1; }
		tv.tv_sec=0; tv.tv_usec = clk0.sel(nj);
		int r = select(nfd, &rset, 0, 0, &tv), cf = pt_chld_flg; 
		if (cf) pt_chld_act(), cf &= 3;
		if (r>0){ FOR_SLC FD_ISSET(k,&rset) && !((1<<i)&cf) && (clk0.ev(48+i), sl_cmd[i].read_f());
			  BVFOR_JM(midi_bv) if (FD_ISSET(midi_fd[j], &rset)) clk0.ev('M'), midi_input(j);  }
		else if (r<0) { perror("select"); }
		if (gui2.pending()) clk0.ev('G'), gui2.flush_all();
		if ((nj=jobq.nj())&&clk0.j0()){ while (clk0.j1(jobq.run())); jobq.upd_gui(0),nj=jobq.purge(); }
		snd0.c_play();
		if ((r=glob_flg&GLF_SILENCE) && !((glob_flg^=r)&GLF_FSTATE) && (r=CFG_ASV_MIN.i) &&
			   clk0.tcond(&asv_ts, r*60000)) Node::save_batch(Node::root(), "/", 0);
	}}

int main(int ac, char** av) { ini(pt_init(ac,av,PFD(0),PFD(2))); sel_loop(); }

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

#define INI_LIST pt_init(), hi(), nd0_init(), cfg_init(), nz_init(), calc_init(), graph_init(), nd_init(), \
		 mx_init(), wrap_init(), track_init(), midi_init()
int glob_flg = 0;
int debug_flags = 0;
int sample_rate = 44100;
int killer_fd = -1;
double sample_length = 1.0/44100.0, natural_bpm = 65.625, natural_bpm_r = 1.0/65.625;
char mostly_zero[0x8080];
double junkbuf[4096];
char save_file_name[1024];
const char *tmp_dir,    *usr_dir,    *wrk_dir,    *hsh_dir, *autosave_name;
int         tmp_dir_len, usr_dir_len, wrk_dir_len, hsh_dir_len;

GuiStub gui2;
JobQ jobq;
ASnd snd0;

#define N_SLCMD 4
static CmdBuf sl_cmd[N_SLCMD];
static int asv_ts[2];

void hi() { log("linear filter lab %d.%02d\n%s\n%s\n%s", v_major, v_minor,
	    "Copyright (C) 2014-2015 Marton Laszlo Toth","This is free software with ABSOLUTELY NO WARRANTY.",
	    "see the file \"COPYING\" for details"); }

static void ini0() {
	const char * dfs = getenv("LF_DEBUGFLG"); if (dfs) debug_flags = atoi_h(dfs);
	imp4097[0] = 1.0; vstring_set(v_major, v_minor);
	extern void INI_LIST; INI_LIST;
	CmdBuf::st_init();
	jobq.init();
	if (CFG_INI_ORDER.i) gui2.start();
}

static void rf(char ** nms, int n) {
	Clock clk; memcpy(save_file_name, "unnamed.lf", 11);
	for (int i=0, rwix=n-1; i<n; i++) {
		if (!nms[i][0] || (nms[i][0]=='/' && !nms[i][1] && (rwix=i+1,1))) continue;
		if (i==rwix) strncpy(save_file_name, nms[i], 1023);
		clk.reset(); log("reading %s... (r%c)", nms[i], 111+8*(i>=rwix));
		int ec = CmdBuf::read_batch(nms[i], NOF_BATCHF|(i<rwix)); if (ec<0) gui_errq_add(ec);
		log("%s: %s (%d)", nms[i], ec<0 ? err_str(ec) : "ok", ec<0 ? ec : clk.get());
	}}

static void ini1() {
	if (!CFG_INI_ORDER.i) gui2.start();
	sl_cmd[0].init(pt_cp_i2m, NOF_ERRMSG);
	sl_cmd[1].init(gui2.outpipe(), NOF_GUI, 126);
	sl_cmd[2].init(-'C', NOF_ERRMSG);
	sl_cmd[3].init(-'A', NOF_ERRMSG);
	gui2.root_expand(); gui2.tree_expand(1,static_cast<ADirNode*>(ANode::lookup_n_q(1)));
	                    gui2.tree_expand(1,static_cast<ADirNode*>(ANode::lookup_n_q(2)));
	sl_cmd[2].sn("v50\n", 4);
	log("###tpipe=%d", gui2.tpipe());
	if (getenv("LF_HITHERE")) snd0.w(-1);
	snd0.cfg(gui2.tpipe(), 0); snd0.start();
	snd0.cond_clk(asv_ts, 1);
}

#define FOR_SLC for (int k,i=0; i<N_SLCMD; i++) if ((k=sl_cmd[i].fd()) >= 0)
static void sel_loop() {
	fd_set rset; struct timeval tv; int nj = 1;
	while (1) {
		int nfd = 0; FD_ZERO(&rset); FOR_SLC { FD_SET(k, &rset); if (k>=nfd) nfd = k+1; }
		BVFOR_JM(midi_bv) { int k=midi_fd[j];  FD_SET(k, &rset); if (k>=nfd) nfd = k+1; }
		tv.tv_sec=0; tv.tv_usec = snd0.time4sel(nj);
		int r = select(nfd, &rset, 0, 0, &tv), cf = pt_chld_flg; 
		if (cf) pt_chld_act(), cf &= 3;
		if (r>0){ FOR_SLC FD_ISSET(k,&rset) && !((1<<i)&cf) && (snd0.mark(48+i), sl_cmd[i].read_f());
			  BVFOR_JM(midi_bv) if (FD_ISSET(midi_fd[j], &rset)) snd0.mark('M'), midi_input(j);  }
		else if (r<0) { perror("select"); }
		if (gui2.pending()) snd0.mark('G'), gui2.flush_all();
		if ((nj=jobq.nj())){ while (snd0.time4job()&&jobq.run());  jobq.upd_gui(0), nj=jobq.purge(); }
		snd0.c_play();
		if ((glob_flg&GLF_SILENCE) && (glob_flg&=~GLF_SILENCE, r=CFG_ASV_MIN.i) &&
				snd0.cond_clk(asv_ts, r*60000000)) Node::save_batch(Node::root(), "/", 0);
	}}

int main(int ac, char** av) { ini0(); rf(av+1, ac-1); ini1(); sel_loop(); }

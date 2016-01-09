#include <signal.h>
#include <glob.h>

#include "util.h"
#include "glob.h"
#include "guistub.h"
#include "mx.h"
#include "pt.h"
#include "asnd.h"

void gui_closewin(int x) { gui2.closewin(x); }
void gui_errq_add(int x) { gui2.errq_add(x); }

void GuiStub::errq_add(int ec) {
	if (!ec) return; ec = (ec==EEE_ERRNO) ? errno&0xffffff : ec&0xffffff;
	for (int v, i=0; i<m_errq_n; i++) if (((v=m_errq_v[i])&0xffffff)==ec && v<0x63000000)
		return (void) (m_errq_v[i] += 0x1000000);
	if (m_errq_n==32) log("errq_add: bug???"), memmove(m_errq_v, m_errq_v+1, 124), --m_errq_n;
	m_errq_v[m_errq_n++] = 0x1000000 + ec;
}

void GuiStub::errq_cfl() {
	if (!m_errq_n || !(snd0.cond_clk(m_errq_t0, 250000))) return;
	memcpy(m_bufp, "\tC47$E>?", 8), m_bufp+=7;
	for (int i=0; i<m_errq_n; i++) hex8(m_errq_v[i]);   m_errq_n = 0; }

int GuiStub::gui_dead(int pid, int stat, int td) {
	int err = !WIFEXITED(stat) || WEXITSTATUS(stat), tlim = 1+err;
	if (td<tlim) log("gui exited again in < %d sec., bye", td), bye(1);
	if (gui2.m_pid!=pid) return log("BUG: unexp. gui_dead() (%d!=%d)", gui2.m_pid, pid), 0;
	ANode::wi_clear();
	gui2.m_pid = 0; gui2.start(); gui2.root_expand();
	if (err) gui2.errq_add(PTE_GUICRASH);
	return gui2.m_pid;
}

int GuiStub::start() {
        if (m_pid) return log("gui2 already started: pid %d", m_pid), PTE_WTF;
	static const char * path = 0;
	if (!path && !(path=getenv("LF_GUI"))) path = "./lf.gui";
	int pfi, pfo, pft;
	if ((m_pid = launch(path, "!><u>", &pfi, &pfo, &pft, (char*)0))<0) log("FATAL: %s\n", path), bye(1);
	pt_reg(PT_GUI, m_pid, &gui_dead);
	set_fd(&m_inpipe, pfi); set_fd(&m_outpipe, pfo); set_fd(&m_tpipe, pft);
	log("m_tpipe = %d", m_tpipe);
	clear(); pf("\tW7$.Vtv%d.%02d\tv%d.%02d", v_major, v_minor, v_major, v_minor); 
	savename(); flush();
	return 0;
}

void GuiStub::stop() {
        if (!m_pid) return;
        close(m_inpipe); close(m_outpipe);
        kill(m_pid, 9);
        m_pid = 0;
}

int GuiStub::flush() {
	int r, i = 0, nx = 0, l = m_bufp - m_buf0 + 1; if (l<2) return 0; else *m_bufp = 10;
	for (int c,j=0; j<l; j++) if ((unsigned int)((c=m_buf0[j])-32)>94u && ((c-9)&254)) ++nx,m_buf0[j]=63;
	if (nx) log("BUG: gui2/flush: %d invalid characters, replaced with '?'");
	if ((debug_flags & DFLG_GUI2)) log_sn("to-gui2:", m_buf0, l);
	do {
		r = write(m_inpipe, m_buf0 + i, l);
		if (r<=0) { log("ERROR: gui2/flush: %s", r ? strerror(errno) : "0 bytes written!"); break; } 
		else if (r>l) { log("ERROR: gui2/flush: %d/%d written!", r, l); break; } 
	} while ((i += r, l -= r));
	clear();
	if (l) {
		log("ERROR: GUI does not seem to work, type 's<filename>' to save, 'q' to quit (with autosave)");
		pt_con_op(-1); return -1;
	} 
	return 1;
}

void GuiStub::hexs(unsigned const char * s, int k) {
	for (int i=0; i<k; i++) m_bufp[0] = hexc1(*s>>4), m_bufp[1]=hexc1(*s&15), m_bufp+=2, s++; }

void GuiStub::pf(const char * fmt, ...) {
	va_list ap; va_start(ap, fmt);
	m_bufp += vsprintf(m_bufp, fmt, ap); va_end(ap);
}

void GuiStub::t_pf(const char * fmt, ...) {
	va_list ap; va_start(ap, fmt); t0();
	m_bufp += vsprintf(m_bufp, fmt, ap); va_end(ap);
}

void GuiStub::cre(int oid, int ty, const char *s) { 
	setwin(oid, ty); lx0('C'); if (s) sz(s);
	int i4 = oid & 15; if (i4==7) return; // TODO: globwin
	ANode * nd = ANode::lookup_n(oid >> 4); 
	if (!nd) return log("gui2/cre: nd_id 0x%x not found", oid >> 4);
	nd -> winflg_or(1 << i4);
}

void GuiStub::t2_sel(int lr, ANode *nd) {
	m_bufp[0] = 9; m_bufp[1] = '+'; m_bufp[2] = 48+(lr&1); m_bufp[3] = nd->cl_id();
	m_bufp += 4; hex5(nd->id()); c1(36); sn(Node::rgb(nd), 6);
	int k = nd -> get_path_uf(m_bufp, 256);
	if (k) m_bufp += k; else *(m_bufp++) = '.';
}

int GuiStub::tree_expand(int lr, ADirNode * dir) {
	if (!dir) return 0;
	ADirNode * up = dynamic_cast<ADirNode*> (dir->up());
	lr += 49; c2(9, 'T');
	hex5(up?up->id():0); c2(lr, 36);
	hex5(dir->id()); c2(lr|48, 36);
	m_bufp += dir->gui_list(m_bufp, 0);
	dir->winflg_or(2+2*lr); invd(); return 0;
}

void GuiStub::tree_sel(int lr, ANode * nd) {
	tree_force(lr, dynamic_cast<ADirNode*>(nd->up())); invd();
	pf("\tN%x%c*%x", nd->up()?nd->up()->id():0, lr+49, nd->id());
}

void GuiStub::clip_box_1(ClipNode * cl, int i) {
	BoxGen * bx = cl -> bx_j(i);
	if (bx) c1(i_to_b32(i)), bxmini(bx); else c4(i_to_b32(i), 48, 48, 48); }

void GuiStub::clip_box(ClipNode * cl, int i, int j, int sel) {
	setwin(16*cl->id()+3,'K'); w0(); c1('K');
	if (sel>=0) c2('+', i_to_b32(sel));
	if (i>=0) clip_box_1(cl, i);
	if (j>=0) clip_box_1(cl, j);
}

void GuiStub::clip_flg(int id, int fid, int v01) {
	setwin(16*id+3,'K'); wupd_i(fid, v01); }

void GuiStub::node_name(int i, ANode * nd) {
	ADirNode * dir = dynamic_cast<ADirNode*>(nd->up()); 
	if (!dir) return log("node_name: no updir");
	c2(9, 'N'); hex5(dir->id()); c2(i+49, '+'+nd->is_dir());
	hex5(nd->id()); c1(36); m_bufp += nd->get_name(m_bufp);
}

void GuiStub::node_rm(int i, ANode * nd) {
	ADirNode * dir = dynamic_cast<ADirNode*>(nd->up()); 
	if (!dir) return log("node_name: no updir"); else invd();
	c2(9, 'N'); hex5(dir->id()); c2(i+49, '-'); hex5(nd->id());
}

void GuiStub::acv_open(int j) { cre(ACV_WIN(j),'A'); hex8(j); }
void GuiStub::flush_all() { cfl(); slr_upd(); errq_cfl(); flush(); }

void GuiStub::j_upd(int wwt, int st, int wwix) {
	w0(); c1(wwt); if (wwix>=0) hexn(wwix, 2);
	c1('s'); c2(48+(st>>6), 48+(st&63)); }

void GuiStub::own_title(int flg) {
	if (m_cwi<0) return log("gui2/own_title: cwi<0");
	ANode * nd = ANode::lookup_n_q(m_cwi>>4);
	if (!nd->cl_id()) return log("gui2/own_title: cwi lookup failed");
	if (flg&1) { w0(); sn("_00t!", 5); m_bufp += nd->title_arg(m_bufp); }
}

void GuiStub::ref_title(int wwt, ANode * nd, int wwix, const char * defstr) {
	wupd_c0(wwt, 't', wwix); if (nd) { m_bufp += nd->title_arg(m_bufp); return; }
	sn("ttt666", 6); if (!defstr) { sn("(none)", 6); return; }
	sn("(no ", 4); sz(defstr); sn(" box)", 5);
}

void GuiStub::savename() { const char *p, *s = save_file_name;
	for (p=s; *p; p++) if (*p=='/') s=p+1;
	setwin( 7,'.'); wupd_s('N',","); sz(s); setwin(23,'/'); wupd_s('1',save_file_name); }


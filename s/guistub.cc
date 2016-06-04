#include <signal.h>
#include <glob.h>

#include "util.h"
#include "glob.h"
#include "guistub.h"
#include "mx.h"
#include "pt.h"
#include "asnd.h"
#include "cfgtab.inc"

void gui_closewin(int x) { gui2.closewin(x); }
void gui_errq_add(int x, const char *s) { gui2.errq_add(x, s); }

int gui_acv_op(int j, int op) { if (op<0) op = (0x73fe>>(4*CFG_AO_ACTION.i)) & 15;
	return op==0xe ? (gui2.cre(ACV_WIN(j),'A'), gui2.hex8(j), 0) : pt_acv_op(j, op|256, 0, 0); }

void gui_sliderwin(int oid, int n, const double * lbl, const unsigned char * v0) {
	gui2.cre(oid, 'J'); gui2.c1(hexc1(n));
	for (int i=0; i<=n; i++) gui2.hdbl(lbl[i]);
	for (int i=0; i< n; i++) gui2.c2(hexc1(v0[i]>>4), hexc1(v0[i]&15));
}

void GuiStub::errq_add(int ec, const char *s) {
	if (ec && s) log("%s: %s", s, err_str(ec));
	if (!ec || ec==EEE_NOEFF || ec==RTE_IWCLOSE) return; 
	ec = (ec==EEE_ERRNO) ? errno&0xffffff : ec&0xffffff;
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
	static const char * path = 0; if (!path && !(path=getenv("LF_GUI"))) path = "./lf.gui";
	int pfi, pfo, pft;
	if ((m_pid = launch(path, "!><u>", &pfi, &pfo, &pft, (char*)0))<0) log("FATAL: %s\n", path), bye(1);
	pt_reg(PT_GUI, m_pid, &gui_dead);
	m_gf0 = glob_flg;
	set_fd(&m_inpipe, pfi); set_fd(&m_outpipe, pfo); set_fd(&m_tpipe, pft);
	clear(); pf("\tW7$.Vtv%d.%02d\tv%d.%02d", v_major, v_minor, v_major, v_minor); 
	set_tlog(); savename(); vol(); flush();
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
		pt_con_op("-1"); return -1;
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
	m_bufp += 4; hex5(nd->id()); c1(36); sn(nd->rgb(), 6);
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

void GuiStub::flush_all() { 
	cfl(); slr_upd(); errq_cfl();
	int x = glob_flg^m_gf0; m_gf0 = glob_flg;
	if (x & GLF_FSTATE) savename();
	if (m_gnaq_n && (snd0.total_played()-m_gnaq_t)>(long long)(44*CFG_GUI_TRKUP.i)) {
		for (int i=0, n=m_gnaq_n; i<n; i++) trk_w_gna(m_gnaq_id[i]);
		m_gnaq_t = snd0.total_played(); m_gnaq_n = 0; }
	flush();
}

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

void GuiStub::set_tlog() { c4(9, 't', 'c', 48+16*!!CFG_TLOG_AUTO.i + CFG_TLOG_BACKUP.i); }
void GuiStub::vol() { setwin(7,'.'); wupd_i2('v', snd0.vol() - 12); }

void GuiStub::mcfg_ud(int wch, cfg_ent * pc, const char *sdef, int ldef) {
	wupd_c0(wch, 't');
	if (pc->i) sn("%%%ttt", 6), sn(pc->s, pc->i); else sn("ppp666(", 7), sn(sdef, ldef), c1(')'); }

void GuiStub::mcfg_win(int flg) {
	if (flg&     1) cre(0x57, 'F'); else setwin(0x57, 'F');
	if (flg&     2) wupd_s('7', ","), sz(*save_file_name?save_file_name:"(no name)");
	if (flg&     4) wupd_i1('s', CFG_SV_EXEC.i);
	if (flg&     8) wupd_i2('a', CFG_ASV_MIN.i);
	if (flg&    16) wupd_i1('t', CFG_TLOG_AUTO.i);
	if (flg&    32) wupd_i2('S', CFG_SV_BACKUP.i);
	if (flg&    64) wupd_i2('A', CFG_ASV_BACKUP.i);
	if (flg&   128) wupd_i2('T', CFG_TLOG_BACKUP.i);
	if (flg&   256) mcfg_ud('k', &CFG_AO_DIR,  tmp_dir, tmp_dir_len);
	if (flg&   512) mcfg_ud('w', &CFG_WAV_DIR, hsh_dir, hsh_dir_len);
	if (flg&  1024) wupd_ls('K', CFG_AO_ACTION.i);
	if (flg&  2048) wupd_i ('L', CFG_AO_TLIM.i);
	if (flg&  8192) wupd_i1('d', CFG_DEVEL.i);
	if (flg& 16384) wupd_s ('x', CFG_XTERM.s);
	if (flg& 32768) wupd_i1('C', CFG_AUTOCON.i);
}

void GuiStub::savename() { 
	char buf[64]; const char *p, *sb, *rgb = (glob_flg&GLF_EMPTY)?"%%%FFF":0, *sp = save_file_name;
	int rf = glob_flg & GLF_RECOVER;
	if (*sp) { if (!rgb) rgb="%%%%^^"; for (p=sb=sp; *p; p++) if (*p=='/') sb=p+1; if (rf) sb="BUG!"; }
	else if (rf) { memcpy(buf,"(recover:0)",12); buf[9]=sp[2]; if (!rgb) rgb="zzz%%h"; sb=sp=buf; }
	else { sp = sb = rgb ? "(empty)" : (rgb = "zz%z%%", "(unnamed)"); }
	setwin( 7,'.'); wupd_c0('N','t'); sn(rgb, 6); sz(sb); 
	setwin(23,'/'); wupd_c0('1','t'); sn(rgb, 6); sz(sp);  mcfg_win(2);
}

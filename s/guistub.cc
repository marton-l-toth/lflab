#include <signal.h>
#include <glob.h>

#include "util.h"
#include "glob.h"
#include "guistub.h"
#include "mx.h"
#include "pt.h"
#include "cfgtab.inc"
#include "asnd.h"
#include "midi.h"
#include "util2.h"

void gui_closewin(int x) { gui2.closewin(x); }
void gui_errq_add (int x, 	 const char *s) { gui2.errq_add (x,    s); }
void gui_errq_add2(int x, int y, const char *s) { gui2.errq_add2(x, y, s); }
void gui_midi(int flg) { gui2.midi(flg); }
void gui_tlog(int i0, int n) { gui2.c2(9, 'c'); gui2.hex8(i0); gui2.hex8(n); }
int  gui_dead(int ef) { return gui2.gui_dead(ef); }

int gui_acv_op(int j, int op) { if (op<0) op = (0x73fe>>(4*CFG_AO_ACTION.i)) & 15;
	return op==0xe ? (gui2.cre(ACV_WIN(j),'A'), gui2.hex8(j), 0) : pt_acv_op(j, op|256, 0, 0, 0); }

void gui_sliderwin(int oid, int n, const double * lbl, const unsigned char * v0) {
	gui2.cre(oid, 'J'); gui2.c1(hexc1(n));
	for (int i=0; i<=n; i++) gui2.hdbl(lbl[i]);
	for (int i=0; i< n; i++) gui2.hex2(v0[i]);
}

void GuiStub::errq_add(int ec0, const char *s) {
	int ec1 = (ec0==EEE_ERRNO) ? errno : ec0, ec2 = ec1 & 0xffffff;
	clk0.ev2('E', ec1);
	if (ec0 && s) log("%s: %s", s, ec0!=EEE_ERRNO ? err_str(ec0) : strerror(ec1));
	if (!ec0 || ec0==EEE_NOEFF || ec0==RTE_IWCLOSE) return; 
	for (int v, i=0; i<m_errq_n; i++) if (((v=m_errq_v[i])&0xffffff)==ec2 && v<0x63000000)
		return (void) (m_errq_v[i] += 0x1000000);
	if (m_errq_n==32) log("errq_add: bug???"), memmove(m_errq_v, m_errq_v+1, 124), --m_errq_n;
	m_errq_v[m_errq_n++] = 0x1000000 + ec2;
}

void GuiStub::errq_cfl() {
	if (!m_errq_n || !(clk0.tcond(&m_errq_t0, 250))) return;
	memcpy(m_bufp, "\tC47$E>?", 8), m_bufp+=7;
	for (int i=0; i<m_errq_n; i++) hex8(m_errq_v[i]);   m_errq_n = 0; }

int GuiStub::gui_dead(int err) {
	errtemp_cond("gui/exit"); ANode::wi_clear();
	gui2.m_pid = 0; gui2.start(); gui2.root_expand();
	if (err) gui2.errq_add(PTE_GUICRASH);
	return gui2.m_pid;
}

int GuiStub::start(int * pfd) {
        if (m_pid) return log("gui2 already started: pid %d", m_pid), PTE_WTF;
	if (pfd) m_pfd = pfd; else if (!m_pfd) bug("gui start without pfd");
	int pfi;  p_close(m_pfd);
	cfg_export(&CFG_GUI_SLFOC); cfg_export(&CFG_GUI_QCPU); cfg_export(&CFG_GUI_F_YN);
	if ((m_pid = launch(QENV('g'), "!><u", &pfi, m_pfd, (char*)0)) < 0)
		log("FATAL: %s\n", QENV('g')), bye(1);
	m_gf0 = glob_flg;
	set_fd(&m_inpipe, pfi, 0); IFDBGX(PT) log("gui_start:  in:%d out:%d", m_inpipe, *m_pfd);
	clear(); pf("\tW7$.Vtv%d.%02d\tv%d.%02d", v_major, v_minor, v_major, v_minor);
	clk0.tcond(&m_errq_t0);
	savename(); vol(); xapp_bv(); flush();
	return 0;
}

void GuiStub::stop() {
        if (!m_pid) return;
        close(m_inpipe); p_close(m_pfd);
        kill(m_pid, 9);
        m_pid = 0;
}

int GuiStub::flush() {
	int r, i=0, nx=0, uec=0, l0 = m_bufp - m_buf0 + 1, l = l0; if (l<2) return 0; else *m_bufp = 10;
	for (int c,j=0; j<l; j++) if ((unsigned int)((c=m_buf0[j])-32)>94u && ((c-9)&254)) 
		++nx, uec = m_buf0[j], m_buf0[j] = 63;
	if ( (nx && (log("BUG: gui2/flush: %d invalid chars (last:0x%x), replaced with '?'", nx, uec),1)) ||
	     ((debug_flags & DFLG_GUI2) && (l!=19 || m_buf0[1]!='c')) ) log_sn(">G>", m_buf0, l, 1);
	do {
		r = write(m_inpipe, m_buf0 + i, l);
		if (r<=0) { log("ERROR: gui2/flush: %s", r ? strerror(errno) : "0 bytes written!"); break; } 
		else if (r>l) { log("ERROR: gui2/flush: %d/%d written!", r, l); break; } 
	} while ((i += r, l -= r));
	clear(); return l ? (errtemp_cond("write to GUI"), -1) : 0;
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

void GuiStub::t0() { cfl(); m_bufp += (m_cwi==m_lwi) ? (*(int*)m_bufp = m_th, 4) : lx0('U'); }
void GuiStub::w0() { cfl(); m_bufp += (m_cwi==m_lwi) ? (*(int*)m_bufp = m_wh, 4) : lx0('W'); }
void GuiStub::setwin(int oid, int ty) { m_cwi=oid; m_cwt=ty; ty<<=24; m_th=ty+0x245509; m_wh=ty+0x245709; }
int  GuiStub::lx0(int c) { c2(9, c); hexn(m_lwi=m_cwi, 6); c2(36, m_cwt); return 0; }
void GuiStub::closewin(int id) { c2(9,'Z'); hexn(id, 6); c2(36, 63); }
void GuiStub::wupd(int wwt, int wwix) { w0(); c1(wwt); if (wwix>=0) hex2(wwix); }
void GuiStub::wupd_s(int wwt, const char *s, int wwix) 		{ wupd(wwt,wwix); c1('t'); sz(s); }
void GuiStub::wupd_cs(int wwt, int c, const char *s, int wwix)  { wupd(wwt,wwix); c1(c); sz(s); }
void GuiStub::wupd_ct(int wwt, const char *s, int wwix)		{ wupd(wwt,wwix); c2('t',','); sz(s); }
void GuiStub::wupd_0(int wwt, const char *s, int wwix)		{ wupd(wwt,wwix); sz(s); }
void GuiStub::wupd_c0(int wwt, int c, int wwix) { wupd(wwt,wwix); c1(c); }
void GuiStub::wupd_i(int wwt, int x, int wwix)  { wupd(wwt,wwix); c1('x'); if(x&~0xfffff) hexn(x,8);
                                                                      	   else           hex5(x); }
void GuiStub::wupd_si(int wwt, int x, int wwix) { if (x<0) wupd_i(wwt, -x, wwix), c1('-');
                                                  else     wupd_i(wwt,  x, wwix); }
void GuiStub::wupd_i1(int wwt, int x, int wwix)   { wupd(wwt,wwix); c2('x', hexc1(x)); }
void GuiStub::wupd_c48(int wwt, int x, int wwix)  { wupd(wwt,wwix); c2('c', x+48); }
void GuiStub::wupd_ls(int wwt, int x, int wwix)   { wupd(wwt,wwix); c2('+', x+48); }
void GuiStub::wupd_i2(int wwt, int x, int wwix)   { wupd(wwt,wwix); c1('x'); hex2(x); }
void GuiStub::wupd_d(int wwt, double x, int wwix) { wupd(wwt,wwix); c1('@'); hdbl(x); }

void GuiStub::t2_sel(int lr, ANode *nd) {
	m_bufp[0] = 9; m_bufp[1] = '+'; m_bufp[2] = 48+(lr&1); m_bufp[3] = nd->cl_id();
	m_bufp += 4; hex5(nd->id()); c1(36); sn(nd->rgb(), 6);
	int k = nd -> get_path_uf(m_bufp, 256);
	if (k) m_bufp += k; else *(m_bufp++) = '.';
}

int GuiStub::tree_expand(int lr, ADirNode * dir) {
	if (!dir) return 0;
	ANode *up0 = dir->up();
	ADirNode * up = (up0 && up0->isADir()) ? static_cast<ADirNode*>(up0) : 0;
	c2(9, 'T');    hex5(up?up->id():0);
	c2(lr+49, 36); hex5(dir->id());     c2(lr+49, 36);
	m_bufp += dir->gui_list(m_bufp, 0);
	dir->winflg_or(2+2*lr); invd(); return 0;
}

void GuiStub::tree_force(int lr, ANode * dir0) {
	ERR_CAST(ADir, dir, dir0, ); if (dir->winflg(2+2*lr)) return;
	tree_force(lr, dir->up()); tree_expand(lr, dir); 
}

void GuiStub::tree_sel(int lr, ANode * nd) {
	tree_force(lr, nd->up()); invd();
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
	setwin(16*id+3,'K'); wupd_i1(fid, v01); }

void GuiStub::node_name(int i, ANode * nd) {
	ERR_CAST(ADir, dir, nd->up(), log("node_name: no updir"));
	c2(9, 'N'); hex5(dir->id()); c2(i+49, '+'+nd->isADir());
	hex5(nd->id()); c1(36); m_bufp += nd->get_name(m_bufp);
}

void GuiStub::node_rm(int i, ANode * nd) {
	ERR_CAST(ADir, dir, nd->up(), log("node_rm: no updir"));
	invd(); c2(9, 'N'); hex5(dir->id()); c2(i+49, '-'); hex5(nd->id());
}

int GuiStub::flush_all() { 
	cfl(); slr_upd(); errq_cfl();
	int x = glob_flg^m_gf0; m_gf0 = glob_flg;
	if (x & GLF_FSTATE) savename();
	if (m_gnaq_n && (snd0.total_played()-m_gnaq_t)>(long long)(44*CFG_GUI_TRKUP.i)) {
		for (int i=0, n=m_gnaq_n; i<n; i++) trk_w_gna(m_gnaq_id[i]);
		m_gnaq_t = snd0.total_played(); m_gnaq_n = 0; }
	return flush();
}

static int k_to_3p(char *q, int x) {
        if (x<=9999) return (x<=0) ? (q[0]=45+3*!x, 1) :
                ((x>999) ? (x/=10, q[0]=x/100+48, q[1]=46) : (q[0]=46, q[1]=x/100+48), d99(q+2,x%100), 4);
        if (x<=99999) return d99(q,x/1000), q[2]=46, q[3]=(x%1000)/100+48, 4;
        return (x<=999999) ? (d999(q, x/1000), 3) : (q[0]=77, 1);
}

void GuiStub::wupd_k3p(int wwt, int x, int wwix) { char s[8]; wupd(wwt,wwix); c1('t'); sn(s, k_to_3p(s, x)); }

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

void GuiStub::vol() { setwin(7,'.'); wupd_i2('v', snd0.vol() - 12); }

void GuiStub::midi(int flg) {
	if (flg&0x40000000) cre(0x27, 'M'); else setwin(0x27, 'M');
	if (flg&31) for (int i=(flg>>5)&31, i1=i+(flg&31); i<i1; i++)
			t0(), c1(i_to_b32(i)), m_bufp += midi_w_ln(m_bufp, i);
	if (flg&1024) t0(), c1('x'), m_bufp += midi_w_slg(m_bufp, 0);
	if (flg&2048) t0(), c1('y'), m_bufp += midi_w_slg(m_bufp, 1);
	if (flg&4096) t0(), c1('z'), hex4(midi_w_flg());
}

void GuiStub::fcfg_ud(int wch, cfg_ent * pc, const char *sdef, int ldef) {
	wupd_c0(wch, 't');
	if (pc->i) sn("%%%ttt", 6), sn(pc->s, pc->i); else sn("ppp666(", 7), sn(sdef, ldef), c1(')'); }

void GuiStub::fcfg_ex(int k) { wupd_s(k+48, (&CFG_XTERM+k)->s); }
void GuiStub::ocfg_draw() { cre(0x37, 'k'); for (const char *s = CFG_IVSUM; *s; s++) ocfg_l(*s); }
void GuiStub::ocfg_l(int c) { wupd_i(c, cfg_tab[c-48].i); }
void GuiStub::xapp_bv() { for (int k,i=0; i<N_XAPP; i++) if ((k=pt_xapp_bv[i])) c3(9, 'x', 48+i), hex8(k); }
void GuiStub::grc(int i, int j, double v) { wupd_c0('g', '#'); b32n(i, 2); c1(i_to_b32(j)); hdbl(v); }

void GuiStub::fcfg_draw() {
	cre(0x57, 'F');
	fcfg_ud('k', &CFG_AO_DIR,  QENV('t'), QENVL('t'));
	fcfg_ud('w', &CFG_WAV_DIR, QENV('h'), QENVL('h'));
	for (int i=0; i<N_XAPP; i++) wupd_s(48+i, (&CFG_XTERM+i)->s); }

void GuiStub::savename() { 
	char buf[64]; const char *p, *sb, *rgb = (glob_flg&GLF_EMPTY)?"%%%FFF":0, *sp = save_file_name;
	int rf = glob_flg & GLF_RECOVER, sf = glob_flg & GLF_SAVED;
	if (*sp) { if (!rgb) rgb="zz%%^^\0 %%%%^^"+8*!!sf;
		   for (p=sb=sp; *p; p++) if (*p=='/') sb=p+1; if (rf) sb="BUG!"; }
	else if (rf) { memcpy(buf,"(recover:0)",12); buf[9]=sp[2]; if (!rgb) rgb="zzz%%h"; sb=sp=buf; }
	else { sp = sb = rgb ? "(empty)" : (rgb = "zz%z%%\0 %zzz%%"+8*!!sf, "(unnamed)"); }
	setwin( 7,'.'); wupd_c0('N','t'); sn(rgb, 6); sz(sb); 
	setwin(23,'/'); wupd_c0('1','t'); sn(rgb, 6); sz(sp);
}

void GuiStub::gn_start(int i, BoxGen * bx, int ni, int no) {
	wupd_c0('g', 'n'); c4(i_to_b32((i+1)>>5), i_to_b32((i+1)&31), i_to_b32(ni), i_to_b32(no));
	if (bx) sn(bx->v_rgb(), 6), nname(bx->node()), c1(36); else sn("zz%z%%!BUG!$", 12); }

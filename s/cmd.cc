#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

#include <cstring>

#include "pt.h"
#include "util.h"
#include "util2.h"
#include "cmd.h"
#include "glob.h"
#include "pzrf.h"
#include "guistub.h"
#include "mx.h"
#include "asnd.h"
#include "midi.h"
#include "cfgtab.inc"

#define CMD_NODE(T) T##Node* nd = dynamic_cast<T##Node*>(p->m_c_node); if (!nd) return GCE_EX##T

#define CF_NODE 0x100
#define CF_TOK 0x200
#define CF_ARG1 0x400
#define CF_ARG2 0x800
#define CF_TOK1 0x600
#define CF_TOK12 0xe00
#define CF_NODE1 0x500
#define CF_NODE12 0xd00

class CmdTab {
	public:
		typedef int (dfun_t) (CmdBuf * p);
		typedef int (*fun_t) (CmdBuf * p);
		struct ent_t { int c; fun_t f; };
		static inline ent_t * ent(int c) { return m0_tab + m0_chtab[c & 127]; }
		static void init() { m0_chtab[0] = 1;
			for (int j,i=1; (j=m0_tab[i].c&127); i++) m0_chtab[j] = i; }
	protected:
		static ent_t m0_tab[];
		static char m0_chtab[128];

                static int c_invalid(CmdBuf * p) { return GCE_WHAT; }
                static int c_nop(CmdBuf * p) { return 0; }
                static int c_bye(CmdBuf * p) { bye((p->m_c_a0[0]&31)<<8); return GCE_WTF; }
                static int c_stfu(CmdBuf * p) { mx_clear(0); return 0; }
                static int c_debug(CmdBuf * p) { if (*p->m_c_a0=='?')
			log("debug_flg=%d (%s)", debug_flags, DFLG_HELP);
			else debug_flags = atoi_h(p->m_c_a0); return 0; }
                static int c_gui2(CmdBuf * p) { gui2.brk(); gui2.pf("%s\n", p->m_c_a0); return 0; }
		static int c_kfw(CmdBuf * p) { CMD_NODE(Clip); return nd->cmd(p); }
		static int c_dsc(CmdBuf * p) { CMD_NODE(ABox); return nd->cmd_H(p); }
		static int c_info(CmdBuf * p) { char *s = p->tok(), k = s ? *s&7 : 7; 
						return k ? (pt_con_op(0), p->m_c_node->debug(k), 0) : 0; }
		static int c_source(CmdBuf * p) { return p->fdok(CmdBuf::read_batch(p->m_c_a0,
								  p->m_c_nof|NOF_BATCHF), '<'); }
                static int c_save(CmdBuf * p) { return p->fdok(Node::save_batch(Node::root(), p->m_c_a0,
						          p->m_c_nof&NOF_FORCE),'>'&-(p->m_c_a0&&*p->m_c_a0));}
                static int c_sv2(CmdBuf * p) { CMD_NODE(ADir); return (!nd->id()) ? NDE_SLROOT :
			p->fdok(Node::save_batch(nd,p->m_c_a1, p->m_c_nof&NOF_FORCE), 'L'); }
		static int c_iofw(CmdBuf * p) { if(*p->m_c_a0=='f')fflush(stderr); return pt_iocmd(p->m_c_a0);}
		static int c_wkfw(CmdBuf * p) { char*s = p->m_c_a0; int r, l = strlen(s);  s[l] = 10;
						r = pt_wrk_cmd(s, l+1); s[l] = 0; return r; }
		static int c_snd(CmdBuf * p) { int k = *p->m_c_a0 - 48;
					       return k ? GCE_PARSE : snd0.cmd(p->m_c_a0+1); }
		static int c_midi(CmdBuf *p) { return midi_cmd(p->m_c_a0); }
		static dfun_t c_cre, c_vol, c_job, c_cfg, c_lib, c_misc, c_nadm, c_tree, c_wav, c_report,
			      c_efw, c_wfw, c_xfw, c_win, c_pfx, c_cont, c_live, c_stat, c_closewin, c_rpl;
};

char * CmdTok::t(char *s) {
	if (s) {
		for (m_s0=s; *s==36; s++);
		if (!*s) return m_p=s, (char*)0;
	} else if (!(s=m_p) || !*s) { return 0; }
	char *p = s+1; while (*p && *p!=36) ++p;
	while (*p==36) *(p++) = 0;
	m_p = p; return s;
}

char * CmdTok::un() {
	if (!m_s0) return 0;
	for (char *p=m_s0; p<m_p; p++) if (!*p) *p=36;
	return m_s0;
}

void CmdBuf::st_init() { CmdTab::init(); }

void CmdBuf::init(int fd, int nof, int px, const char * inm, int bs, int rs) {
	if (fd>=0) m_fd = fd; else if (fd<-31) m_fd = open(inm = tpipe_name(-fd), O_RDWR|O_APPEND);
	m_nof0 = nof; m_prefix = px; m_rsiz = rs; m_buf = (char*)malloc(m_bsiz = bs); 
	m_lineno = m_rpos = m_cpos = m_try = 0;
	memset(m_nd_var, 0, 8*sizeof(int));
	m_iname = inm ? strdup(inm) : 0;
	delete(m_cont); m_cont = 0;
}

int CmdBuf::read_batch(const char * name, int nof1) {
	Clock clk; clk.reset();
	int ec, fd = open(name, O_RDONLY);
	if (fd<0) return GCE_FOPEN;
	if (nof1&1) { if (!(glob_flg&GLF_EMPTY)) return NDE_EXPVIRG;  Node::lib_start(); }
	else if(!*save_file_name) {
		if ((ec=is_asv_name(name))) glob_flg|=GLF_RECOVER, save_file_name[2]=ec;
		else strncpy(save_file_name, name, 1023), glob_flg &= ~GLF_RECOVER;
		if (!(glob_flg & GLF_INI1)) gui2.savename(); }
	CmdBuf cb; cb.init(fd, nof1&(NOF_FLAGS&~NOF_FGUI), 0, name);
	do ec = cb.read_f(); while (ec>=0);
	if (nof1&1) Node::lib_end();
	if (ec==GCE_EOF) ec = (EEE_SUM &-!!cb.m_errcnt);
	close(fd); if (!ec) log("read_batch OK, t=%d", clk.get());    return ec;
}

int CmdBuf::bprep(int siz) {
        if (m_rpos + siz < m_bsiz) return siz;
        if (!m_cpos) return m_bsiz - m_rpos - 1;
        if (m_rpos > m_cpos) memmove(m_buf, m_buf+m_cpos, m_rpos-m_cpos);
        m_rpos -= m_cpos; m_cpos = 0;
        return min_i(siz, m_bsiz - m_rpos - 1);
}

int CmdBuf::read_f() {
        int r, len = bprep(m_rsiz);   if (len<1) return GCE_FULLBUF;
        if ((r = read(m_fd, m_buf+m_rpos, len))>0) return m_try=0, chunk(r); else ++m_try;
	if (m_nof0&NOF_EXPEOF) return p_close(&m_fd), 
		r ? (log("read_f(%s): %s", m_iname, strerror(errno)), GCE_FREAD) : GCE_EOF;
	log("read_f(): %s(%d): %s (%dx)", m_iname?m_iname:"<unnamed>", m_fd,
		                          r ? strerror(errno) : "EOF", m_try);
	if (!m_fd) return m_fd=-1, pt_con_op("-4");
	if (!m_iname || m_try>5) return p_close(&m_fd), r ? GCE_FREAD : GCE_EOF;
	return (close(m_fd), m_fd=open(m_iname,O_RDWR))<0 ? EEE_ERRNO : 0;
}

int CmdBuf::vpf(const char * fmt, va_list ap) {
        int len = bprep(m_rsiz); if (len<1) return GCE_FULLBUF;
	int r = vsnprintf(m_buf+m_rpos, len, fmt, ap);
	int ec = chunk(r); return ec < 0 ? ec : r;
}

int CmdBuf::sn(const char * s, int n) {
        int len = bprep(n); if (len<n) return GCE_FULLBUF;
        memcpy(m_buf+m_rpos, s, n); return chunk(n);
}

int CmdBuf::chunk(int len) {
        char *s = m_buf+m_rpos, *c0 = m_buf+m_cpos; if (!*s) return 0;
        int ec, r = 1;  s[len] = 0;
        if (s==c0 && *s<11) goto sep;
tok:    while (*s>10) ++s;
        if (*s==10) ++m_lineno, *(s++) = 0; else goto done;
        if (m_prefix && m_prefix != *(c0++))
                log("cmd: exp '%c'..., got \"%s\"", m_prefix, c0-1), ec = GCE_NOPFX;
        else ++r, ec = line(c0);
        if (ec<0 && (m_c_nof&NOF_THROW)) return ec;
sep:    while(*s==10) s++, m_lineno++;
        if (*(c0=s)>10) goto tok;
done:   if (s==c0) m_cpos = m_rpos = 0; else m_cpos = c0 - m_buf, m_rpos = s - m_buf;
        return r;
}

int CmdBuf::line(char * s) {
	if ( (debug_flags&DFLG_GUICMD) && is_gui()) log("gui: %s", s);
	m_c_nof = m_nof0;
	int ec, flg; CmdTab::ent_t * pe;
start: 	pe = CmdTab::ent(*s); ec = 0; flg = pe->c;
	if (flg&CF_NODE) {
		char *nm = tok(s); 
		if (!(m_c_node = nm[1] ? lookup(nm+1) : curnode())) {
			ec = nm[1] ? GCE_LOOKUP : GCE_ZNODE; goto err; }
	} else if (m_c_node=0, flg&CF_TOK) { m_c_a0 = tok(s) + 1; }
	  else { m_c_a0 = s+1; m_c_tok.rs(); goto run; }
	if (!(flg&CF_ARG1)) goto run;
	if (!(m_c_a1 = tok())) { flg &= ~(CF_ARG1|CF_ARG2); ec = GCE_ZARG1; goto err; }
	if ((flg&CF_ARG2) && !(m_c_a2 = tok())) { flg &= ~CF_ARG2; ec = GCE_ZARG2; goto err; }
run:    if ((ec = (*pe->f)(this)) >= 0) return ec;
        if (ec==GCE_PREFIX) { s = m_c_a0; goto start; }
err:    if (m_c_nof & NOF_ERRMSG) show_error(ec);
	if (ec == EEE_NOEFF) return 0; else ++m_errcnt;
	if (m_c_nof&NOF_FGUI) gui2.errq_add(ec);
	return ec;
}

int CmdBuf::fdok(int ec, int ty) {
	if ( (ec|(ty-1))>=0 && is_gui() ) gui2.c4(9, 'f', ty, 'Z');     return ec; }

void CmdBuf::show_error(int ec) {
	// TODO: NOF_FGUI, NDE_CONFIRM
	char * s0 = untok(); if (!s0) s0 = m_c_a0 - 1;
	log("%s: %d: \"%s\": %s (%d)", m_iname, m_lineno, s0, err_str(ec), ec);
}

int CmdBuf::cf_p2i(CmdBuf *p, char *q, int l) {
	q[l-1] = 0; ANode * nd = p->lookup(q+1); if (!nd) return GCE_RPCONV;
	return hx5(q, nd->id()); }

int CmdBuf::cf_i2p(CmdBuf *p, char *q, int l) {
	q[l-1] = 0; int upf = q[1]=='^';
	ANode *nd2, *nd = Node::lookup_n(atoi_h(q+1+upf)); if (!nd) return GCE_RPCONV;
	if (upf && (nd2=nd->up())) nd = nd2;
	return l = nd->get_path(q+1,999), q[l+1] = '`', l+2; }

int CmdBuf::rpl_cp(char *to, const char *s, CmdBuf::conv_fun f) {
	int k, r; if (!memcmp(s,"^<",2)) return r = strlen(s), memcpy(to, s, r), r;
	char *q0 = 0, *q = to;
	while(*s) { if ((*(q++) = *(s++)) != '`') continue;
		    if (!q0) q0 = q-1;
		    else if ((r = (*f)(this, q0, k=q-q0)) < 0) return r;
		    else q += r-k, q0 = 0; }
	return *q = 0, q - to;
}

////////////////////////////////////////////////////////////////////////////////////////////////

static void log_load(int n) { for (int i=0; i<n; i++) log("log load test ------------------------- #%d", i); }

static int ocfg_iv(int j, const char *s) {
	if ((unsigned int)(j-1) >= (unsigned int)(sizeof(cfg_tab)/sizeof(cfg_ent)-1)) return GCE_UCFG;
	cfg_ent * q = cfg_tab + j;
	if (!(intv_cmd_cfg(q, s, -1))) return 0;
	gui2.setwin(0x37,'k'); gui2.ocfg_l(j+48); return 0;
}

static int fcfg_dir(cfg_ent * q, char *s, int dq) {
	int l = strlen(s+1); if (l>254) return EEE_LONGSTR;
	cfg_setstr(q, s+1); 
	gui2.setwin(0x57, 'F'); gui2.fcfg_ud(*s, q, QENV(dq), QENVL(dq));
	gui2.c4(9, 'f', *s, 'Z'); return 0;
}

static int fcfg_exe(char *s) {
	int k = s[1]-48; if ((unsigned int)k >= (unsigned int)N_XAPP) return GCE_UCFG;
	int ec = cfg_setstr(&CFG_XTERM+k, s+2), l = (&CFG_XTERM+k)->i;
	s[l+2] = 10; pt_iocmd_sn(s, l+3); s[l+2] = 0;
	if ((*s&32) && (ec>=0)) gui2.setwin(0x57, 'F'), gui2.fcfg_ex(k);
	return ec;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int CmdTab::c_cont(CmdBuf * p) {
	if (!p->m_cont) return GCE_NOCONT;
	int ec = p->m_cont->line(p->m_c_a0); 
	if (ec<1) delete(p->m_cont), p->m_cont = 0; 
	return ec;
}

int CmdTab::c_cre(CmdBuf * p) { 
	ANode * nd; const char * s;
	int ty = *p->m_c_a1;
	if (ty==63) { if (!(s=p->tok())) return GCE_PARSE; else ty = *s; }
	if (!p->m_c_a1[1]) return NDE_ZNAME;
	int ec = Node::mk(&nd, p->m_c_node, p->m_c_a1+1, ty, p->m_c_nof|NOF_PARSE);
	if (ec>=0) p->set_curnode(nd);
	else if (ec!=NDE_NDUP || (ty|32)!='d') return p->set_curnode(0), ec;
	else if (s=p->m_c_a1+1, !(nd=p->m_c_node->sn(&s)) || nd->cl_id()!='D') return ec;
	else if (p->set_curnode(nd), ec=0, !(p->cnof()&NOF_OVRRD)) 
			log("warning: double mkdir \"%s\" ignored", p->m_c_a1+1);
	if (!(p->m_c_a1 = p->tok())) return 0;
	while (ec>=0 && p->m_c_a1) ec = p->curnode_ccmd(), p->m_c_a1 = p->tok();
	return ec;
}

int CmdTab::c_nadm(CmdBuf * p) {
	char *name, *arg = p->m_c_a1;
	ANode *trg, *nd = p->m_c_node;
	int ec, nof = p->m_c_nof;
	ADirNode *dir = dynamic_cast<ADirNode*> (nd);
	switch(arg?*arg:0) {
		case 'm': name = arg+1; ec = Node::parse_target(p, &name, &trg);
			  return ec<0 ? ec : Node::move(nd, trg, name, nof | NOF_PARSE);
		case 'c': name = arg+1; ec = Node::parse_target(p, &name, &trg);
			  return ec<0 ? ec : Node::copy(nd, trg, name, nof | NOF_PARSE);
		case 'd': return Node::del(nd, nof);
		case 'p': if (!dir) return GCE_EXADir;
			  for (int i=1;i<5;i++) if ((arg[i]&56)!=48) return GCE_EXOCT;
			  dir -> set_perm(8*(arg[1]&7)+(arg[2]&7), 8*(arg[3]&7)+(arg[4]&7));
			  return 0;
		case 'L': Node::set_lr(0, nd); return 0;
		case 'R': Node::set_lr(1, nd); return 0;
		case 'l': gui2.tree_sel(0, nd); return 0;
		case 'r': gui2.tree_sel(1, nd); return 0;
		case '+': return *p->m_nd_var = nd->id();
		case 'V': return nd->cl_id()=='w' ? (p->m_cont = wrap_avreader(nd->box0(), p->before(0,7)), 0) 
			  		          : GCE_EXWR;
		case 'W': return nd->draw_window(16);
		case '=': if (!(ec = arg[1] & 7)) return EEE_NOEFF;
			  return p->m_nd_var[ec] = nd->id(); 
		case 'F': return Node::lib_cfg(nd);
		case '?': return Node::obj_help(nd->cl_id()+256*(arg[1]=='?'));
		case 'I': return pt_con_op(0), p->m_c_node->debug(arg[1]&7), 0;
		case 'D': if ((ec=nd->show_dsc(1))>=0 || !(arg[1]&1)) return ec;
			  return (ec=nd->draw_window(16))>=0 ? ec : Node::obj_help(nd->cl_id());
		case 0: return GCE_ZARG0;
		default: return GCE_UNADM;
	}}

int CmdTab::c_tree(CmdBuf * p) {
	char * arg = p->m_c_a1;
	int ec, c0 = *arg, lr = c0>96 && (c0-=32, 1);
	switch(c0) {
		case 'W': if ((ec = p->m_c_node->draw_window(16))!=NDE_NOWIN) return ec;
		case 'E': ADirNode * dir; dir = dynamic_cast<ADirNode*>(p->m_c_node);
			  return (dir) ? gui2.tree_expand(lr,dir) : GCE_EXADir;
		default:  return GCE_UTREE;
	}}

int CmdTab::c_xfw(CmdBuf * p) { BoxGen * bx = p->m_c_node->box0(); return bx ? bx->cmd(p) : GCE_EXABox; }

int CmdTab::c_wfw(CmdBuf * p) {
	if (!(p->m_c_a1 = p->tok())) return p->m_c_node->draw_window(16);
	BoxGen * bx = p->m_c_node->box0(); return bx ? bx->gui_cmd(p) : GCE_EXABox; }

int CmdTab::c_efw(CmdBuf * p) { 
	CMD_NODE(ABox); int ec;
	if (*p->m_c_a1=='C') return (p->m_cont = nd->dsc_reader(p->m_c_a1[1]-48)) ? 0 : NDE_NOGUI;
	do { ec = nd->ui_cmd(p); } while (ec>=0 && (p->m_c_a1 = p->tok()));
	return ec;
}

int CmdTab::c_win(CmdBuf * p) {
	char * s = p->tok(); return p->m_c_node->draw_window(s?b32_to_i(*s):0); }

int CmdTab::c_vol(CmdBuf * p) {
	int skipgui = (*p->m_c_a0=='G');
        const char *s = p->m_c_a0 + skipgui;
	if (!*s) return GCE_ZARG0;
	int x = ivlim((s[1]&80) ? hex2(s) : hxd2i(*s), 0, 99);
	snd0.set_vol(x+12); if (!skipgui) gui2.vol();
	return 0;
}

int CmdTab::c_stat(CmdBuf * p) { 
	log("total_played: %g", sample_length * (double)snd0.total_played());
	return 0; // TODO: more debug
}

int CmdTab::c_cfg(CmdBuf * p) {
	char *s = p->m_c_a0; switch(*s) {
		case '@': return ocfg_iv(s[1]-48, s+2);
		case '>': return cfg_write(1);
		case 'k': return fcfg_dir(&CFG_AO_DIR,  s, 't');
		case 'w': return fcfg_dir(&CFG_WAV_DIR, s, 'h');
		case 'x': case 'X': return fcfg_exe(s);
		case 'V': return gui2.ocfg_draw(), 0;
		case 'W': return gui2.fcfg_draw(), 0;
		default:  return GCE_UCFG;
	}}

int CmdTab::c_report(CmdBuf * p) { switch(*p->m_c_a0) {
		case 'g': pt_calc_xbv(); gui2.xapp_bv();
			  glob_flg |= GLF_GUIOK; log("gui says hi"); return 0;
		case 'S': return qstat_op(p->m_c_a0[1]);
		case 'p': return snd0.pump_rprt(p->m_c_a0+1);
		default: return GCE_UREPORT;
	}}

int CmdTab::c_wav(CmdBuf * p) {
	int id = atoi_h(p->m_c_a0), op = atoi_h(p->m_c_a1);
	const char *s1 = 0, *s2 = 0;
	if ((op|4)==4){ if ((s1=p->tok())) { while (*s1==32) ++s1; if (!*s1) s1 = 0; } 
		 	if ((s2=p->tok())) { while (*s2==32) ++s2; if (!*s2) s2 = 0; } 
			op += s1 ? (1+!!s2) : 2*(s2 && (s1="0")); }
	log("c_wav 0x%x 0x%x '%s' '%s'", id, op, s1?s1:"NULL", s2?s2:"NULL");
	return pt_acv_op(id, op, s1, s2);
}

int CmdTab::c_misc(CmdBuf * p) {
	int i,j; char *s = p->m_c_a0; switch(*s) {
		case 'F': pzrf_show_last(); return 0;
		case ':': p->m_nof0 = p->m_c_nof; return 0;
		case 'w': i = s[1]&1 ? snd0.hcp_start(209203*CFG_AO_RECLIM.i-818) : snd0.hcp_end();
			  gui2.setwin(7,'.'); gui2.wupd_i1('W',!!snd0.hcp()); return i;
		case 'c': return pt_con_op(s+1); 
		case 'L': log_load(atoi(s+1)); return 0;
		case 'l': return pt_show_lic();
		case 'V': if (sscanf(s+1, "%d.%d", &i, &j)!=2) return GCE_PARSE;
		          if (i!=v_major || j!=v_minor) log("save file is from version %d.%d", i, j);
			  p->m_sv_M = i; p->m_sv_m = j; return 0;
		case 'K': return cfg_write(1);
		case 'M': return mx_debug(s+1);
		case 'm': return nd_mem_debug();
		case 'p': return pt_debug(), 0;
		case 's': return i=(int)lround(1e3*atof(s+1)), u_sleep(((s[1]-43)&~2) ? i:i+(clk0.t()/1000)),0;
		case 'G': switch(j=s[1], i=j?atoi_h(s+2):0, j) {
				  case '-': i=~i; case '&': i &= glob_flg; break;
				  case '+':       case '|': i |= glob_flg; break; default:i=glob_flg; break; }
			  return j=glob_flg, glob_flg=i, log("glob_flg: 0x%x -> 0x%x", j, i), 0;
		case 'E': return gui2.errq_add(s[1] ? -atoi_h(s+1) : -32, s[1]?"_E":0), 0;
		case 'W': return ANode::wi_debug(), 0;
		case 'R': if ( !!(glob_flg&GLF_RECORD) == (j = s[1]&1) ) return EEE_NOEFF;
		          glob_flg ^= GLF_RECORD; gui2.pf("\tdf%x", j<<11);
			  log("cmd/event recording st%sed", "opp\0art"+4*j);
			  if (j) log("cmd_rec:0000 ^_Sn%s", s[2]?s+2:"63"),  qstat_cfg('N', s+2);
			  else   fflush(stderr), i='f', pt_iocmd((char*)&i), qstat_cfg('-', 0);   return 0;
		case 'S': switch(s[1]) { case 'n': return qstat_cfg('n', s+2);
				  	 case 'c': return qstat_chk0(&p->m_cont, s+2);
					 default : return GCE_PARSE; }
		case 'D': return (j=s[1]) ? (strncpy(DEBUG_UTXT(hxd2i(j)), s+2, 63), 0) : GCE_PARSE;
		case '>': return log("_> %s", s+1), 0;
		case '@': return clk0.wrk(s[1]);
		default: return GCE_UMISC;
	}}

int CmdTab::c_closewin(CmdBuf * p) {
	char * s = p->m_c_a0;
	if (!(*s & 80)) return GCE_WRONGWT;
	int ty = hxd2i(*(s++)); if (ty==7) return 0; // TODO
	while (1) { if (*s==36) ++s; if (!*s) return 0;
		    ANode * nd = Node::lookup_n(atoi_h_pp(&s));
		    if (nd) nd -> close_window(ty);
	}}

int CmdTab::c_job(CmdBuf * p) {
	char *s = p->m_c_a0, *a1 = p->tok();
	ANode *nd = 0; int n = 4096;
	if (*s=='D') return jobq.debug(), 0; else if (!a1) return GCE_ZARG1;
	switch(*s) {
		case 'X': return jobq.cmd(0, atoi_h(s+1), a1);
		case '_': return (n |= b32_to_i(s[1])) < 0 ? JQE_PARSE : jobq.cmd(0, n, a1);
		case  0 : if (!(nd=p->curnode())) return GCE_ZNODE; else break;
		default:  if (!(nd=p->lookup(s))) return GCE_LOOKUP; else break;
	}
	return (n = *a1) ? jobq.cmd(nd, (n&7)-(n&8), a1+1) : JQE_PARSE;
}

int CmdTab::c_lib(CmdBuf *p) {
	char *s = p->m_c_a0;
	switch(*s) {
		case 'W': return (glob_flg&GLF_EMPTY) ? (gui2.sn("\tflW",4), 0) : NDE_EXPVIRG;
		case '.': return p->fdok(CmdBuf::read_batch(s+1, p->m_c_nof|NOF_BATCHF|1), 'l');
		default:  return GCE_ULIB; }}

int CmdTab::c_pfx(CmdBuf * p) {  // --RC TEFG
	static const int pfxtab[8] = {0,0,NOF_OVRRD,NOF_YES,NOF_THROW,NOF_STRICT,NOF_FORCE,NOF_FGUI};
	int x, y;
   	if ( !((x = *(p->m_c_a0++)) & 64) || !(y = pfxtab[x&7]) ) return GCE_UPREFIX;
	if (x&32) p->m_c_nof &= ~y; else p->m_c_nof |= y;
	return GCE_PREFIX;
}

int CmdTab::c_live(CmdBuf *p) {
	const char * s = p->m_c_a0;
	for (int i=0; i<7; i++) if (!s[i]) return GCE_PARSE;
	log("c_live: id:0x%04x, ix:0x%02x, v:0x%02x", qh4rs(s), hex2(s+4), s[7] ? hex2(s+6) : hxd2i(s[6]));
	mx_l_op(qh4rs(s), hex2(s+4), s[7] ? hex2(s+6) : hxd2i(s[6]));
	return 0;
}

int CmdTab::c_rpl(CmdBuf *p) {
	static long long rec_ts = 0ll;
	const char * s = p->m_c_a0;
	if ((*s|2) != 'R') return (*s=='_') ? (log("cmd_rec:%s", s+1), 0) : GCE_PARSE;
	int l, rf = *s & 2;
	char buf[1024];
	if(rf){ if ((l=p->rpl_cp(buf, s+1, CmdBuf::cf_i2p)) < 0) return l;
		long long ttp = snd0.total_played();
		int dif = rec_ts ? (int)(ttp - rec_ts) : 0;  rec_ts = ttp;
		dif = (dif<2888999) ? ((dif*743)>>15) : 65535;
		log("cmd_rec:%04x %s", dif, buf); return 0;
	}else { if ((l=p->rpl_cp(buf+1, s+1, CmdBuf::cf_p2i)) < 0) return l; else buf[0] = 9;
		gui2.sn(buf, l+1); return 0; }}

char CmdTab::m0_chtab[128];
CmdTab::ent_t CmdTab::m0_tab[] = {
{'?',c_invalid}, {'#',c_nop}, {'q',c_bye}, {'r',c_stfu}, {'v',c_vol}, {'^',c_gui2}, {'c', c_cfg}, {'l',c_lib},
{'/',c_stat}, {'C'|CF_NODE1,c_cre}, {'N'|CF_NODE1,c_nadm}, {'I'|CF_NODE,c_info}, {'S'|CF_NODE1,c_sv2},
{'K'|CF_NODE,c_kfw}, {'D'|CF_NODE1,c_tree}, {'X'|CF_NODE1,c_xfw}, {'J'|CF_TOK,c_job}, {':',c_pfx}, {'Q',c_rpl},
{'W'|CF_NODE,c_wfw}, {'E'|CF_NODE1,c_efw}, {'M'|CF_NODE,c_win}, {'s',c_save}, {'*', c_iofw}, {'L', c_live},
{'d',c_debug}, {'_',c_misc}, {'<',c_cont}, {'.', c_source}, {'x', c_closewin}, {'A', c_snd}, {'R', c_report},
{'w'|CF_TOK1,c_wav}, {'H'|CF_NODE, c_dsc}, {'m', c_midi}, {'%', c_wkfw}, {0,0} };

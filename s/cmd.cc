#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

#include <cstring>

#include "util.h"
#include "util2.h"
#include "cmd.h"
#include "glob.h"
#include "pzrf.h"
#include "gp.h"
#include "guistub.h"
#include "mx.h"
#include "asnd.h"
#include "pt.h"
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
                static int c_bye(CmdBuf * p) { bye((p->m_c_a0[0]==33)<<8); return GCE_WTF; }
                static int c_stfu(CmdBuf * p) { mx_clear(0); return 0; }
                static int c_debug(CmdBuf * p) { if (*p->m_c_a0=='?')
			log("debug_flg=%d (%s)", debug_flags, DFLG_HELP);
			else debug_flags = atoi_h(p->m_c_a0); return 0; }
                static int c_gui2(CmdBuf * p) { gui2.brk(); gui2.pf("%s\n", p->m_c_a0); return 0; }
		static int c_kfw(CmdBuf * p) { CMD_NODE(Clip); return nd->cmd(p); }
		static int c_info(CmdBuf * p) { char * s = p->tok(); p->m_c_node->debug(s ? *s&7 : 7); return 0; }
		static int c_source(CmdBuf * p) { return CmdBuf::read_batch(p->m_c_a0, p->m_c_nof); }
                static int c_save(CmdBuf * p) { return Node::save_batch(Node::root(), p->m_c_a0, p->m_c_nof&NOF_FORCE); }
                static int c_sv2(CmdBuf * p) { CMD_NODE(ADir); return Node::save_batch(nd, p->m_c_a1, p->m_c_nof&NOF_FORCE); }
		static int c_iofw(CmdBuf * p) { return pt_iocmd(p->m_c_a0); }
		static int c_snd(CmdBuf * p) { int k = *p->m_c_a0-48; return k ? GCE_PARSE : snd0.cmd(p->m_c_a0+1); }
		static dfun_t c_cre, c_vol, c_misc, c_nadm, c_closewin, c_tree, c_stat, c_cfg,
			      c_efw, c_wfw, c_xfw, c_win, c_job, c_pfx, c_cont, c_live;
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
	m_fd = fd>-32 ? fd : open(inm = tpipe_name(-fd), O_RDWR|O_APPEND);
	m_nof0 = nof; m_prefix = px; m_rsiz = rs; m_buf = (char*)malloc(m_bsiz = bs); 
	m_lineno = m_rpos = m_cpos = 0;
	memset(m_nd_var, 0, 8*sizeof(int));
	m_iname = inm ? strdup(inm) : 0;
	delete(m_cont); m_cont = 0;
}

int CmdBuf::read_batch(const char * name, int nof1) {
	int ec, fd = open(name, O_RDONLY);
	if (fd<0) return GCE_FOPEN;
	if (nof1&1) Node::lib_start();
	CmdBuf cb; cb.init(fd, nof1&NOF_FLAGS);
	do ec = cb.read_f(); while (ec>=0);
	if (nof1&1) Node::lib_end();
	if (ec==GCE_EOF) ec = 0;
	close(fd); return ec ? ec : (EEE_SUM &-!!cb.m_errcnt);
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
        return ((r = read(m_fd, m_buf+m_rpos, len))>0) ? chunk(r)
		:( log("read_f(): %s: %s (%cx)", m_fd?(m_iname?m_iname:"(unnamed)"):"(stdin)",
		                                 r ? strerror(errno) : "EOF", ++m_try+48),
		   (m_fd ? (close(m_fd), (m_iname&&m_try<5&&(m_fd=open(m_iname,O_RDWR))>=0))
			 : (m_try<5 || (close(0),0))) ? 0 : (m_fd=-1, r ? GCE_FREAD:GCE_EOF) );}

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
	if (is_gui()) gui2.errq_add(ec);
	return ec;
}

void CmdBuf::show_error(int ec) {
	// TODO: NOF_FGUI, NDE_CONFIRM
	char * s0 = untok(); if (!s0) s0 = m_c_a0 - 1;
	log("%d: \"%s\": %s (%d)", m_lineno, s0, err_str(ec), ec);
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
	if (ty==63) ty = ((s=p->tok())) ? *s : 0; 
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
		case 'V': return nd->cl_id()=='w' ? (p->m_cont = wrap_avreader(nd->box0()), 0) : GCE_EXWR;
		case 'W': return nd->draw_window(16);
		case '=': if (!(ec = arg[1] & 7)) return EEE_NOEFF;
			  return p->m_nd_var[ec] = nd->id(); 
		case 'F': return Node::lib_cfg(nd);
		case '?': return Node::obj_help(nd->cl_id()+256*(arg[1]=='?'));
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
	do { ec = nd->ui_cmd(p); } while (ec>=0 && (p->m_c_a1 = p->tok()));
	return ec;
}

int CmdTab::c_win(CmdBuf * p) {
	char * s = p->tok(); return p->m_c_node->draw_window(s?b32_to_i(*s):0); }

int CmdTab::c_vol(CmdBuf * p) {
	int skipgui = (*p->m_c_a0=='G');
        const char *s = p->m_c_a0 + skipgui;
	if (!*s) return GCE_ZARG0;
	int x = (s[1]&80) ? hex2(s) : hxd2i(*s);
	if ((unsigned int)x > 99u) x = 0;
	snd0.set_vol(x+12);
	if (!skipgui) gui2.setwin(7, '.'), gui2.wupd_i2('v', x);
	return 0;
}

int CmdTab::c_stat(CmdBuf * p) { 
	// TODO: mx_debug
	log("total_played: %g", sample_length * (double)snd0.total_played());
	unsigned int buf[16384];
	int n = snd0.stat(buf, 16384);
	char nm[16];
	nm [ sprintf(nm, "t%08x.bin", (unsigned int)time(0)) ] = 0;
	int fd = creat(nm, 0644);
	if (fd<0) return EEE_ERRNO;
	write(fd, buf, n*sizeof(unsigned int));
	close(fd);
	return 0;
}

static void log_load(int n) { for (int i=0; i<n; i++) log("log load test ------------------------- #%d", i); }

int CmdTab::c_cfg(CmdBuf * p) {
	const char *s = p->m_c_a0; int f = 0; switch(*(s++)) {
		case '>': return cfg_write();
		case 's': CFG_SV_EXEC.i = *s&1; f = 4; break;
		case 'a': intv_cmd(&CFG_ASV_MIN.i, s, 0, 35, 0x33330501); f = 8; break;
		case 't': CFG_TLOG_AUTO.i = *s&1; f = 16; break;
		case 'S': intv_cmd(&CFG_SV_BACKUP.i,   s, 0, 9); f =  32; break;
		case 'A': intv_cmd(&CFG_ASV_BACKUP.i,  s, 0, 9); f =  64; break;
		case 'T': intv_cmd(&CFG_TLOG_BACKUP.i, s, 0, 9); f = 128; break;
		case 'k': return (f=strlen(s))<255 ? (cfg_setstr(&CFG_AO_DIR, s),0) : EEE_LONGSTR;
		case 'd': CFG_DEVEL.i = *s&1; f = 512; break;
		case 'W': f = -1; break;
		default:  return GCE_UCFG;
	}
	return gui2.mcfg_win(f), 0;
}

int CmdTab::c_misc(CmdBuf * p) {
	int i,j; switch(*p->m_c_a0) {
		case 'g': Gnuplot::sg() -> restart(); return 0;
		case 'F': pzrf_show_last(); return 0;
		case ':': p->m_nof0 = p->m_c_nof; return 0;
		case 'w': i = p->m_c_a0[1]&1 ? snd0.hcp_start(20903400) : snd0.hcp_end();
			  gui2.setwin(7,'.'); gui2.wupd_i1('W',!!snd0.hcp()); return i;
		case 'W': return sscanf(p->m_c_a0+1, "%x:%x", &i, &j)==2 ? pt_acv_op(i, j) : GCE_PARSE;
		case 'c': return pt_con_op(atoi_h(p->m_c_a0+1)); // TODO
		case 'L': log_load(atoi(p->m_c_a0+1)); return 0;
		case 'l': return pt_show_lic();
		case 'V': if (sscanf(p->m_c_a0+1, "%d.%d", &i, &j)!=2) return GCE_PARSE;
		          if (i!=v_major || j!=v_minor) log("save file is from version %d.%d", i, j);
			  p->m_sv_M = i; p->m_sv_m = j; return 0;
		case 'K': return cfg_write();
		case 'M': return mx_debug(p->m_c_a0+1);
		case '_': return midi_cmd(p->m_c_a0+1);
		case 's': return u_sleep((int)lround(1000.0*atof(p->m_c_a0+1))), 0;
		default: return GCE_UMISC;
	}}

int CmdTab::c_closewin(CmdBuf * p) {
	char * s = p->m_c_a0;
	if (!(*s & 80)) return GCE_WRONGWT;
	int ty = hxd2i(*(s++)), msk = ~(1<<ty);
	while (ty!=7) {
		if (*s==36) ++s;
		if (!*s) return 0;
		int id = 0; while (*s & 80) id = 16*id + hxd2i(*(s++));
		ANode * nd = Node::lookup_n(id);
		if (nd) nd -> winflg_and(msk);
	}
	return 0; // TODO (?)
}

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

char CmdTab::m0_chtab[128];
CmdTab::ent_t CmdTab::m0_tab[] = {
{'?',c_invalid}, {'#',c_nop}, {'q',c_bye}, {'r',c_stfu}, {'v',c_vol}, {'^',c_gui2}, {'c', c_cfg},
{'/',c_stat}, {'C'|CF_NODE1,c_cre}, {'N'|CF_NODE1,c_nadm}, {'I'|CF_NODE,c_info}, {'S'|CF_NODE1,c_sv2},
{'K'|CF_NODE,c_kfw}, {'D'|CF_NODE1,c_tree}, {'X'|CF_NODE1,c_xfw}, {'J'|CF_TOK,c_job}, {':',c_pfx},
{'W'|CF_NODE,c_wfw}, {'E'|CF_NODE1,c_efw}, {'M'|CF_NODE,c_win}, {'s',c_save}, {'*', c_iofw}, {'L', c_live},
{'d',c_debug}, {'_',c_misc}, {'<',c_cont}, {'.', c_source}, {'x', c_closewin}, {'A', c_snd}, {0,0} };


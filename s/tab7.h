#ifndef __qwe_tab7_h__
#define __qwe_tab7_h__

#define TAB7_DECL(NM, TY) \
        typedef struct { int x; TY f; } NM##_ent_t; \
        static char m0_##NM##_ctab[128]; \
        static NM##_ent_t m0_##NM##_tab[]; \
        static inline NM##_ent_t * NM##_ent(int c) { return m0_##NM##_tab + m0_##NM##_ctab[c&127]; } \
        static inline int NM##_ent_i(int c) { return m0_##NM##_ctab[c&127]; } \
        static void NM##_init()

#define TAB7_DEF(CL, NM) \
void CL :: NM##_init() { \
        for (int k, i=0; (k = m0_##NM##_tab[i].x&127); i++) m0_##NM##_ctab[k] = i; } \
char CL :: m0_##NM##_ctab[128]; \
CL :: NM##_ent_t CL :: m0_##NM##_tab[] =

#define BXCMD_DECL(CL) \
        typedef int (*cmd_t) (CL *, const char *, CmdBuf *); \
        typedef int (dcmd_t) (CL *, const char *, CmdBuf *); \
	TAB7_DECL(cmd, cmd_t); \
	virtual int cmd(CmdBuf* cb); \
        static dcmd_t

#define BXCMD_DEF(CL) int CL :: cmd(CmdBuf* cb) { \
	char *s = cb->a1(); \
	cmd_ent_t * ep = cmd_ent(*s); \
	int ec, x = ep->x; if (x&8192) return BXE_UCMD; \
	if (!(x&256) && !cb->cperm(DF_EDBOX)) return NDE_PERM; \
	if ((ec = (*ep->f)(this, s, cb))>0) { \
		if ((ec&1) && m_node->winflg(2048)) box_window(); \
		if (ec&2) m_node->nio_change(); \
	} return ec; } \
	TAB7_DEF(CL, cmd)

#define BXCMD_H(CL, NM) int CL::c_##NM(CL *p, const char * s, CmdBuf * cb)

#endif // __qwe_tab7_h__

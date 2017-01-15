#include "util.h"
#include "job.h"
#include "guistub.h"
#include "glob.h"
#include "cfgtab.inc"
#include "asnd.h"

#define JF_LU(E,X) int j = lu(nd, ix); return j<0 ? (E) : (X)
Job* JobQ::job (ANode * nd, int ix)	   { JF_LU(0,	 m_ent[j].p); }
int  JobQ::jst (ANode * nd, int ix) 	   { JF_LU(1016, m_ent[j].st); }
int  JobQ::wake(ANode * nd, int ix, int k) { JF_LU(j,    m_ent[k].st &= ~(1024<<(k&3))); }
int  JobQ::kill(ANode * nd, int ix, int e) { JF_LU(j,    kill_j(j, e)); }

void JobQ::init() { 
	m_nj = 0; m_ent_bv = 0u; m_last_upd = -999999; m_upd_t = 11025;
	for (int i=0; i<32; i++) m_ent[i].jid = 0xfe0 + i;
	m_jlfun[0] = &asnd_mkjob; for (int i=1; i<16; i++) m_jlfun[i] = &jl_dummy;
}

int JobQ::launch(ANode * nd, int ix, char * arg) {
	int i = __builtin_ffs(~m_ent_bv) - 1;
	if ((i | (31-m_nj)) < 0) return JQE_FULL;
	ent_t * p = m_ent + i; m_upd_flg = 1;
	int ji0 = p->jid, ji = (ji0+32) & 4095, i4 = ix&15;
	p->p   = 0; p->st = 0; p->xst = 9999; p->jid = ji; p->i4f = i4;
	int ec = nd ? (p->nid = nd->id(), nd->start_job_2(p, arg))
		    : (p->nid = -1	, (*m_jlfun[i4]) (p, arg));
	if (ec<0) return p->jid = ji0, ec;
	int j    = m_nj++,
	    prio = (p->plttwwii >> 23) & 224;
	m_ent_bv |= 1u << i;
	while (j>0 && m_px[j-1] > prio+31) m_px[j] = m_px[j-1], --j; 
	m_px[j] = (unsigned char) (prio + (ji&31));
	return ji;
}


int JobQ::kill_j(int j, int ec) { ent_t *p = m_ent+j; j = (p->p && (p->st&1023)<1005); m_upd_flg = 1;
				  return j ? (p->p->abort(), p->st = ec&1023, 0) : EEE_NOEFF; }

int JobQ::kill5(int j) {
	if (j<0) return 0; else j &= 31;
	ent_t *p = m_ent+j; if (p->p) delete(p->p), p->p = 0;
	if (p->nid>=0) ANode::lookup_n_q(p->nid)->set_job_result(p->i4f, p->st), p->nid = -1;
	p->jid |= 4096; return 1;
}

int JobQ::lu(ANode * nd, int ix) {
	int j, k;
	if (!nd) return j = ix&31, k = m_ent[j].jid, (ix>4095 ? k<4096 : k==ix) ? j : JQE_LU;
	if (!nd->winflg(WF_JOBF)) return JQE_LU;
	return j = nd->winflg(WF_JOB5) >> 24, (ix<0 || (m_ent[j].i4f&7)==ix) ? j : JQE_LU;
}

int JobQ::cmd(ANode * nd, int ix, char * arg) {
	if (*arg=='+') return nd ? nd->start_job(ix, arg+1) : JQE_PARSE;
	int j = lu(nd, ix);
	if (*arg=='-') return j<0 ? JQE_LU : kill_j(j, JQE_KILL); 
	if (j<0) return j;
	Job * jp = m_ent[j].p; if (!jp) return JQE_ZOMBIE;
	int ec = jp->cmd(arg); if (ec<32767) return ec;
	if (((m_ent[j].st = ec&32767) & 1023) > 1004) m_upd_flg = 1;
	return ec;
}

bool JobQ::run() {
	for (int j,i=0; i<m_nj; i++) {
		ent_t * p = m_ent + (j = m_px[i]&31);
		if (p->st>1004) continue;
		int ec = p->p->run1();
		if (ec<0) return p->st = ec<-19 ? 1019 : ec&1023, m_upd_flg = 1;
		if ((ec&1023) > 1004)  return    p->st = ec&1023, m_upd_flg = 1;
		return p->st = ec & 32767, 1;
	}
	return 0;
}

void JobQ::upd_gui(int force) {
	if (!(force | m_upd_flg) && snd0.total_played() < m_last_upd+m_upd_t) return;
	m_last_upd = snd0.total_played(); m_upd_flg = 0;
	for (int i=0; i<m_nj; i++) {
		ent_t * p = m_ent + (m_px[i] & 31);
		if (force || p->xst!=p->st) upd_gui_1p(p);
	}}

void JobQ::upd_gui_1p(JobQ::ent_t * p) {
	int k = p->plttwwii, ni = p->nid;
	IFDBGX(JQ) log("upd_gui_1p: plttwwii=0x%x, nid=0x%x, st=0x%x", k, ni, p->st);
	if (!(k&0xff00) || ni<0) return;
	int kl = (k>>24)&15, wi = k&255, oid = ni<8 ? 16*kl+23 : 16*ni+kl;
	if (!Node::chkwin(oid)) return;
	gui2.setwin(oid, (k>>16)&127);
	gui2.j_upd((k>>8)&127, (p->xst=p->st), wi==255 ? -1 : wi);
}

int JobQ::purge() {
	int j, i1=0, naj = 0;
	for (int x,i0=0; i0<m_nj; i0++){ ent_t * p = m_ent + (j = m_px[i0] & 31);
					 if (((x=p->st)&1023)<1005) m_px[i1++] = m_px[i0], naj+=(x<1005);
					 else kill5(j), m_ent_bv &= ~(1u << j);	}
	return m_nj = i1, naj;
}

void JobQ::debug() {
	log("jobq debug: t=%lld, last=%lld, maxdiff=%d", snd0.total_played(), m_last_upd, m_upd_t);
	for (int i=0; i<31; i++) 
		log("%02d: (%c) jid 0x%x, nid 0x%x, st:%d, xst:%d",
				i, 45-2*!!(m_ent_bv&(1u<<i)), m_ent[i].jid,
				m_ent[i].nid, m_ent[i].st, m_ent[i].xst);
}

#include "util.h"
#include "job.h"
#include "guistub.h"
#include "glob.h"
#include "asnd.h"

Job * JobQ::job(ANode * nd, int ix) { int j = lu(nd,ix); return j<0 ? 0  : m_ent[j].p; }
int   JobQ::jst(ANode * nd, int ix) { int j = lu(nd,ix); return j<0 ?1016: m_ent[j].st; }

void JobQ::init() { 
	m_nj = 0; m_ent_bv = 0u; m_last_upd = -999999; m_upd_t = 11025;
	for (int i=0; i<32; i++) m_ent[i].jid = 0xfe0 + i; }

int JobQ::launch(ANode * nd, int ix, char * arg) {
	int i = __builtin_ffs(~m_ent_bv) - 1;
	if ((i | (31-m_nj)) < 0) return JQE_FULL;
	ent_t * p = m_ent + i;
	int ji0 = p->jid, ji = (ji0+32) & 4095;
	p->p   = 0;  p->nid = nd->id(); p->st = 0; p->xst = 9999;
	p->jid = ji; p->i4f = ix & 15; p->plttwwii = 0x60000000;
	int ec = nd->start_job_2(p, arg);
	if (ec<0) return p->jid = ji0, ec;
	int j    = m_nj++,
	    prio = (p->plttwwii >> 23) & 224;
	m_ent_bv |= 1u << i;
	while (j>0 && m_px[j-1] > prio+31) m_px[j] = m_px[j-1], --j; 
	m_px[j] = (unsigned char) (prio + (ji&31));
	return ji;
}

int JobQ::kill5(int j, int ec) {
	if (j<0) return 0; else j &= 31;
	ent_t * p = m_ent + j;
	if (p->p && (p->st&1023)<1005) p->p->abort(), p->st = ec&1023; 
	if (p->p) delete(p->p), p->p = 0;
	if (p->nid>=0) ANode::lookup_n_q(p->nid)->set_job_result(p->i4f, p->st),
		       jobq.upd_gui_1p(m_ent), p->nid = -1;
	p->jid |= 4096;
	return 1;
}

int JobQ::lu(ANode * nd, int ix) {
	int j, k;
	if (!nd) return j = ix&31, k = m_ent[j].jid, (ix&4096 ? k<4096 : k==ix) ? j : JQE_LU;
	if (!nd->winflg(WF_JOBF)) return JQE_LU;
	return j = nd->winflg(WF_JOB5) >> 24, (ix<0 || (m_ent[j].i4f&7)==ix) ? j : JQE_LU;
}

int JobQ::cmd(ANode * nd, int ix, char * arg) {
	if (*arg=='+') return nd ? nd->start_job(ix, arg+1) : JQE_PARSE;
	int j = lu(nd, ix);
	if (*arg=='-') return kill5(j); 
	if (j<0) return j;
	Job * jp = m_ent[j].p; if (!jp) return JQE_ZOMBIE;
	int ec = jp->cmd(arg); if (ec<32767) return ec;
	if (((m_ent[j].st = ec&32767) & 1023) > 1004) kill5(j);
	return ec;
}

bool JobQ::run() {
	for (int j,i=0; i<m_nj; i++) {
		ent_t * p = m_ent + (j = m_px[i]&31);
		if (p->st>1004) continue;
		int ec = p->p->run1();
		if (ec<0) return p->st = ec<-19 ? 1019 : ec&1023, kill5(j), 1;
		if ((ec&1023) > 1004)  return    p->st = ec&1023, kill5(j), 1;
		return p->st = ec & 32767, 1;
	}
	return 0;
}

bool JobQ::need_upd() {
	if (snd0.total_played() >= m_last_upd + m_upd_t) return true;
	int j;
	for (unsigned int m = m_ent_bv; (j=__builtin_ffs(m)-1) >= 0; m &= ~(1u<<j)) {
		if (m_ent[j].xst==m_ent[j].st) continue;
		if (m_ent[j].xst==9999 || m_ent[j].st > 1004) return true;
	}
	return false; 
}

void JobQ::upd_gui(bool force) {
	if (!force && !need_upd()) return;
	m_last_upd = snd0.total_played();
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
	int j, i1=0; for (int i0=0; i0<m_nj; i0++) {
		ent_t * p = m_ent + (j = m_px[i0] & 31);
		if ((p->st&1023)<1005) m_px[i1++] = m_px[i0];
		else kill5(j), m_ent_bv &= ~(1u << j);
	}
	return (m_nj = i1);
}

void JobQ::debug() {
	log("jobq debug: t=%lld, last=%lld, maxdiff=%d", snd0.total_played(), m_last_upd, m_upd_t);
	for (int i=0; i<31; i++) 
		log("%02d: (%c) jid 0x%x, nid 0x%x, st:%d, xst:%d",
				i, 45-2*!!(m_ent_bv&(1u<<i)), m_ent[i].jid,
				m_ent[i].nid, m_ent[i].st, m_ent[i].xst);

}

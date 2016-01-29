#ifndef __qwe_box_h__
#define __qwe_box_h__

#include "node.h"
#include "box0.h"

#define BXW_GET     sthg * bxw_rawptr =      m_node->wdat_raw(); if (!bxw_rawptr) return BXE_NOWDAT
#define BXW_GETV(S) sthg * bxw_rawptr =      m_node->wdat_raw(); if (!bxw_rawptr) return log("BUG: %s: no wdat", (S))
#define BXW_GETP(P) sthg * bxw_rawptr = (P)->m_node->wdat_raw(); if (!bxw_rawptr) return BXE_NOWDAT
#define BXSV2_HEAD if (!sv) return BXE_NOARG; int ec, r = 1; AOBuf * f = sv->out

class CmdBuf;
class BoxGen {
	public:
		friend void boxg_init();
		inline static ABoxNode* node0(const BoxGen *p) { return p?p->m_node:0; }
		BoxGen() : m_model(0), m_node(0) {}
		BoxGen(ABoxNode * nd) : m_model(0), m_node(nd) {}
		virtual ~BoxGen();
		virtual int n_in() const = 0;
		virtual int n_out() const = 0;
		virtual int in_mask() { return 0; }
		virtual int cmd(CmdBuf* cbf) { return BXE_NOCMD; }
		virtual int gui_cmd(CmdBuf* cbf) { return BXE_NOGCMD; }
		virtual const char * cl_name() = 0;
		virtual void upd_subbox(BoxGen * sub) {}
		virtual int cond_del() { return 0; }
		virtual int save2(SvArg * sv) { return 1; } // nonemp:save2(0)==BXE_NOARG
		virtual bool io_alias_perm() const = 0;
		virtual void spec_debug();
		virtual void notify_nio(BoxGen * bx) {}
		virtual void box_window();
		virtual int aux_window() { return BXE_NOAWIN; }
		virtual int extra_window() { return BXE_NOEWIN; }
		virtual int save_sob(SvArg *p) { p->st = 3; return 1; }
		virtual void wdat_cons(sthg * p) {}
		virtual void wdat_del(sthg * p) {}
		virtual int start_job_3(JobQ::ent_t * ent, char * arg) { return JQE_UNDEF; }
		virtual int mini(char * to) { memcpy(to, "B!U!G!a:zz%z%%", 14); return 14; }
		virtual int df_ui_ix() const { return 0; }
		virtual void mxc_notify(int ky, int f) { log("BUG: mxc_notify undefined"); }

		inline int sob_rw_arg() const { return (glob_flg&GLF_LIBMODE) + m_node->id(); }
		inline ABoxNode* node() const { return m_node; }
		inline int id() { return m_node->id(); }
		inline int w_oid(int k = 11) { return 16 * m_node->id() + k; }
		inline int get_ionm(char *to, int io, int j) { return m_node->get_ionm(to, io, j); }
		void set_node(ABoxNode* p) { m_node = p; }
		int set_title(char * to, int wid);
		int titlecmd_arg(char * to, int wid);
		BoxModel * model0() const { return m_model; }
		BoxModel * model() { if (!m_model) set_model(); return m_model; } // \/ TODO: fix model-ref
		inline void unset_model() { if (m_model) { /*BoxModel::unref(m_model);*/ m_model = 0; m_node->unset_model_rec(); }}
		inline bool unset_model_1() { return m_model && (m_model=0, 1); }
		BoxInst * mk_box() { if (!m_model) set_model();
			return m_model ? m_model -> mk_box() : 0; }
		void doc_window(int id4 = 13);
		inline int wnfl(int m = 2048) { return m_node->winflg(m); }
	protected:
		TAB7_DECL(iolbl, const char*);
		virtual void set_model() = 0; // only called when (m_model == 0)
		int set_boxp(BoxGen ** pp, BoxGen * to);
		BoxModel * m_model;
		ABoxNode * m_node;
};

class PrimBoxGen: public BoxGen {
	public:
		PrimBoxGen(BoxModel * m, int ni, int no, const char * cl);
		PrimBoxGen(qmb_arg_t qa, int k, int ni, int no, const char * cl);
		virtual void set_model() { bug("prim.box/set_model"); }
                virtual int n_in() const;
                virtual int n_out() const;
		virtual bool io_alias_perm() const;
                virtual const char * cl_name() { return m_cl; }
		virtual void box_window();
	protected:
		char m_ni, m_no, m_cl[6],
		     m_mdl_spc[24]; // TODO
};

extern BoxGen * box_bookmark[8];

#endif // __qwe_box_h__

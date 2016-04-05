#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/mman.h>

#define QWE_UTILC_DEF
#include "uc0.h"
#include "uc1.h"
#define QWE_DEFINE_ERRTAB
#include "errtab.inc"

///////////////// types/globals/decls ////////////////////////////////////////

#define TWF_UNDEAD 1
#define TWF_BOX 2
#define TWF_CLIP 4
#define TWF_OBJ 6
#define TWF_YSIZE 8
#define TWF_XSIZE 16
#define TWF_XYSIZE (TWF_YSIZE|TWF_XSIZE)
#define TWF_ABOVE 32

#define WF_RESIZE 1
#define WF_SURF 2
#define WF_CONT 4
#define WF_BIGDA1 8 
#define WF_BIGDA2 16
#define WF_BIGDA3 24 
#define WF_FULLSURF 32
//#define WF_RAWCLK 64
#define WF_CLIPL8R 128
#define WF_DTOR 256
#define WF_KEYEV 512
#define WF_XM1EV 1024
#define WF_EV_SHIFT 9

#define OI_ID 255
#define OI_NEW 1024

#define DF_NEXP 1
#define DF_NCOLL 2
#define DF_WCLOSE 4
#define DF_OIDEL 8
#define DF_WRAP 16
#define DF_BOXCONF 32
#define DF_TRK 64
#define DF_CWLU 128
#define DF_GRAPH 256
#define DF_WIDG 512
#define DF_MENU 1024
#define DF_REC  2048

#define RGB_C(X) (.0117647058823529 * (double)((X)-37))
#define RGB_C3(X,Y) cairo_set_source_rgb(X, RGB_C((Y)[0]), RGB_C((Y)[1]), RGB_C((Y)[2]));

#define EVMASK_DEF (GDK_BUTTON_PRESS_MASK | GDK_EXPOSURE_MASK)
#define EVMASK_KEY (GDK_KEY_PRESS_MASK|GDK_KEY_RELEASE_MASK)
#define EVMASK_XM1 (GDK_BUTTON_RELEASE_MASK|GDK_BUTTON1_MOTION_MASK)

struct _topwin;
struct _ww_t;

typedef void (tw_skel_fun_t) (struct _topwin * tw, char * arg); 
typedef void (tw_cmd_fun_t)  (struct _topwin * tw, char * arg);
typedef void (ww_skel_fun_t) (struct _ww_t * ww, const char **pp);
typedef int  (ww_get_fun_t)  (void * to, struct _ww_t * ww, int ty);
typedef void (ww_cmd_fun_t)  (struct _ww_t * ww, const char * arg);
typedef void (ww_clk_fun_t)  (struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev);
typedef GtkWidget * (vbox_line_fun_t) (struct _ww_t * ww, int ix);
#define FTYD(n) typedef n##_fun_t * n##_fun
FTYD(tw_skel); FTYD(tw_cmd); FTYD(ww_skel); FTYD(ww_get); FTYD(ww_cmd); FTYD(vbox_line); FTYD(ww_clk);

typedef struct _tw_cl {
	int ch;
	int flg;
	tw_skel_fun skel;
	tw_cmd_fun cmd;
} tw_cl;

typedef struct _ww_cl {
	int ch;
	ww_skel_fun skel;
	ww_get_fun get;
	ww_cmd_fun cmd;
	ww_clk_fun clk;
	int flg;
} ww_cl;

typedef struct _obidtab {
	int otab[256];
	unsigned char ftab[256];
	int n, n0, f0;
} obidtab;

typedef struct _topwin {
	int id;
	int state; // 0:none 1:skel 2:closecmd
	char cmdpref[20]; int cmdp_len;
	char title[64];
	tw_cl * cl;
	GtkWidget * w;
	int ns_base, ns_rec1;
	struct _ww_t *pp_sub[32];
	chtab sub_ch;
	sthg arg[128];
	void * etc;
	int anon_w;
	int ix4;
	int vb_c_ix;
} topwin;

typedef struct _ww_t {
	ww_cl * cl;
	GtkWidget * w;
	struct _topwin * top;
	int ix;
	char cmd[16];
	sthg arg[5];
	void * etc;
} ww_t;

#define DA_SURF(x) ((x)->arg[0].p)

typedef struct _ob_dir {
	int id, flg;
	topwin * clip;
	GtkTreeIter nd[2];
} ob_dir;

typedef struct _ob_box {
	int id, flg;
	topwin * tw[8];
} ob_box;


static tw_skel_fun_t mwin_skel, t2win_skel, clip_skel, wrap_skel, tgrid_skel, graph_skel, acfg_skel, mcfg_skel,
		     calc_skel, pz_skel, gconf_skel, doc_skel, ttrk_skel, err_skel, a20_skel, in01_skel,
		     itb_skel, tkcf_skel;
static tw_cmd_fun_t wrap_cmd, tgrid_cmd, gconf_cmd, doc_cmd, ttrk_cmd, err_cmd, tkcf_cmd;
static ww_skel_fun_t pv_skel, button_skel, entry_skel, scale_skel, daclip_skel, dasep_skel, dakcf_skel,
                     dacnt_skel, dacntvs_skel, dalbl_skel, daclb_skel, dawr1_skel, daprg_skel,
                     dagrid_skel, dagraph_skel, dapz_skel, vbox_skel, dabmp_skel, datrk_skel;
static ww_cmd_fun_t pv_cmd, entry_cmd, daclip_cmd, dalbl_cmd, daclb_cmd, dawr1_cmd, dacnt_cmd, dasep_cmd,
                    daprg_cmd, dagrid_cmd, dagraph_cmd, dapz_cmd, vbox_cmd, dabmp_cmd, datrk_cmd, dakcf_cmd;
static ww_get_fun_t pv_get, entry_get, dagraph_get, daclip_get;
static ww_clk_fun_t debug_clk, daclip_clk, dlmenu_clk, dlyn_clk, dlbtn_clk, dacnt_clk, dacntvs_clk, dakcf_clk,
		    daclb_clk, dawr1_clk, daprg_clk, dagrid_clk, dagraph_clk, dabmp_clk, datrk_clk;
static vbox_line_fun_t wrap_vbl_i, wrap_vbl_t, calc_vbl, gconf_vbl, doc_vbl, err_vbl, clip_vbl;

static tw_cl tw_cltab[] = { {'?',0,NULL,NULL}, 
	{'.', TWF_UNDEAD|TWF_ABOVE, mwin_skel, NULL },
	{'/', TWF_UNDEAD, t2win_skel, NULL },
	{'K', TWF_XYSIZE, clip_skel, NULL },
	{'w', TWF_YSIZE , wrap_skel, wrap_cmd },
	{'#', 0         , tgrid_skel, tgrid_cmd },
	{'g', 0         , graph_skel, NULL },
	{'P', 0         , pz_skel, NULL },
	{'c', TWF_YSIZE , calc_skel, NULL },
	{'i', 0         , itb_skel, NULL },
	{'t', 0         , ttrk_skel, ttrk_cmd },
	{'C', 0         , gconf_skel, gconf_cmd },
	{'D', 0         , doc_skel, doc_cmd },
	{'E', TWF_YSIZE , err_skel, err_cmd },
	{'A', 0         , a20_skel, NULL },
	{'S', 0         , acfg_skel, NULL },
	{'F', 0         , mcfg_skel, NULL },
	{'J', 0         , in01_skel, NULL },
	{'k', 0         , tkcf_skel, tkcf_cmd },
	{ 0 , 0, 0, NULL } };

static ww_cl ww_cltab[] = { {'?', pv_skel, pv_get, pv_cmd, debug_clk, 0 },
	{'b', button_skel, NULL, NULL, NULL, 0 },
	{'e', entry_skel, entry_get, entry_cmd, NULL, 0 },
	{'s', scale_skel, NULL, NULL, NULL, 0 },
	{':', vbox_skel, NULL, vbox_cmd, NULL, 0 },
	{'K', daclip_skel, daclip_get, daclip_cmd, daclip_clk, WF_BIGDA3 | WF_KEYEV | WF_RESIZE },
	{'M', dalbl_skel, NULL, dalbl_cmd, dlmenu_clk, WF_RESIZE },
	{'8', dacnt_skel, NULL, dacnt_cmd, dacnt_clk, WF_RESIZE },
	{'!', dacntvs_skel, NULL, dacnt_cmd, dacntvs_clk, WF_RESIZE | WF_CONT },
	{'Y', dalbl_skel, NULL, dalbl_cmd, dlyn_clk, WF_RESIZE },
	{'L', dalbl_skel, NULL, dalbl_cmd, debug_clk, WF_RESIZE },
	{'E', dalbl_skel, NULL, dalbl_cmd, debug_clk, WF_RESIZE },
	{'B', dalbl_skel, NULL, dalbl_cmd, dlbtn_clk, WF_RESIZE },
	{'C', daclb_skel, NULL, daclb_cmd, daclb_clk, WF_RESIZE },
	{'P', dapz_skel,  NULL, dapz_cmd, debug_clk, WF_BIGDA3 },
	{'1', dawr1_skel, NULL, dawr1_cmd, dawr1_clk, WF_RESIZE },
	{'2', dabmp_skel, NULL, dabmp_cmd, dabmp_clk, WF_RESIZE },
	{'_', dasep_skel, NULL, dasep_cmd, debug_clk, WF_RESIZE },
	{'%', daprg_skel, NULL, daprg_cmd, daprg_clk, WF_RESIZE },
	{'*', dakcf_skel, NULL, dakcf_cmd, dakcf_clk, WF_RESIZE|WF_BIGDA1|WF_KEYEV|WF_XM1EV},
	{'#', dagrid_skel, NULL, dagrid_cmd, dagrid_clk, WF_RESIZE|WF_BIGDA1|WF_DTOR|WF_KEYEV|WF_XM1EV},
	{'+', dagrid_skel, NULL, dagrid_cmd, dagrid_clk, WF_RESIZE|WF_BIGDA1|WF_DTOR|WF_XM1EV},
	{'t', datrk_skel,  NULL, datrk_cmd, datrk_clk, WF_RESIZE|WF_BIGDA2|WF_DTOR|WF_KEYEV|WF_XM1EV},
	{'g', dagraph_skel, dagraph_get, dagraph_cmd, dagraph_clk, WF_SURF | WF_CONT | WF_FULLSURF },
	{ 0, NULL, NULL, NULL } };

static chtab tw_clch, ww_clch;

static obidtab  oi_box, oi_dir, oi_etc;
static ob_dir ot_dir[256];
static ob_box ot_box[256];
static topwin ot_etc[256];

static int expand_cmd_flg = 1;

static int conf_lbh = 18, conf_lbfh = 15, conf_lbfs = 15, conf_lbfh_s = 12, conf_lbfs_s = 12;
static int conf_nomwin = 0;
static int conf_portwid;

static int tlog_c_onq = 0, tlog_c_bk = 0;
static int dflg = 0;
static const char * dflg_s = "1:node_expand 2:node_collapse 4:closewin 8:oi_del 10:wrap 20:boxconf\n"
			     "40:track 80:lookuperr-0x[56]7 100:graph 200:widg 400:menu 800:rec";
// general
#define MYPRINTF(NM, L)             			\
void NM(const char * fmt, ...) {      			 \
        va_list ap; va_start(ap, fmt); 			  \
        int n = vsnprintf(pf_buf+1, 1022, fmt, ap);        \
        va_end(ap); if (n<=0) return; if (n>1022) n = 1022; \
	pf_buf[n+1]='\n'; write(1+L, pf_buf+L, n+2-L); }

static char pf_buf[1024];
MYPRINTF(CMD,0)
MYPRINTF(LOG,1)

char * bigblk(int n) {
	        return (char*)mmap(0, n, PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE, -1, 0); }

char * alloc_64k() {
        static char *cur=0, *lim = 0;
        static int bs = 262144;
        if (cur==lim) {
                cur = bigblk(bs); lim = cur+bs;
                if (bs<(1<<24)) bs += bs;
        }
        char * p = cur; cur += 65536; return p;
}

static char * freeptr_3k = NULL;
char * mk_3k_blk() { char * p = alloc_64k(); int i; 
		     for (i=0; i<61440; i+=3072) *(char**)(p+i) = p+i+3072; *(char**)(p+i)=NULL; return p;}
char * alloc3k() { char * p = freeptr_3k ? freeptr_3k : (freeptr_3k = mk_3k_blk());
		   freeptr_3k = *(char**)p; return p; }
void free3k(void *p) { *(char**)p = freeptr_3k; freeptr_3k = (char*)p; }
	
static int get_tok(char * to, int n, const char ** pp, int delim) {
	int cls = (delim&128) ? 0 : '}'; delim &= 127;
	int i,z=1; if (n<=0) z=0, n=-n;
	for (i=0; i<n-1 && **pp && **pp!=delim && **pp!=cls; i++, (*pp)++) to[i] = **pp;
	if (delim && **pp==delim) ++(*pp);
	if (z) to[i] = 0; return i;
}

static int get_dec(const char ** pp, int delim) {
	int x,a=0; while (x=**pp-'0',0<=x&&x<=9) a=10*a+x, ++*pp;
	if (delim && delim==**pp) ++*pp;
	return a;
}

static void write_tlog();
static void bye(int op) {
	LOG("bye(0x%x)", op );
	if (tlog_c_onq) write_tlog();
	write(2, "bye\n&CMD0\n", 10);
	if (op<32) exit(op);
	CMD("q%c", 32+(op&31));
	while(1) sleep(60); // patiently wait to be killed
}

static gboolean fifo_bye (GIOChannel *src, GIOCondition condition, gpointer data) {
	int k = (char*)data - pf_buf;
	return LOG("input pipe %d: %s", k>>1, "HUP\0ERR"+4*(k&1)), bye(1), TRUE; }

static void hello2(GtkWidget *w, gpointer p) { LOG ("Hello bazmeg"); }

static gboolean not_so_easy(GtkWidget *widget, GdkEvent *event, gpointer data) {
    LOG ("delete event ignored for main window. hahahaha!");
    return TRUE;
}

static void mw_bye(GtkWidget *w, gpointer p) { bye(1); }
#define OOPS gtk_button_new_with_label("oops");

static void pv_skel(struct _ww_t * ww, const char **pp) { LOG("ERROR: pv_skel called!"); }
static void pv_cmd(struct _ww_t * ww, const char * arg) { LOG("ERROR: pv_cmd called!"); }
static int pv_get (void * to, struct _ww_t * ww, int ty) { LOG("ERROR: pv_get called!"); return 0; }

int smallnan(char * to, long long xl) {
	union { int i[16]; char c[64]; } uu;
	int dsc = nan_unpk(uu.c, uu.i, xl, 0); if (dsc<0) return dsc;
	int j, n = dsc & 255, bits = (dsc>>8)&255, sg = dsc >> 16;
	if (n==254) return memcpy(to, "-nan"+1-sg, 3+sg), 3+sg;
	if (n==253) return memcpy(to, "-nan?"+1-sg, 4+sg), 4+sg;
	if (bits==1) return memcpy(to, "{...}", 5), 5;
	if (bits<8 && n<7) return to[0]=39, memcpy(to+1, uu.c, n), n+1;
	*to = '['; j = 1; if (n>9) n-=10, to[j++]='1';
	to[j] = 48+n; to[j+1] = 'x'; to[j+2] = ']'; return j+3;
}

int smallnum(char * to, double v) {
        const char * pfx="afpnum kMGTPE";
        const char * fmt1[2] = {"%.5g", "%.4g"};
        const char * fmt2[3] = {"%.3e", "%.2e", "%.1e"};
	long long ll; memcpy(&ll, &v, 8); int r = smallnan(to, ll); if (r>=0) return r;
        if (v!=v) { memcpy(to, "NaN", 3); return 3; }
        int i, j, k, sg = 0; if (v<0.0) v=-v, to[sg++]='-';
        if (v<1e-301) { to[sg] = 48; return sg+1; }
        char buf[16];
        if (v>9e20 || v<1.1e-18) {
                k = sprintf(buf, fmt2[sg+1], v);
                for (i=0, j=sg; buf[i]; i++) {
                        if ((buf[i]|33)=='e') continue;
                        else if (buf[i]=='+') to[j++] = '!';
                        else if (buf[i]=='-') to[j++] = ';';
                        else to[j++] = buf[i];
                }
                return j; 
        }
        if (v>=(sg ? .0995 : .09995) &&  v < (sg ? 9999.5 : 99999.5)) return sg+sprintf(to+sg, fmt1[sg], v);
        k = sprintf(buf, fmt2[sg], v);
        for (i=0; i<k; i++) if ((buf[i]|33)=='e') break;
        if (i==k) { memcpy(to+sg, buf, i); buf[i+sg] = '?'; return i+sg+1; }
        int xp = atoi(buf+i+1) + 18;
        buf[i]='$'; buf[i+1]=0;
        if (dflg & DF_GRAPH) LOG("xp=%d, buf=\"%s\"", xp, buf);
        switch(xp % 3) {
                case 0: break;
                case 1:
                        if (buf[1]=='.') buf[1] = buf[2], buf[2] = '.';
                        else buf[1]=48, i++;
                        break;
                case 2:
                        if (buf[1]!='.') buf[1] = buf[2] = 48, i=3;
                        else if (k<4) buf[1] = buf[2], buf[2] = 48, i=3;
                        else buf[1]=buf[2], buf[2]=buf[3], buf[3]='.';
                        break;
        }
        memcpy(to+sg, buf, i); to[sg+i] = pfx[xp/3];
        return sg+i+1;
}

static char * hxdoub_str(char * to, const char * s, int prec) {
	double v = hx2doub(s);
	long long xl; memcpy(&xl, &v, 8);
	static char buf[256]; if (!to) to = buf;
	int l1 = nan_2str(to, xl); if (l1>=0) return to[l1]=0, /*LOG("nan_str: '%s'", to),*/ to;
	char fmt[8], *p=fmt; *(p++)='%'; *(p++)='.';
	if (prec<2) prec=2; else if (prec>15) prec=15;
	if (prec>9) prec-=10, *(p++)=49;
	*(p++) = prec+48; *(p++) = 'g'; *p = 0;
	sprintf(to, fmt, v); return to;
}

static char * hxdoub_lbl(const char * s) { // q.mark remove, min.len=3
	char * r = hxdoub_str(NULL, s, 15); if (*r=='"') ++r;
	int l = strlen(r); if (r[l-1]=='"') r[--l] = 0;
	if (l<3) { do r[l++] = ' '; while(l<3);  r[3] = 0; }
	return r;
}

///////////////// font metrics /////////////////////////////////////////////////

static char font_h[20], font_o[20], font_dig[20], font_b32[20], font_w[1920];
#define FONT_HEIG(x) (font_h[(x)-6])
#define FONT_OFFS(x) (font_o[(x)-6])
#define FONT_NWID(x) (font_dig[(x)-6])
#define FONT_VWID(x) (font_b32[(x)-6])
#define FONT_WTAB(x) (font_w + 96 * ((x)-6) - 32)

static int tx_len(int fs, const char * s) {
	char * p = font_w + 96 * (fs - 6) - 32;
	int r; for (r=0; *s; s++) if (*s>31 && *s<127) r+=p[(int)*s];
	return r;
}

static int get_fontsize(int h, int w, const char * s) {
	int i = h - 6;
	if (i<0) return 6; else if (i>19) i=19;
	while (i && font_h[i]>h) --i;
	while (i && s && tx_len(i+6, s)>w) --i;
	return i+6;
}

static void txtm_init() {
	cairo_surface_t * surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, 1024, 64);
	int cs = cairo_surface_status(surf);
	if (cs) { LOG("test surf: %s", cairo_status_to_string(cs)); bye(1); }
	cairo_t * cr = cairo_create(surf);
	cairo_text_extents_t ext;
	char buf[2]; buf[1] = 0; 
	int i,j,k; double y, ym, yp; //char c1,c2; TODO: env
	for (i=0; i<20; i++) {
		ym=yp=0.0;
		double xt = 0.0;
		for (j=0; j<95; j++) {
			buf[0] = j+32; 
			cairo_set_font_size(cr, (double)(i+6));
			cairo_move_to(cr, 1.0, 32.0);
			cairo_text_extents(cr, buf, &ext);
			if ((y=ext.y_bearing)<ym) ym = y; //c1=j+32;
			if ((y+=ext.height)>yp) yp = y; //c2=j+32;
			font_w[96*i+j] = (int)lround(ext.x_advance);
			xt += ext.x_advance;
		}
		// LOG("ext(%d): h=%g o=%g xt=%g '%c%c'", i+6, yp-ym, -ym, xt, c1, c2); TODO: env
		font_h[i] = (int)lround(yp-ym);
		font_o[i] = -(int)lround(ym);
		font_dig[i] = font_w[96*i+13];
		for (j=16; j<26; j++) if (font_dig[i]<font_w[k=96*i+j]) font_dig[i]=font_w[k];
		font_b32[i] = font_dig[i];
		for (j=65; j<87; j++) if (font_b32[i]<font_w[k=96*i+j]) font_b32[i]=font_w[k];
	}
}

///////////////// obidtab //////////////////////////////////////////////////////

static void obidtab_ini(obidtab *p, int n) {
	p->n = p->n0 = p->f0 = n; 
	p->ftab[255] = 255;
	int i,j; for (i=n; i<255; i++) p->ftab[i] = i+1;
	for (i=j=0; i<n; i++, j+=257) p->otab[i] = j;
}

static void obidtab_debug(obidtab *p) {
	int i, j, k, n = p->n;
	LOG("obidtab: n=%d, n0=%d", n, p->n0);
	fprintf(stderr,"otab:"); for (i=0; i<n; i++) { fprintf(stderr, " %x", p->otab[i]); }; 
	fprintf(stderr, "\n!free:");
	for (i=p->f0,k=0; i>=0; j=p->ftab[i], i=(i==j)?-1:j,k++) {
		if (k<10) fprintf(stderr, " %x", i); 
		else if (k==10) fprintf(stderr, " ...");
	}
	fprintf(stderr, " -- %d/%d\n", k, 256 - n); fflush(stderr);
}

static int obidtab_alloc(obidtab *p) {
	int r = p->f0; if (r<0) return r;
	int k = p->ftab[r]; p->f0 = (k==r) ? -2 : k;
	return r;
}

static void obidtab_free(obidtab *p, int i) {
	i &= 255;
	if (p->f0 < 0) p->f0 = (unsigned char)i, p->ftab[i] = p->f0;
	else p->ftab[i] = (unsigned char)(p->f0), p->f0 = i;
}

static int obidtab_cut(obidtab *p, int ix) {
	if (ix<0 || ix>=p->n) {
		LOG("oi_cut: invalid idx %x", ix); return 0; }
	if (dflg & DF_OIDEL) LOG("oi_cut: %d 0x%x", ix, p->otab[ix]);
	int k = (--p->n - ix);
	obidtab_free(p, p->otab[ix] & 255);
	if (k) memmove(p->otab+ix, p->otab+ix+1, k*sizeof(int));
	return 1;
}

static int obidtab_lookup(obidtab *p, int id, int flg) {   // 1:force 2:del
	int id8 = (id & 0xfffff0)<<4, id1 = id >> 4;
	if (id1 < p->n0) return id1;
	int n = p->n, *ot = p->otab;
	int lo=0, hi=n-1;
	while (lo<=hi) {
		int md = (lo+hi)>>1;
		int d = id8 - (0xfffff00 & ot[md]);
		if (!d) return (flg&2) ? obidtab_cut(p,md) : (ot[md] & 255);
		d<0 ? (hi = md-1) : (lo = md+1);
	}
	if (!(flg&1)) return -1;
	int oi = obidtab_alloc(p); if (oi<0) return oi;
	if (lo<n) memmove(ot+lo+1, ot+lo, sizeof(int)*(n-lo));
	ot[lo] = id8 | oi;
	++p->n; 
	return oi + OI_NEW;
}

static void ob_dir_ini(ob_dir * p, int id) { p->id = id; p->flg = 0; p->clip = NULL; }
static void ob_box_ini(ob_box * p, int id) { p->id = id; p->flg = 0; }

static int dirtab_del(ob_dir * p, int flg) {
	return p && p->id && !((p->flg) &= ~flg) && obidtab_lookup(&oi_dir, p->id, 2); }

static int boxtab_del(ob_dir * p, int flg) {
	return p && p->id && !((p->flg) &= ~flg) && obidtab_lookup(&oi_dir, p->id, 2); }

static void tw_remove(topwin *tw) {
	int id = tw->id, ty=id&15;
	if (ty&8) {
		int k = obidtab_lookup(&oi_box, id, 0);
		if (k<0) { LOG("tw_remove: %x: boxtab lookup failed", id); return; }
		ob_box * ob = ot_box + (k & OI_ID);
		if (!(ob->flg &= ~(1<<(ty-8)))) obidtab_lookup(&oi_box, id, 2);
	} else if (ty==7) {
		int k = obidtab_lookup(&oi_etc, id, 0);
		if (k<0) { LOG("tw_remove: %x: globwin lookup failed", id); return; }
		obidtab_lookup(&oi_etc, id, 2);
	} else if (ty==3) {
		int k = obidtab_lookup(&oi_dir, id, 0);
		if (k<0) { LOG("tw_remove: %x: dirtab lookup failed", id); return; }
		ob_dir * od = ot_dir + (k & OI_ID);
		if (!(od->flg &= ~4)) obidtab_lookup(&oi_dir, id, 2);
	}
}

///////////////// top window / widget (gen) //////////////////////////////////

GtkWidget * ww_widg(ww_t * ww) {
	if (ww->cl->flg & WF_CONT) return GTK_WIDGET(ww->arg[4].p);
	else return GTK_WIDGET(ww->w);
}

static void ww_free(ww_t * ww) {
	if (ww->etc) (ww->cl->flg&WF_DTOR) ? (*ww->cl->cmd)(ww, "~") : free(ww->etc);
	if ((ww->cl->flg & WF_SURF) && DA_SURF(ww)) cairo_surface_destroy(DA_SURF(ww));
	memset(ww, 0, sizeof(ww_t));
}

static void tw_free(topwin *tw) {
	int i,j;
	for (i=0; i<32; i++) {
		ww_t * p = tw->pp_sub[i]; if (!p) continue;
		for (j=0; j<32; j++) if (p[j].cl) ww_free(p+j);
		free(p);
	}
	if ((tw->id&15) != 7) free(tw); else memset(tw, 0, sizeof(topwin));
}

static void cltab_init() {
	int i,j; 
	for (i=0; (j=tw_cltab[i].ch); i++) 
		if (chtab_force(&tw_clch,j) != i) LOG("cltab_ini: ixt error");
	for (i=0; (j=ww_cltab[i].ch); i++) 
		if (chtab_force(&ww_clch,j) != i) LOG("cltab_ini: ixw error");
}

static void tw_bye(GtkWidget *w, gpointer p) {
	topwin * tw = (topwin*)p;
	if (dflg & DF_WCLOSE) LOG("tw_bye: 0x%x, cl=0x%x '%c'", tw->id, tw->cl->ch, tw->cl->ch);
	if (tw->state!=2) {
		int j = tw->id, j4 = j&15, j20 = j>>4, cl = tw->cl->ch;
		if (dflg&DF_REC) (j4==7) ? CMD("QRZ%x7$%c", j20, cl) : CMD("QRZ`%x`%c$%c", j20, hexc1(j4), cl);
		if (cl=='J') CMD("L%04xFF01", j20&65535); else CMD("x%c%x", hexc1(j4), j20);
	}
	tw_remove(tw); tw_free(tw);
}

static gboolean box_conf_2(GtkWidget *w, int c, gpointer p) {
	if (dflg & DF_BOXCONF) { char ch = c; write(2, &ch, 1); }
	topwin * tw = (topwin*) p;
	if (!tw->arg[0].p) return TRUE;
	GtkRequisition rq; gtk_widget_size_request (GTK_WIDGET(tw->arg[0].p),&rq);
	tw->arg[2].s[0] = rq.width; tw->arg[2].s[1] = rq.height;
	if (dflg & DF_BOXCONF) LOG("box config:wi=%d, he=%d", rq.width, rq.height);
	int cw, ch, flg = tw->cl->flg, xf = flg&TWF_XSIZE, yf = flg&TWF_YSIZE;
	gtk_window_get_size(GTK_WINDOW(tw->w), &cw, &ch);
	if ( ( ((cw-rq.width ) & (xf ? -1 : INT_MIN)) && (cw = rq.width,  1)) |
	     ( ((ch-rq.height) & (yf ? -1 : INT_MIN)) && (ch = rq.height, 1)) )
		gtk_window_resize(GTK_WINDOW(tw->w), cw, ch);
	return TRUE;
}

static gboolean box_conf(GtkWidget *w, GdkEventConfigure * ev, gpointer p) {
	return box_conf_2(w, '#', p); }
static gboolean box_xps(GtkWidget *w, GdkEventExpose * ev, gpointer p) {
	return box_conf_2(w, '&', p); }

static void topwin_skel(topwin * p, int cl_ch, char * arg) {
	int nf = (!p->state); if (!nf) goto st1;
	p->cl = tw_cltab + chtab_get(&tw_clch, cl_ch);
	if (!p->cl) { LOG("topwin_skel: unknown tw class 0x%x '%c'", cl_ch, cl_ch); return; }
	if (!p->cl->skel) { LOG("topwin_skel: undef skel fun: 0x%x '%c'", cl_ch, cl_ch); return; }
	p->w = gtk_window_new (GTK_WINDOW_TOPLEVEL);
	p->etc = NULL;
	if (p->cl->flg & TWF_UNDEAD) {
		g_signal_connect (p->w, "delete-event", G_CALLBACK(not_so_easy), NULL);
		g_signal_connect (p->w, "destroy", G_CALLBACK(mw_bye), NULL);
	} else {
		g_signal_connect (p->w, "destroy", G_CALLBACK(tw_bye), (gpointer)p);
	}
	gtk_container_set_border_width (GTK_CONTAINER (p->w), 0);
	chtab_ini(&p->sub_ch,111);
	memcpy(p->title, "...", 4);
st1:	(*p->cl->skel)(p,arg); p->state = 1;
	gtk_window_set_title(GTK_WINDOW(p->w), p->title);
	if (nf) {
		GtkWidget * w = (GtkWidget*)(p->arg[0].p);
		gtk_container_add (GTK_CONTAINER (p->w), w);
		gtk_widget_show(w);
		gtk_widget_show(p->w);	
		if (p->cl->flg & TWF_ABOVE) 
			gdk_window_set_keep_above(gtk_widget_get_window(GTK_WIDGET(p->w)), TRUE);
		if (p->cl->flg & TWF_XYSIZE) 
			g_signal_connect(w, "configure-event", G_CALLBACK(box_conf), (gpointer)p),
			g_signal_connect(w, "expose-event", G_CALLBACK(box_xps), (gpointer)p);
	}
}

static topwin * topwin_mk(int id, int cl_ch, char * arg) {
	int ty = id&15; topwin * tw = NULL;
	if (ty&8) {
		int ty2 = ty-8, k = obidtab_lookup(&oi_box, id, 1);
		if (k<0) { LOG("topwin_mk: lookup(f) failed"); return NULL; }
		ob_box * ob = ot_box + (k & OI_ID);
		if (k&OI_NEW) { ob_box_ini(ob, id&0xfffff0); }
		else if (ob->flg & (1<<ty2)) { topwin_skel(ob->tw[ty2], cl_ch, arg);
					       return ob->tw[ty2]; }
		tw = ob->tw[ty2] = calloc(sizeof(topwin),1);
		ob->flg |= 1<<ty2; goto node;
	} else if (ty==3) {
		if (cl_ch!='K') { LOG("topwin_mk: ty==3, cl!='K'(0x%x)",cl_ch); return NULL; }
		int k = obidtab_lookup(&oi_dir, id, 1);
		if (k<0) { LOG("topwin_mk: lookup(f) failed"); return NULL; }
		ob_dir * od = ot_dir + (k & OI_ID);
		if (k&OI_NEW) { ob_dir_ini(od, id&0xfffff0); }
		else if (od->flg&4) { topwin_skel(od->clip, 'K', arg); return od->clip; }
		tw = od->clip = calloc(sizeof(topwin),1);
		od->flg |= 4; goto node;
	} else if (ty==7) {
		int j = obidtab_lookup(&oi_etc, id, 1);
		//LOG("topwin 0x%x found: 0x%x %d", id, j, j&OI_ID);
		tw = ot_etc + (j&OI_ID);
		if (j&OI_NEW) tw->state = 0;
		tw->id = id; tw->cmdp_len = 0;
		topwin_skel(tw, cl_ch, arg);
		gtk_widget_show(tw->w);
		return tw;
	} else {
		LOG("topwin_mk: invalid id 0x%x", id); return NULL;
	}
node:	tw->id = id;
	topwin_skel(tw, cl_ch, arg);
	tw->cmdpref[0]='#';
	int j = hx5(tw->cmdpref+1, (id>>4)) + 1;
	tw->cmdpref[j++] = 36; tw->cmdp_len = j;
	gtk_widget_show(tw->w);
	return tw;
}

static topwin * tw_lookup(int id) {
	int ty = id & 15;
	if (ty&8) {
		int j = obidtab_lookup(&oi_box, id, 0);
		if (j<0) return NULL;
		return (ot_box[j].flg & (1<<(ty-8))) ? ot_box[j].tw[ty-8] : NULL;
	} else if (ty==3) {
		int j = obidtab_lookup(&oi_dir, id, 0);
		if (j<0) return NULL;
		return (ot_dir[j].flg&4) ? ot_dir[j].clip : NULL;
	} else if (ty==7) {
		int j = obidtab_lookup(&oi_etc, id, 0);
		return j<0 ? NULL : ot_etc + j;
	} else {
		LOG("tw_lookup: invalid id 0x%x", id); return NULL;
	}
}

static ww_t * widg_p(topwin * tw, int ix) {
	if (ix&~1023) return NULL;
	int ix0 = ix >> 5, ix1 = ix & 31;
	ww_t ** pp = tw->pp_sub + ix0;
	if (!*pp) *pp = calloc(sizeof(ww_t), 32);
	(*pp)[ix1].ix = ix;
	return *pp + ix1;
}

static int widg_lookup(topwin * tw, const char ** pp, int force) {
	if(!**pp) return -1;
	if(**pp=='_') {
		++*pp; if (force) return (--tw->anon_w)&1023;
		int k = hex2(*pp); *pp+=2; return 1023 - k;
	}
	if (force&2) {
		int i = chtab_get(&tw->sub_ch, **pp);
		if (i<99) { LOG("widg 0x%x already ex. (%c)", **pp, **pp); return -1; }
	}
	int ix = force ? chtab_force(&tw->sub_ch, *((*pp)++))
	               : chtab_get(&tw->sub_ch, *((*pp)++));
	return (ix>99) ? -1 : ix;
}

static ww_t * widg_lookup_p(topwin * tw, const char ** pp) { return widg_p(tw, widg_lookup(tw, pp, 0)); }
static ww_t * widg_lookup_ps(topwin * tw, const char * p) { return widg_p(tw, widg_lookup(tw, &p, 0)); }
static int widg_lookup_ci(topwin * tw, int c, int i) {
	char s[4]; const char *p = s; 
	s[0] = c; s[1] = hexc1(i>>4); s[2] = hexc1(i&15), s[3]=0; 
	return widg_lookup(tw, &p, 0); }
static ww_t * widg_lookup_pci(topwin * tw, int c, int i) { return widg_p(tw, widg_lookup_ci(tw, c, i)); }

static ww_t * widg_new(topwin * tw, const char ** pp) {
	if(!**pp) return NULL; 
	ww_cl * cl = ww_cltab + chtab_get(&ww_clch, *((*pp)++));
	if (cl->ch=='?') return NULL;
	int ix = tw->vb_c_ix ? tw->vb_c_ix + (*(*pp)++) - 48
			     : widg_lookup(tw, pp, 3);
	ww_t * p = widg_p(tw, ix);
	if (!p) { LOG("widget lookup(%d,f) failed", ix); return p; }
	p->cl = cl; p->top = tw; p->ix = ix; p->w = NULL; p->etc = NULL;
	(*cl->skel) (p, pp);
	if (!p->w) p->w = OOPS;
	gtk_widget_show(p->w);
	return p;
}

static void tw_close(topwin * tw) { if(tw->state)tw->state=2; gtk_widget_destroy(tw->w); }
static inline int tw_defcmd(topwin * tw, char* to) {
	int l = tw->cmdp_len; memcpy(to, tw->cmdpref, l); return l; }

static void cmd1(char * str);
static void fw_cmd(topwin * tw, const char *s1, const char *s2) {
	int l1 = s1 ? strlen(s1) : 0,
	    l2 = s2 ? strlen(s2) : 0;
	char buf[l1+l2+1]; if (l1) memcpy(buf   , s1, l1);
			   if (l2) memcpy(buf+l1, s2, l2); buf[l1+l2] = 0;
	ww_t * ww;
	if (!tw) cmd1(buf);
	else if (tw->cl->cmd) (*tw->cl->cmd)(tw, buf);
	else s1=buf, ww=widg_lookup_p(tw, &s1), (ww->cl->cmd)(ww, s1);
}

static int cmd_esc(ww_t * ww, const char *s1, const char *s2) { switch(*s1) {
	case '>': return fw_cmd(ww->top, s1+1, s2), 1;
	case '!': return fw_cmd(NULL,    s1+1, s2), 1;
	case '/': return gtk_window_present(GTK_WINDOW(ot_etc[1].w)), 1;
	case '?': if (s1[1]) return CMD("W#2.%s", s1+1), 1;
		  switch(ww->top->cl->ch) {
			  case '.': s2 = "main";   break;
			  case '/': s2 = "tree";   break;
			  case 'F': s2 = "config"; break;
			  case 'E': s2 = "errors"; break;
			  case 'S': s2 = "audio";  break;
			  case 'A': s2 = "auconv"; break;
			  default: return LOG("builtin help not found"), 0; }
		  return CMD("W#2.win.%s", s2), 1;
	default: return 0; 
}}

static int widg_defcmd(ww_t * ww, const char * arg) {
	char cmd[4096];
	const char * s = ww->cmd;
	topwin * tw = ww->top;

	if (*s==36) return cmd_esc(ww, s+1, arg) ? 0 : (LOG("invalid esc-cmd 0x%x(%c)", s[1], s[1]), 0);
	cmd[0] =  '~';
	int mode=0;
	if (*s) { cmd[1] = (*s=='~') ? (mode=1, s+=2, s[-1]) : *(s++); }
	else {  if (*arg==36) { switch(*++arg) {
			case '.': mode = -1; ++arg; break;
			case '~': mode = 1; ++arg; break;
			default: return cmd_esc(ww, arg, NULL) ? 0 : 
				 (LOG("invalid esc-cmd 0x%x(%c)", s[1], s[1]), 0);
		}}
		cmd[1] = *(arg++); if (!*arg && !tw->cmdp_len) return cmd[2]='\n', write(1, cmd, 3);
	}
	int k, i = (mode<1) ? 2 + tw_defcmd(tw, cmd+2) : 2;
	if (mode<0) cmd[i-1] = '.';
	mode = 0;
	while (1) {
		if (i>3840) break;
		if (!*s) { if (mode || ww->cl->ch != 'M') break; else mode=1, s=arg; }
		if (*s=='%') {
			if (!*++s) break;
			ww_t * ww2; switch(*s) {
				case '@': i += tw_defcmd(tw, cmd+i) - 1; ++s; continue;
				case '.': memcpy(cmd+i, ww->arg[4].c, 6); i+=6; s++; continue;
				case '$': kill(getppid(), hxd2i(s[1])); return 0;
				case '-': case '+':
					  ww2 = widg_p(tw, ww->ix + (44-*s)*hxd2i(s[1]));
					  if (ww2) i += (*ww2->cl->get)(cmd+i, ww2, 's'); else cmd[i++]='?';
					  s += 2; continue;
				default: break;
			}
			if ((k=chtab_get(&tw->sub_ch, *(s++)))<99) {
				ww_t * ww2 = widg_p(tw, k);
				i += (ww2->cl->get)(cmd+i, ww2, 's'); 
				continue;
			}}
		cmd[i++] = *(s++);
	}
	if (!mode) { while (*arg) cmd[i++] = *(arg++); }
	cmd[i] = '\n'; return write(1, cmd, i+1);
}

GtkWidget * t2_widg(int t);

GtkWidget * parse_w(topwin * tw, const char **s) {
	GtkWidget * w; ww_t * ww; int cls;
	char lb[64];
	int xpf = 0;
	if ((**s&120)==48) xpf = *((*s)++) & 7;
	switch(**s) {
		case '(': cls=')'; w = gtk_hbox_new(0,0); break;
		case '[': cls=']'; w = gtk_vbox_new(0,0); break;
		case 'l': ++*s; get_tok(lb, 64, s, '$'); return gtk_label_new(lb);
		case 'h': ++*s; return gtk_button_new_with_label("he");
		case 'g': ++*s; return gtk_button_new_with_label("gecc");
		case 'T': ++*s; return t2_widg(*((*s)++));
		case '{': 
			  ++*s; ww = widg_new(tw, s);
			  if (**s=='}') ++*s; else LOG("parse_w: exp. '}', got '%c'",**s);
			  return ww ? ww_widg(ww) : OOPS;
		default: LOG("parse_w: unexp. '%c'", *((*s)++)); return OOPS;
	}
	++*s;
	while(1) {
		int trg = 0;
		if (!**s) break;
		if (**s==cls) { ++*s; break; }
		if ((**s&120)==48) xpf = **s & 7, ++*s;
		if (**s=='<') trg = (*s)[1]-48, *s += 2;
		GtkWidget * w2 = parse_w(tw, s);
		gtk_box_pack_start(GTK_BOX(w),w2,xpf&1,(xpf&2)>>1,0);
		if (xpf&4) tw->arg[tw->ix4++].p = w2; else gtk_widget_show(w2);
		if (trg) tw->arg[trg].p = w2;
		if (**s==36) { ++*s; }
		else if (**s==cls) { ++*s; break; }
	}
	return w;
}

GtkWidget * parse_w_s(topwin * tw, const char *s) { return parse_w(tw, &s); if (*s) LOG("WARN: parse_w_s: extra chars"); }
#define TW_TOPH "(3{C_300$1$eeeeee333333...}0{__}{M_+$|+0})"

static int ww_rev_lu(char * to, ww_t * ww);
static gboolean ww_debug(ww_t * ww, int btn) {
	topwin * tw = ww->top;
	char buf[1024]; buf[ww_rev_lu(buf, ww)] = 0;
	LOG("ww_debug: twid=0x%x, wwix=0x%x, rlu=\"%s\"", tw->id, ww->ix, buf);
	return TRUE;
}

static int rec_qrv(char *to, ww_t * ww) {
	int tid = ww->top->id, id4 = tid&15;
	char * q = to+4; memcpy(to, "~QRV", 4);
	q += (id4==7) ? sprintf(q, "%06x", tid) : sprintf(q, "`%05x`%x", tid>>4, id4);
	*q = 36; q[1] = ww->top->cl->ch; q += 2+ww_rev_lu(q+2, ww);
	return q - to;
}

static void rec_vpf(ww_t * ww, const char * fmt, ...) {
	char buf[2048], *q = buf + rec_qrv(buf, ww);
	va_list ap; va_start(ap, fmt);
	q += vsnprintf(q, 999, fmt, ap);
	va_end(ap); *(q++) = 10; write(1, buf, q-buf);
}

///////////////// simple widgets /////////////////////////////////////////////

#define ENT_SL(x) ((x)->arg[1].i[0])
#define ENT_SA(x) ((x)->arg[1].i[1])

void button_click(GtkButton *btn, gpointer p) { widg_defcmd((ww_t*)p, ""); }

static void button_skel(struct _ww_t * ww, const char **pp) {
	char lb[32]; get_tok(lb, 32, pp, 36);
	ww -> w = gtk_button_new_with_label(lb);
	get_tok(ww->cmd, 16, pp, 0);
	g_signal_connect (ww->w, "clicked", G_CALLBACK (button_click), (gpointer)ww);
}

static void entry_chg(GtkEditable * ed, gpointer p) {
	ww_t * ww = (ww_t*)p;
	int flg = !!(*ww->cmd) + 2*!!(dflg & DF_REC);  if (!flg) return;
	char * s = gtk_editable_get_chars(ed, 0, -1);
	if (!ENT_SL(ww) || memcmp(ww->etc, s, ENT_SL(ww))) {
		if (flg&2) rec_vpf(ww, "e%s", s); if (flg&1) widg_defcmd(ww, s); ENT_SL(ww) = 0; }
	free(s);
}

static void entry_set(ww_t * ww, const char * s) {
	int l = ENT_SL(ww) = strlen(s) + 1;
	if (l>ENT_SA(ww)) {
		if (ENT_SA(ww)) free(ww->etc); else ENT_SA(ww) = 16;
		while (l>ENT_SA(ww)) ENT_SA(ww) *= 2;
		ww->etc = malloc(ENT_SA(ww));
	}
	memcpy(ww->etc, s, l);
	gtk_entry_set_text(GTK_ENTRY(ww->w), s);
}

static void entry_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) { LOG("entry_cmd: arg==0, cl->ch=0x%x(%c)", ww->cl->ch, ww->cl->ch); return; }
	switch(*arg) {
		case 't': entry_set(ww, arg+1); return;
		case '@': entry_set(ww, hxdoub_str(NULL, arg+1, 15)); return;
		default: LOG("entry: invalid cmd 0x%x(%c)", *arg, *arg); return;
	}
}

static int entry_get (void * to, struct _ww_t * ww, int ty) {
	int l; const char *s = gtk_entry_get_text(GTK_ENTRY(ww->w));
	switch(ty) {
		case 'i': return atoi(s);
		case 'd': *(double*)to = atof(s); return 0xffff;
		case 's': if ((l=strlen(s)) > 255) l = 255; break;
		case 'S': if ((l=strlen(s)) > 1023) l = 1023; break;
		default: return 0xfffff;
	}
	memcpy(to, s, l);
	return l;
}

static void entry_act(GtkEntry *entry, gpointer p) { ww_debug((ww_t*)p, 1); }
static void entry_skel(struct _ww_t * ww, const char **pp) {
	int x,w=0; while (x=**pp-'0',0<=x&&x<=9) w=10*w+x, ++*pp;
	GtkEntryBuffer * eb = gtk_entry_buffer_new("",0);
	if (**pp==36) ++*pp; get_tok(ww->cmd, 16, pp, 0);
	ww -> w = gtk_entry_new_with_buffer(eb);
	gtk_entry_set_width_chars(GTK_ENTRY(ww->w), w?w:10);
	gtk_entry_set_has_frame(GTK_ENTRY(ww->w), 0);
	ENT_SL(ww) = ENT_SA(ww) = 0;
	g_signal_connect (ww->w, "changed",  G_CALLBACK (entry_chg), (gpointer)ww);
	g_signal_connect (ww->w, "activate", G_CALLBACK (entry_act), (gpointer)ww);
}

///////////////// file/dir chooser dialog ////////////////////////////////////

typedef struct { int ch, ty; const char *title, *b0, *b1, *c0, *c1; GtkWidget * pg; } choo_t;
static void choo_cmd(choo_t *q, const char *s), local_error(int ec);

static chtab choo_ch;
static choo_t choo_tab[] = {{'?', 0, NULL, NULL, NULL, NULL, NULL, NULL},
	{'<'|256, GTK_FILE_CHOOSER_ACTION_OPEN, "load",     "load", "see autosaves", ".", "$A", NULL},
	{'l'|256, GTK_FILE_CHOOSER_ACTION_OPEN, "load lib", "load", NULL, "l.", "", NULL},
	{'>'|256, GTK_FILE_CHOOSER_ACTION_SAVE, "save as",  "save", NULL, "s",  "", NULL},
	{'L'|384, GTK_FILE_CHOOSER_ACTION_SAVE, "save lib", "save (R)", "save (L)", "S@R$", "S@L$", NULL},
	{'k'|512, GTK_FILE_CHOOSER_ACTION_CREATE_FOLDER,"a20(tmp) dir", "select","default", "ck","#ck",NULL},
	{'w'|512, GTK_FILE_CHOOSER_ACTION_CREATE_FOLDER,"WAV/FLAC dir", "select","default", "cw","#cw",NULL},
	{0, 0, NULL, NULL, NULL, NULL, NULL}};

static void choo_init() {
	int i,j; for (i=0; (j=choo_tab[i].ch); i++) 
		if (chtab_force(&choo_ch,j&127) != i) LOG("choo_init: ixw error"); }

static void choo_ucmd(const char * cmd, char * uri) {
	if (!uri) return local_error(EEE_FDNOSEL), LOG("choo_ucmd: missing uri");
	if (memcmp(uri, "file://", 7)) return LOG("choo_ucmd: unexp uri \"%s\"", uri);
	int clen = strlen(cmd), ulen = strlen(uri+7), cpos = 7 - clen;
	if (cpos<1) return LOG("choo_ucmd: cmd too long/me too lazy");
	uri[cpos-1] = '~'; memcpy(uri+cpos, cmd, clen); uri[7+ulen] = 10;
	write(1, uri+cpos-1, ulen+clen+2); g_free(uri);
}

static void choo_bye(GtkWidget *w, gpointer p) { ((choo_t*)p)->pg = NULL; LOG("choo_bye: %c", ((choo_t*)p)->ch&127);}
static void choo_resp(GtkDialog* choo, gint resp, void * p) {
	GtkFileChooser * fc = GTK_FILE_CHOOSER(choo);
	choo_t * q = (choo_t *)p;
	LOG("choo choo resp %d \"%s\"", resp, choo ?  gtk_file_chooser_get_uri(fc) : "BUG"); 
	switch(resp) {
		case GTK_RESPONSE_HELP:   return CMD("W#2.win.fd.%s",(q->ch&512) ? "set directory" : q->title);
		case GTK_RESPONSE_ACCEPT: return choo_ucmd(q->c0, gtk_file_chooser_get_uri(fc));
		case GTK_RESPONSE_REJECT: case GTK_RESPONSE_OK:
			switch(*q->c1) {
				case 0:   return LOG("choo_resp: unexp rej for 0x%x'%c'", q->ch&127,q->ch&127);
				case '$': if (*q->c1=='$') return choo_cmd(q, q->c1+1);
				case '#': return CMD("%s\n", q->c1+1);
				default:  return choo_ucmd(q->c1, gtk_file_chooser_get_uri(fc));     }
		case GTK_RESPONSE_CANCEL: return gtk_widget_destroy(q->pg); //, (void)(q->pg = 0);
		default: return LOG("choo_resp: unknown resp %d", resp);
	}}

#define GTBPAIR(X) GTK_STOCK_
#define CHOO_ARGL q->title, NULL, q->ty, GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL, GTK_STOCK_HELP,GTK_RESPONSE_HELP
static void choo_open(choo_t * q) {
	GtkFileFilter * ffilt = NULL;
	q->pg = q->b1 ? gtk_file_chooser_dialog_new(CHOO_ARGL, q->b1, (q->ch&128) ?
					GTK_RESPONSE_OK : GTK_RESPONSE_REJECT, q->b0,GTK_RESPONSE_ACCEPT,NULL)
		      : gtk_file_chooser_dialog_new(CHOO_ARGL, 		       q->b0,GTK_RESPONSE_ACCEPT,NULL);
	GtkFileChooser * fc = GTK_FILE_CHOOSER(q->pg);
	if (q->ch&256) { 
		ffilt = gtk_file_filter_new(); gtk_file_filter_set_name(ffilt,"lflab saves");
		gtk_file_filter_add_pattern(ffilt, "*.lf"); gtk_file_chooser_add_filter(fc, ffilt);
		ffilt = gtk_file_filter_new(); gtk_file_filter_set_name(ffilt,"all files");
		gtk_file_filter_add_pattern(ffilt, "*"); gtk_file_chooser_add_filter(fc, ffilt); }
	g_signal_connect (q->pg, "response", G_CALLBACK(choo_resp), (gpointer)q);
	g_signal_connect (q->pg, "destroy", G_CALLBACK(choo_bye), (gpointer)q);
	gtk_file_chooser_set_do_overwrite_confirmation(fc, TRUE);
	gtk_widget_show(GTK_WIDGET(fc));
}

static void choo_cmd(choo_t *q, const char *s) {
	if (!q->title) return LOG("choo_cmd: unknown choo 0x%x/%c", q->ch, q->ch&127);
	GtkFileChooser * fc = q->pg ? GTK_FILE_CHOOSER(q->pg) : NULL;
	switch(*s) {
		case 'W': if (fc) gtk_window_present(GTK_WINDOW(q->pg)); else choo_open(q);  return;
		case 'A': if (q->pg) gtk_file_chooser_set_current_folder(fc, getenv("LF_USERDIR"));
			  else LOG("choo_cmd/A: window not open");			     return;
		case 'Z': return (q->pg) ? gtk_widget_destroy(q->pg)
			  		 : LOG("choo_cmd/Z: nothing happens");
		default: return LOG("choo_cmd: unknown cmd 0x%x/%c", *s, *s);
	}}

///////////////// track div. TODO: fix prod //////////////////////////////////

typedef struct { char i,d,s,f; } tdiv_idsf_t;
typedef struct { unsigned int d; unsigned short v; unsigned char m,j; } tdiv_gd_t;
#define TDIV_N_SD 27
#define TDIV_N_D 30
const char tdiv_d[30] = {8,9,10,12,14,15,16,18,20,21,24,28,30,32,35,36,40,42,45,48,56,60,63,64,70,72,80,84,90,
        96};
const char tdiv_sd[27] = {1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,28,30,32,35,36,40,42,45,48};
const unsigned int tdiv_sdbv[30] = {0xb,0x5,0x13,0x2f,0x43,0x15,0x8b,0x127,0x21b,0x45,0x4af,0x84b,0x1237,
        0x208b,0x51,0x452f,0x829b,0x10867,0x1115,0x224af,0x408cb,0x8963f,0x10145,0x10208b,0x200a53,0x4245af,
        0x80a29b,0x1050c6f,0x2085337,0x41224af};
const tdiv_idsf_t tdiv_idsf[189] = {{0,0,0,0},{0,8,1,1},{0,8,2,0},{0,8,4,0},{1,9,1,1},{1,9,3,0},{2,10,1,1},
        {2,10,2,0},{2,10,5,0},{3,12,1,1},{3,12,2,0},{3,12,3,0},{3,12,4,0},{3,12,6,0},{4,14,1,1},{4,14,2,0},
        {4,14,7,0},{5,15,1,1},{5,15,3,0},{5,15,5,0},{6,16,1,1},{6,16,2,0},{6,16,4,0},{6,16,8,0},{7,18,1,1},
        {7,18,2,0},{7,18,3,0},{7,18,6,0},{7,18,9,0},{8,20,1,1},{8,20,2,0},{8,20,4,0},{8,20,5,0},{8,20,10,0},
        {9,21,1,1},{9,21,3,0},{9,21,7,0},{10,24,1,1},{10,24,2,0},{10,24,3,0},{10,24,4,0},{10,24,6,0},
        {10,24,8,0},{10,24,12,0},{11,28,1,1},{11,28,2,0},{11,28,4,0},{11,28,7,0},{11,28,14,0},{12,30,1,1},
        {12,30,2,0},{12,30,3,0},{12,30,5,0},{12,30,6,0},{12,30,10,0},{12,30,15,0},{13,32,1,1},{13,32,2,0},
        {13,32,4,0},{13,32,8,0},{13,32,16,0},{14,35,1,1},{14,35,5,0},{14,35,7,0},{15,36,1,1},{15,36,2,0},
        {15,36,3,0},{15,36,4,0},{15,36,6,0},{15,36,9,0},{15,36,12,0},{15,36,18,0},{16,40,1,1},{16,40,2,0},
        {16,40,4,0},{16,40,5,0},{16,40,8,0},{16,40,10,0},{16,40,20,0},{17,42,1,1},{17,42,2,0},{17,42,3,0},
        {17,42,6,0},{17,42,7,0},{17,42,14,0},{17,42,21,0},{18,45,1,1},{18,45,3,0},{18,45,5,0},{18,45,9,0},
        {18,45,15,0},{19,48,1,1},{19,48,2,0},{19,48,3,0},{19,48,4,0},{19,48,6,0},{19,48,8,0},{19,48,12,0},
        {19,48,16,0},{19,48,24,0},{20,56,1,1},{20,56,2,0},{20,56,4,0},{20,56,7,0},{20,56,8,0},{20,56,14,0},
        {20,56,28,0},{21,60,1,1},{21,60,2,0},{21,60,3,0},{21,60,4,0},{21,60,5,0},{21,60,6,0},{21,60,10,0},
        {21,60,12,0},{21,60,15,0},{21,60,20,0},{21,60,30,0},{22,63,1,1},{22,63,3,0},{22,63,7,0},{22,63,9,0},
        {22,63,21,0},{23,64,1,1},{23,64,2,0},{23,64,4,0},{23,64,8,0},{23,64,16,0},{23,64,32,0},{24,70,1,1},
        {24,70,2,0},{24,70,5,0},{24,70,7,0},{24,70,10,0},{24,70,14,0},{24,70,35,0},{25,72,1,1},{25,72,2,0},
        {25,72,3,0},{25,72,4,0},{25,72,6,0},{25,72,8,0},{25,72,9,0},{25,72,12,0},{25,72,18,0},{25,72,24,0},
        {25,72,36,0},{26,80,1,1},{26,80,2,0},{26,80,4,0},{26,80,5,0},{26,80,8,0},{26,80,10,0},{26,80,16,0},
        {26,80,20,0},{26,80,40,0},{27,84,1,1},{27,84,2,0},{27,84,3,0},{27,84,4,0},{27,84,6,0},{27,84,7,0},
        {27,84,12,0},{27,84,14,0},{27,84,21,0},{27,84,28,0},{27,84,42,0},{28,90,1,1},{28,90,2,0},{28,90,3,0},
        {28,90,5,0},{28,90,6,0},{28,90,9,0},{28,90,10,0},{28,90,15,0},{28,90,18,0},{28,90,30,0},{28,90,45,0},
        {29,96,1,1},{29,96,2,0},{29,96,3,0},{29,96,4,0},{29,96,6,0},{29,96,8,0},{29,96,12,0},{29,96,16,0},
        {29,96,24,0},{29,96,32,0},{29,96,48,0}};
unsigned const char tdiv_dvix[30] = {1,4,6,9,14,17,20,24,29,34,37,44,49,56,61,64,72,79,86,91,100,107,118,123,
        129,136,147,156,167,178};
#define TDIV_MENU_LN \
{0,31,0,3,1,"def/08/09/10/12/14/15/16/18/20/21/24/28/30/32/35/36/40/42/45/48/56/60/63/64/70/72/80/84/90/96",\
        "0123456789:;<=>?@ABCDEFGHIJKLMN"},\
{0,27,0,3,1,":01:02:03:04:05:06:07:08:09:10:12:14:15:16:18:20:21:24:28:30:32:35:36:40:42:45:48",\
        "0123456789:;<=>?@ABCDEFGHIJ"},
#define TDIV_N_PPB 45
const tdiv_gd_t tdiv_gd[96] = {{0x0,1,0,0},{0x0,2,0,0},{0x0,4,0,0},{0x1,8,15,0},{0x41,16,31,0},
        {0x2041,32,31,0},{0x802041,64,30,7},{0x0,128,0,18},{0x0,5,0,0},{0x4,10,11,0},{0x104,20,27,0},
        {0x10105,40,91,0},{0x4010145,80,91,10},{0x4012145,160,83,21},{0x4812145,320,80,32},{0x0,640,0,43},
        {0x0,7,0,0},{0x10,14,7,0},{0x810,28,23,0},{0x100811,56,151,4},{0x100851,112,151,15},
        {0x102851,224,147,26},{0x902851,448,144,37},{0x0,896,0,0},{0x4000,35,3,0},{0x1004014,70,211,8},
        {0x1004914,140,211,19},{0x1114915,280,209,29},{0x5114955,560,208,40},{0x0,1120,0,0},{0x0,2240,0,0},
        {0x0,4480,0,0},{0x0,3,0,0},{0x0,6,0,0},{0x8,12,15,0},{0x409,24,63,0},{0x80449,48,63,3},
        {0x20082449,96,63,13},{0x20882449,192,50,24},{0x0,384,0,35},{0x20,15,11,0},{0x1024,30,59,0},
        {0x20112c,60,123,5},{0x21152d,120,115,16},{0x429156d,240,115,27},{0x2429356d,480,112,38},{0x0,960,0,0},
        {0x0,1920,0,0},{0x200,21,7,0},{0x20210,42,55,1},{0x8020a18,84,183,11},{0x8120e19,168,179,22},
        {0x81a0e59,336,177,33},{0x281a2e59,672,176,44},{0x0,1344,0,0},{0x0,2688,0,0},{0x4220,105,227,14},
        {0x1025234,210,243,25},{0x9225b3c,420,240,36},{0x0,840,0,0},{0x0,1680,0,0},{0x0,3360,0,0},
        {0x0,6720,0,0},{0x0,13440,0,0},{0x2,9,13,0},{0x82,18,29,0},{0x808a,36,61,0},{0x200848b,72,61,9},
        {0x20884cb,144,53,20},{0x2208a4cb,288,49,30},{0x2288a4cb,576,48,41},{0x0,1152,0,0},{0x40022,45,105,2},
        {0x100410a6,90,121,12},{0x102491ae,180,113,23},{0x122595af,360,113,34},{0x162d95ef,720,112,45},
        {0x0,1440,0,0},{0x0,2880,0,0},{0x0,5760,0,0},{0x400202,63,165,6},{0x420292,126,181,17},
        {0x8428a9a,252,177,28},{0xa528e9b,504,176,39},{0x0,1008,0,0},{0x0,2016,0,0},{0x0,4032,0,0},
        {0x0,8064,0,0},{0x444222,315,225,31},{0x114652b6,630,240,42},{0x0,1260,0,0},{0x0,2520,0,0},
        {0x0,5040,0,0},{0x0,10080,0,0},{0x0,20160,0,0},{0x0,40320,0,0}};
const char tdiv_ppb_ix[47] = {0,49,72,36,19,42,80,6,25,67,12,50,73,37,56,20,43,81,7,26,68,13,51,74,38,57,21,44,
        82,27,69,88,14,52,75,39,58,22,45,83,28,70,89,15,53,76,0};

///////////////// menu / shortlist ///////////////////////////////////////////

typedef struct _noderef_t {
	int id;
	char ty, len, nm[21], rsrv[5];
} noderef_t;

noderef_t lsr_nodes[96];
int n_lsr[3];

typedef struct _menu_t {
	int ty, ni, widf;
	int lbll, cmdl;
	const char *lbl, *cmd;
} menu_t;

menu_t menutab[] = { {'?',0,0,0,0,NULL,NULL},
{'0',4, 0,7,2,"one daya menu will behere   ","########"},
{'W',8, 1,2,1,"x y s1s2s3s4s5s6", "xy123456"},
{0,  9, 1,3,1,"concn1linlg sq cu 1/x/sq/cu", "'\"-lqchQC"},
{0,  5, 1,1,1,".=*/A",".=*/A"},
{0,  9, 1,2,1, "/-/L/Q/e/g/s/E/G/S", "-LQegsEGS"},
{0,  9, 1,2,1, "\\A\\L\\Q\\e\\g\\s\\E\\G\\S", "ALQegsEGS"},
{0,  32,5,3,1, "2:22:33:22:44:23:32:55:22:63:44:36:22:77:23:55:32:84:48:22:93:66:39:24:55:4"
	       "3:77:33:84:66:48:35:5", "0182@93H4:AP5X;I6B`7<QhCJ=Y>DRaK"},
{0,  32,5,3,1, "C  2:33:22:44:23:32:55:22:63:44:36:22:77:23:55:32:84:48:22:93:66:39:24:55:4"
	       "3:77:33:84:66:48:35:5", "p182@93H4:AP5X;I6B`7<QhCJ=Y>DRaK"},
{0,  13,7,3,1, "C  C# D  D# E  F  F# G  G# A  A# H  3:4", "pqrstuvwxyz{:" /*"}"*/ },
{'.',18,0,12,4, "filter disp.audio configmain config console     error list  ------------flush log   "
		"write tlog  save config ------------exit(autosv)restart(asv)restart GUI ------------"
		"exit w/o a/srestart-noASSIGABRT     SIGKILL     ",
		"_F  A0W cW  _c-1_E  ####*f  $!t c>  ####$!q0$!q2$!q ####$!q1$!q3#%$6#%$9"},
{'/',7, 0, 7,1, "folder clipbrdinstr. graph  iter.f.calc   track  ","DCwgict"},
{0,  9, 0,8,5,"load    save    save as --------load libsave lib--------exit+AS rstrt+AS",
	          "$!f<Ws    $!f>W#####lW   $!fLW#####$!q0 $!q2 "},
{'#',5, 0, 8, 3, "[focus] config  help    redraw  wav/flac", "$>*W$ N??M9 XAW"},
{0,  4, 1, 4, 1, "clikholdtggluniq", "0123"},
{'_',0, 0, 0,0, "(none)","0"},
{'+',6, 0, 7,2, "help   info1  info2  info*  GUI cfgdelete ", "N?I1I2I3EWNd"},
{0, 24, 0, 8,13, ">@L/src >@R/src >clip/s.+anon/s.>@L/filt>@R/filt>clip/f.+anon/f.anon#00 anon#01 anon#02 "
"anon#03 anon#04 anon#05 anon#06 anon#07 anon#08 anon#09 anon#10 anon#11 anon#12 anon#13 anon#14 anon#15 ",
"$~C@L.wr$wb%@$~C@R.wr$wb%@$~C@C._$wb%@ Cw_$wb%@     "
"$~C@L.wr$w>%@$~C@R.wr$w>%@$~C@C._$w>%@ Cw_$w>%@     $.Ww0        $.Ww1        $.Ww2        $.Ww3        $.Ww4        $.Ww5        "
"$.Ww6        $.Ww7        $.Ww8        $.Ww9        $.WwA        $.WwB        $.WwC        $.WwD        $.WwE        $.WwF        "},
{'T', 9, 0, 8, 1, "[toggle]cut     copy    paste   [grid]  [config]copy selmove selnew...  ", "txcv96CMN"},
TDIV_MENU_LN
{0, 8,0,2,1,"*2*3*5*7/2/3/5/7","01234567"},
{0, 32,0,6,1,"play  stop  play1 play2 play3 play4 play6 play8 play12play16play22play30loop1 loop2 loop3 loop4 loop5 loop6 "
"loop7 loop8 loop9 loop10loop11loop12loop13loop14loop15loop16loop17loop18loop19loop20","103579=AIQ]m2468:<>@BDFHJLNPRTVX"},
{'K', 9, 0,5,6, "help info del/snew  -----cleardel  -----WAV/s", "N?    I3    Kd    Cw%K  ##    KZ    Nd    ######$.X*$AW"},
{'k', 4, 1,4,1, "ask keepwav flac", "0123"},
{'S', 1, 0,1,1, "??", "##" },
{0,   3, 1,5,1, "availdelayretBS", "012"},
{0,   7, 0,4,4, "lr  rl  lrc clr lrcclrlrzzlr", "lr  rl  lrc clr lrcclrlrzzlr"},
{'i', 9, 1,3,1, "conask/cu/sq1/xloglinsq cu", "012345678"},
{'g', 4, 0,9,2, "[shuffle]rgb:sel  inlbl:selinlbl:all", "s UrUiUI"},
{-1,0,0,0,0,NULL,NULL} };

static char * menu_txt = NULL;
static int n_menu_t = 0;

char lr_dirname[48];

static chtab menu_ctab;
static ww_t * cmenu_w = NULL;
static menu_t * cmenu_mt = NULL;
static GtkWidget * cmenu_gw = NULL;
static char cmenu_dcmd[1024], cmenu_dlbl[1024];
static int cmenu_dic;
static unsigned int cmenu_msk;

void lsr_upd(int k, const char * s) {
	noderef_t * np = lsr_nodes + 32 * k;
	int i, j, c;
	for (i=0; i<32; i++, np++) { 
		if (*s==36) ++s; if (!*s) break;
		np->ty = *(s++);
		j=0; while (*s&80) j = 16*j+hxd2i(*s), s++; np->id = j;
		for (j=0,++s; c=s[j], c && c-36 && j<20; j++) np->nm[j] = c;
		s += j; np->nm[j] = 0; np->len = (char)j;
	}
	n_lsr[k] = i;
}

void lsr_debug() {
	int i, j, n; for (i=0; i<3; i++) {
		LOG("lsr_debug: n(%c) = %d", "LSR"[i], (n = n_lsr[i]));
		for (j=0; j<n; j++) {
			noderef_t * p = lsr_nodes + 32*i + j;
			LOG("%02d: %05x %c %02d %s", 
				j, p->id, p->ty, p->len, p->nm);
		}}}

static int str_unpack(char *to, const char * s, int ni, int l) {
	char * p = to;
	int i,j,k;
	for (i=0; i<ni; i++) {
		for (j=0; j<l; j++) p[j]=s[j];
		p[j]=0; k=j;
		for (--j; j>0 && p[j]==32; j--) p[j] = 0;
		p += (k+1); s += k;
	}
	return p - to;
}

static int menutab_lu(int t, int j) {
	int i = chtab_get(&menu_ctab, t); if (!i) {
		if (t-63) LOG("menutab_lu: not found: 0x%x (%c)", t, t); return 1; }
	i = (menutab[i].ty>>16) + j; 
	if (i<1 || i>=n_menu_t) {
		LOG("menutab_lu: invalid idx %d", i); return 1; }
	return i;
}

static int tsc_mi, tsc_x, tsc_y, tsc_id, tsc_md,  // trk
	   tsc_drag_id, tsc_drag_xd, tsc_drag_x, tsc_drag_y;
static ww_t * tsc_ww;
static void tsc_act(ww_t * ww, int ix);

static void menutab_init() {
	int i, j, w, w2, k = 0, n = 0;
	chtab_ini(&menu_ctab, 0);
	menu_t * p = menutab; 
	while (p->ty>=0) k += p->ni * (p->lbll+2+p->cmdl), p++, n++;
	n_menu_t = n; p = menutab;
	char * q = menu_txt = malloc(k);
	k = 0;
	for (i=0; i<n; i++, p++) {
		if (p->ty) j=chtab_force(&menu_ctab, k=p->ty), menutab[j].ty += 65536*i; 
		else p->ty = (k += 256);
		j = str_unpack(q, p->lbl, p->ni, p->lbll); p->lbl = q; q += j; ++p->lbll;
		j = str_unpack(q, p->cmd, p->ni, p->cmdl); p->cmd = q; q += j; ++p->cmdl;  w2 = 0;
		if (p->widf&1) for (j=0; j<p->ni; j++) if ((w=tx_len(conf_lbfs,p->lbl+j*p->lbll)) > w2) w2=w;
		p->widf <<= 10; p->widf |= w2;
	}
	tsc_mi = menutab_lu('T', 0);
}

static int menu_findlbl(const char *s) {
	int i, n = cmenu_mt->ni;
	for (i=0; i<n; i++) {
		if (strcmp(s, cmenu_mt->lbl + i*cmenu_mt->lbll)) continue;
		if (!((1u<<i)&cmenu_msk)) return i; 
		LOG("menu_findlbl(%s) : skipping %d (hidden)", s, i); }
	for (i=0; i<cmenu_dic; i++) if (!strcmp(s, cmenu_dlbl+32*i)) return i;
	return -1;
}

static const char * menu_lbl(int ix) {
	if (ix < cmenu_mt->ni) return cmenu_mt->lbl + ix*cmenu_mt->lbll;
	return ((ix-=cmenu_mt->ni) < cmenu_dic) ? cmenu_dlbl + 32*ix : "???"; }

static void menu_act(GtkMenuItem *it, gpointer p) {
	size_t ix = (char*)p - cmenu_dcmd;
	if (it && (dflg & DF_REC)) CMD("QRm:%s", menu_lbl(ix));
	if (ix<0 || ix>31) LOG("menu_act: ptr:%p dcmd:%p diff:%d", (char*)p, cmenu_dcmd, ix);
	else if (cmenu_w->cl->ch=='t') tsc_act(cmenu_w, ix);
	else if (ix < cmenu_mt->ni) widg_defcmd(cmenu_w, cmenu_mt->cmd + ix*cmenu_mt->cmdl);
	else if ((ix-=cmenu_mt->ni) < cmenu_dic) widg_defcmd(cmenu_w, cmenu_dcmd + 32*ix);
	else LOG("menu_act: invalid ix2=%d", ix);
}

static void menu_del() { if (cmenu_gw) gtk_widget_destroy(cmenu_gw); cmenu_gw = NULL; }

static void menu_act_s(const char *s) {
	if (!cmenu_w) return LOG("menu_act_s: no active popup menu");
	int j = menu_findlbl(s); if (j<0) return LOG("menu_act_s: \"%s\" not found", s);
	LOG("menu_act_s: found: %d", j); menu_act(NULL, cmenu_dcmd + j); menu_del();
}

static int add_dyn(const char * nm, const char * cmd, int x) {
	int i = cmenu_dic; if (i==32) return -1; else cmenu_dic++;
	GtkWidget * it2 = gtk_menu_item_new_with_label(nm);
	strncpy(cmenu_dlbl + 32*i, nm, 32);
	gtk_menu_shell_append(GTK_MENU_SHELL(cmenu_gw), it2);
	gtk_widget_show(it2);
	g_signal_connect(it2, "activate", G_CALLBACK(menu_act), cmenu_dcmd + i);
	char *q0 = cmenu_dcmd + 32*i, *q = q0, *qlim = q0 + 31;
	while (*cmd && q<qlim) *(q++) = *(cmd++);
	if (x>=0 && q+5<qlim) q += hx5(q, x);
	*q = 0; return i;
}

static void lsr_popup(int flg) {
	if (dflg & DF_MENU) LOG("lsr_popup: flg = 0x%x", flg);
	int k = flg & 3, zf = (flg>>2)&1;
	if (k==3) { LOG("lsr_popup: invalid ix: 3"); return; }
	noderef_t * np = lsr_nodes + 32 * k;
	if (zf) add_dyn("(none)", "0", -1);
	int i, n = n_lsr[k];
	char nm[24]; nm[1] = 60+k; nm[2] = 32;
	for (i=0; i<n; i++) {
		nm[0] = np[i].ty; memcpy(nm+3, np[i].nm, np[i].len+1);
		if (add_dyn(nm, "#", np[i].id)<0) break;
	}
}

static void pcm_popup() {
	unsigned char ls[16]; 
	char buf[12]; memcpy(buf, "plughw:x,y\0",12);
	int i, n = find_dev(ls, 0, 15);
	add_dyn("default", "default", -1);
	for (i=0; i<n; i++) buf[7] = 48+(ls[i]>>4), buf[9] = 48+(ls[i]&15), add_dyn(buf+4, buf+4, -1);
	for (i=0; i<n; i++) buf[7] = 48+(ls[i]>>4), buf[9] = 48+(ls[i]&15), add_dyn(buf  , buf  , -1);
}

static void popup2(ww_t * ww, int tid, unsigned int msk, int btn, GdkEventButton * ev) {
	menu_del();
	if (tid<1 || tid>=n_menu_t) { LOG("popup2: invalid tid %d", tid); return; }
	cmenu_w = ww; cmenu_mt = menutab + tid; cmenu_dic = 0; cmenu_msk = msk;
	int i, n = cmenu_mt->ni;
	const char * lbl = cmenu_mt->lbl;
	if (*lbl=='[' && btn==1) return (void) widg_defcmd(ww, cmenu_mt->cmd);
	cmenu_gw = gtk_menu_new(); gtk_widget_set_name(cmenu_gw, "lflabPU");
	if (n) { for (i=0; i<n; i++, lbl += cmenu_mt->lbll) {
		if (msk & (1u<<i)) continue;
		GtkWidget * it2 = memcmp(lbl, "----", 4) ? gtk_menu_item_new_with_label(lbl) 
							 : gtk_separator_menu_item_new();
		gtk_menu_shell_append(GTK_MENU_SHELL(cmenu_gw), it2);
		gtk_widget_show(it2);
		g_signal_connect(it2, "activate", G_CALLBACK(menu_act), cmenu_dcmd + i);
	}} else { switch(msk&15) {
		case 1: pcm_popup(); break;
		case 3: lsr_popup((int)(msk>>4u)); break;
		default: LOG("undef _menu %d", (int)(msk&15)); break;
	}}
	gtk_menu_popup(GTK_MENU(cmenu_gw), NULL, NULL, NULL, NULL, btn, ev ? ev->time : 0);
	gtk_widget_show(cmenu_gw);
}

///////////////// lazy vbox //////////////////////////////////////////////////

#define VB_LMAX(x) ((x)->arg[1].c[0])
#define VB_LCRE(x) ((x)->arg[1].c[1])
#define VB_LBV(x) ((x)->arg[4].i[0])
#define VB_WPL(x) ((x)->arg[1].c[3])
#define VB_WBASE(x) ((x)->arg[2].i[0])
#define VB_ARG(x) ((x)->arg[2].i[1])
#define VB_LINE(x, i) ( ((GtkWidget**)(x)->etc)[i] )

static ww_t * wlu_any_pp(topwin * tw, const char ** pp) { 
	ww_t * ww = widg_p(tw, widg_lookup(tw, pp, 0));
	if (!ww) return NULL;
	if (!ww->cl) return LOG("widget %p(%x:%x) has zero class", ww, tw->id, ww->ix), NULL;
	if (ww->cl->ch!=':') return ww;
	if (**pp=='.') return ++*pp, ww;
	int k = hex2(*pp); *pp += 2;
	//LOG("wlu_any: Cxx k=%d(0x%x), WBASE=%d", k, k, VB_WBASE(ww));
	return widg_p(tw, VB_WBASE(ww) + k);
}

static ww_t * wlu_any_s(topwin * tw, const char * p) { return wlu_any_pp(tw, &p); }

static int ww_rev_lu(char * to, ww_t * ww) {
	topwin * tw = ww->top;
	const char * q = tw->sub_ch.tab;
	int i, ix = ww->ix;
	for (i=0; i<96; i++) {
		int j = q[i];
		if (j==ix) return *to=i+32, 1;
		if (j>99) continue;
		ww_t * w2 = widg_p(tw, j); if (w2->cl->ch != ':') continue;
		int ns = VB_WPL(w2) * VB_LMAX(w2), k = ix - VB_WBASE(w2);
		if (k>=0 && k<ns) return *to=i+32, to[1]=hexc1(k>>4), to[2]=hexc1(k&15), 3;
	}
	return ((ix = 1023-ix) < 256) ? (*to='_', to[1]=hexc1(ix>>4), to[2]=hexc1(ix&15), 3) : 0;
}

static GtkWidget* vbox_mkline(struct _ww_t * ww, int ix) {
	topwin * tw = ww->top;
	if (tw->vb_c_ix) { LOG("ERROR: lazy vboxes cannot be nested"); return OOPS; }
	tw->vb_c_ix = VB_WBASE(ww) + ix * VB_WPL(ww);
	GtkWidget * ln = (*(vbox_line_fun)ww->arg[3].p) (ww, ix);
	tw->vb_c_ix = 0;
	gtk_box_pack_start(GTK_BOX(ww->w),ln,0,0,0);
	return VB_LINE(ww, ix) = ln;
}

static void vbox_show_bv(struct _ww_t * ww, int bv) {
	int i, msk, bv0 = VB_LBV(ww), n=VB_LMAX(ww); 
	if (bv==bv0) return; else VB_LBV(ww) = bv;
	for (i=0, msk=1; i<n; i++, msk*=2) {
		if (!(msk & (bv^bv0))) continue;
		if (msk&bv0) { gtk_widget_hide(VB_LINE(ww, i)); continue; }
		while (VB_LCRE(ww)<i+1) vbox_mkline(ww, VB_LCRE(ww)++);
		gtk_widget_show(VB_LINE(ww, i));
	}
	gtk_container_check_resize(GTK_CONTAINER(ww->top->arg[0].p));
}

static void upd_flgvec(topwin * tw, const char* s, int n, int oldbv, int newbv);
static void vbox_cmd(struct _ww_t * ww, const char * arg) {
	if (dflg&DF_WIDG) LOG("vbox_cmd=\"%s\"", arg);
	int k=0; switch(*arg) {
		case 'W': upd_flgvec(ww->top, "G", 4, VB_LBV(ww), k = 1<<(arg[1]-48)); break;
		case '.': k = (1<<(arg[1]-48)); break;
		case '^': k = VB_LBV(ww) ^ (1<<(arg[1]-48)); break;
		case '+': k = (1<<(arg[1]-48)) - 1; break;
		case '*': ++arg; while(*arg&80) k = 16*k + hxd2i(*(arg++));
			  break;
		default: LOG("vbox: unknown cmd 0x%x(%c)",*arg,*arg); return;
	}
	vbox_show_bv(ww, k);
}

static void vbox_skel(struct _ww_t * ww, const char **pp) {
	int ty = *((*pp)++); switch(ty) {
		case 'w': ww->arg[3].p = &wrap_vbl_i; break;
		case 'W': ww->arg[3].p = &wrap_vbl_t; break;
		case 'c': ww->arg[3].p = &calc_vbl; break;
		case 'g': ww->arg[3].p = &gconf_vbl; break;
		case 'd': ww->arg[3].p = &doc_vbl; break;
		case 'e': ww->arg[3].p = &err_vbl; break;
		case 'K': ww->arg[3].p = &clip_vbl; break;
		default: LOG("vbox: unknown subcl 0x%x(%c)",ty,ty); return;
	}
	VB_LMAX(ww) = *((*pp)++) - 48;
	VB_WPL(ww) = *((*pp)++) - 48;
	int n = VB_LMAX(ww) * VB_WPL(ww); if (n>256) {
		LOG("vbox: n(%d) > 256", n); return; }
	int k = 0; while (**pp != '}') k = 16*k + hxd2i(*((*pp)++));
	VB_ARG(ww) = k;
	VB_LCRE(ww) = VB_LBV(ww) = 0;
	VB_WBASE(ww) = (ww->top->anon_w -= n) & 1023;
//	LOG("vbox: n=%d, WBASE=%d", n, VB_WBASE(ww));
	ww->etc = calloc(n, sizeof(GtkWidget*));
	ww->w = gtk_vbox_new(0, 0);
}

///////////////// drawing area (gen) /////////////////////////////////////////

#define DA_W(x) ((x)->arg[1].s[0])
#define DA_H(x) ((x)->arg[1].s[1])
#define DA_FLG(x) ((x)->arg[1].s[2])
#define DA_XY01 short * a4 = ww->arg[4].s; int x0 = a4[0], y0 = a4[1], x1 = a4[2], y1 = a4[3]

static gboolean da_conf(GtkWidget *w, GdkEventConfigure * ev, gpointer p) {
	ww_t * ww = (ww_t*)p;
	GtkAllocation alc; gtk_widget_get_allocation(ww->w, &alc);
	DA_W(ww) = alc.width; DA_H(ww) = alc.height;
	if (!(ww->cl->flg & WF_SURF)) return TRUE;
	cairo_surface_t * surf = (cairo_surface_t*)DA_SURF(ww);
	if (surf) {
		if (!(ww->cl->flg&WF_RESIZE)) return TRUE;
		if (cairo_surface_get_reference_count(surf)) cairo_surface_destroy(surf);
	}
	int wi = alc.width, he = alc.height;
	if (ww->cl->flg & WF_FULLSURF) {
		int hw = (*ww->cl->get)(NULL, ww, 'B'), sw = hw&0xffff, sh = hw>>16;
		if (wi<sw || he<sh) return TRUE; else wi=sw, he=sh;
	}
	DA_SURF(ww) = gdk_window_create_similar_surface(gtk_widget_get_window(w),
			CAIRO_CONTENT_COLOR, wi, he);
	(*ww->cl->cmd)(ww, NULL);
	return TRUE;
}

#define LAZYSURF_HEAD(W,H) static cairo_surface_t * p = 0; if (p) return p; \
	p = gdk_window_create_similar_surface(gtk_widget_get_window(ww->w), CAIRO_CONTENT_COLOR, W, H); \
	cairo_t * cr2 = cairo_create(p)

static cairo_surface_t * default_lzs(ww_t * ww) {
	LAZYSURF_HEAD(100, 100); cairo_set_source_rgb(cr2, 1.0, 0.0, 0.0); cairo_paint(cr2);
	cairo_destroy(cr2); return p; }

static cairo_surface_t *daclip_lzs(ww_t * ww), *dapz_lzs(ww_t * ww);
static cairo_surface_t * lazysurf(ww_t * ww) { int cl = ww->cl->ch; switch(cl) {
	case 'K': return daclip_lzs(ww);
	case 'P': return dapz_lzs(ww);
	default: LOG("undef lazysurf, ix 0x%x, class 0x%x '%c'", ww->ix, cl, cl); return default_lzs(ww);
}}

#define DA_DREV "expose-event"
static int bigda_counter = 0;
static gboolean da_draw(GtkWidget *w, GdkEventExpose *ev, gpointer p) {
	ww_t * ww = (ww_t*)p;
	cairo_t * cr = gdk_cairo_create(gtk_widget_get_window(w));
	gdk_cairo_region (cr, ev->region); 
	if (!(ww->cl->flg&WF_CLIPL8R)) cairo_clip (cr);
	int k; if ((k = ww->cl->flg & WF_BIGDA3)) {
		if (k==WF_BIGDA1) cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
		else if (k==WF_BIGDA2) cairo_set_source_rgb(cr, 0.2, 0.2, 0.2);
		else cairo_set_source_surface(cr, lazysurf(ww), 0.0, 0.0);
		cairo_paint(cr);
		int i, nrec; GdkRectangle * prec;
		gdk_region_get_rectangles(ev->region, &prec, &nrec);
		for (i=0; i<nrec; i++) {
			bigda_counter = i;
			short x0 = ww->arg[4].s[0] = (short)prec[i].x; 
			short y0 = ww->arg[4].s[1] = (short)prec[i].y; 
			ww->arg[4].s[2] = x0 + (short)prec[i].width - 1;
			ww->arg[4].s[3] = y0 + (short)prec[i].height - 1;
			ww->arg[0].p = cr; (*ww->cl->cmd)(ww, NULL);
		}
		g_free(prec);
	} else if (ww->cl->flg & WF_SURF) {
		double x0 = 0.0, y0 = 0.0;
		if (ww->cl->flg & WF_FULLSURF) {
			cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
			cairo_paint(cr);
			int hw = (*ww->cl->get)(NULL, ww, 'B'), sw = hw&0xffff, sh = hw>>16;
			int xi = DA_W(ww)-sw, yi = DA_H(ww)-sh; 
			if(xi>1) x0 = (double)(xi>>1); if(yi>1) y0 = (double)(yi>>1);
		}
		cairo_surface_t * sf = ((ww_t*)p)->arg[0].p;
		cairo_set_source_surface (cr, sf, x0, y0);
		cairo_paint (cr); 
	} else {
		ww->arg[0].p = cr; (*ww->cl->cmd)(ww, NULL);
	}
	cairo_destroy(cr); return FALSE;
}

static void debug_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	LOG("click: button %d, cx=%d, cy=%d -- nothing happens", b9, cx, cy); }

static gboolean da_click(GtkWidget *w, GdkEventButton * ev, gpointer p) {
	if (ev->type != GDK_BUTTON_PRESS) return TRUE;
	ww_t * ww = (ww_t*)p;
	int btn = ev->button + 3*!!(ev->state&GDK_SHIFT_MASK) + 6*!!(ev->state&GDK_CONTROL_MASK),
	    x = (int)lround(ev->x), y = (int)lround(ev->y);
	if (dflg&DF_REC) { if (!(ww->cl->flg&(WF_BIGDA3|WF_SURF))) rec_vpf(ww, "c%c", 48+btn);
			   else if (ww->cl->ch=='K') rec_vpf(ww, "c%c,%x", 48+btn, (y<<12)|x); }
	if (btn>9) return ww_debug(ww, btn);
	(ww->cl->clk) (ww, btn, (int)lround(ev->x), (int)lround(ev->y), ev);
	return TRUE;
}

static gboolean da_key(GtkWidget *w, GdkEventKey * ev, gpointer p) {
	static unsigned int glob_keyflg[4];
	ww_t * ww = (ww_t*)p;
	int kpf = (ev->type==GDK_KEY_PRESS), kc = ev->hardware_keycode & 127;
	unsigned int kmsk = 1u<<(kc&31), *pkf = glob_keyflg + (kc>>5);
	if (!kpf) (ev->type!=GDK_KEY_RELEASE) ? LOG("da_key: unknown type 0x%x", ev->type)
		     : (*pkf &= ~kmsk, (ww->cl->clk)(ww, -2, kc, ev->keyval, NULL));
	else if (!(*pkf & kmsk)) *pkf |= kmsk,  (ww->cl->clk)(ww, -1, kc, ev->keyval, NULL);
	return TRUE;
}
				
static void da_skel(struct _ww_t * ww, int wid, int heig) {
	static const int evmask[4] = { EVMASK_DEF, EVMASK_DEF|EVMASK_KEY, EVMASK_DEF|EVMASK_XM1, 
					EVMASK_DEF|EVMASK_KEY|EVMASK_XM1 };
	gtk_widget_set_size_request(ww->w = gtk_drawing_area_new(), wid, heig);
	ww->arg[0].p = NULL;
	g_signal_connect (ww->w, "configure-event", G_CALLBACK (da_conf), (gpointer)ww);
	g_signal_connect (ww->w, DA_DREV, G_CALLBACK (da_draw), (gpointer)ww);
	int flg = ww->cl->flg; if (!flg) return; 
	gtk_widget_add_events(ww->w, evmask[(flg>>WF_EV_SHIFT)&3]);
	g_signal_connect (ww->w, "button-press-event", G_CALLBACK(da_click), (gpointer)ww);
	if (flg&WF_KEYEV) g_signal_connect (ww->w, "key-press-event",   G_CALLBACK(da_key), (gpointer)ww),
			  g_signal_connect (ww->w, "key-release-event", G_CALLBACK(da_key), (gpointer)ww),
			  gtk_widget_set_can_focus(ww->w, TRUE),
			  gtk_window_set_focus(GTK_WINDOW(ww->top->w), ww->w);

}

static void da_fullre(struct _ww_t * ww) {
	if (!ww->w) return;
	struct _GdkRectangle rct; rct.x = rct.y = 0;
	rct.width = DA_W(ww); rct.height = DA_H(ww);
	if (DA_H(ww)<=0 || DA_W(ww)<=0) return;
	gdk_window_invalidate_rect(gtk_widget_get_window(ww->w), &rct, 0);
}

static void fullsurf_invd(ww_t * ww, int x0, int y0, int x1, int y1) {
	int xy00 = (*ww->cl->get)(NULL, ww, 'B');
	int w = x1-x0+3, h = y1-y0+3;
	x0 += ((DA_W(ww)-(xy00&0xffff))>>1) - 1; 
	y0 += ((DA_H(ww)-(xy00>>16))>>1) - 1; 
	if (x0<0) x0 = 0;
	if (y0<0) y0 = 0; 
	if (w>(DA_W(ww)-x0)) w = DA_W(ww)-x0;
	if (h>(DA_H(ww)-y0)) h = DA_H(ww)-y0;
	struct _GdkRectangle rct; rct.x = x0; rct.y = y0; rct.width = w; rct.height = h;
	gdk_window_invalidate_rect(gtk_widget_get_window(ww->w), &rct, 0);
}

///////////////// drawing area (aux) /////////////////////////////////////////

static void dasep_skel(struct _ww_t * ww, const char **pp) {
	int k = **pp<'}' ? *((*pp)++) - 48 : 1;
	da_skel(ww, k, k); }

static void dasep_cmd(struct _ww_t * ww, const char * arg) {
	cairo_t * cr = (cairo_t*)ww->arg[0].p;
	cairo_set_source_rgb(cr, .5, .5, .5); cairo_paint(cr); }

static void da_frame(struct _ww_t * ww, cairo_t * cr) {
	cairo_set_line_width (cr, 1.0);
	cairo_set_source_rgb (cr, .6667, .6667, .6667);
	cairo_rectangle(cr, 0.5, 0.5, (double)(DA_W(ww)-1), (double)(DA_H(ww)-1));
	cairo_stroke(cr); }

static inline void cr_line(cairo_t * cr, double x0, double y0, double x1, double y1) {
	cairo_move_to(cr, x0, y0); cairo_line_to(cr, x1, y1); }

static void get_rgb8(double * to, int hue8, double c) {
	c *= (1.0 / 1120.0);
	int h7 = hue8&127, flg = -((hue8&128)>>7), sg = 2*flg + 1;
	int r,g,b;
	if      (h7<56) r=1120&~flg, g=(1120&flg)+sg*20*h7, b=1120&flg;
	else if (h7<96) r=(1120&~flg)-sg*28*(h7-56), g=1120&~flg, b=1120&flg;
	else            r=1120&flg, g=1120&~flg, b=(1120&flg)+sg*35*(h7-97);
	to[0] = c*(double)r; to[1] = c*(double)g; to[2] = c*(double)b;
}

static void please_wait(cairo_t * cr2) {
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE);
	cairo_set_source_rgb (cr2, 1.0, 1.0, .2);
	cairo_set_line_width (cr2, 1.0);
	cairo_move_to(cr2, 7.5,  2.5); cairo_line_to(cr2, 25.5, 29.5);
	cairo_line_to(cr2, 7.5, 29.5); cairo_line_to(cr2, 25.5, 2.5);
	cairo_line_to(cr2, 7.5, 2.5); cairo_stroke(cr2);
	return;
}

///////////////// drawing area (bmp/cpu-meter) + tlog ////////////////////////

static unsigned char cpu_bmp[160];
static cairo_surface_t * cpu_surf;
static int cpu_austat = 1;
static unsigned int tlog_dat[0x101000];
static int tlog_ix_c = 0, tlog_ix_i = 0, tlog_wr = 0;

static void dabmp_draw(ww_t * ww, cairo_t * cr2) {
	static const double fg[2]={.2, 1.0}, bg[2]={.1,.3};
	static int state = 0;
	unsigned int * q = (unsigned int*)cpu_bmp;
	if (!state && q[2]!=0x80000100) state=1, CMD("Rg");
	int rf = !!(q[9]&0xf8000000), gf = !(q[3]&0xf8000000);
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE);
	if (cpu_austat) cairo_set_source_rgb (cr2, bg[rf], bg[gf], *bg); 
	else cairo_set_source_rgb(cr2, 0.0, 0.0, 0.0);
	cairo_paint(cr2);
	da_frame(ww, cr2);
	if (cpu_austat) cairo_set_source_rgb (cr2, fg[rf], fg[gf], .2);
	else cairo_set_source_rgb(cr2, 0.5, 0.5, 0.5);
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE);
	cairo_mask_surface(cr2, cpu_surf, 1.0, 1.0);
}

static void dabmp_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	widg_defcmd(ww, b9==9 ? "A0Y" : "#hello"); }

static void dabmp_skel(struct _ww_t * ww, const char **pp) {
	unsigned int * q = (unsigned int *) cpu_bmp;
	int i; for (i=0; i<37; i++) q[i] = ((i-2)&3) ? 0 : 0x80000100;
	q[37] = -1; q[38] = 0;
	cpu_surf = cairo_image_surface_create_for_data(cpu_bmp, CAIRO_FORMAT_A1, 32, 38, 4);
	da_skel(ww, 34, 40); }

static void dabmp_cmd(struct _ww_t * ww, const char * arg) { return dabmp_draw(ww, (cairo_t*)ww->arg[0].p); }

static void dabmp_upd(struct _ww_t * ww, const char * dat, int nf) {
	if (!ww && !(ww = widg_lookup_ps(ot_etc, "2"))) return LOG("default bitmap not found");
	unsigned int *q = (unsigned int *)(cpu_bmp);
	int k, i, j, n = nf&63; 
	for (j=0; j<n; j++) {
		if ((k=dat[j])<0) k = 37; else if (k>37) k = 0; else k = 37-k;
		for (i=0; i<k;  i++) q[i]>>=1;
		for (i=k; i<38; i++) q[i]>>=1, q[i]|=0x80000000;
		if (++q[38]==25) { q[38] = 0; for (i=2; i<37; i+=4) q[i]^=0x8000000; }
	}
	if (nf&64) da_fullre(ww);
}

static void write_tlog() {
	static const char *fn = 0; if (!fn && !(fn=getenv("LF_TLOG"))) fn = "lf.tlog";
	backup(fn, tlog_c_bk);
	LOG("write_tlog, expected file size: %d bytes", tlog_wr ? (1<<22) - 4*!!(tlog_ix_c&3) : 4*tlog_ix_i);
	int n, k, fd = creat(fn, 0644); if (fd<0) return perror(fn);
	if (tlog_wr && (n=(1<<22)-(k=(tlog_ix_c+3)&~3))) write(fd, (char*)tlog_dat + k, n);
	if (tlog_ix_i) write(fd, tlog_dat, 4*tlog_ix_i);
	close(fd);
}

static gboolean tpipe_in (GIOChannel *src, GIOCondition cond, gpointer data) {
	static int maxv = 0, cnt = 0, unexp = 0;
	if (cond!=G_IO_IN) return LOG("ioc=%d", cond), FALSE;
	int r = read( g_io_channel_unix_get_fd  (src), (char*)tlog_dat + tlog_ix_c, 4096);
	if (r<1) return LOG("tpipe errno: r=%d, %d/%s", r, errno, strerror(errno)), TRUE;
	int k, i, j=0, m=0, nxc = tlog_ix_c + r, nxi = nxc>>2;
	unsigned int * q = tlog_dat;  char bbuf[32];
	for (i=tlog_ix_i; i<nxi; i++) {
		int c = q[i]>>18; if ((c|1)=='q') cpu_austat = 'q'-c; else continue;
		if (!i || (k=q[i-1]>>18)<4097) {
			if (++unexp<1001) LOG("unexp 'p' in tlog (cnt:%d)%s", unexp, 
					unexp==1000?" -- this is the last message":"");
			if (unexp>999888777) LOG("this is impossible."), kill(getppid(), 9);
			continue; 
		}
		int v = (51 * ((q[i]&262143)+(q[i-1]&262143)) / (k-4096)) >> 5;
		if (v>maxv) maxv = v;
		if (++cnt>3) cnt = 0, bbuf[j] = min_i(maxv, 99), maxv = 0, j = j<31 ? j+1 : (++m,0);
	}
	if (j|m) { if (m) LOG("tipe_in: lots of input, m=%d", m), dabmp_upd(NULL, bbuf+j, 32-j);
		   if (j) dabmp_upd(NULL, bbuf, j+64); }
	tlog_ix_c = nxc & 0x3fffff; tlog_ix_i = nxi & 0xfffff;
	if (nxc>0x400000) memcpy((char*)q, (char*)(q+(1<<20)), tlog_ix_c), tlog_wr = 1;
	return TRUE;
}

//////// mini-wrap(18x44)  arg: nmxy12v:rgbrgb /////////////////////////

static void da_wr18(cairo_t * cr2, int xi, int yi, const char * s0, int xclip) {
	//LOG("da_wr18: x=%d, y=%d, clp=%d", xi, yi, xclip);
	if (!xclip) return;
	double x = (double)xi, y = (double)yi, y2, rgb[6];
	char txt[9]; const char * s = s0;
	if (s[2]=='0') {
		cairo_set_source_rgb (cr2, .2, .2, .2);
		cairo_rectangle(cr2, x, y, 18.0, 44.0);
		cairo_fill(cr2);
		return;
	}
	int i,j; for (i=j=0; i<3; i++) {
		if (!(txt[j++] = *(s++))) goto err;
		if (!(txt[j++] = *(s++))) goto err;
		txt[j++] = 0;
	}
	if (!*s) goto err;
	double v = 0.8 * (double)aA2i(*(s++));
	if ((*s++) != ':') goto err;
	if (xclip<18) cairo_save(cr2), cairo_rectangle(cr2, x, y, (double)xclip, 44.0), cairo_clip(cr2);
	for (i=0; i<6; i++) rgb[i] = RGB_C(*(s++));
	cairo_set_source_rgb (cr2, rgb[3], rgb[4], rgb[5]);
	cairo_rectangle(cr2, x+1.0, y+1.0, 16.0, 42.0);
	cairo_fill(cr2);
	cairo_set_source_rgb (cr2, rgb[0], rgb[1], rgb[2]);
	cairo_select_font_face (cr2, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
	cairo_set_font_size(cr2, 10.0);
	for (i=0, y2=y; i<3; i++) {
		cairo_move_to(cr2, x+2.0, y2+12.0); y2 += 13.0;
		cairo_show_text(cr2, txt+3*i);
	}
	cairo_set_line_width (cr2, 1.0);
	cairo_move_to(cr2, x+16.5, y+42.5);
	cairo_line_to(cr2, x+16.5, y+41.5-v);
	cairo_stroke(cr2);
	cairo_set_source_rgb (cr2, .5, .5, .5); 
	cairo_rectangle(cr2, x+.5, y+.5, 17.0, 43.0);
	cairo_stroke(cr2);
	if (xclip<18) cairo_restore(cr2);
	return;
err:    LOG("wr18: invalid format \"%s\"", s0); if (!*s0) abort(); // TODO: rm
}

///////////////// drawing area (text) ////////////////////////////////////////

static int tx_box(cairo_t * cr, int x, int y, int w, int h, int h0, const char * s) {
	int fs = get_fontsize(h0&255, w, s);
	if (FONT_HEIG(fs) < (h0>>8)) return 0;
	int x0 = (x&65536) ? 3 : ((w - tx_len(fs, s)) >> 1); 
	int y0 = ((h - FONT_HEIG(fs) + 1) >> 1) + FONT_OFFS(fs);
	cairo_set_font_size (cr, (double)fs);
	cairo_move_to(cr, (double)((x&65535)+x0), (double)(y+y0));
	cairo_show_text(cr, s); return 1;
}

static int tx_split(int fs, const char *s, int wid, int force) {
	int i, k, spc = wid;
	const char * wt = FONT_WTAB(fs);
	for (i=0; s[i] && spc >= wt[(int)s[i]]; i++) spc -= wt[(int)s[i]];
	k=i; spc = wid;
	for (; s[i] && spc >= wt[(int)s[i]]; i++) spc -= wt[(int)s[i]];
	if (!s[i]) return k;
	if (!force) return 0;
	spc = wid - 3*wt['.']; while (s[i]) i++;
	for (--i; i>k && spc >= wt[(int)s[i]]; i--) spc -= wt[(int)s[i]];
	return k + 256*i;
}

static void tx_box_split(cairo_t * cr, int x, int y, int w, int h, int h0, char * s) {
	if (tx_box(cr, x, y, w, h, h0 + 256*(h0-2), s)) return;
	if (2*h0+6 > h) h0 = h/2 - 3;
	int fs, k = 0, k2 = 0, df = 0;
	char sv[4];
	for (fs = get_fontsize(h0, 0, NULL); fs>11; fs--)
		if (( k2 = k = tx_split(fs, s, w, 0 ))) break;
//	LOG("split1: k=%d k2=%d", k, k2);
	if (!k) k = tx_split(fs, s, w, 1);
	if (k<256) k2 = k; else k2 = k>>8, k &= 255;
//	LOG("split2: k=%d k2=%d", k, k2);
	if (!s[k]) return (void) tx_box(cr, x, y, w, h, h0, s);
	sv[0] = s[k]; s[k] = 0;
	tx_box(cr, x+2, y+2, w-4, h/2-3, fs, s);
	s[k] = sv[0];
	if (k2>=3&&k2>k) {
		k2 -= 3; df = 1;
		memcpy(sv, s+k2, 4); memcpy(s+k2, "...", 3);
	}
	tx_box(cr, x+2, y+h/2+1, w-4, h/2-3, fs, s+k2);
	if (df) memcpy(s+k2, sv, 4);
}

static int tx_box_split_c(cairo_t * cr, int x, int y, int w, int h, int h0, const char * s) {
	int l = strlen(s) + 1; char s2[l]; memcpy(s2, s, l);
	tx_box_split(cr, x, y, w, h, h0, s2); return 1; }

#define DACLR(w,i)  RGB_C((w)->arg[2].c[i])
#define DACLR4(w,i) RGB_C((w)->arg[4].c[i])
#define DABOOL(w) ( (w)->arg[3].c[0] )
#define DALBL_TXT(x) ((x)->arg[2].c)

static const char * daerr_txt(int i, int j);
static void dalbl_draw(ww_t * ww, cairo_t * cr2) {
	double fr, fg, fb; const char *s = 0; int lal = 0;
	if (ww->cl->ch=='E') {
		int *p2 = ww->arg[2].i; lal = !!*p2; s = lal ? daerr_txt(p2[0], p2[1]) : "---";
		if (*p2&0x80000000) cairo_set_source_rgb(cr2, fr=1.0, 0.0, 0.0);
		else cairo_set_source_rgb(cr2, 0.7, 0.7, 0.7), fr = lal ? 0.0 : 0.333;
		fg = fb = fr; goto dr;
	}
	s = DALBL_TXT(ww);
	if (ww->arg[4].c[6]) {
		RGB_C3(cr2, ww->arg[4].c+3);
		fr = DACLR4(ww,0), fg = DACLR4(ww,1), fb = DACLR4(ww,2);
	} else { switch(ww->cl->ch) {
		case 'B': if (*s=='X' && !s[1]) 
				  cairo_set_source_rgb(cr2, .2, .1, .1), fr=1.0, fg=fb=.5; 
			  else    cairo_set_source_rgb(cr2, .1, .2, .2), fr=.5, fg=fb=1.0;
			  break;
		case 'L': cairo_set_source_rgb (cr2, .2, .2, .2);
			  fr = fg = fb = .9; break;
		case 'Y': if (ww->arg[3].c[0]) 
				  cairo_set_source_rgb (cr2, 0.0, 0.0, 1.0), fr=fg=1.0, fb=0.0;
			  else    cairo_set_source_rgb (cr2, 0.0, 0.0, 0.33333), fr=fg=0.5, fb=0.0;
			  break;
		case 'M': cairo_set_source_rgb (cr2, .2, .2, .1), fr = fg = 1.0, fb = .5; break;
		default:  LOG("dalbl_draw: invalid class 0x%x(%c)", ww->cl->ch, ww->cl->ch); return;
	}}
dr:	cairo_paint(cr2);
	cairo_set_source_rgb (cr2, fr, fg, fb);
	tx_box(cr2, 2+(lal<<16), 2, DA_W(ww)-4, DA_H(ww)-4, conf_lbfs, s);
	da_frame(ww, cr2);
}

static void daclb_draw(ww_t * ww, cairo_t * cr2) {
	char *q, *s = ww->etc, buf[32];
	if (!memcmp(s, "@@@", 4)) { s=buf; strcpy(s,"n/m/upuiAhnbjm/dpn"); for (q=s; *q; q++) --*q; } 
	RGB_C3(cr2, ww->arg[2].c+3); cairo_paint(cr2); RGB_C3(cr2, ww->arg[2].c);
	if (DABOOL(ww)) tx_box_split_c (cr2, 1, 1, DA_W(ww)-1, DA_H(ww)-1, conf_lbfs, s);
	else tx_box (cr2, 2, 2, DA_W(ww)-4, DA_H(ww)-4, conf_lbfs, s);
	if (s==buf) memset(s, 0, 32);
	da_frame(ww, cr2);
}

#define DLM_MT(x)  ((x)->arg[3].i[0])
#define DLM_MSK(x) (*(unsigned int *)&((x)->arg[3].i[1]))

static void dlmenu_ilb(ww_t * ww, int j) {
	int k, i = DLM_MT(ww), f = menutab[i].widf;
	char * s = DALBL_TXT(ww);
	if (f&4096) { if (j<64) { if (f&2048) --i, --DLM_MT(ww);
				  s[0]=50+(j>>3), s[1]=':', s[2]=50+(j&7), s[3]=0; 
				  return da_fullre(ww); }
		      else      { if (!(f&2048)) ++i, ++DLM_MT(ww); }}
	k = menutab[i].lbll, memcpy(s, menutab[i].lbl+(j&31)*k, k);
	return da_fullre(ww);
}

static void dalbl_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return dalbl_draw(ww, (cairo_t*)ww->arg[0].p);
	int k, ch = ww->cl->ch; char *s = DALBL_TXT(ww);
	switch(arg[0]) {
		case '+': case 'c': if (ch!='M') goto err; else return dlmenu_ilb(ww, arg[1]-48);
		case 'C': s = ww->arg[4].c; s[6] = 1; ++arg;
			  for (k=0; arg[k] && k<6; k++) s[k] = arg[k];
			  if (k<6 || *(arg+=k)!=':') goto draw;
		case 't': strncpy(DALBL_TXT(ww), arg+1, 7); ww->arg[2].c[7] = 0; goto draw;
		case 'x':
		case 's': if (ch!='Y') goto err;
			  k = arg[1]-'0';
			  if (k&~1) { LOG("dalbl: Y/s: invalid value 0x%x", arg[1]); return; }
			  if (DABOOL(ww) != k) DABOOL(ww) = k; else return;
			  goto draw;
		case 'S':
			  ++arg; get_tok(ww->arg[2].c, 8, &arg, '$'+128);
			  get_tok(ww->cmd, 16, &arg, 128); goto draw;
		case 'M':
			  if (ch!='M') goto err;
			  DLM_MSK(ww) = atoi_h(arg+1); return;
		case ':': s[1]=0; memcpy(ww->arg[4].c, (s[0]=arg[1])=='+' ? "zzz%z%\001":"zz%z%%\01", 8);
			  goto draw;
		default:
			  LOG("dalbl: invalid cmd 0x%x",arg[0]); return;
	}
draw:   da_fullre(ww); return;
err:	LOG("dalbl: invalid cmd 0x%x for class 0x%x (%c/%c)", *arg, ch, *arg, ch);
}

static void dalbl_mw18(ww_t * ww, const char * s) {
	char *p = DALBL_TXT(ww), *q = ww->arg[4].c;
	if (s) memcpy(p, s, 6), p[7] = 0, memcpy(q, s+8, 6), q[6] = 1;
	else   memcpy(p, "no sel.", 8),   memcpy(q, "qqq---\01", 7);
	da_fullre(ww); }

static void upd_flgvec(topwin * tw, const char* s, int n, int of, int nf) {
	int ix = widg_lookup(tw, &s, 0);
	if (ix<0) { LOG("upd_flgvec: widget not found"); return; }
	int i, msk=1, sg=1; if (n<0) n = -n, sg = -1;
	for (i=0; i<n; i++, ix+=sg, msk<<=1) {
		if (!(msk & (of^nf))) continue;
		ww_t * ww = widg_p(tw, ix);
		DABOOL(ww) = !!(msk&nf); da_fullre(ww);
	}
}

static void dlbtn_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	if (b9==1) widg_defcmd(ww, ""); }

static void dlyn_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	if (b9==1) widg_defcmd(ww, "1\0""0"+2*ww->arg[3].c[0]); }

static void dlmenu_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
        if ((b9|2)==3) popup2(ww, DLM_MT(ww), DLM_MSK(ww), b9, ev); }

static void dalbl_skel(struct _ww_t * ww, const char **pp) {
        int ch = ww->cl->ch, h2 = 0, wmul = 5, mt = 0, mm = 0;
	if (ch=='E') return da_skel(ww, 50*FONT_NWID(conf_lbfs), conf_lbh);
	if (**pp=='|') h2 = (conf_lbh>>2), ++*pp;
	else if (**pp=='`') wmul = 2, ++*pp;
	else if (**pp==' ') h2 = conf_lbfh_s - conf_lbfh, ++*pp;
        if (ch=='M' && **pp==36) ww->arg[2].c[0] = 0, ++*pp;
        else get_tok(ww->arg[2].c, 8, pp, '$');
        get_tok(ww->cmd, 16, pp, '|');
        if (ch=='M') {
		if (**pp == '_') wmul = 3;
		mt = menutab_lu(**pp, (*pp)[1] - 48); *pp += 2;
		if (**pp!='}') for (; 47<**pp&&**pp<103; ++*pp)
			mm = 16*mm + hxd2i(**pp);
		DLM_MT(ww) = mt; DLM_MSK(ww) = mm; ww->arg[4].c[6] = 0;
	}
        int w = menutab[mt].widf & 1023, wid = 6 + (w ? w : tx_len(conf_lbfs, ww->arg[2].c));
        if (ch!='L' && 4*wid < wmul*conf_lbh) wid = (wmul*conf_lbh)>>2;
        else if (wid==6) wid = 2;
        da_skel(ww, wid, conf_lbh + h2);
}

static void daclb_set(struct _ww_t * ww, const char **pp, int flg) {
	char buf[256]; int force = 0;
	if (**pp=='!') force=1, ++*pp;
	if (!memcmp(*pp, "==> ", 4)) {
		memcpy(ww->arg[2].c, "Azz666_", 6); ww->arg[3].i[0] = -1;
	} else if (**pp=='('&&!(flg&8)) { 
		++*pp; int x = 0; while (**pp&80) x = 16*x+hxd2i(**pp), ++*pp;
		ww->arg[3].i[0] = x; if (**pp==')') ++*pp;
	} else { ww->arg[3].i[0] = 0; }
	if (**pp==',') ++*pp; else if (!(flg&2)) memcpy(ww->arg[2].c, *pp, 6), *pp += 6;
	int l = get_tok(buf, 256, pp, ((flg&1)<<7)+36);
	if (l) ++l; else buf[0]=32, buf[1]=0, l=2;
	ww->etc = realloc(ww->etc, l); memcpy(ww->etc, buf, l);
	if (ww->arg[3].i[0]<0) {
		int l2 = 0; while (buf[l2] && memcmp(buf+l2, " -- ", 4)) ++l2; ww->arg[3].i[0] = -l2; }
	if (force) {
		if (buf[0]!='.') { LOG("daclb_set: path begins with 0x%x", buf[0]); return; }
		int j; for (j=l-1; buf[j]!='.'; j--) ;
		memcpy(ww->top->title, buf+j+1, l-j-1);
		gtk_window_set_title(GTK_WINDOW(ww->top->w), ww->top->title);
	}
	da_fullre(ww);
}

static void daclb_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	int k = ww->arg[3].i[0]; if (!k) return;
	if(k<0){char * q = (char*)ww->etc; int c = q[-k]; q[-k] = 0;
		CMD("W%s%s", (q[4]=='.')?"":".!b.?.", q + 4); q[-k] = c; return; }
	int tid = ww->top->id, f = (k==(tid>>4));
	char buf[16]; buf[0] = '~'; buf[1] = 87-10*f; buf[2] = '#';
	int l = 3 + hx5(buf+3, k);
	if (f) buf[l] = 36, buf[l+1] = hexc1(tid&15), l += 2;
	buf[l] = 10; write(1, buf, l+1);
}

static void daclb_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return daclb_draw(ww, (cairo_t*)ww->arg[0].p);
	if (arg[0]!='s' && arg[0]!='t') { LOG("daclb: invalid cmd 0x%x",arg[0]); return; }
	const char * s = arg + 1; daclb_set(ww, &s, 1);
}

static void daclb_skel(struct _ww_t * ww, const char **pp) {
	int f2 = DABOOL(ww) = (**pp=='/' && ++*pp);
	int w = get_dec(pp, 36), w2 = 0;
	int h = get_dec(pp, 36);
	if (h<2) w2 = h, h = f2 ? 2*conf_lbh-6 : conf_lbh;
	memcpy(ww->arg[2].c, "ppp666_", 8);
	daclb_set(ww, pp, 0);
	if (!w || w2) {
		w2 = tx_len(conf_lbfh, ww->etc);
		if (!w || w2<w) w = w2;
	}
	da_skel(ww, w, h);
}

///////////////// progress bar / button //////////////////////////////////////
#define DAPRG_VAL(x) ((x)->arg[3].i[0]) // 1023:def 1022:ok 1021:stop 1020:error

static void daprg_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	char buf[2]; buf[0] = b9+42; buf[1] = 0; widg_defcmd(ww, buf); }

static void daprg_skel(struct _ww_t * ww, const char **pp) {
	get_tok(DALBL_TXT(ww), 8, pp, '$');
	int wid = max_i(3*FONT_NWID(conf_lbfs) + tx_len(conf_lbfs_s, ".%"),
			tx_len(conf_lbfs_s, DALBL_TXT(ww)));
	DAPRG_VAL(ww) = 1016;
	da_skel(ww, 6 + wid, conf_lbh);
	get_tok(ww->cmd, 16, pp, 0);
}

static void daprg_draw(ww_t * ww, cairo_t * cr2) {
	int i, diag=0, v = DAPRG_VAL(ww), v2;
	char *s=DALBL_TXT(ww), buf[6];
	cairo_set_source_rgb (cr2, .1, .2, .2); cairo_paint(cr2);
	switch(v) {
		case 1016: break;
		case 1017: diag = 2; break;
		case 1018: diag = 066; break;
		case 1019: diag = 044; break;
		case 1000: *(s=buf) = '-'; s[1] = 0; break;
		case 1001: *(s=buf) = '/'; s[1] = 0; break;
		case 1002: *(s=buf) = '|'; s[1] = 0; break;
		case 1003: *(s=buf) = '\\';s[1] = 0; break;
		case 1004: *(s=buf) = '*'; s[1] = 0; break;
		default:
			if (v>1004) { diag = 044; break; }
			if (v<0) v = 333;
			s = buf;
			if (v<100) s[0]=32, v2=v; else s[0]=48+v/100, v2=v%100;
			s[1] = 48+v2/10; s[2]='.'; s[3] = 48+v2%10; s[4]='%'; s[5] = 0;
			cairo_set_source_rgb(cr2, .25, .5, .5);
			cairo_rectangle(cr2, 0.0, 0.0, .001*(double)(DA_W(ww)*v), (double)DA_H(ww));
			cairo_fill(cr2);
			break;
	}
	for (i=0; i<2; i++) {
		int d = (diag>>(3*i)) & 7; if (!d) continue;
		LOG("diag %d %d", diag, i);
		cairo_set_source_rgb(cr2, .2+.1*(double)(d&4), .2+.2*(double)(d&2), .2+.4*(double)(d&1));
		cairo_set_line_width(cr2, 3.0);
		cairo_move_to(cr2, (double) (DA_W(ww) & (i-1)), 1.0);
		cairo_line_to(cr2, (double) (DA_W(ww) &  (-i)), (double)DA_H(ww)-1.0);
		cairo_stroke(cr2);
	}
	if (v>999) cairo_set_source_rgb(cr2, .5, 1.0, 1.0);
	else cairo_set_source_rgb(cr2, 1.0, 1.0, .5);
	cairo_set_font_size(cr2, conf_lbfs);
	int x = (DA_W(ww) - tx_len(conf_lbfs, s)) >> 1;
	int y = (DA_H(ww) - FONT_HEIG(conf_lbfs) + 1) >> 1;
	cairo_move_to(cr2, (double)x, (double)(y+FONT_OFFS(conf_lbfs)));
	cairo_show_text(cr2, s);
	da_frame(ww, cr2);
}

static void daprg_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return daprg_draw(ww, (cairo_t*)ww->arg[0].p);
	if (arg[0]=='t') { strncpy(ww->arg[2].c, arg+1, 7); ww->arg[2].c[7] = 0; }
	else if (arg[0]!='s') { LOG("daclb: invalid cmd 0x%x",arg[0]); return; }
	else { DAPRG_VAL(ww) = 64*arg[1] + arg[2] - 3120; }
	da_fullre(ww);
}

///////////////// counter (+/-) + vscale /////////////////////////////////////
// arg: f:zpad ff0:lbwid 1000:vert/n 2000:vert/aA 4000:mvol+- 8000:cont ff0000:x8 ff000000:x8max
#define DACNT_ARG(x) ((x)->arg[2].i[1])
#define DACNT_CONT(x) ((x)->arg[2].i[1]&32768)
#define DACNT_X(x) ( DACNT_CONT(x) ? (DACNT_ARG(x)>>16)&255 : (x)->arg[2].i[0] )
#define DACNT_LIM(x) ( (DACNT_ARG(x)>>24) & 255 )
#define DACNT_LBL(x) ( DACNT_CONT(x) ? (x)->arg[2].c : (x)->arg[3].c ) 
#define DACNT_VNFS(z) ((z)<2 ? conf_lbfs : conf_lbfs_s - ((z)>2))
#define DACNT_SLIDER(x) (GTK_RANGE((x)->arg[3].p))
#define DACNT_DOT(x) ((x)->arg[4].c[0])
#define DACNT_DOTC(x) ( DACNT_CONT(x) ? 0 : (x)->arg[4].c[0] )

static void dacnt_act(ww_t * ww, int x) {
	char buf[16], *s = ww->cmd; int j;
	if (*s!=36) { buf[hx5(buf, x)] = 0; widg_defcmd(ww, buf); return; }
	switch(s[1]) {
		case 'c': ww->top->arg[s[2]-48].c[s[3]&7] = x; return;
		case 'C': for (j=s[2]&7, s+=3; *s;) 
				  (ww=widg_lookup_p(ww->top, (const char**)(&s)))->arg[4].c[j] = x + 37, 
					  da_fullre(ww); 
			  return;
		default:  LOG("dacnt_act: unknown esc 0x%x(%c)", s[1], s[1]); return;
	}}

static void dacnt_set_x(ww_t * ww, int x, int flg) { // 255:lim 256:da_fullre 512:sl_upd 1024:cmd 2048:cond_r
	if (flg&1024) dacnt_act(ww, x);
	if (flg& 255) { DACNT_ARG(ww) &= 0xffffff, DACNT_ARG(ww) |= (flg&255)<<24;
			if (x==-1234567890) x = min_i(DACNT_X(ww), flg&255); }
	if (DACNT_CONT(ww)) {
		DACNT_ARG(ww) &= 0xff00ffff; DACNT_ARG(ww) |= (x&255)<<16;
		if (flg&512) gtk_range_set_value (DACNT_SLIDER(ww), (double)x / (double)DACNT_LIM(ww));
	} else if (ww->arg[2].i[0]!=x) {
			ww->arg[2].i[0] = x; flg |= (flg>>3);
	}
   	if (flg&256) da_fullre(ww);
}

static const char * fmts[] = {"%d","%01d","%02d","%03d","%04d","%05d","%06d","%07d",
                              "%X","%01X","%02X","%03X","%04X","%05X","%06X","%07X"};

static void dacnt_draw(ww_t * ww, cairo_t * cr2) {
	cairo_set_source_rgb (cr2, .1, .2, .1);
	cairo_paint (cr2);
	char buf[16]; buf[0] = 0;
	int zpad = DACNT_ARG(ww) & 15;
	int lwid = (DACNT_ARG(ww)>>4) & 255;
	int vmode = (DACNT_ARG(ww) >> 12) & 15;
	int nflg = -(~(vmode>>3) | (vmode&1));
	int val = DACNT_X(ww);
	if (vmode&4) {
		if (val==0)   buf[0]='0', buf[1] = 0;
		else if (val==100) buf[0]='1', buf[1] = 0;
		else buf[0]='.', buf[1]=48+val/10, buf[2]=48+val%10, buf[3] = 0;
	} else if (nflg) {
		int k = sprintf(buf, fmts[zpad], val);
		if (DACNT_DOTC(ww)==1) buf[k]=buf[k-1], buf[k-1]='.', buf[k+1]=0;
		else if (DACNT_DOTC(ww)) buf[k]=buf[k-1], buf[k-1]=buf[k-2], buf[k-2]='.', buf[k+1]=0;
		else buf[k]=0;
	}
	buf[14] = i2aA(val), buf[15] = 0;
	int xl, xn, xc=0, yl, yn, yc=0, nfs = conf_lbfs;
	if (vmode) {
		nfs = DACNT_VNFS(zpad);
		int hn = (FONT_HEIG(nfs)+1) &- (vmode&1);
		int hc = (FONT_HEIG(conf_lbfs)+1) &- ((vmode>>1)&1);
		xl = 3; 
		xn = (DA_W(ww) - tx_len(nfs, buf)) >> 1,
		xc = (DA_W(ww) - tx_len(conf_lbfs, buf+14)) >> 1;
		yl = (DA_H(ww) - FONT_HEIG(conf_lbfs) - hn - hc + 1) >> 1;
		yn = yl + 1 + FONT_HEIG(conf_lbfs_s);
		yc = yn + hn;
	} else {
		int cwid = lwid + 1 + tx_len(nfs, buf);
		xl = (DA_W(ww) - cwid) >> 1; if (xl>3) xl = 3;
		xn = xl + lwid + 1;
		yn = (DA_H(ww) - FONT_HEIG(conf_lbfs) + 1) >> 1;
		yl = yn + FONT_HEIG(conf_lbfs) - FONT_HEIG(conf_lbfs_s);
	}
	cairo_set_source_rgb (cr2, .8, .8, .8);
	cairo_set_font_size(cr2, (double)conf_lbfs_s);
	cairo_move_to(cr2, (double)xl, (double)(yl+FONT_OFFS(conf_lbfs_s)));
	cairo_show_text(cr2, DACNT_LBL(ww));
	cairo_set_source_rgb (cr2, .0, 1.0, .0);
	if (nflg) {
		cairo_set_font_size(cr2, (double)nfs);
		cairo_move_to(cr2, (double)(xn), (double)(yn+FONT_OFFS(nfs)));
		cairo_show_text(cr2, buf);
	}
	if (vmode&2) {
		cairo_set_font_size(cr2, conf_lbfs);
		cairo_move_to(cr2, (double)xc, (double)(yc+FONT_OFFS(conf_lbfs)));
		cairo_show_text(cr2, buf+14);
	}
	da_frame(ww, cr2);
}

static void dacnt_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return dacnt_draw(ww, (cairo_t*)ww->arg[0].p);
	int k,x;
	switch(*arg) { 
		case 't': 
			++arg; get_tok(DACNT_LBL(ww), 8, (const char**)&arg, 128); 
			da_fullre(ww); return;
		case 'c': case 'n': dacnt_set_x(ww, arg[0]+arg[1]-147, 768); return;
		case 'x':
			  x = 0; while ((k=*++arg)&80) x = 16*x + hxd2i(k);
			  dacnt_set_x(ww, (*arg=='-') ? -x : x, 768); return;
		case 'C': x = arg[1]-48; k = arg[2]-48; arg += 3; break;
		case 'X': x = hex2(arg+1); k = hex2(arg+3); arg += 5; break;
		default: LOG("dacnt: unknown cmd 0x%x (%c)", *arg, *arg); return;
	}
	if (*arg) get_tok(DACNT_LBL(ww), 8, (const char**)&arg, 128);
	dacnt_set_x(ww, x, 768+k);
}

#define	CNT_CLK int sh = 3*b9, btn = (03213213210>>sh)&7, mod = (02221110000>>sh)&7
static void dacnt_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	CNT_CLK; char buf[4]; buf[0] = 42+btn; buf[1] = 48+mod; buf[2] = 0;
	widg_defcmd(ww, buf);
	return;
}

static void dacntvs_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	CNT_CLK, sg = 2 - btn; if (!sg) return;
	int x0 = DACNT_X(ww), x = x0, k = DACNT_LIM(ww);
	switch(mod) { case 0 : x +=   sg; break;
		      case 1 : x += 8*sg; break;
		      default: x = sg<0 ? 0 : k; break; }
	if (x<0) x = 0; else if (x>k) x = k;
	if (x==x0) return;
	dacnt_set_x(ww, x, 0x700);
	return;
}

gboolean scale_chg(GtkRange * w, gpointer p) {
	ww_t * ww = (ww_t*) p;
	int k = DACNT_LIM(ww);
	double v0 = gtk_range_get_value(w);
	double vx = v0*(double)k; int x = (int)lround(vx);
	double v = round(vx) / (double)k;
	if (fabs(v-v0) > 1e-4) gtk_range_set_value(GTK_RANGE(w), v);
	if (x==DACNT_X(ww)) return TRUE;
	dacnt_set_x(ww, x, 0x500);
	return TRUE;
}

static void scale_skel(struct _ww_t * ww, const char **pp) {
	int k = 17, ch = **pp; if (ch>49&&ch<124) k = ch-48, ++*pp;
	if (k==73) k=101; else if (k==74) k = 256;
	LOG("scale: step=%g chr=0x%x(%c)", 1.0/(double)(k-1), ch, ch);
      	ww -> w = gtk_vscale_new_with_range(0.0, 1.0, 1.0/(double)(k-1));
	gtk_scale_set_draw_value(GTK_SCALE(ww->w), FALSE);
//	get_tok(ww->cmd, 16, pp, 0);
  	g_signal_connect (ww->w, "value-changed", G_CALLBACK (scale_chg), (gpointer)ww);
}

static void dacnt_skel(struct _ww_t * ww, const char **pp) {
	int cont = (ww->cl->ch=='!'); DACNT_ARG(ww) = 0x8000*cont;
	get_tok(DACNT_LBL(ww), 8, pp, '$');
	int zf = 0, wid = 0, lbwid, nlen;
	switch(**pp) {
		case '=': zf = 1;
		case '-': wid = font_w[96 * (conf_lbfs - 6) + '-' - 32]; ++*pp; break;
		case '0': zf = 1; ++*pp; break;
		case '.': zf = 1; DACNT_DOT(ww) = 1; wid = font_w[96*(conf_lbfs - 6)+'.'-32]; ++*pp; break;
		case ':': zf = 1; DACNT_DOT(ww) = 2; wid = font_w[96*(conf_lbfs - 6)+'.'-32]; ++*pp; break;
		default: break;
	}
	wid += FONT_NWID(conf_lbfs) * ((nlen = hxd2i(**pp))&7); ++*pp;
	dacnt_set_x(ww, (nlen>7) ? (int)(0xabcdef01u>>(64-4*nlen)) : 1<<(3*nlen), 0);
	lbwid = tx_len(conf_lbfs_s, DACNT_LBL(ww));
	DACNT_ARG(ww) |= 16 * lbwid + (15 & nlen & -zf);
	da_skel(ww, 9 + wid + lbwid, conf_lbh);
	get_tok(ww->cmd, 16, pp, '|');
}
	
static void dacntvs_skel(struct _ww_t * ww, const char **pp) {
	DACNT_ARG(ww) = 0x8000;
	get_tok(DACNT_LBL(ww), 4, pp, '$');
	int ty = **pp & 7, k = hex2(*pp+1); *pp += 3;
	get_tok(ww->cmd, 16, pp, 0);
	int zpad = (ty&4) ? 2 : (1 + (k>9) + (k>99));
	int nfs = (ty&1) ? DACNT_VNFS(zpad) : 0;
	int wid2, wid = tx_len(conf_lbfs_s, DACNT_LBL(ww));
	int heig = FONT_HEIG(conf_lbfs_s) + 4;
	(ty&1) && (heig += 1+FONT_HEIG(nfs)) && (wid2 = FONT_NWID(nfs)*zpad) > wid && (wid = wid2);
	(ty&2) && (heig += 1+FONT_HEIG(conf_lbfs)) && (wid2=tx_len(conf_lbfs,"W"))>wid && (wid = wid2);
	DACNT_ARG(ww) |= zpad | (wid<<4) | (ty<<12) | (k<<24);
	da_skel(ww, wid+4+(ty&4), heig);
	double step = 1.0 / (double) k;
	//LOG("scale: ty:%c step=%g w:%d h:%d", 48+ty, step, wid, heig);
      	GtkWidget * scl = gtk_vscale_new_with_range(0.0, 1.0, step);
	gtk_range_set_inverted (GTK_RANGE(scl), TRUE);
	gtk_scale_set_draw_value(GTK_SCALE(scl), FALSE);
  	g_signal_connect (scl, "value-changed", G_CALLBACK (scale_chg), (gpointer)ww);
	GtkWidget * bx = gtk_vbox_new(0,0);
	gtk_box_pack_start(GTK_BOX(bx),scl,1,1,0);
	gtk_box_pack_start(GTK_BOX(bx),ww->w,0,0,0);
	gtk_widget_show(scl); gtk_widget_show(bx);
	ww->arg[3].p = scl; ww->arg[4].p = bx;
}

///////////////// wrap/main //////////////////////////////////////////////////

#define WRAP_TABIX(x) ((x)->arg[100].i[0])
#define WRAP_TABCONT(x, i) GTK_WIDGET((x)->arg[101+(i)].p)

static const char * err_box = "error!v:zz%z%\0\0";

static void dawr1_draw(ww_t * ww, cairo_t * cr2) {
	cairo_set_source_rgb (cr2, .2, .2, .2);
	cairo_paint(cr2);
	da_frame(ww, cr2);
	char buf[16]; memcpy(buf, ww->arg[3].c, 8); memcpy(buf+8, ww->arg[4].c, 8);
	da_wr18(cr2, (DA_W(ww)-18)>>1, (DA_H(ww)-44)>>1, buf, 18);
}

static void dawr1_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	int x = 0x305058 + (b9<<16); widg_defcmd(ww, (const char*)&x); }

static void dawr1_skel(struct _ww_t * ww, const char **pp) {
	memcpy(ww->arg[3].c, err_box, 8); memcpy(ww->arg[4].c, err_box+8, 8);
	da_skel(ww, 40, 60); }

static void dawr1_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return dawr1_draw(ww, (cairo_t*)ww->arg[0].p);
	if (*arg!='s') { LOG("dawr1: invalid cmd"); return; }
	memcpy(ww->arg[3].c, arg+1, 8); memcpy(ww->arg[4].c, arg+9, 6); ww->arg[4].c[6] = 0;
	da_fullre(ww);
}

GtkWidget * wrap_vbl_i (struct _ww_t * ww, int ix) {
	static const char * head[3] = {
		"({M0[$Xb|_03}{M1S$Xb|_013}{M2]$Xb|_023}3{C3280$1$eeeeee333333(...src...)}0{B4<>$XW5})",
		"({M0[$X>|_03}{M1F$X>|_053}{M2]$X>|_023}3{C3280$1$eeeeee333333(...filter...)}0{B4<>$XW6})",
		"(3{C1280$1$ttt666vol./envlp./LR}0{M2$XU|W3}{M3$XD|W4}{B4<>$XW4})"};
	static const char s0[] = "({L0QW00}3{e19$Xi[]=}0"
	"{M2$Xi[]s|W0}{M3$Xi[]f|W1}3{e49$Xi[]<}0"
	"{L7...}3{e59$Xi[]>}0{M6.$Xi[]F|W2})";
	topwin * tw = ww->top;
	if (!ix) return parse_w_s(tw, head[VB_ARG(ww)&3]);
	GtkWidget * rw = parse_w_s(tw, s0);
	int i, j, k = VB_ARG(ww);
	int ix2 = VB_WBASE(ww) + 8*ix;
	if (k&2) j = ix + (ix>6 ? 89 : 31),
		memcpy(DALBL_TXT(widg_p(tw, ix2)), "vol1\0vol2\0up\0  wait\0dn/l\0LR\0  Fwt\0 Flim"+5*ix-5, 5);
	else j = ix - 1 + (65 &- (k&1));
	char h2[2]; h2[0] = hexc1(j>>4), h2[1] = hexc1(j&15);
	for (i=1; i<7; i++) memcpy(widg_p(tw, ix2+i)->cmd+2, h2, 2);
	return rw;
}

GtkWidget * wrap_vbl_t (struct _ww_t * ww, int ix) { static const char * str[4] = { 
"({80x|$2X#d0}{M8$X#R0|W6}{81y|$2X#d1}{M9$X#R1|W5}{821|$2X#d2}{832|$2X#d3}{843|$2X#d4}{854|$2X#d5}{865|$2X#d6}{876|$2X#d7})",
"({Y0stereo$XA2}{Y1from-to$XA-}3{B2write$XAW})", "({Y0mute$XT_}{L1t:}3{e28$XTT})",
"({M7+$XV|?1}{80x|$1XVx}{81y|$1XVy}{821|$1XVz}{832|$1XVw}{85t|$:3XVt}{84r|$1XVr}3{%6calc$J1})"};
return parse_w_s(ww->top, str[ix]); }

static void wrap_cmd (struct _topwin * tw, char * arg) {
	int i, ix = hex2(arg), flg = hex2(arg+2), sh = (ix>>2)&24, j = (ix+(0x7000101>>sh)) & 31;
	ww_t *p, *ww = widg_lookup_pci(tw, (0x455a4553 >> sh) & 127, -1);
	int wi = VB_WBASE(ww) + j * VB_WPL(ww);
	const char * s = arg + 4;
	for (i=0; i<4; i++) if (flg&(1<<i))
		get_tok((p = widg_p(tw, wi+((0x6320>>(4*i))&7)))->arg[2].c, 8, &s, 164), da_fullre(p);
	for (i=0; i<3; i++) if (flg&(32<<i))
		entry_set(widg_p(tw,wi+((0x541>>(4*i))&7)), hxdoub_str(NULL, s, 15)), s += 16;
}

static void wrap_skel (struct _topwin * tw, char * arg) {
	const char * str = 
	"[" TW_TOPH
	"([{Mm$Xm|#1}{B_stp$X.1}{B_kill$X.2}]{1w}"
	"3[3{+#XG.}0{__}]0"
	"3[({L_}[(3{et6$Xt1}0{L_...}3{eT6$Xt2}0{L_s}{B_plot(t)$XPT})"
	        "(3{ef6$Xt4}0{L_...}3{eF6$Xt8}0{L_Hz}{B_plot(F)$XPF})])"
	        "(3()0{YG[#]$XWt0}{YWwav$XWt1}{Y:tlim$XWt2}{YAa.v$XWt3})])"
		"{:YW5:0}{:Ew982}{:SwO80}{:ZwN81}"
		"]";
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, str);
	const char *s = arg;
	if (*s) daclb_set(widg_lookup_ps(tw, "."), &s, 1);
	if (tw->state) return;
	if (tw->id & 1) tw->ix4 = 101, WRAP_TABIX(tw) = 0;
}

///////////////// wrap/grid //////////////////////////////////////////////////

#define GRID_SX(x)  ((x)->arg[2].c[0])
#define GRID_SY(x)  ((x)->arg[2].c[1])
#define GRID_X(x)   ((x)->arg[2].c[2])
#define GRID_Y(x)   ((x)->arg[2].c[3])
#define GRID_XRY(x) ((x)->arg[2].c[4])
#define GRID_YRY(x) ((x)->arg[2].c[5])
#define GRID_CX(x)  ((x)->arg[2].c[6])
#define GRID_CY(x)  ((x)->arg[2].c[7])
#define GRID_X0(x) ( (DA_W(x)-GRID_SX(x)*GRID_X(x)-1)>>1 )
#define GRID_Y0(x) ( (DA_H(x)-GRID_SY(x)*GRID_Y(x)-1)>>1 )
#define GRID_UX(x)  ((x)->arg[3].c[6])
#define GRID_UY(x)  ((x)->arg[3].c[7])
#define GRID_EX(x)  ((x)->arg[3].c[2])
#define GRID_EY(x)  ((x)->arg[3].c[3])
#define GRID_KX(x)  ((x)->arg[3].c[4])
#define GRID_KY(x)  ((x)->arg[3].c[5])
#define GRID_PM(x)  ((x)->arg[3].c[1])
#define GRID_WH1(W) int wid1 = GRID_SX(W) = (DA_W(W)-5) / GRID_X(W), \
			hgt1 = GRID_SY(W) = (DA_H(W)-5) / GRID_Y(W); \
		    if (wid1<1 || hgt1<1) return; \
		    double w1d = (double)wid1, h1d = (double)hgt1; \
		    int x00 = GRID_X0(ww), y00 = GRID_Y0(ww)
#define GRID_KEY(W,X,Y) (((gr_keytab_t*)(W)->etc)->xy2k[51*(Y)+(X)])
		    
typedef struct { unsigned char xy2k[2604]; unsigned short k2xy[128]; } gr_keytab_t;
static const char * grid_defaults = "\x19\x11\x0a\x0a\034\034\0\0";
static const char rawkeys[] = "?????????E1234567890-=BTqwertyuiop[]RCasdfghjkl;'`L\\zxcvbnm,./R"
			      "*A_UFFFFFFFFFF??????????????????????????????????????????????????";

static int grid_lpos(int * to, int z0, int z1, int ofs, int step, int skip) {
	int z, n = 0; 	z0 -= ofs, z1 -= ofs;
	for (z=((z0+step-1)/step)*step; z<=z1; z+=step) if (!skip || z%skip) to[n++] = z+ofs;
	return n; }

static void grid_vl(cairo_t * cr2, const int * px, int n, int y0, int y1) {
	if (!n) return; int i; double d[2]; d[0] = (double)(y1-y0+1); d[1] = -d[0]; 
	cairo_move_to(cr2, (double)(*px)+.5, (double)y0); cairo_rel_line_to(cr2, 0.0, d[0]);
	for (i=1; i<n; i++) cairo_rel_move_to(cr2, (double)(px[i]-px[i-1]), 0.0),
			    cairo_rel_line_to(cr2, 0.0, d[i&1]);
	cairo_stroke(cr2); }

static void grid_hl(cairo_t * cr2, const int * py, int n, int x0, int x1) {
	if (!n) return; int i; double d[2]; d[0] = (double)(x1-x0+1); d[1] = -d[0]; 
	cairo_move_to(cr2, (double)x0, (double)(*py)+.5); cairo_rel_line_to(cr2, d[0], 0.0);
	for (i=1; i<n; i++) cairo_rel_move_to(cr2, 0.0, (double)(py[i]-py[i-1])),
			    cairo_rel_line_to(cr2, d[i&1], 0.0);
	cairo_stroke(cr2); }

#define CONDK1(J) if ((cx=ww->arg[J].c[6])>=x0s && cx<=x1s && (cy=ww->arg[J].c[7])>=y0s && cy<=y1s) \
			cx*=wid1, cx+=x00, cy*=hgt1, cy+=y00, 
#define GRSQ_LOOP(M, X) for (cy=y00+hgt1*(i=y0s),j0=x00+wid1*x0s; i<=y1s; i++, cy+=hgt1) { \
			        unsigned char * q = &GRID_KEY(ww, x0s, i); \
			        for (cx=j0, j=x0s; j<=x1s; j++, cx+=wid1, q++) if (*q&M) (X); }
static void dagrid_draw(ww_t * ww, cairo_t * cr2) {
//	LOG("gr-draw: %d,%d - %d,%d", ww->arg[4].s[0], ww->arg[4].s[1], ww->arg[4].s[2], ww->arg[4].s[3]);
	if (!GRID_X(ww) || !GRID_Y(ww)) return (void)  LOG("BUG: grid %dx%d", GRID_X(ww), GRID_Y(ww)); 
	GRID_WH1(ww); short * a4 = ww->arg[4].s;
	int x0 = max_i(a4[0], x00), x1 = min_i(a4[2], x00+GRID_X(ww)*wid1),
	    y0 = max_i(a4[1], y00), y1 = min_i(a4[3], y00+GRID_Y(ww)*hgt1),
	    x0s = (x0-x00) / wid1, x1s = (x1-x00-1) / wid1,
	    y0s = (y0-y00) / hgt1, y1s = (y1-y00-1) / hgt1,
	    y_r = (GRID_YRY(ww)&7)+2, y_y = y_r * ((GRID_YRY(ww)>>3)+2), lxy[52],
	    cx, cy, i, j, j0, x_ry = GRID_XRY(ww), x_r = 0, x_y = 0, pm;
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE); cairo_set_line_width(cr2, 1.0);
	if (x_ry<64) { x_r = (x_ry&7)+2, x_y = x_r*((x_ry>>3)+2), pm = 0; }
	else {  pm = (x_ry&63)+1; cairo_set_source_rgb(cr2,.333,.333,.333);
		for (i=x0s; i<=x1s; i++) if ((1<<((i-1+pm)%12))&0xab5)
			cairo_rectangle(cr2, x00+i*wid1, y0, wid1, y1-y0+1), cairo_fill(cr2);
	}
	cairo_set_source_rgb(cr2,0.1,0.2,1.0);
	GRSQ_LOOP(128, (cairo_rectangle(cr2, cx, cy, wid1, hgt1), cairo_fill(cr2)))
	int ex = GRID_EX(ww), ey = GRID_EY(ww);
	if (ex<=x1s && ex>=x0s && ey<=y1s && ey>=y0s) 
		cairo_set_source_rgb(cr2,0.5, 0.0, .8), 
		cairo_rectangle(cr2, x00+wid1*ex, y00+hgt1*ey, wid1, hgt1), cairo_fill(cr2);
	if (wid1>5 && hgt1>6) { char b[2]; b[1] = 0; int h0 = min_i(conf_lbfs, hgt1-2);
		cairo_set_source_rgb(cr2, 0.2, 1.0, 0.2);
		GRSQ_LOOP(127, (*b=rawkeys[*q&127], tx_box(cr2, cx, cy, wid1, hgt1, h0, b))) }
	CONDK1(2) cairo_set_source_rgb(cr2, 0.3, 1.0, 1.0),
	          cairo_move_to(cr2, (double)cx+.5, (double)cy+.5), cairo_rel_line_to(cr2, w1d, h1d),
		  cairo_rel_move_to(cr2, -w1d, 0.0), cairo_rel_line_to(cr2, w1d, -h1d), cairo_stroke(cr2);
	CONDK1(3) cairo_set_source_rgb(cr2, 0.3, 1.0, 1.0),
		  cairo_move_to(cr2, (double)cx+.5, (double)cy+.5*h1d+.5), cairo_rel_line_to(cr2,w1d,0.0),
		  cairo_rel_move_to(cr2,-.5*w1d,-.5*h1d), cairo_rel_line_to(cr2,0.0,h1d), cairo_stroke(cr2);
	cairo_set_source_rgb(cr2, .5, .5, .5);
	if (!pm) grid_vl(cr2, lxy, grid_lpos(lxy, x0, x1, x00, wid1, x_r*wid1), y0, y1);
	grid_hl(cr2, lxy, grid_lpos(lxy, y0, y1, y00, hgt1, y_r*hgt1), x0, x1);
	cairo_set_source_rgb(cr2, 1.0, .0, .0);
	if (pm) grid_vl(cr2, lxy, grid_lpos(lxy, x0, x1, x00+wid1*((18-pm)%12),  12*wid1, 0), y0, y1);
	else    grid_vl(cr2, lxy, grid_lpos(lxy, x0, x1, x00, x_r*wid1, x_y*wid1), y0, y1);
	grid_hl(cr2, lxy, grid_lpos(lxy, y0, y1, y00, y_r*hgt1, y_y*hgt1), x0, x1);
	cairo_set_source_rgb(cr2, 1.0, 1.0, .0);
	if (pm) grid_vl(cr2, lxy, grid_lpos(lxy, x0, x1, x00+wid1*((13-pm)%12), 12*wid1, 0), y0, y1);
	else	grid_vl(cr2, lxy, grid_lpos(lxy, x0, x1, x00, x_y*wid1, 0), y0, y1);
	grid_hl(cr2, lxy, grid_lpos(lxy, y0, y1, y00, y_y*hgt1, 0), x0, x1);
}

static void dagrid_inv(ww_t * ww, int x, int y) {
	if (!ww->w) return;
	GRID_WH1(ww);
	struct _GdkRectangle rct;
	rct.x = x00 + wid1*x+1; rct.width  = wid1 - 1;
	rct.y = y00 + hgt1*y+1; rct.height = hgt1 - 1;
	gdk_window_invalidate_rect(gtk_widget_get_window(ww->w), &rct, 0);
}

static void dagrid_set_bu_xy(ww_t * ww, char * q, int x, int y) {
	if (x<0) x = 127, y = 0;
	if (x==q[0] && y==q[1]) return; 
	if ( q[0]   <127) dagrid_inv(ww, q[0], q[1]);
	if ((q[0]=x)<127) dagrid_inv(ww, q[0], q[1]=y);
}

static void dagrid_set_k7(ww_t * ww, int x, int y, int v) {
	if ((GRID_KEY(ww, x, y)&128) != v) GRID_KEY(ww, x, y) ^= 128, dagrid_inv(ww, x, y); }

static void dagrid_cmd(struct _ww_t * ww, const char * s) {
	if (!s) return dagrid_draw(ww, (cairo_t*)ww->arg[0].p);
	int i, j, k;
	while(1) switch(*s) {
		case 0: return;
		case '~': free3k(ww->etc); return;
		case 'g':
			for (i=0; i<6; i++) if(!(k=s[i+1])) return; else if ((k-=48)>=0) ww->arg[2].c[i+2] = k;
			da_fullre(ww); s+=7; continue;
		case 'b': case 'u':
			if (!s[1] || !s[2]) return; 
			dagrid_set_bu_xy(ww, ww->arg[2+(*s&1)].c+6, s[1]-48, s[2]-48); s+=3; continue;
		case 'k': case 'K':
			if (!s[1] || !s[2]) return;
			dagrid_set_k7(ww, s[1]-48, s[2]-48, 4*(*s&32)); s+=3; continue;
		case 'L':
			for (++s; (i=s[0]-48)>0 && (k=s[1]-48)>0; s+=2) GRID_KEY(ww, k, i) |= 128;
			s += (*s=='.'); continue;
		case 'B':
			memset(ww->etc, 0, 3072);
			for (++s; s[0]>47 && s[1]>47 && s[2]>47 && s[3]>47; s+=4)
				i = hex2(s) & 127, k = 64*s[2]+s[3] - 3120, j = 51*s[2]+s[3]-2547,
				((gr_keytab_t*)ww->etc)->xy2k[j] = i, ((gr_keytab_t*)ww->etc)->k2xy[i] = k;
			s += (*s=='.'); continue;
		case 'P': GRID_PM(ww) = hxd2i(s[1]); da_fullre(ww); s+=2; continue;
		default: LOG("dagrid: invalid cmd 0x%x(%c)", *s, *s); return;
	}}

#define GR_MKXY int x=(cx-GRID_X0(ww)-1)/GRID_SX(ww); if(x<0)x=0; else if(x>=GRID_X(ww)) x=GRID_X(ww)-1;\
		int y=(cy-GRID_Y0(ww)-1)/GRID_SY(ww); if(y<0)y=0; else if(y>=GRID_Y(ww)) y=GRID_Y(ww)-1

static void dagrid_kcmd(ww_t * ww, int b9, int x, int y) {
	if ((dflg & DF_REC) && !(b9&32)) rec_vpf(ww, "k%05x", (b9<<16)|(y<<8)|x);
	char buf[12], *q = ww->top->arg[9].c; buf[0] = 48|b9; buf[1] = 48+x; buf[2] = 48+y;
	int i; for (i=0; i<6; i++) buf[3+i] = q[i]+48; buf[9] = 0;
	widg_defcmd(ww, buf);
}

static void grid_setkey(struct _ww_t * ww, int x, int y, int kc) {
	char buf[64]; buf[0]='~'; buf[1]='X'; 
	int l = 2 + tw_defcmd(ww->top,  buf+2); buf[l++] = 'k';
	gr_keytab_t * kt = (gr_keytab_t*)(ww->etc);
	int k64 = 64*(y+1)+x, k51 = 51*y + x, oldkc = kt->xy2k[k51], oldxy = kt->k2xy[kc];
	if (kc==oldkc) return;
	if (oldkc) kt->k2xy[oldkc] = 0, buf[l] = 48+(oldkc>>4), buf[l+1] = hexc1(oldkc&15),
		buf[l+2] = buf[l+3] = 48, l += 4;
	if (oldxy) { int oldy = (oldxy>>6)-1, oldx = oldxy & 63;
		kt->xy2k[51*oldy+oldx] = 0; dagrid_inv(ww, oldx, oldy); }
	kt->xy2k[k51] = kc; kt->k2xy[kc] = k64; dagrid_inv(ww, x, y);
	if (kc && x<GRID_X(ww)-1) dagrid_inv(ww, ++GRID_EX(ww), y);
	buf[l] = 48+(kc>>4), buf[l+1] = hexc1(kc&15), buf[l+2] = y+49, buf[l+3] = x+48, l+=4;
	buf[l] = 10; write(1, buf, l+1);
}

static void dagrid_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	if (b9<0) {
		if (GRID_EX(ww)==99) {
			if (cy==GDK_KEY_Escape) { if (b9==-1) GRID_EX(ww)=GRID_EY(ww)=0, dagrid_inv(ww,0,0); 
						  return; }
			int i; gr_keytab_t * kt = (gr_keytab_t*)(ww->etc);
			if ((i = kt->k2xy[cx])) dagrid_kcmd(ww, b9&1, i&63, (i>>6)-1);
		} else if (b9==-1) {
			if (cy==GDK_KEY_Escape) dagrid_inv(ww,GRID_EX(ww),GRID_EY(ww)), GRID_EX(ww)=99;
			else grid_setkey(ww,GRID_EX(ww),GRID_EY(ww),cx);
		}
		return;
	}
	GR_MKXY;
	if (ww->cl->flg & WF_KEYEV) gtk_window_set_focus(GTK_WINDOW(ww->top->w), ww->w);
	if (GRID_EX(ww)!=99) { switch(b9) {
		case 1: dagrid_inv(ww, GRID_EX(ww)  , GRID_EY(ww)  );
			dagrid_inv(ww, GRID_EX(ww)=x, GRID_EY(ww)=y); return;
		case 2: grid_setkey(ww, x, y, 0); return;
		default: LOG("edit keys: L-clk:select M-clk:remove ESC:leave edit mode"); return;
	}}
	if (b9==1) GRID_KX(ww) = x, GRID_KY(ww) = y;
	dagrid_kcmd(ww, b9, x, y);
}

static gboolean dagrid_rel(GtkWidget *w, GdkEventButton * ev, gpointer p) {
	ww_t * ww = (ww_t*)p;
	if (GRID_EX(ww)==99 && ev->button == 1 && GRID_KX(ww)!=99)
		dagrid_kcmd(ww, 0, GRID_KX(ww), GRID_KY(ww)), GRID_KX(ww) = 99;
	return TRUE;
}

static gboolean dagrid_mv(GtkWidget *w, GdkEventMotion * ev, gpointer p) {
	ww_t * ww = (ww_t*)p; if (GRID_KX(ww)==99) return TRUE;
	int cx = (int)lround(ev->x), cy = (int)lround(ev->y); GR_MKXY;
	if (GRID_KX(ww)!=x || GRID_KY(ww)!=y)
		dagrid_kcmd(ww, 0, GRID_KX(ww),     GRID_KY(ww)),
		dagrid_kcmd(ww, 1, GRID_KX(ww) = x, GRID_KY(ww) = y);
	return TRUE;
}

static void dagrid_skel(struct _ww_t * ww, const char **pp) {
	memcpy(ww->arg[2].c, grid_defaults, 8);
	get_tok(ww->cmd, 16, pp, 0);
	ww->etc = alloc3k(); memset(ww->etc, 0, 3072);
	da_skel(ww, 57, 57); GRID_EX(ww) = GRID_KX(ww) = 99;
	g_signal_connect (ww->w, "button-release-event", G_CALLBACK(dagrid_rel), (gpointer)ww);
	g_signal_connect (ww->w, "motion-notify-event", G_CALLBACK(dagrid_mv), (gpointer)ww);
}

static void tgrid_cmd (struct _topwin * tw, char * arg) {
	int i, j, k, kd; ww_t * ww; char * q;
	while(*arg) { switch(*arg) {
		case '+':
			if (!arg[1]||!arg[2]||!arg[3]) { LOG("tgrid/+: missing flg(hx3)"); return; }
			k = 256*hxd2i(arg[1]) + hex2(arg+2);
			kd = k & ~tw->arg[8].i[0]; tw->arg[8].i[0] = k;
			for (i=0; i<12; i++) {
				*arg = 65 + (i>>1) + 32 * (i&1);
				ww = widg_lookup_ps(tw, arg);
				GtkWidget * w = GTK_WIDGET(ww->arg[4].p);
				if (!(k & (1<<i))) gtk_widget_hide(w);
				else if (gtk_widget_show(w), kd & (1<<i))
					dacnt_set_x(ww, tw->arg[9].c[i>>1], 768);
			}
			(*((k&0xaaa)?&gtk_widget_show:&gtk_widget_hide))(GTK_WIDGET(tw->arg[3].p));
			(*((k&0x555)?&gtk_widget_show:&gtk_widget_hide))(GTK_WIDGET(tw->arg[4].p));
			arg += 4; continue;
		case '!':
			for (i=0; i<6; i++) {
				if (!(j=arg[2*i+1]) || !(k=arg[2*i+2])) return;
				if ((j-=48)<0) j = -1234567890; if ((k-=48)<2) k = 0;
				tw->arg[9].c[i] = j;
				*arg = 65+i; dacnt_set_x(widg_lookup_ps(tw, arg), j, 768+k);
				*arg = 97+i; dacnt_set_x(widg_lookup_ps(tw, arg), j, 768+k);
			}
			arg += 13; continue;
		case 'M':
			ww = widg_lookup_ps(tw, arg); 
			q = ww->arg[2].c; q[0] = arg[1]; q[1] = arg[2]; q[2] = 0;
			q = ww->arg[4].c; memcpy(q, arg+3, 6); q[6] = 1;
			da_fullre(ww); arg += 9; continue;
		case '#':
			ww = widg_lookup_ps(tw, arg); (*ww->cl->cmd)(ww, arg+1); return;
		case 'T':
			tw->title[0] = '#';
			for (i=1; i<20 && arg[i]&&arg[i]!=36; i++) tw->title[i] = arg[i]; tw->title[i] = 0;
			if (!arg[i]) return;  arg += i+1; continue;
		case '*':
			ww = widg_lookup_ps(tw, "#");
			gtk_window_set_focus(GTK_WINDOW(tw->w), ww->w); ++arg; continue;
		case 0: return;
		default: LOG("tgrid: unknown cmd 0x%x(%c)", *arg, *arg); return;
	}}}

static void tgrid_skel (struct _topwin * tw, char * arg) { const char * str = 
  "([7({!as1$332$c90}{!bs2$332$c91}{!cs3$332$c92}{!ds4$332$c93}{!es5$332$c94}{!fs6$332$c95})"
    "0{MMW.W$|#0}" //"0{MM+$G$|3$_3|set values (r-click)$_2|copy2clipb.(m-click)$?|help}"
    "7({!As1$332$c90}{!Bs2$332$c91}{!Cs3$332$c92}{!Ds4$332$c93}{!Es5$332$c94}{!Fs6$332$c95})]"
     "3{##XG})";
	if (tw->state) { if (arg) tgrid_cmd(tw, arg); return; }
	tw->ix4 = 3;
	tw->arg[0].p = parse_w(tw, &str);
	gtk_widget_show(GTK_WIDGET(tw->arg[3].p));
	gtk_widget_show(GTK_WIDGET(tw->arg[4].p));
	if (arg) tgrid_cmd(tw, arg);
}

///////////////// keycfg /////////////////////////////////////////////////////

#define KCF_IX(x) ((x)->arg[2].c[0])
#define KCF_X(x)  ((x)->arg[2].c[1])
#define KCF_Y(x)  ((x)->arg[2].c[2])

static void dakcf_draw(ww_t * ww, cairo_t * cr2) {
        short * a4 = ww->arg[4].s;
        int x0 = a4[0], x1 = a4[2],
            y0 = a4[1], y1 = a4[3];
	int i, x, y, ix = KCF_IX(ww), wid = KCF_X(ww), heig = KCF_Y(ww);
	int h1 = conf_lbh, w1 = (5*h1+2)>>2, wtot = wid*w1, htot = heig*h1;
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE); cairo_set_line_width(cr2, 2.0);
	cairo_set_source_rgb(cr2, .05*(ix&4), .1*(ix&2), .2*(ix&1)); cairo_paint(cr2);
	cairo_set_source_rgb(cr2, .8, .8, .8);
	for (i=0; i<=wid; i++) if (x = 1+i*w1, x>=x0-1, x<=x1+1)
		cairo_move_to(cr2, (double)x, 0.0), cairo_rel_line_to(cr2, 0.0, (double)(htot+2));
	for (i=0; i<=heig;i++) if (y = 1+i*h1, y>=y0-1, y<=y1+1)
		cairo_move_to(cr2, 0.0, (double)y), cairo_rel_line_to(cr2, (double)(wtot+2), 0.0);
	cairo_stroke(cr2);
}

static void dakcf_cmd(struct _ww_t * ww, const char * s) {
	if (!s) return dakcf_draw(ww, (cairo_t*)ww->arg[0].p);
	switch(*s) {
		case 's': KCF_X(ww) = b32_to_i(s[1]), KCF_Y(ww) = b32_to_i(s[2]); da_fullre(ww); return;
		default: return LOG("dakcf:unknown cmd 0x%x/%c", *s, *s);
	}}

static void dakcf_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	LOG("kcf click!!"); }

static void dakcf_skel(struct _ww_t * ww, const char **pp) {
	int i = (*pp)[0] - 48; if ((unsigned int)i > 7u) goto err; else KCF_IX(ww) = i;
	int x = (*pp)[1] - 48; if ((unsigned int)i > 30u) goto err; else KCF_X(ww)  = x;
	int y = (*pp)[2] - 48; if ((unsigned int)i > 30u) goto err; else KCF_Y(ww)  = y;
	da_skel(ww, ((5*conf_lbh+2)>>2)*x+2, conf_lbh*y+2); *pp += 3;
err:    LOG("dakcf_skel: wrong arg");
}

static void tkcf_cmd (struct _topwin * tw, char * arg) {
	LOG("tkcf: no commands yet"); }

static void tkcf_skel (struct _topwin * tw, char * arg) { 
	const char * str = "[{*0044}]";
	if (!tw->state) tw->arg[0].p = parse_w(tw, &str);
  	if (arg) tkcf_cmd(tw, arg);
}

///////////////// track //////////////////////////////////////////////////////

#define TBXF_CRE 0x1000000
#define TBXF_RNG 0x2000000
#define TBXF_DEL 0x4000000
#define TBXF_CLK 0x8000000
#define TCF_PEND 1
#define TCF_OOM  2

#define TC_GD(X)  (tdiv_gd + (X)->gwfr[0])
#define TC_PPB(X) (tdiv_gd[(X)->gwfr[1]].v)
#define TC_UPP(X) (tdiv_gd[95-(X)->gwfr[1]].v)
typedef struct { int id, x15y4; unsigned short nx; char s[14]; } trk_24_b;
typedef union  { unsigned short s[12]; trk_24_b b;  } trk_24;

typedef struct {
        ww_t * ww;
        int x0, y0, w_nb, w_nl, dl_px, dl_py, dl_xcur, dl_tz,
            x0r, x1r, y0r, y1r, x0a, x1a, y0a, y1a, cb0, cb1, cl0, cl1;
	trk_24 * dl_bcur;
	cairo_t * dl_cr;
        unsigned int rbv, ybv[16];
        trk_24 * b3k[290];
        int b3k_n, t24_f, flg, nid, selx, sely;
        unsigned short blk[128];
	unsigned char gwfr[4], div[256];
} tc_t;

static int tc_mk_b3k(tc_t * tc) {
        int i, k, n = tc->b3k_n;
        if (n > 289) return LOG("trk 0x%x out of mem!!!", tc->nid), 0; else ++tc->b3k_n, k = 128*n;
        trk_24 * p = tc->b3k[n] = (trk_24*)alloc3k();
        for (i=0; i<127; i++) p[i].b.nx = k+i+1; p[i].b.nx = 0;
        return k + !n;
}

static inline int le40320(int x, int y) { return x<=y && (x&7)<=(y&7) && !(x&~y&24); }
static inline int tc_effdiv(tc_t * tc, int i) { return i=tc->div[i], i?i:tc->div[0]; }
static inline trk_24 * tc_p24(tc_t * tc, int ix) { return tc->b3k[ix>>7]+(ix&127); }

static int tc_f24(tc_t * tc, int ix) {
	if (!ix) return 0; trk_24 * p = tc_p24(tc, ix); int k = p->b.nx;
	p->b.nx = tc->t24_f; memcpy(p->b.s, "BUG", 4);  return k; }

static trk_24 * tc_a24(tc_t * tc, unsigned short * to) {
	trk_24 * r = tc_p24(tc, ((*to=tc->t24_f)||(*to=tc_mk_b3k(tc))||(tc->flg|=TCF_OOM),*to));
	if (*to) tc->t24_f = r->b.nx;     return r; }

static void tc_vw_cmd(tc_t * tc) { CMD("X#%x$V%c%c$%x$%x",
		tc->nid, tc->y0a+48, tc->y1a+48, tc->x0a, tc->x1a); }

static void tc_clear(tc_t * tc, int re) {
	trk_24 sv, *q;
	if (re) { memcpy(&sv, tc_p24(tc, 1), sizeof(trk_24));
		  int i; for (i=0; i<tc->b3k_n; i++) free3k(tc->b3k[i]); }
        tc->b3k_n = tc->t24_f = tc->flg = 0; tc->x1r = tc->x1a = -1;
	unsigned short qw=0; q = tc_a24(tc, &qw); if (qw!=1) LOG("BUG: tc_clear: %d!=1", qw);
	if (re) memcpy(q, &sv, sizeof(trk_24)); else q->b.id = -1; 
}

#define TC_BLK(P,X,Y) ((P)->blk[((X)&15)+16*((Y)&7)])
#define TC_BLKZ(P,X,Y) memset(tc_a24((P),(P)->blk+(((X)&15)+16*((Y)&7))), 0, 16)

static tc_t* tc_mk(ww_t*ww) {
	tc_t*r=(tc_t*)alloc3k(); tc_clear(r,0); r->nid=(r->ww=ww)->top->id>>4; 
	r->x0a = 0; r->x1a = -1; r->selx = 0; r->sely = 1; 
	return r;
}

static int tc_ddef(int gd) {
        int i,k;
        unsigned int bv = tdiv_gd[gd].d;
        for (i=0; i<8; i++) if (bv & (1u << (k = (0xe9452130>>(4*i)) & 15))) goto ok;
        LOG("BUG!!! ddef: gd=%d v=%d\n", gd, tdiv_gd[gd].v); return 1;
ok:     return tdiv_dvix[k] + 1;
}

static void tc_mk_xy_a(tc_t * tc) {
        int i, j, i1 = tc->y1a = tc->y1r, j1 = tc->x1a = tc->x1r;
	if (j1<0) { if (dflg & DF_TRK) LOG("tc_mk_xy_a: x1r = %d", j1); return; }
        for (j=tc->x0a=tc->x0r; j<=j1; j++) {
                tc->ybv[j&15]=0u; for (i=tc->y0a=tc->y0r; i<=i1; i++) TC_BLKZ(tc,j,i); }
	tc_vw_cmd(tc);
}

static void tc_del_b64(tc_t * tc, int ix) {
	if (!ix) { LOG("BUG: tc_del_b64: ix=0"); return; }
        int i, j, k;   trk_24 *q, *p = tc_p24(tc, ix);
        for (i=0; i<8; i++) if ((k=p->s[i])) for (q=tc_p24(tc, k), j=0; j<8; j++)
                for (k=q->s[j]; k; k = tc_f24(tc, k)); }

static void tc_add_ns(tc_t * tc, int ns) {
        int j, y = ns ? ++tc->y1a : --tc->y0a;
        unsigned int msk = ~(1u << y);
        for (j = tc->x0a; j <= tc->x1a; j++) tc->ybv[j&15] &= msk, TC_BLKZ(tc, j, y); }

static void tc_drop_ns(tc_t * tc, int ns) {
        int j, y = ns ? tc->y1a-- : tc->y0a++;
        for (j = tc->x0a; j <= tc->x1a; j++) tc_del_b64(tc, TC_BLK(tc, j, y)); }

static void tc_add_we(tc_t * tc, int we) {
        int i, x = we ? ++tc->x1a : --tc->x0a;  tc->ybv[x&15] = 0u;
        for (i = tc->y0a; i <= tc->y1a; i++) TC_BLKZ(tc, x, i); }

static void tc_drop_we(tc_t * tc, int we) {
        int i, x = we ? tc->x1a-- : tc->x0a++;
        for (i = tc->y0a; i <= tc->y1a; i++) tc_del_b64(tc, TC_BLK(tc, x, i)); }

static void tc_adj_xy_a(tc_t * tc) {
	//LOG("tc_adj_xy_a xr:%d:%d yr:%d:%d xa:%d:%d ya:%d:%d", tc->x0r, tc->x1r, tc->y0r, tc->y1r, tc->x0a, tc->x1a, tc->y0a, tc->y1a);
        if (tc->x1a<0) return tc_mk_xy_a(tc);
        if (tc->x1a<tc->x0r || tc->x0a>tc->x1r || tc->y1a<tc->y0r || tc->y0a>tc->y1r)
                return tc_clear(tc, 1), tc_mk_xy_a(tc);
        int dw = tc->x0a - tc->x0r, de = tc->x1r - tc->x1a,
            dn = tc->y0a - tc->y0r, ds = tc->y1r - tc->y1a;
	if (!(((dn+1)|(dw+1)|(de+1)|(ds+1))&(-2))) return;
        while (dw<-1) tc_drop_we(tc, 0),++dw;  while (de<-1) tc_drop_we(tc, 1),++de;
        while (dn<-1) tc_drop_ns(tc, 0),++dn;  while (ds<-1) tc_drop_ns(tc, 1),++ds;
        while (dw>0) { --dw; if (tc->x1a-tc->x0a>14) tc_drop_we(tc,1);  tc_add_we(tc,0); }
        while (de>0) { --de; if (tc->x1a-tc->x0a>14) tc_drop_we(tc,0);  tc_add_we(tc,1); }
        while (dn>0) { --dn; if (tc->y1a-tc->y0a> 6) tc_drop_ns(tc,1);  tc_add_ns(tc,0); }
        while (ds>0) { --ds; if (tc->y1a-tc->y0a> 6) tc_drop_ns(tc,0);  tc_add_ns(tc,1); }
	tc_vw_cmd(tc);
}

static void trk_d1b2(tc_t * tc, trk_24 * bx, int px, int wi) {
	if (bx) /*LOG("d1b2: px=%d, dl_px=%d, '%.7s'",  px, tc->dl_px, bx->b.s),*/ da_wr18(tc->dl_cr, tc->dl_px+px, tc->dl_py, wi<0?"B!U!G!a:zz%z%%":bx->b.s, wi<0?18:wi);}

static void trk_d1b(tc_t * tc, trk_24 * bx, int x16) {
	//LOG("trk_d1b: bx=%p", bx);
	if (!bx) return trk_d1b2(tc, tc->dl_bcur, tc->dl_xcur, 18); 
	int px = ((bx->b.x15y4>>4) + x16 - tc->dl_tz) / TC_UPP(tc);
	int wi = px - tc->dl_xcur;
	if (wi<3) { if (!wi) tc->dl_bcur = bx; if (wi<0) LOG("WTF:%d", wi), trk_d1b2(tc, bx, px, -1); return; }
	if (tc->dl_bcur) trk_d1b2(tc, tc->dl_bcur, tc->dl_xcur, wi);
	tc->dl_bcur = bx; tc->dl_xcur = px;
}

#define BXCL(L) ((dflg&DF_TRK)?LOG("tc_bx_cmd/%s: x15y4=0x%x, idr_flg=0x%x s14=\"%s\"", \
				    #L, x15y4, idr_flg, s14):(void)0)
static int tc_bx_cmd(tc_t * tc, int ib, int jb, int idr_flg, const char * s14) {
	if (dflg&DF_TRK) LOG("tc_bx_cmd: ib=0x%x jb=0x%x, idr_flg=0x%x, s14=%.14s",
				ib, jb, idr_flg, s14?s14:"(null pointer)");
        int y64 = ib >> 7, x64 = jb >> 18;
        if (y64<tc->y0a || y64>tc->y1a || x64<tc->x0a || x64>tc->x1a) return -1;
        int y8 = (ib >> 4) & 7, x8 = (jb >> 15) & 7, cf = idr_flg & TBXF_CRE;
	if (cf && !*s14) LOG("tc_bx_cmd: wtf???"), abort();
        unsigned short *b8 = tc_p24(tc, TC_BLK(tc, x64, y64))->s + x8;
        if (!*b8) { if (cf) memset(tc_a24(tc, b8)->s, 0, 16); else return 0; }
        trk_24 *cur = 0;
        int r, x16, lim, ci = 0, x15y4 = 16*(jb&32767) + (ib&15);
        unsigned short *b1 = tc_p24(tc, *b8)->s + y8;
        for(;;){if (!(ci = *b1)) goto eol;
                cur = tc_p24(tc, ci);
                int d19 = x15y4 - cur->b.x15y4;
                if (!d19) goto bingo;  if (d19<0) goto ins;
                b1 = &cur->b.nx; }
bingo:  BXCL(bingo); if (cf) return memcpy(cur->b.s, s14, 14), 1;
        if (idr_flg & TBXF_RNG) goto rng;
        return (idr_flg & TBXF_DEL) ? (*b1 = cur->b.nx, tc_f24(tc, ci), 1) : -2;
eol:    BXCL(eol); if (!cf) return 0; else goto fill;
ins:    BXCL(ins); if (idr_flg & TBXF_RNG) goto rng;   if (!cf) return 0;
fill:   BXCL(fill); (cur = tc_a24(tc, b1))->b.nx = ci; memcpy(cur->b.s, s14, 14);
        cur->b.id = idr_flg&0xfffff; cur->b.x15y4 = x15y4; return 1;
rng:    BXCL(rng); if ((lim = x15y4 + 16*(idr_flg&65535)) <= cur->b.x15y4) return 0;
        if (idr_flg & TBXF_CLK) {
                do r = cur->b.id; while (cur->b.nx && (cur=tc_p24(tc, cur->b.nx))->b.x15y4 < lim); return r;
        } else {
		x16 = jb&0xffff8000;
		do trk_d1b(tc,cur,x16); while (cur->b.nx && (cur=tc_p24(tc, cur->b.nx))->b.x15y4 < lim);
		return 0;
	}}

static void trk_bxline(ww_t * ww, cairo_t * cr2, int by, int x0, int x1) {
	int k;  tc_t * tc = (tc_t*) ww->etc;
	/*LOG("trk_bxline: y=%d, x=%d..%d", by, x0, x1);*/
	tc->dl_xcur = -9999; tc->dl_bcur = NULL; tc->dl_tz = 40320 * tc->x0; tc->dl_cr = cr2;
	if (++x1 <= x0) return;
	if (x0&32767) tc_bx_cmd(tc, by, x0, 32768|TBXF_RNG, NULL), x0 &= 0xffff8000, x0 += 32768;
	while (x0 < (x1 & 0xffff8000)) tc_bx_cmd(tc, by, x0, 32768|TBXF_RNG, NULL),  x0 += 32768;
	if ((k=x1&32767)) tc_bx_cmd(tc, by, x0, k|TBXF_RNG, NULL);
	trk_d1b(tc, NULL, 0);
}

static void arrow(cairo_t * cr2, int x, int y, int d) {
	cairo_set_source_rgb(cr2, .1, .2, .2); 
	cairo_rectangle(cr2, (double)(x+1), (double)(y+1), 18.0, 18.0); cairo_fill(cr2);
	cairo_set_source_rgb (cr2, .6667, .6667, .6667);
	cairo_rectangle(cr2, (double)(x+.5), (double)(y+.5), 19.0, 19.0); cairo_stroke(cr2);
	static const char xydat[6] = { 2, 10, 18, 2, 18, 18 };
	char xy[6]; int i,k;
	for (i=0; i<6; i++) k = xydat[i^(d&1)], xy[i] = (d&2)?20-k:k;
	cairo_set_source_rgb(cr2, .5, 1.0, 1.0);
	cairo_move_to(cr2, (double)(x+xy[0]), (double)(y+xy[1]));
	cairo_line_to(cr2, (double)(x+xy[2]), (double)(y+xy[3]));
	cairo_line_to(cr2, (double)(x+xy[4]), (double)(y+xy[5])); cairo_fill(cr2);
	cairo_set_source_rgb(cr2, .8, .0, .0); xy[0] = d>>2; xy[1] = 0;
	tx_box(cr2, x+3, y+3, 14, 14, conf_lbfs, xy);
}

typedef void (*trk_dfun) (ww_t*, cairo_t*, int *, int *);

static void trk_ns(ww_t * ww, cairo_t * cr2, int * blim, int * bclip) {
	tc_t * tc = (tc_t*) ww->etc;
	cairo_set_source_rgb (cr2, 1.0, 1.0, 1.0);
	int i, j, c, ppb = TC_PPB(tc), flg = (blim[2] > 40), nb = tc->w_nb;
	for (i=0; i<nb; i++) cairo_move_to(cr2, (double)(i*ppb)+60.5, (double)blim[2]+.5),
			     cairo_rel_line_to(cr2, 0.0, 19.5);
	char buf[8]; cairo_stroke(cr2); cairo_set_source_rgb (cr2, .8, .8, .8);
	for (i=0; i<nb; i++) buf[sprintf(buf, "%d", tc->x0+i)] = 0, 
			     tx_box(cr2, i*ppb+61, blim[2]+1, ppb-2, 18, conf_lbfs, buf);
	int xul = bclip[1] + 1, xll = bclip[0] - 20;
	for (i=0; i<3; i++) if (j=20*i, xll<j&&j<xul) arrow(cr2, j, blim[2], 4*"IXL"[i]+2*flg+1);
	if (flg) { for (i=0; i<4; i++) { c = 4*"MCXI"[i];
		if (c=4*"MCXI"[i], j=20*i+60, xll<j&&j<xul) arrow(cr2, j, blim[2], c);
		if (j = DA_W(ww)+40-j,        xll<j&&j<xul) arrow(cr2, j, blim[2], c+2); }}
	else if (xll-30<(j=DA_W(ww)-50)) {
		cairo_set_source_rgb(cr2, .5, 1.0, 1.0);
		cairo_rectangle(cr2, (double)(j+1), 1.0, 49.0, 19.0); cairo_fill(cr2);
		cairo_set_source_rgb (cr2, .6667, .6667, .6667);
		cairo_rectangle(cr2, (double)(j+.5), .5, 50.0, 20.0); cairo_stroke(cr2);
		cairo_set_source_rgb(cr2, .0, .0, .0); tx_box(cr2, j+3, 3, 46, 16, conf_lbfs, "clip");
	}
}

static void trk_e(ww_t * ww, cairo_t * cr2, int * blim, int * bclip) { 
	tc_t * tc = (tc_t*) ww->etc;
	cairo_set_source_rgb (cr2, .8, .8, .8); cairo_paint(cr2);
	cairo_set_source_rgb (cr2, .4, .4, .4);
	int i, x = blim[0], y = blim[2], hh = blim[3]-y+1;
	char buf[2]; buf[1] = 0;
	double x0 = (double)x+.5;
	for (i=tc->cl0; i<=tc->cl1; i++) cairo_move_to(cr2, x0, (double)(48*i+y)+.5), 
					 cairo_rel_line_to(cr2, 11.0, 0.0);
	cairo_stroke(cr2); cairo_set_source_rgb (cr2, .0, .0, .0);
	for (i=tc->cl0; i<=tc->cl1; i++) 
		buf[0] = hexc1((i+tc->y0)>>4), tx_box(cr2, x+1, 48*i+y+6 , 9, 16, conf_lbfs_s, buf),
		buf[0] = hexc1((i+tc->y0)&15), tx_box(cr2, x+1, 48*i+y+24, 9, 16, conf_lbfs_s, buf);
	cairo_set_source_rgba (cr2, 0.0, 0.0, 1.0, 0.5);
	double sb0 = (double)hh * (1.0/256.0) * (double)tc->y0,
	       sbh = (double)hh * (1.0/12288.0) * (double)(blim[3]-blim[2]+1);
	cairo_rectangle(cr2, x0, (double)y + sb0, 12, sbh); cairo_fill(cr2);
}

#define QR65 (1e-5*(double)(random()&65535))
static void trk_w(ww_t * ww, cairo_t * cr2, int * blim, int * bclip) {
	tc_t * tc = (tc_t*) ww->etc;
	int i,j,k,x; char buf[4];
	cairo_set_source_rgb (cr2, .8, .8, .8);
	for (i=tc->cl0, buf[2] = 0; i<=tc->cl1; i++) 
		cairo_move_to(cr2, .5, (double)(48*i+blim[2]) + 47.5), cairo_rel_line_to(cr2, 59.0, 0.0),
		cairo_rel_move_to(cr2, 0.0, -24.0), cairo_rel_line_to(cr2, -25.0, 0.0);
	cairo_stroke(cr2);
	for (i=tc->cl0, buf[1] = 0; i<=tc->cl1; i++) 
		k = i+tc->y0, j = 48*i + blim[2],
		buf[0] = hexc1((k>>4)&15), tx_box(cr2, 1, j+3 , 10, 20, conf_lbfs, buf),
		buf[0] = hexc1(k&15)     , tx_box(cr2, 1, j+24, 10, 20, conf_lbfs, buf);
	cairo_set_source_rgb(cr2, 1.0, 1.0, 0.5);
	for (i=tc->cl0, j = i+tc->y0, buf[3] = 0; i<=tc->cl1; i++, j++)
		k = (tc->div[j]), *buf=k?'/':'(',  k = k?k:tc->div[0], x = 48*i+blim[2],
		          d99(buf+1, tdiv_idsf[k].d), tx_box(cr2, 33, x+ 2, 25, 22, conf_lbfs, buf),
		*buf=':', d59(buf+1, tdiv_idsf[k].s), tx_box(cr2, 33, x+24, 25, 22, conf_lbfs, buf);
	cairo_set_source_rgb(cr2, 1.0, 1.0, 0.5);
}

static void trk_c(ww_t * ww, cairo_t * cr2, int * blim, int * bclip) { 
	tc_t * tc = (tc_t*) ww->etc;
	cairo_set_source_rgb (cr2, .8, .8, .8);
	int i, j, k, k0, i0, ppb = TC_PPB(tc), xd;
	double xx0 = (double)bclip[0], wid = (double)(bclip[1] - bclip[0] + 1),
	       	   yy0 = (double)bclip[2], heig = (double)(bclip[3] - bclip[2] + 1);
	cairo_set_source_rgb (cr2, .3, .3, .3);
	for (i = i0 = tc->cl0, j = i+tc->y0, k=k0=0; k>=0; i++,j++) {
		k = (i>tc->cl1) ? -1 : tc->div[j], k = k?k:tc->div[0];
		if (k==k0) continue; if (!k0) {k0=k; continue; }
		int d7 = tdiv_idsf[k0].d, s6 = tdiv_idsf[k0].s, inc1 = ppb/d7, inc2 = ppb/s6, j1, j2, j3;
		if (!s6 || d7==ppb) goto c2;
		double dy = (double)(blim[2]+i0*48) + .5, dh = (double)((i-i0)*48);
		for (j1 = tc->cb0, xd = blim[0]+ppb*j1; j1 <= tc->cb1; j1++, xd+=ppb)
			for (j2=0; j2<ppb; j2+=inc2) for (j3=inc1; j3<inc2; j3+=inc1)
				cairo_move_to(cr2, (double)(xd+j2+j3)+.5, dy), cairo_rel_line_to(cr2, 0.0, dh);
	    c2: k0 = k; i0 = i;
	}
	cairo_stroke(cr2); cairo_set_source_rgb (cr2, .6, .6, .6);
	for (i = i0 = tc->cl0, j = i+tc->y0, k=k0=0; k>=0; i++,j++) {
		k = (i>tc->cl1) ? -1 : tc->div[j], k = k?k:tc->div[0];
		if (k==k0) continue; if (!k0) {k0=k; continue; }
		int s6 = tdiv_idsf[k0].s;
		if (s6==ppb) s6 = (s6&1) ? (s6%3?s6/5:s6/3) : (s6>>1);
		int inc1 = ppb/s6, j1, j2;
		double dy = (double)(blim[2]+i0*48) + .5, dh = (double)((i-i0)*48);
		for (j1 = tc->cb0, xd = blim[0]+ppb*j1; j1 <= tc->cb1; j1++, xd+=ppb)
			for (j2=0; j2<ppb; j2+=inc1) 
				cairo_move_to(cr2, (double)(xd+j2)+.5, dy), cairo_rel_line_to(cr2, 0.0, dh);
	        k0 = k; i0 = i;
	}
	cairo_stroke(cr2); cairo_set_source_rgb (cr2, .8, .8, .8);
	for (i=tc->cl0; i<=tc->cl1; i++)  k = blim[2]+48*i,
		cairo_move_to(cr2, xx0, (double)k +  2.5), cairo_rel_line_to(cr2, wid, 0.0),
		cairo_move_to(cr2, xx0, (double)k + 45.5), cairo_rel_line_to(cr2, wid, 0.0);
	cairo_stroke(cr2); cairo_set_source_rgb (cr2, 1.0, 1.0, 1.0);
	for (i=tc->cb0; i<=tc->cb1; i++) 
		cairo_move_to(cr2, (double)(blim[0]+ppb*i)+.5, yy0), cairo_rel_line_to(cr2, 0.0, heig);
	cairo_stroke(cr2);
	int upp = TC_UPP(tc), tz = 40320 * tc->x0, t0 = max_i(0, tz + upp * (bclip[0]-78)),
	    t1 = tz + upp * (bclip[1] + 1), syw = (tc->sely>>4) - tc->y0;
	//LOG("cl0 = %d, cl1 = %d, t0 = %d, t1 = %d", tc->cl0, tc->cl1, t0, t1);
	tc->dl_px = *blim; tc->dl_py = blim[2]+48*tc->cl0+2;
	for (i=tc->cl0; i<=tc->cl1; i++, tc->dl_py+=48) trk_bxline(ww, cr2, 16*(i+tc->y0), t0, t1);
	//LOG("trkc/sel: t0s1: %d %d %d  l0s1: %d %d %d", t0,tc->selx,t1, tc->cl0, syw, tc->cl1);
	if (t0<=tc->selx && tc->selx<=t1 && tc->cl0<=syw && syw<=tc->cl1) {
		trk_24 * q = tc_p24(tc, 1);
		syw *= 48; syw += blim[2];
		int sxw = (tc->selx - tz) / upp + blim[0];
		if (q->b.id) da_wr18(cr2, sxw, syw+2, q->b.s,18);
		cairo_set_source_rgb (cr2, 1.0, 0.0, 0.0);
		cairo_move_to(cr2, (double)sxw+.5, (double)syw+2.5);
		cairo_rel_line_to(cr2,  17.0, 0.0); cairo_rel_line_to(cr2, 0.0,  43.0);
		cairo_rel_line_to(cr2, -17.0, 0.0); cairo_rel_line_to(cr2, 0.0, -43.0); cairo_stroke(cr2);
		cairo_set_source_rgb (cr2, 1.0, 1.0, 0.0);
		cairo_move_to(cr2, (double)sxw-.5, (double)syw+1.5);
		cairo_rel_line_to(cr2,  19.0, 0.0); cairo_rel_line_to(cr2, 0.0,  45.0);
		cairo_rel_line_to(cr2, -19.0, 0.0); cairo_rel_line_to(cr2, 0.0, -45.0); cairo_stroke(cr2);
	}
}

static void trk_5(trk_dfun fun, ww_t * ww, cairo_t * cr2, int * blim) {
	int clip[4]; short * cl = ww->arg[4].s;
	if ( (clip[0] = max_i(blim[0], cl[0])) > (clip[1] = min_i(blim[1], cl[2])) ||
	     (clip[2] = max_i(blim[2], cl[1])) > (clip[3] = min_i(blim[3], cl[3])) ) return;
	cairo_save(cr2); 
	cairo_rectangle(cr2, (double)clip[0], (double)clip[2], 
		(double)(clip[1]-clip[0]+1), (double)(clip[3]-clip[2]+1));
	cairo_clip(cr2); (*fun)(ww, cr2, blim, clip); cairo_restore(cr2);
}

static void trk_reqdr(tc_t * tc, int x, unsigned int ym) {
	char buf[24]; buf[0] = 'X'; buf[1] = 'R'; *(int*)(buf+2) = qh4(x);
	buf[6] = 36; *(int*)(buf+7) = qh4((int)(ym>>16)); *(int*)(buf+11) = qh4(ym&65535); buf[15] = 0;
	widg_defcmd(tc->ww, buf); tc->flg &= TCF_PEND;
	if (dflg&DF_TRK) LOG("trk 0x%x: req x=%d, msk=0x%x", tc->nid, x, ym); }

static int trk_cond_req(tc_t * tc) {
	int i; unsigned int m2;
	for (i=tc->x0r; i<=tc->x1r; i++) if ((m2 = ~tc->ybv[i&15] & tc->rbv)) return trk_reqdr(tc, i, m2), 1;
	return 0; }

static int trk_draw_pre(ww_t * ww) {
        tc_t * tc = (tc_t*) ww->etc;
        int r = 0, k;
        if (tc->x0<0) r=1, tc->x0=0;
        else if (tc->x0 > (k = 52002-tc->w_nb)) r=1, tc->x0 = k;
        if (tc->y0<1) r=1, tc->y0=1;
        else if (tc->y0 > (k = 257 - tc->w_nl)) r=1, tc->y0 = k;
        if (r && bigda_counter) { if (!tc->flg&TCF_PEND) da_fullre(ww);  return 1; }
        int x0 = 40320*tc->x0, x1 = x0 + 40320*tc->w_nb + 262143,
            y0 = tc->y0, y1 = y0 + tc->w_nl + 7;
        tc->x0r = x0>>18; tc->x1r = x1>>18; tc->y0r = y0>>3; tc->y1r = y1>>3;
	tc->rbv = 0xffffffffu<<tc->y0r; if (tc->y1r!=31) tc->rbv &= ((1u<<tc->y1r)-1u);
	tc_adj_xy_a(tc); return (tc->flg&TCF_PEND) ? 1 : trk_cond_req(tc);
}

static void datrk_draw(ww_t * ww, cairo_t * cr2) {
        tc_t * tc = (tc_t*) ww->etc;
        int wx = DA_W(ww), wy = DA_H(ww);
        tc->w_nl = (wy+8)/48, tc->w_nb = (wx-70)/TC_PPB(tc) + 1;
        if (trk_draw_pre(ww)) return;
	int cbz = 51999 - tc->x0, clz = 255 - tc->y0;
        tc->cb0 = max_i((ww->arg[4].s[0]-60)/TC_PPB(tc), 0);
        tc->cb1 = min_i((ww->arg[4].s[2]-71)/TC_PPB(tc)+1, cbz);
        tc->cl0 = max_i((ww->arg[4].s[1]-20)/48, 0);
        tc->cl1 = min_i((ww->arg[4].s[3]+28)/48, clz);
        cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE); cairo_set_line_width(cr2, 1.0);
        int bxlim[4];
        bxlim[0] = bxlim[2] = 0; bxlim[1] = wx-1; bxlim[3] = 19; trk_5(trk_ns, ww, cr2, bxlim);
        bxlim[2] = wy-20; bxlim[3] = wy-1;                       trk_5(trk_ns, ww, cr2, bxlim);
        bxlim[2] = 20; bxlim[3] = wy-21; bxlim[1] = 59;          trk_5(trk_w , ww, cr2, bxlim);
        bxlim[0] = wx-12; bxlim[1] = wx-1;                       trk_5(trk_e , ww, cr2, bxlim);
        bxlim[0] = 60;    bxlim[1] = wx-13;                      trk_5(trk_c , ww, cr2, bxlim);
}

static void trk_mv(ww_t * ww, int dx, int dy) { 
        tc_t * tc = (tc_t*) ww->etc;
	if (dflg&DF_TRK) LOG("trk_mv: x:%d y:%d", dx, dy); 
	tc->x0 += dx; tc->y0 += dy; da_fullre(ww);
}

static void trk_ns_playm(ww_t * ww, int b9, int cx, GdkEventButton * ev) {
	tc_t * tc = (tc_t*) ww->etc; tsc_x = tc->x0 + (cx - 60) / TC_PPB(tc); 
	tsc_ww = ww; tsc_md = 3; popup2(ww, tsc_mi+4, 0, b9, ev); }

static void trkclk_n(ww_t * ww, int b9, int cx, GdkEventButton * ev) {
	tc_t * tc;
	if (cx<60) trk_mv(ww, 0, -((0621201>>((cx/20)*6))&63));
	else if (DA_W(ww)-cx<50) tsc_drag_id = tsc_drag_x = tsc_drag_y = -1, tc = (tc_t*)ww->etc,
		tsc_ww = ww, tsc_drag_xd = (cx -= 60),
		CMD("X#%x$Cv0$%x0$%x", ww->top->id>>4, tc->y0, 40320*tc->x0 + TC_UPP(tc)*(cx-9));
	else trk_ns_playm(ww, b9, cx, ev);
}

static void trkclk_s(ww_t * ww, int b9, int cx, GdkEventButton * ev) {
	static const int mcxi[4] = { 1000, 100, 10, 1 };
	if (cx<140) return ((cx/=20)<3) ? trk_mv(ww, 0, (0621201>>(cx*6))&63)
				        : trk_mv(ww, -mcxi[cx-3], 0);
	if (DA_W(ww)-cx<80) return        trk_mv(ww, mcxi[(DA_W(ww)-cx)/20], 0);
	trk_ns_playm(ww, b9, cx, ev);
}

static void trkclk_e(ww_t * ww, int b9, int cy) { LOG("trkclk_e: b%c y:%d", 48+b9, cy); }

static void tsc_cc(int c) {
	CMD("X#%x$C%c%x$%x$%x", tsc_ww->top->id>>4, c, tsc_id, tsc_y, tsc_x);
	if (dflg&DF_TRK) LOG("tsc_cmd:  id:0x%x, x:%d y:%d c:0x%x '%c'", tsc_id, tsc_x, tsc_y, c, c);
}

static void tsc_act(ww_t * ww, int ix) { int buf[3];
	if (ww != tsc_ww) { LOG("BUG: tsc_act: %p!=%p", ww, tsc_ww); return; }
	switch(tsc_md) { case 0: tsc_cc(menutab[tsc_mi].cmd[2*ix]); return;
			 case 1: tsc_x = ix; ttrk_cmd(ww->top, "1"); return;
			 case 2: tsc_x = ix; ttrk_cmd(ww->top, "2"); return;
			 case 3: buf[0] = (menutab[tsc_mi+4].cmd[2*ix]<<24)+0x705800; buf[1] = qh4(tsc_x);
				 buf[2] = 0; widg_defcmd(ww, (char*)buf + 1); return;
			 default: LOG("BUG: tsc_act: md=%d", tsc_md); return; }}

static void trkclk_w(ww_t * ww, int b9, int cx, int cy, GdkEventButton *ev) {
	if (dflg&DF_TRK) LOG("trkclk_w: b%c x:%d y:%d", 48+b9, cx, cy);
	if (cx<20) return LOG("trkclk_w: nothing happens");
        tc_t * tc = (tc_t*) ww->etc;   tsc_y = cy/48+tc->y0; tsc_ww = ww;
	if (cy%48<24) tsc_md=1, popup2(ww, tsc_mi+1, DLM_MSK(widg_lookup_ps(ww->top,"d"))&~1u, b9, ev);
	else tsc_md=2, popup2(ww, tsc_mi+2, 0x7ffffff^tdiv_sdbv[(int)tdiv_idsf[tc_effdiv(tc,tsc_y)].i], b9,ev);
}

static void trkclk_c(ww_t * ww, int b9, int cx, int cy, GdkEventButton *ev) { 
        tc_t * tc = (tc_t*) ww->etc;
	int k, y1 = (cy/48+tc->y0), y = 16*y1;
	//if (tc->gwfr[2]&1) k = TC_PPB(tc)/tdiv_idsf[tc_effdiv(tc,y1)].d, cx = k*((2*cx/k+1)/2);
	int r, tz = 40320*tc->x0, upp = TC_UPP(tc), t0 = max_i(0,tz+upp*(cx-18)), t1 = tz+upp*cx;
	if (!((t0^t1)&0xffff8000)) r = tc_bx_cmd(tc, y, t0, (t1-t0)|TBXF_RNG|TBXF_CLK, 0);
	else if (!(r = tc_bx_cmd(tc, y, t1&0xffff8000, (t1&32767)|TBXF_RNG|TBXF_CLK, 0)))
		   r = tc_bx_cmd(tc, y, t0, (32768-(t0&32767))|TBXF_RNG|TBXF_CLK, 0);
	if (!r && (tc->gwfr[2]&1)) k = TC_PPB(tc)/tdiv_idsf[tc_effdiv(tc,y1)].d,
			           cx = k*((2*cx/k+1)/2), t1 = tz+upp*cx;
	tsc_md = 0; tsc_x = t1; tsc_y = y; tsc_id = r; tsc_ww = ww;
	if (dflg&DF_TRK) LOG("trkclk_c: b%c x:%d(%d) y:%d(%d) --> 0x%x", 48+b9, cx, t1, cy, y, r);
	if (b9==3) return popup2(ww, tsc_mi, r ? 1 : 0x37, b9, ev);
	tsc_cc(48+b9); 
	if (r && b9==1) tsc_drag_id = tsc_drag_x = tsc_drag_y = -1, tsc_drag_xd = cx;
}

static gboolean trk_rel(GtkWidget *w, GdkEventButton * ev, gpointer p) { tsc_drag_id=0; return TRUE; }

static gboolean trk_drag(GtkWidget *w, GdkEventMotion * ev, gpointer p) {
	if (tsc_drag_id <= 0) return TRUE;
	if (tsc_ww != (ww_t*)p) return LOG("trk_drag: wigdet mismatch"), tsc_drag_id = 0, TRUE;
	tc_t * tc = (tc_t*)tsc_ww->etc;
	int cx = ivlim((int)lround(ev->x)-60-tsc_drag_xd, 0, DA_W(tsc_ww)-78), k,
	    cy = ivlim((int)lround(ev->y)-20, 		  0, DA_H(tsc_ww)-41), ty = cy/48+tc->y0;
	if (tc->gwfr[2]&1) k = TC_PPB(tc)/tdiv_idsf[tc_effdiv(tc,ty)].d, cx = k*((2*cx/k+1)/2);
	int tx = 40320*tc->x0+TC_UPP(tc)*cx;
	if (!((tx-tsc_drag_x)|(ty-tsc_drag_y))) return TRUE; else tsc_drag_x = tx, tsc_drag_y = ty;
	if (dflg&DF_TRK) LOG("drag (%d,%d) --> (%x:%x) (x0=%d, y0=%d, upp=%d)", 
			            cx,cy,      ty,tx, tc->x0, tc->y0,TC_UPP(tc));
	CMD("X#%x$CM%x$%x0$%x", tsc_ww->top->id>>4, tsc_drag_id, ty, tx);
	return TRUE;
}

static void datrk_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	tsc_drag_id = 0;
	if (b9<0) { if (b9!=-1) return;  // TODO: 'Q' (?)
		    LOG("trk_key: 0x%x 0x%x", cx, cy);
		    int k = b32_to_i(cy); if (k>=0) return CMD("K@p$1%c", k+48);
		    if ((cy&=255)==8) cy = 'X';
		    else if (cy<32 || cy>126) return LOG("trk: unknown key %d", cy);
		    k = (cy<<16) + 0x4b58; widg_defcmd(ww, (char*)&k); return;  }
	int wx = DA_W(ww), wy = DA_H(ww);
	if (cy<20)    return trkclk_n(ww, b9, cx, ev); 
	if (cy>wy-21) return trkclk_s(ww, b9, cx, ev); else cy -= 20;
	if (cx<60) return trkclk_w(ww, b9, cx, cy, ev);
	if (cx>wx-13) return trkclk_e(ww, b9, cy);
	trkclk_c(ww, b9, cx-60, cy, ev);
}

static void tsr_op(ww_t * ww, int y, int x);
static void trk_sel(ww_t * ww, int i, int j, int id, const char * s14) {
        tc_t * tc = (tc_t*) ww->etc;   topwin * tw = ww->top;  int f;
	if (i<0) f = 3, i = tc->sely, j = tc->selx; 
	else tsr_op(ww,tc->sely,tc->selx), f=(tc->selx!=j&&(tc->selx=j,1))|(3*(tc->sely!=i&&(tc->sely=i,1))),
	     tsr_op(ww,tc->sely,tc->selx);
	if (tsc_drag_id==-1) tsc_drag_id = id, tsc_drag_xd -= (j-40320*tc->x0)/TC_UPP(tc);
	int k, d = tdiv_idsf[tc_effdiv(tc, i>>4)].d, d2 = 40320/d;
	if (f&2) dacnt_set_x(widg_lookup_ps(tw,"L"), i, 256), d99(DACNT_LBL(widg_lookup_ps(tw,"D")), d);
	if(f&1){int a[4]; k=TC_UPP(tc), a[0]=j/40320; j%=40320; a[1]=j/d2; j%=d2; a[2]=j/k; a[3]=1000*(j%k)/k;
		for (k=0;k<4;k++) dacnt_set_x(widg_lookup_ps(tw, "B\0D\0P\0U" + 2*k), a[k], 256);        }
	trk_24 * q = tc_p24(tc, 1);
	if (id) q->b.id = id, memcpy(q->b.s, s14, 14), dalbl_mw18(widg_lookup_ps(tw, "S"), s14);
	else if (q->b.id) q->b.id = 0, dalbl_mw18(widg_lookup_ps(tw, "S"), 0);
}

static void trk_send_gcf(ww_t * ww, unsigned int * bv9) {
        tc_t * tc = (tc_t*) ww->etc;
	int buf[300], i, n = 0, id = tc->nid, *q = buf+3;
	buf[0] = 0x587e0000; buf[1] = qh4(id>>8) - 13;
	buf[2] = hexc1((id>>4)&15) + 256*hexc1(id&15) + 0x67240000;
	for (i=0; i<4; i++) if (bv9[8]&(1u<<i)) q[n++] = qh4((i<<8)  + tc->gwfr[i]) + 47;
	for (i=0; i<8; i++) { BVFOR_JMC(bv9[i]) q[n++] = qh4((i<<13) + (j<<8) + tc->div[32*i+j]); }
	if (n) q[n] = 10, write(1, (char*)buf+2, 4*n+11);
}

#define TSR_START (-1)
#define TSR_DROP  (-2)
#define TSR_ALL   (-3)
#define TSR_END   (-4)
static void tsr_op(ww_t * ww, int y, int x) {
	static ww_t * cww = 0;
	static unsigned int flg, cbv[2], gbv[9];
	static short cx0[64], cx1[64];
	int i;
	if (y==TSR_START) { if (cww) LOG("BUG: TSR_START: %p != 0", cww); cww = ww; cbv[0]=cbv[1]=flg=0u;
			    for (i=0; i<9; i++) gbv[i] = 0u; return; }
	if (ww != cww) { LOG("BUG: unexp. tsr_op(%p,%d,%d): cur=%p", ww, y, x, cww); return; }
        tc_t * tc = (tc_t*) ww->etc;
	if (y>=0) {
		if(x<-1) { if (y>4095) return (void)(gbv[8] |= 1u<<(y&3));
			   else gbv[(y>>9)&7] |= 1u<<((y>>4)&31); }
		if ((tc->flg&TCF_PEND) || (flg&1u)) return; else y>>=4;
		int wid = DA_W(ww), heig = DA_H(ww), yb = y-tc->y0, yp = 48*yb + 20;
		if (yp<0 || yp>heig) return;
		if (x<0) { if ((tc->sely>>4)==y) trk_sel(ww, -1, 0, 0, NULL);
			   cbv[yb>>5] |= 1u<<(yb&31), cx0[yb] = 0, cx1[yb]=wid-12; return; }
		int xp0 = (x - 40320*tc->x0) / TC_UPP(tc) + 58, xp1 = min_i(wid-10, xp0+22);
		if (xp0<60) xp0 = 60;    if (xp0>xp1) return;
		unsigned int msk = 1u<<(yb&31), *pbv = cbv + (yb>>5);
		if (!(*pbv&msk)) return (void)(*pbv|=msk, cx0[yb] = xp0, cx1[yb] = xp1);
		if (xp0<cx0[yb]) cx0[yb] = xp0;
		if (xp1>cx1[yb]) cx1[yb] = xp1;    return;
	}
	switch(y) {
		case TSR_DROP: trk_send_gcf(ww, gbv);  cww = 0; return; 
		case TSR_END:  break;
		case TSR_ALL:  flg|=1u; return;
		default:       LOG("BUG: invalid tsr_op(%p,%d,%d) cur=%p", ww, y, x, cww); return;
	}
	trk_send_gcf(ww, gbv);
	cww = 0; if (tc->flg&TCF_PEND || !(ww->w)) return;
	if (flg&1u) return da_fullre(ww);
	unsigned int j,m;
	for (i=0; i<2; i++) { for (m=cbv[i]; j=__builtin_ffs(m), j-- > 0; m &= ~(1u<<j)) {
		int yi = 32*i+j; struct _GdkRectangle rct; 
		rct.x = cx0[yi]; rct.width = cx1[yi] - cx0[yi] + 1;
		rct.y = 48*yi + 20; rct.height = 48;
		gdk_window_invalidate_rect(gtk_widget_get_window(ww->w), &rct, 0);
	}}}

static void datrk_cmd(struct _ww_t * ww, const char * s) {
	if (!s) return datrk_draw(ww, (cairo_t*)ww->arg[0].p);
	if (*s=='>') ++s; else tsr_op(ww, TSR_START, 0);
        tc_t * tc = (tc_t*) ww->etc;
	int i,hx[4],r, x, y;
	while (1) { switch(*s) {
		case 0:   goto done;
		case '.': for (i=tc->x0r; i<=tc->x1r; i++) tc->ybv[i&15]|=tc->rbv; 
			  tc->flg &= ~TCF_PEND; da_fullre(ww); s++; continue;
		case '+': for (i=0; i<4; i++) hx[i] = qh4rs(s+1+4*i); 
			  x = (hx[2]<<16)+hx[3]; y = hx[0]>>4;  tsr_op(ww, y, x);
			  r=tc_bx_cmd(tc, y, x, ((hx[0]&15)<<16)|hx[1]|TBXF_CRE, s+17);
			  if (r<0) LOG("trk/+: ret=%d", r);
			  s += 31; continue;
		case '-': for (i=0; i<4; i++) hx[i] = qh4rs(s+1+4*i);
			  x = (hx[2]<<16)+hx[3]; y = hx[0]>>4;  tsr_op(ww, y, x);
			  r=tc_bx_cmd(tc, y, x, ((hx[0]&255)<<16)|hx[1]|TBXF_DEL, s+17);
			  if (r<0) LOG("trk/-: ret=%d", r);
			  s += 17; continue;
		case '*': for (i=0; i<3; i++) hx[i] = qh4rs(s+1+4*i);
			  LOG("trk/*: %x %x %x", hx[0], hx[1], hx[2]);
			  if (*hx<tc->x0a||*hx>tc->x1a) { LOG("datrk/*: x(%d) out of range", *hx); goto done;}
			  tc->ybv[*hx&15] |= ((unsigned int)hx[1]<<16)+hx[2]; tc->flg &= ~TCF_PEND;
			  da_fullre(ww); tsr_op(ww, TSR_DROP, 0); return;
		case ',': widg_defcmd(ww, "XR,"); return;
		default:  LOG("datrk_cmd: unknow ch 0x%x(%c)", *s, *s); goto done;
	}}
done:	tsr_op(ww, TSR_END, 0);
}

static void trk_upd_gds(topwin *tw, int i, int v, unsigned int m) {
	ww_t * ww2 = widg_lookup_ps(tw, "M\0d\0s"+2*i);
	int k; char * s = DALBL_TXT(ww2); 
	s[0] = i ? "?/:"[i] : (k=v/100, v-=100*k, k+48-16*!k);
	s[1] = 48+v/10; s[2] = 48+v%10; s[3] = 0;  DLM_MSK(ww2) = m;  da_fullre(ww2);
}

static void trk_upd_wgd(ww_t * ww, int flg) { // 1:div 2:gd 4:wid
        tc_t * tc = (tc_t*) ww->etc;   topwin * tw = ww->top;
	const tdiv_idsf_t * q;
	if (flg&4) dacnt_set_x(widg_lookup_ps(tw, "W"), TC_PPB(tc), 256);
	if (flg&2) trk_upd_gds(tw, 0, TC_GD(tc)->v, 255 ^ TC_GD(tc)->m),
		 DLM_MSK(widg_lookup_ps(tw, "d")) = (2 * ~TC_GD(tc)->d) + 1;
	if (flg&1) q = tdiv_idsf+tc->div[0], trk_upd_gds(tw, 1, q->d, ~(2*TC_GD(tc)->d)),
					     trk_upd_gds(tw, 2, q->s, ~tdiv_sdbv[(int)q->i]);
}

static void tc_set_d1(tc_t * tc, int i, int d) {
	if (!i) return tc->div[0] = d, trk_upd_wgd(tc->ww, 1);
	int k = (tc->div[i]==d);  tc->div[i] = d; tsr_op(tc->ww, 16*i, k-2); }

static void tc_upd255(tc_t * tc, int flg) {  
	int i; unsigned char *p = tc->div;
	if (!(flg&2)) {
		if (flg&1) for (i=1; i<256; i++) if (!p[i]) tc_set_d1(tc, i, 0);     return; }
	unsigned int m = TC_GD(tc)->d;
	if(flg&1) for (i=1; i<256; i++) if (!p[i] || !((1u<<tdiv_idsf[p[i]].i)&m)) tc_set_d1(tc, i, 0);
	else	  for (i=1; i<256; i++) if ( p[i] && !((1u<<tdiv_idsf[p[i]].i)&m)) tc_set_d1(tc, i, 0);
}

static void tc_set_ppb(tc_t * tc, int v) { if (dflg&DF_TRK)  LOG("tc_set_ppb %d", v);
	if (v==tc->gwfr[1]) return;
	if (!le40320(tc->gwfr[0], v)) return LOG("BUG: tc_set_ppb: v=%d, g=%d", v, tc->gwfr[0]);
	tc->gwfr[1] = v;  tsr_op(tc->ww, TSR_ALL, 0);  trk_upd_wgd(tc->ww, 4);
}

static int tc_find_ppb(tc_t * tc, int i, int d) { int x, b = tc->gwfr[0]; LOG("tc_find_ppb %d %d", i, d);
	for (; (x=tdiv_ppb_ix[i]); i+=d) if (le40320(b,x)) return tc_set_ppb(tc, x), 1;   return 0; }

static void tc_set_gd(tc_t * tc, int gd) {
	if (gd==tc->gwfr[0]) return;
	unsigned int dbv0 = TC_GD(tc)->d, dbv = tdiv_gd[tc->gwfr[0]=gd].d;
	int j, ppb = tc->gwfr[1]; 
	le40320(gd, ppb) || tc_find_ppb(tc,j=tdiv_gd[ppb].j,1) || tc_find_ppb(tc,j,-1);
	trk_upd_wgd(tc->ww, 2);
	if (dbv0 & ~dbv) tc_upd255(tc, (dbv&(1u<<tc->div[0])) ? 2 : (tc->div[0]=tc_ddef(gd), 3) );
}

static void tc_dcmd(tc_t * tc, int i, int v, int cf) {
	if (dflg&DF_TRK) LOG("tc_dcmd: %d %d(%d) %d", i, v, tc->div[i], cf);
	if (tc->div[i] == v) return;
	tc_set_d1(tc, i, v); if (!i) tc_upd255(tc, 1);
}

static void datrk_skel(struct _ww_t * ww, const char **pp) {
	tc_t * tc = tc_mk(ww); ww->etc = tc; memcpy(tc->gwfr, TRK_DEF_GWFR, 4);
	int i;  tc->x0 = 0; tc->y0 = 1;
	tc->div[0] = 3; for (i=1; i<256; i++) tc->div[i] = 0;
	da_skel(ww, 425, 232);
	g_signal_connect (ww->w, "button-release-event", G_CALLBACK(trk_rel), (gpointer)ww);
	g_signal_connect (ww->w, "motion-notify-event", G_CALLBACK(trk_drag), (gpointer)ww);
}

static void tc_adj_gd(tc_t * tc, int k) {
	static const signed char dif[8] = {1, 32, 8, 16, -1, -32, -8, -16};
	if (TC_GD(tc)->m & 1<<(k&=7)) tc_set_gd(tc, tc->gwfr[0]+dif[k]); else LOG("tc_adj_gd: out of rng"); }

static void tc_adj_w(tc_t * tc, int k) {
	int w0 = tc->gwfr[1]; switch(k&3) {
		case 0: tc_find_ppb(tc,		      1,  1); return;
		case 1: tc_find_ppb(tc, tdiv_gd[w0].j-1, -1); return;
		case 2: tc_find_ppb(tc,	     TDIV_N_PPB, -1); return;
		case 3: tc_find_ppb(tc, tdiv_gd[w0].j+1,  1); return;
	}}

static void tc_adj_d(tc_t * tc, int i, int v) {
	if (!v) { if(i) tc_dcmd(tc, i, 0, 1); else LOG("tc_adj_d: 0=0");  return; }
	int d0 = tc->div[i]; if (!d0 || tdiv_idsf[d0].i!=v-1) tc_dcmd(tc, i, tdiv_dvix[v-1]+1, 1); }
	
static void tc_adj_s(tc_t * tc, int i, int v) {
	int d0 = tc_effdiv(tc, i), di = tdiv_idsf[d0].i;
	unsigned int dm = tdiv_sdbv[di], m1 = (1u<<v);
	if (!(dm & m1)) return LOG("tc_adj_s(%d,%d): not allowed", i, v);
	tc_dcmd(tc, i, tdiv_dvix[di] + bitcnt(dm&(m1-1u)), 1);
}

static void ttrk_cmd (struct _topwin * tw, char *s) {
	ww_t *wn, *wt = widg_lookup_ps(tw, "t");
        tc_t * tc = (tc_t*) wt->etc;
	int i, k, a[4]; char buf[8];
	tsr_op(wt, TSR_START, 0);
	while (1) { switch(*s) {
		case 0: goto done;
		case '^': 
			for (i=0,++s; i<4; i++,s+=4) a[i] = qh4rs(s);
			trk_sel(wt, a[0]>>4, (a[2]<<16)+a[3], k = ((a[0]&15)<<16)+a[1], s);
			if (k) s += 14;     continue;
		case '>':
			return (*wt->cl->cmd)(wt, s);
		case '1': tc_adj_d(tc, tsc_y, tsc_x); goto done;
		case '2': tc_adj_s(tc, tsc_y, tsc_x); goto done;
		case 'g': ++s; s += trk_g_parse(s, tc->div, tc->gwfr); 
			  DABOOL(wn=widg_lookup_ps(tw,"A")) = tc->gwfr[2]&1; da_fullre(wn); continue;
		case '*':
			memcpy(buf, "Xmd", 3); memcpy(buf+3, s+1, 2);
			k = TC_PPB(tc) / tdiv_idsf[tc_effdiv(tc, tc->sely>>4)].d;
			buf[5] = hexc1(k>>4); buf[6] = hexc1(k&15); buf[7] = 0; 
			widg_defcmd(wt, buf); goto done;
		case 'G': switch(s[1]) {
			case 'D': tc_adj_d(tc, 0, s[2]-48); goto done;
			case 'S': tc_adj_s(tc, 0, s[2]-48); goto done;
			case 'W': tc_adj_w(tc, s[2]^s[3]); goto done;
			case 'M': tc_adj_gd(tc, s[2]); goto done;
			case 'A': DABOOL(wn=widg_lookup_ps(tw,"A")) = (tc->gwfr[2]^=1)&1;
				  da_fullre(wn); tsr_op(wt, 4098, -2); goto done;
			default: break; }
		default: LOG("ttrk: invalid cmd 0x%x \"%s\"", *s, s); goto done;
	}}
done:	tsr_op(wt, TSR_END, 0);
	
}

static void ttrk_skel (struct _topwin * tw, char * arg) { const char * str = 
	"[" TW_TOPH
	  "({B<`<$Xs<}{MS*QWEQWE$Xc|T0}{B>`>$Xs>}{YAali$$>GA}"
	   "{8LL$.bXml}{8BB$5Xmb}{8D08\\$02$>*}{8Pp$02Xmp}{8Um$03Xmu}3()0"
	   "{Ypplay$Xp}{8bbpm$.4Xb}{YRrec$Xr}{MM999$$>GM|T3}{Md$$>GD|T1}{Ms$$>GS|T2}{8Wp/b$3$>GW})3{tt}]";
	if (!tsc_mi) tsc_mi = menutab_lu('T', 0);
	if (tw->state) { if (arg) ttrk_cmd(tw, arg); return; }
	tw->arg[0].p = parse_w(tw, &str);
	if (arg) ttrk_cmd(tw, arg); 
	trk_upd_wgd(widg_lookup_ps(tw, "t"),  7);
}

///////////////// doc ////////////////////////////////////////////////////////

#define ECHK(v,e) if ((v)<0) return perror(e), -1
#define SCHK(v,s,l) if (memcmp((v),(s),(l))) return fprintf(stderr,"help: '%s' expected\n", (s)), -1

#define HFUN(L,R) inline static int hfun##L(const char *s, int m) {        \
        unsigned int r = 0; while (*s) r = (r<<L) + (r>>R) + *(s++);        \
        (m<65536) && (r^=(r>>16), m<256) && (r^=(r>>8), m<16) && (r^=(r>>4));\
        return m & (int)r; }
HFUN(3,29) 
HFUN(5,27)


static char *help_buf, **help_htab, *help_vpos;
static int help_hmsk;
static char help_notfound[] = "4help$text$not$found";

int help_readbuf(const char * fn) {
        int fd = open(fn, O_RDONLY); ECHK(fd, fn);
        int len = lseek(fd, 0, SEEK_END); ECHK(len, "lseek");
        if (!len) return LOG("empty help file"), -1;
        help_buf = malloc(len); ECHK(lseek(fd, 0, SEEK_SET),"lseek(0)");
        int r = read(fd, help_buf, len);
        if (r!=len) { if (r<0) perror("read"); else LOG("read %d/%d",r,len); return -1; }
        char *p = help_buf, *q = p, *z = p+len-1, *ky, *tx;
        if (*z!=10) return LOG("help: last chr is 0x%x, should be 0xa", *z), -1; else *z=0;
        int i, nw = 0;
	for (i=0; p[i]!=10; i++);  for(++i; p[i]!=10; i++);   help_vpos = p+i+2;
        SCHK(p, "{{{", 3);
        while(q<z) {
                p = q; q += 3; ++nw;
                int nl = 0;
                for (ky = q; *q!='}'; ++q) ;
                SCHK(q,"}}}\n",4); *q = 0; p[2] = 48+q-ky; tx = (q+=4);
                while(*q>31 || (++nl, *q==10 && (*q = memcmp(q+1,"{{{",3) ? 36 : 0))) q++;
                if (*q) return LOG("help: invalid chr 0x%x", *q), -1;
                int l = ++q-p;  p[0] = 48+(l>>6), p[1] = 48+(l&63), tx[-1] = 48+nl;
                //LOG("info: Ttlk* \"%s\", nv* \"%s\"", p, tx-1);
        }
        int hmin=(7*nw+2)>>2, hsiz=1;
        while (hsiz<hmin) hsiz+=hsiz;
        LOG("help: \"%s\": %d help windows, hsiz=%d", fn, nw, hsiz);
        help_htab = calloc(hsiz, sizeof(char*)); help_hmsk = hsiz - 1;
        for (i=0, p=help_buf; i<nw; i++, p += 64*p[0]+p[1]-06060) {
                int h0 = hfun3(p+3, help_hmsk);
                if (help_htab[h0]) { int h1 = 2*hfun5(p+3,help_hmsk)+1; 
				     do h0 += h1, h0 &= help_hmsk; while (help_htab[h0]); }
                help_htab[h0] = p;
                //LOG("info(%d): Ttlk* \"%s\", nv* \"%s\", h0:%d", i, p, p+p[2]-42, h0);
        }
        return 0;
}

char * help_lu(const char * ky) {
        int h0 = hfun3(ky, help_hmsk);
        char * p = help_htab[h0]; if (!p) return help_notfound;
        if (!strcmp(ky, p+3)) return p+p[2]-42;
        int h1 = 2*hfun5(ky, help_hmsk) + 1;
	while (1) {
                p = help_htab[h0+=h1,h0&=help_hmsk];
                if (!p) return help_notfound;
                if (!strcmp(ky, p+3)) return p+p[2]-42;
        }}

GtkWidget * doc_vbl (struct _ww_t * ww, int ix) { return parse_w_s(ww->top, "{C0448$0$%%%ccc}"); }

static void doc_skel (struct _topwin * tw, char * arg) {
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, "[(3{C_300$1$ttt666...}0{__}{M_+$|+0}){:DdS10}]");
   	if (arg) doc_cmd(tw, arg);
}

static void doc_cmd (struct _topwin * tw, char * s) {
	ww_t * w = widg_lookup_ps(tw, "D.");
	if (*s=='!') s = help_lu(s+1);
	int i, i0 = VB_WBASE(w), n = *(s++) - 48;
	if (n>32) LOG("BUG: help txt has %d lines, cut to 32", n), n = 32;
	vbox_show_bv(w, (n==32 ? 0 : (1<<n))-1); const char * ahh = s;
	for (i=0; i<n; i++) daclb_set(widg_p(tw, i0+i), &ahh, 15), ahh+=(*ahh==36);
}

///////////////// gui conf ///////////////////////////////////////////////////

GtkWidget * gconf_vbl (struct _ww_t * ww, int ix) {
	topwin * tw = ww->top;
	GtkWidget * rw = parse_w_s(tw, "({L0i.99}[{e15$E00}{_2}]"
			"{L4o.99}[{e55$E20}{_3}])"); // TODO: def.v, desc
	int ix0 = VB_WBASE(ww) + 8*ix,
	    ten = ix<10 ? 48 : 49 + (ix>19), one = ix + 528 - 10 * ten;
	char * s = DALBL_TXT(widg_p(tw, ix0  )); s[1] = ten; s[2] = one; s[3] = 0;
	       s = DALBL_TXT(widg_p(tw, ix0+4)); s[1] = ten; s[2] = one; s[3] = 0;
	           //DALBL_TXT(widg_p(tw, ix0+6))[0] = ">>>short description 63c/ln>>>"[ix];
	s = widg_p(tw, ix0 + 1)->cmd; s[1] += (ix>>4); s[2] = hexc1(ix&15);
	s = widg_p(tw, ix0 + 5)->cmd; s[1] += (ix>>4); s[2] = hexc1(ix&15);
        return rw;
}

static void gconf_skel (struct _topwin * tw, char * arg) { const char * str = 
  "([3({!rred$155$C0wy}{!ggrn$155$C1wy}{!bblu$155$C2wy})0(3{1w}[3{Byok$EG%.}{Bnno$EG!}])"
    "3({!Rred$155$C3wy}{!Ggrn$155$C4wy}{!Bblu$155$C5wy})]3{:LgN80})";
	if (tw->state) { if (arg) gconf_cmd(tw, arg); return; }
	tw->ix4 = 3;
	tw->arg[0].p = parse_w(tw, &str);
	ww_t * ww = widg_lookup_ps(tw, "w"); memcpy(ww->arg[3].c, "FG^vBGM:", 8);
	ww = widg_lookup_ps(tw, "L."); vbox_cmd(ww, "+N");
	if (arg) gconf_cmd(tw, arg);
}

static void gconf_cmd (struct _topwin * tw, char * s) {
	int i,i0=0; ww_t *w, *w2; char *p; while (*s) { switch (*s) {
		case 'R': for (++s, i=0; i<6; i++) dacnt_set_x(widg_lookup_ps(tw,"rgbRGB"+i), s[i]-37,  0x355);
			  for (i=0; i<3; i++) memcpy(p=(w2=widg_lookup_ps(tw, "wny"+i))->arg[4].c, s, 6),
				  	      p[6] = i, da_fullre(w2);
			  s += 6; break;
		case 48: case 49: case 50: case 51:
			  if (!i0) i0 = VB_WBASE(w = widg_lookup_pci(tw, 'L', -1));
			  i = i0 + 1 + 2*(*s&2) + 8*(16*(s[0]&1)+hxd2i(s[1])); s += 2;
			  p = s; while (*s && *s!=36) ++s; if (*s==36) *(s++) = 0;
			  entry_set(widg_p(tw, i), p);
			  break;
		case 'T':
			  for (i=0, s++; i<59 && *s && *s!=36; i++, s++) tw->title[i] = *s;
			  memcpy(tw->title+i, "(UI)", 5); break;

		case '$': s++; break;
		default: LOG("gconf_cmd: unknown cmd \"%s\"",s); return;
}}}

///////////////// graph //////////////////////////////////////////////////////

typedef struct {
	char nm[24];
	short x0, x1, y0, y01, y02, y1;
	int iw8, ow8;
	float x, y, w, h;
	char * s; // 15*ni+5*no
	int ni, no;
	char rgb[6];
} gr_node;

typedef struct {
	short fr, to;
	short len;
	float * p;
} gr_edge;

typedef struct {
	int seq, nn, ne, flg; // 1:i0 2:v0 4:i1 8:v1
	float bx, by, x0, y0, scx, scy;
	int bxp, byp;
	int nn0, nn1, ne1;
	int sel[3]; // in box out
	gr_node * node;
	gr_edge * edge;
} gr_dat;

#define GR_IL(p,i) ((p)->s + 15*(i))
#define GR_IV(p,i) ((p)->s + 15*(i) + 5)
#define GR_OL(p,i) ((p)->s + 15*(p)->ni + 5*(i))
#define GR_ICOL(p,i) ( ((128+(i)*(p)->iw8)>>8) + (p)->x0 )
#define GR_OCOL(p,i) ( ((128+(i)*(p)->ow8)>>8) + (p)->x0 )

static gr_dat * gr_new(int nn, int ne) {
	int k = sizeof(gr_dat) + nn*sizeof(gr_node);
	char * p = calloc(k + ne*sizeof(gr_edge), 1);
	gr_dat * q = (gr_dat*) p;
	q->node = (gr_node*) (p + sizeof(gr_dat));
	q->edge = (gr_edge*) (p + k);
	q->nn = nn; q->ne = ne;
	q->sel[0] = q->sel[1] = q->sel[2] = -1;
	return q;
}

static void gr_del(gr_dat * p) {
	int i; if (!p) return;
	for (i=0; i<p->nn; i++) free(p->node[i].s);
	for (i=0; i<p->ne; i++) free(p->edge[i].p);
	free(p);
}

static void gr_node_ini(gr_node *p, int ni, int no) {
	p->ni = ni; p->no = no; p->s = malloc(15*ni + 5*no); }


#define DAGRAPH_DAT(x)  (*(gr_dat**)(&((x)->arg[2].p)))
#define DAGRAPH_XDAT(x) (*(gr_dat**)(&((x)->arg[3].p)))

static int dagraph_active(ww_t * ww, int i) {
	gr_dat * p = DAGRAPH_XDAT(ww); i=9*i+3; return p && (p->flg&i)==i; }

static void gr_draw_edge(cairo_t * cr, gr_dat * dat, int i) {
	gr_edge * ed = dat->edge + i;
	double xo = dat->x0, yo = dat->y0, sx = dat->scx, sy = dat->scy;
	float * p = ed->p; if (!p) { LOG("gr_draw_edge: #%d empty", i); return; }
	int j, n = ed->len;
	double rgb[3]; get_rgb8(rgb, 45*i, 1.0);
	cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
	cairo_move_to(cr, xo + sx*p[0], yo + sy*p[1]);
	for (j=1; j<n-2; j+=3) cairo_curve_to(cr, xo+sx*p[2*j  ], yo+sy*p[2*j+1],
					 	  xo+sx*p[2*j+2], yo+sy*p[2*j+3],
						  xo+sx*p[2*j+4], yo+sy*p[2*j+5]);
	for (; j<n; j++) cairo_line_to(cr, xo + sx*p[2*j], yo + sy*p[2*j+1]);
	cairo_stroke(cr);
}

static void gr_node_geom(gr_dat * dat, int i) {
	gr_node * nd = dat->node + i;
	int x0, x1, y0, y1, ni=nd->ni, no=nd->no;
	double xo = dat->x0, yo = dat->y0, sx = dat->scx, sy = dat->scy;
	nd->x0 = x0 = (int)lround(xo + sx * (nd->x - .5*nd->w));
	nd->x1 = x1 = (int)lround(xo + sx * (nd->x + .5*nd->w));
	nd->y0 = y0 = (int)lround(yo + sy * (nd->y + .5*nd->h));
	nd->y1 = y1 = (int)lround(yo + sy * (nd->y - .5*nd->h));
	if (!i) nd->y01 = nd->y02 = y0;
	else if (!no) nd->y01 = nd->y02 = y1;
	else if (!ni) nd->y01 = y0, nd->y02 = (0x800 + 0x555*y0 + 0xaaa*y1) >> 12;
	else nd->y01 = (0x800 + 0x99a*y0 + 0x666*y1) >> 12,
	     nd->y02 = (0x800 + 0x333*y0 + 0xccd*y1) >> 12;
	nd->iw8 = ni ? (256*(x1-x0)) / ni : -1;
	nd->ow8 = no ? (256*(x1-x0)) / no : -1;
}

static void gr_draw_node(cairo_t * cr, gr_dat * dat, int i) {
	gr_node * nd = dat->node + i;
	if (!nd->s) { LOG("gr_draw_node: %d empty", i); return; }
	int j, ni = nd->ni, no = nd->no, x0  = nd->x0,  x1  = nd->x1,
	       y0 = nd->y0, y1 = nd->y1, y01 = nd->y01, y02 = nd->y02,
	       x01 = x0 + 3 + FONT_VWID(conf_lbfs_s);
	cairo_set_source_rgba (cr, .0, .0, .0, .667);
	if (y0<y01) { cairo_rectangle(cr, x0, y0, x1-x0, y01-y0); cairo_fill(cr); }
	if (y02<y1) { cairo_rectangle(cr, x0, y02, x1-x0, y1-y02); cairo_fill(cr); }
	for (j=0; j<ni; j++) {
		if (no) cairo_set_source_rgb (cr, .333, .667, 1.0);
		else cairo_set_source_rgb (cr, 1.0, 1.0, 1.0);
		tx_box(cr, GR_ICOL(nd, j)+2, y0 + 2, (nd->iw8>>8)-3, (y01-y0)/2-3,
				conf_lbfs, GR_IL(nd, j));
		const char * vs = GR_IV(nd, j);
		if (*vs==':') cairo_set_source_rgb (cr, 1.0, 0.0, 0.0);
		else if (*vs=='<') cairo_set_source_rgb (cr, 1.0, 1.0, 0.0);
		else cairo_set_source_rgb (cr, .667, 1.0, .667);
		tx_box(cr, GR_ICOL(nd, j)+2, y0 + (y01-y0)/2+1, (nd->iw8>>8)-3, (y01-y0)/2-3,
				conf_lbfs, vs+(*vs=='<'));
	}
	cairo_set_source_rgb (cr, 1.0, (double)(y01==y02), 0.0);
	for (j=0; j<no; j++) tx_box(cr, GR_OCOL(nd, j)+2, y02+2, (nd->ow8>>8)-3, y1-y02-3,
			                conf_lbfs, GR_OL(nd, j));
	if (y01 < y02) {
		cairo_set_source_rgba (cr, .0, .0, .0, .667);
		cairo_rectangle(cr, x0, y01, x01-x0, y02-y01); cairo_fill(cr);
		cairo_set_source_rgb(cr, RGB_C(nd->rgb[3]), RGB_C(nd->rgb[4]), RGB_C(nd->rgb[5]));
		cairo_rectangle(cr, x01+1.0, y01+1.0, x1-x01-2.0, y02-y01-2.0); cairo_fill(cr);
		char buf[4]; buf[1] = buf[3] = 0;
		cairo_set_font_size(cr, conf_lbfs_s); cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
		buf[0] = i_to_b32((i-1)>>5); buf[2] = i_to_b32((i-1)&31);
		cairo_move_to(cr, x0+3.0, y01 + (double)(3+FONT_OFFS(conf_lbfs_s)));
		cairo_show_text(cr, buf);
		cairo_move_to(cr, x0+3.0, y01 + (double)(23+FONT_OFFS(conf_lbfs_s)));
		cairo_show_text(cr, buf + 2);
		cairo_set_source_rgb(cr, RGB_C(nd->rgb[0]), RGB_C(nd->rgb[1]), RGB_C(nd->rgb[2]));
		tx_box_split(cr, x01, y01, x1-x01, y02-y01, 30, nd->nm);
		cairo_set_source_rgb (cr, .667, .667, .667); cairo_set_line_width(cr, 1.0);
		cr_line(cr, x01+.5, y01, x01+.5, y02); cairo_stroke(cr);
	}
	cairo_set_line_width(cr, 2.0);
	cairo_set_source_rgb (cr, .667, .667, .667);
	for (j=1; j<ni; j++) { double x = GR_ICOL(nd, j); cr_line(cr, x, y0, x, y01); }
	for (j=1; j<no; j++) { double x = GR_OCOL(nd, j); cr_line(cr, x, y02, x, y1); }
	cairo_rectangle(cr, x0, y0, x1-x0, y1-y0);
	if (y0<y01&&y01<y1) { cr_line(cr, x0, y01, x1, y01); }
	if (y0<y02&&y02<y1) { cr_line(cr, x0, y02, x1, y02); }
	cairo_stroke(cr);
}

static void dagraph_set_dc(ww_t * ww, int i, int j, const char * s) {
	if (i<0) return;
	gr_dat * dat = DAGRAPH_DAT(ww); if (!dat) { LOG("gr/dc: no graph"); return; }
	cairo_surface_t * surf = DA_SURF(ww); if (!surf) { LOG("gr/dc: no surface"); return; }
	cairo_t * cr = cairo_create(surf);
	gr_node * nd = dat->node + i;
	if (j>=nd->ni) { LOG("gr/dc: j=%d, nd->ni=%d", j, nd->ni); return; }
	int x0 = GR_ICOL(nd, j)+1, x1 = GR_ICOL(nd, j+1)-1,
	    y1 = nd->y01 - 1, y0 = (nd->y0 + y1) / 2 + 1;
	cairo_set_source_rgb (cr, .1, .1, .1);
	cairo_rectangle(cr, x0, y0, x1-x0, y1-y0); cairo_fill(cr);
	cairo_set_source_rgb (cr, .667, 1.0, .667);
	tx_box(cr, x0+1, y0, (x1-x0)-1, (y1-y0)-1, conf_lbfs, s);
	fullsurf_invd(ww, x0, y0, x1, y1);
}

static void dagraph_sel_2(ww_t * ww, int t, int ix, int ny) {
	if (ix<0) return;
	gr_dat * dat = DAGRAPH_DAT(ww); 
	if (!dat) { LOG("gr_sel_2: no graph"); return; }
	cairo_surface_t * surf = DA_SURF(ww);
	if (!surf) { LOG("gr_sel_2: no surface"); return; }
	cairo_t * cr = cairo_create(surf);
	gr_node * nd = dat->node + (ix>>5);
	int x0, x1, y0, y1,  ix2 = ix & 31;
	cairo_set_line_width(cr, 2.0);
	if (!ny) cairo_set_source_rgb (cr, .667, .667, .667);
	switch (t) {
		case 0:
			x0 = GR_ICOL(nd, ix2); x1 = GR_ICOL(nd, ix2+1);
			y0 = nd->y0; y1 = nd->y01;
			if (ny) cairo_set_source_rgb (cr, 0.0, 0.0, 1.0);
			break;
		case 1:
			x0 = nd->x0; x1 = nd->x1; y0 = nd->y01; y1 = nd->y02;
			if (ny) cairo_set_source_rgb(cr, 1.0, 0.0, 1.0);
			break;
		case 2:
			x0 = GR_OCOL(nd, ix2); x1 = GR_OCOL(nd, ix2+1);
			y0 = nd->y02; y1 = nd->y1;
			if (ny) cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
			break;
		default:
			LOG("gr_sel_2: invalid t=%d", t); return;
	}
	cairo_rectangle(cr, x0, y0, x1-x0, y1-y0); 
	cairo_stroke(cr); cairo_destroy(cr);
	fullsurf_invd(ww, x0, y0, x1, y1);
}

static void dagraph_sel(ww_t * ww, int t, int i, int j) {
	gr_dat * dat = DAGRAPH_DAT(ww); 
	if (!dat) { LOG("gr_sel: no graph"); return; }
	if (i<0 || i>=dat->nn) { LOG("gr_sel: invalid i=%d", i); return; }
	if (j<0 || j>31) { LOG("gr_sel: invalid j=%d", j); return; }
	if (t<0 || t> 2) { LOG("gr_sel: invalid t=%d", t); return; }
	int ix = 32*i + j; if (ix==dat->sel[t]) return;
	int * p = dat->sel;
	dagraph_sel_2(ww, t, p[t], 0);
	if (p[t]<0) {} else if (t==1) {
		if (!((p[1] ^ p[0])&~31)) dagraph_sel_2(ww, 0, p[0], 1);
		if (!((p[1] ^ p[2])&~31)) dagraph_sel_2(ww, 2, p[2], 1);
	} else {
		if (!((p[t] ^ p[1])&~31)) dagraph_sel_2(ww, 1, p[1], 1);
	}
	dagraph_sel_2(ww, t, p[t] = ix, 1);
}

static void dagraph_draw(ww_t * ww) {
	if (dflg&DF_GRAPH) LOG("dagraph_draw %dx%d", DA_W(ww), DA_H(ww));
        cairo_surface_t * surf = DA_SURF(ww);
        cairo_t * cr2 = cairo_create(surf);
        cairo_set_source_rgb (cr2, .3, .3, .3);
        cairo_paint(cr2);
	gr_dat * dat = DAGRAPH_DAT(ww); 
	if (!dat) { please_wait(cr2); return; }
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE);
	cairo_set_line_width(cr2, 1.0);
	cairo_set_source_rgb(cr2, .8, .8, .8);
	cairo_rectangle(cr2, .5, .5, (double)(dat->bxp-1), (double)(dat->byp-1));
	cairo_stroke(cr2);
	int i;
	for (i=0; i<dat->ne; i++) gr_draw_edge(cr2, dat, i);
	for (i=0; i<dat->nn; i++) gr_draw_node(cr2, dat, i);
}

static void dagraph_ini(ww_t * ww, int cd, const char * s) {
	if (dflg&DF_GRAPH) LOG("dagraph_ini(%d,\"%s\")", cd, s);
	gr_dat * dat; int m1 = 3*cd + 1;
	int seq= 0; for(; *s>47 && *s!='_'; s++) seq = 16*seq + hxd2i(*s);
	if (*s=='_') ++s; else goto ex_;
	int nn = 0; for(; *s>47 && *s!='_'; s++) nn = 16*nn + hxd2i(*s);
	if (*s=='_') ++s; else goto ex_;
	int ne = 0; for(; *s>47 && *s!='_'; s++) ne = 16*ne + hxd2i(*s);
	if ((dat = DAGRAPH_XDAT(ww))) {
		if (seq == dat->seq) goto chk;
		if ((seq - dat->seq) & 128) goto obs1;
		gr_del(dat); DAGRAPH_XDAT(ww) = NULL;
	} else if ((dat = DAGRAPH_DAT(ww))) {
		if (seq == dat->seq) goto dup;
		if ((seq - dat->seq) & 128) goto obs2;
	}
	DAGRAPH_XDAT(ww) = dat = gr_new(nn, ne);
	dat->seq = seq; dat->flg = 3 * m1;
	return;
ex_:
	LOG("dagraph_ini: err: exp '_', got 0x%x '%c'", *s, *s); return;
obs1:
	if (dat->flg & (3*m1)) goto dup;
	dat->flg |= m1; 
obs2:
	LOG("WARNING: obsolete dagraph_ini(%d,%c), seq:%d after %d", 
			cd, 49-!DAGRAPH_XDAT(ww), seq, dat->seq);
	return;
chk:
	if (dat->flg & (3*m1)) goto dup;
	if ((nn-dat->nn)|(ne-dat->ne))
		LOG("dagraph_ini: ERROR: g%c(%d,%d) after (%d,%d)", 48+cd, 
				nn, ne, dat->nn, dat->ne);
	dat->flg |= 3*m1;
	return;
dup:
	LOG("dagraph_ini: ERROR: i%c,s%d after i%c,s%d", 48+cd, seq, 48+cd, dat->seq);
}

static void dagraph_d_bbox(ww_t * ww, const char * s) {
	if (!dagraph_active(ww, 1)) return;
	gr_dat * dat = DAGRAPH_XDAT(ww);
	float qw;
	int k = sscanf(s, " %f %f %f", &qw, &dat->bx, &dat->by);
	if (k!=3) LOG("dagraph_d_bbox: error: parsed %d/3, \"%s\"", k, s);
}

static void dagraph_cd_node(ww_t * ww, int cd, const char * s) {
	if (!dagraph_active(ww, cd)) return;
	gr_dat * dat = DAGRAPH_XDAT(ww);
	char buf[4]; const char * ids = s;
	float x=0.0, y=0.0, w=0.0, h=0.0;
	if (cd) {
		int k = sscanf(s, " b%4c %f %f %f %f", buf, &x, &y, &w, &h);
		if (k!=5) LOG("dagraph_d_node: err: parsed %d/5, \"%s\"", k, s);
		ids = buf;
	}
	int ix = 32*b32_to_i(ids[0]) + b32_to_i(ids[1]);
	int ni = b32_to_i(ids[2]), no = b32_to_i(ids[3]);
	if (ix<0 || ix>=dat->nn) { LOG("dagraph_%c_node: invalid ix:%d", 99+cd, ix); return; }
	gr_node * p = dat->node + ix;
	if (p->s) {
		if ((p->ni-ni)|(p->no-no)) 
			LOG("dagraph_%c_node: #%d old: %d,%d new:%d,%d",99+cd, ix,p->ni,p->no,ni,no);
	} else {
		gr_node_ini(p, ni, no);
	}
	if (cd) {
		p->x=x; p->y=y; p->w=w; p->h=h; ++dat->nn1;
		double sx = (double)(ni*conf_portwid+2) / w;
		double sy = (double)(2+conf_lbh*(ni?(no?5:2):(ix?3:1))) / h;
		if (sx > dat->scx) dat->scx = sx;
		if (sy > dat->scy) dat->scy = sy;
	} else {
		int i; memcpy(p->rgb, s+4, 6); // for (i=0; i<6; i++) p->rgb[i] = 16*hxd2i(s[2*i+4]) + hxd2i(s[2*i+5]);
		s += 10; get_tok(p->nm, 24, &s, 164);
		for (i=0; i<ni; i++) {
			char *s2, *q = GR_IV(p, i);
			if (*s!='o') {
				double v = hx2doub(s); s += 16;
				q[smallnum(q,v)] = 0;
			} else if (s[1]=='v'&&s[2]=='u') {
				*(q++) = '<'; *(q++) = ':'; 
				s2 = GR_OL(dat->node, b32_to_i(s[3]));
				while (*s2) *(q++) = *(s2++);
				q[0] = ':'; q[1] = 0; s += 4;
			} else {
				q[0] = q[3] = q[5] = 58; q[6] = 0;
				q[1] = s[1]; q[2] = s[2]; q[4] = s[3]; s += 4;
			}
			get_tok(GR_IL(p,i), 5, &s, 164);
		}
		for (i=0; i<no; i++) get_tok(GR_OL(p,i), 5, &s, 164);
	}
}

static void dagraph_d_edge(ww_t * ww, char * s0) {
	if (!dagraph_active(ww, 1)) return;
	gr_dat * dat = DAGRAPH_XDAT(ww);
	char * toks;
	char * sfr = strtok_r(s0," \t",&toks); if (!sfr) goto tok;
	char * sto = strtok_r(NULL," \t",&toks); if (!sto) goto tok;
	char * sn = strtok_r(NULL," \t",&toks); if (!sn) goto tok;
	int n = atoi(sn);
	if (n<2) { LOG("dagraph_d_edge: err: n=%d", n); return; }
	gr_edge * p = dat->edge + (dat->ne1++);
	p->len = n; p->p = malloc(2*n*sizeof(float));
	int i; for (i=0; i<2*n; i++) {
		char * s = strtok_r(NULL," \t",&toks);
		if (!s) { LOG("dagraph_d_edge: found %d/%d coord.", i, 2*n); return; }
		p->p[i] = atof(s);
	}
	return;
tok:    LOG("dagraph_d_edge: expected: <from> <to> <n>");
}

static void graph_sendcmd(ww_t * ww, int b9, int i, int sn) {
	if ((dflg & DF_REC) && !(b9&32)) rec_vpf(ww, "g%c%x", 48|b9, 256*i+sn);
	char buf[8]; buf[0] = b9|48; dagraph_get(buf+1, ww, 's');  --i;
	buf[3] = i_to_b32((i>>5)&31); buf[4] = i_to_b32(i&31); buf[5] = sn; buf[6] = 0;
	if (dflg&DF_GRAPH) LOG("graph_sendcmd: \"%s\"", buf);
	widg_defcmd(ww, buf);
}

static void dagraph_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	gr_dat * dat = DAGRAPH_DAT(ww);
	if (!dat) { LOG("dagraph_click: no graph"); return; }
	int x = cx - ((DA_W(ww)-(dat->bxp))>>1);
	int y = cy - ((DA_H(ww)-(dat->byp))>>1);
	if (x<0 || y<0 || x>dat->bxp || y>dat->byp) {
		LOG("dagraph_click: click outside graph"); return; }
	if (dflg&DF_GRAPH) LOG("dagraph_click: b%d %d,%d", b9, x, y);
	int i; gr_node * nd;
	for (i = 0, nd = dat->node; i<dat->nn; i++,nd++) 
		if (nd->x0<x && x<nd->x1 && nd->y0<y && y<nd->y1) goto found;
	return LOG("dagraph_click: no box at %d,%d", x, y);
found:; int xr = x-nd->x0, sn = (y<nd->y01) ?      48 + (256*xr)/nd->iw8 
		             : ((y<nd->y02) ? 42 : 80 + (256*xr)/nd->ow8);
	graph_sendcmd(ww, b9, i, sn);
	return;
}

static int dagraph_get (void * to, ww_t * ww, int ty) {
	gr_dat * dat = DAGRAPH_DAT(ww);
	if (!dat) { LOG("dagraph_get: no graph"); return 0; }
	char * q = (char*) to;
	switch (ty) {
		case 'B': return dat ? (dat->byp<<16)+dat->bxp : 0x200020;
		case 's': if (dat) q[0] = hexc1(dat->seq>>4), q[1] = hexc1(dat->seq&15);
			  else q[0] = q[1] = '?';
			  return 2;
		default:  LOG("dagraph_get: error: 's' or 'B' expected, got 0x%x(%c)", ty, ty);
			  return 0;
	}}

static void dagraph_resize(ww_t * ww, int wid, int heig) {
	GtkWidget * ch = gtk_bin_get_child(GTK_BIN(ww->arg[4].p));
	gtk_container_remove(GTK_CONTAINER(ww->arg[4].p), ch);
	da_skel(ww, wid, heig);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(ww->arg[4].p), ww->w);
	gtk_widget_show_all(GTK_WIDGET(ww->arg[4].p));
}

static void dagraph_done(ww_t * ww, int cd) {
	gr_dat * dat = DAGRAPH_XDAT(ww);
	if (!dat) return; 
	dat -> flg &= ~(1<<(2*cd));
	if (dat->flg != 10) return;
	gr_del(DAGRAPH_DAT(ww)); DAGRAPH_DAT(ww) = dat; DAGRAPH_XDAT(ww) = NULL;
	dat->bxp = 12 + (int)(dat->scx * dat->bx);
	dat->byp = 12 + (int)(dat->scy * dat->by);
	dat->x0 = 6.0; dat->y0 = (double)dat->byp - 6.0; 
	dat->scy = -dat->scy;
	if (dflg&DF_GRAPH) LOG("dagraph_done: sx:%g sy:%g", dat->scx, dat->scy);
  	int i=0; for (; i<dat->nn; i++) gr_node_geom(dat, i);
	dagraph_resize(ww, dat->bxp, dat->byp);	
}

static void dagraph_cmd(ww_t * ww, const char * s) {
	if (!s) return dagraph_draw(ww);
	double v; char sbuf[16];
	switch(*s) {
		case 'i': dagraph_ini(ww, 0, s+1); return;
		case 'n': dagraph_cd_node(ww, 0, s+1); return;
		case 'z': dagraph_done(ww, 0); return;
		case 'R': dagraph_resize(ww, 10*s[1]-480, 10*s[2]-480); return;
		case '+': dagraph_sel(ww, s[1]-48, 
					  32*b32_to_i(s[2]) + b32_to_i(s[3]),
					  b32_to_i(s[4]) );  return;
		case '#': v = hx2doub(s+4); sbuf[smallnum(sbuf, v)] = 0;
			  dagraph_set_dc(ww, 32*b32_to_i(s[1]) + b32_to_i(s[2]),
					     b32_to_i(s[3]), sbuf); return;
		case '=': dagraph_set_dc(ww, 32*b32_to_i(s[1]) + b32_to_i(s[2]),
					     b32_to_i(s[3]), s+4); return;
		default: LOG("dagraph: invalid cmd 0x%x(%c)", *s, *s); return;
	}
}

static void dagraph_skel(struct _ww_t * ww, const char **pp) {
	da_skel(ww, 32, 32);
	GtkWidget * sw = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(sw), GTK_WIDGET(ww->w));
	gtk_widget_show(ww->w); gtk_widget_show(sw);
	gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW(sw), GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
	ww->arg[4].p = sw; gtk_widget_show_all(sw);
	get_tok(ww->cmd, 16, pp, 0);
}

static void graph_skel (struct _topwin * tw, char * arg) {
	const char *s = arg, *str = 
	"[" TW_TOPH
	"({8i#in$2X<}{8o#out$2X>}{8f#fb$2X@}"
	"{Maadd...$Xiz$|_?0}{M0[$XB%g|_03}{M1sl$XB%g|_013}{M2]$XB%g|_023}"
	"{B_--->$XG%g+}{B_==$XG%g=}{ex25$XG%gv}3{__}0{BXX$XG%gX}{Ms++$X|g0})"
	"3{ggX}]";
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, str);
	if (s && *s) daclb_set(widg_lookup_ps(tw, "."), &s, 1);
	if (!tw->state) gtk_window_set_default_size(GTK_WINDOW(tw->w), 300, 400);
}

///////////////// pole/zero //////////////////////////////////////////////////

#define DAPZ_N(x) ((x)->arg[2].i[0])
#define DAPZ_A(x) ((x)->arg[2].i[1])

static void dapz_rad(cairo_t * cr2, double ang, double r0, double r1) {
	double si = sin(ang), co = cos(ang);
	cairo_move_to(cr2, 383.5 + r0*si, 383.5 + r0*co);
	cairo_line_to(cr2, 383.5 + r1*si, 383.5 + r1*co); }

static cairo_surface_t * dapz_lzs(ww_t * ww) {
	LAZYSURF_HEAD(767, 767);  int i,j;
	cairo_set_source_rgb (cr2, .2, .2, .2);
	cairo_paint(cr2);
	cairo_set_source_rgb (cr2, .1, .1, .1); 
	cairo_set_line_width (cr2, 1.0);
	for (i=0; i<543; i+=30) for (j=6; j<30; j+=6) cairo_new_sub_path(cr2),
		cairo_arc(cr2, 383.5, 383.5, (double)(i+j), 0.0, 2.0*M_PI);
	double ang1 = M_PI / 32.0;
	for (i=0; i<64; i+=4) {
		double ang = (double)(i+2)*ang1;
		dapz_rad(cr2, ang, 30.0, 545.0);
		dapz_rad(cr2, ang-ang1, 60.0, 545.0);
		dapz_rad(cr2, ang+ang1, 60.0, 545.0);
	}
	cairo_stroke(cr2); cairo_set_source_rgb (cr2, .0, .0, .0); 
	for (i=30; i<543; i+=30) cairo_new_sub_path(cr2), 
		cairo_arc(cr2, 383.5, 383.5, (double)i, 0.0, 2.0*M_PI);
	for (i=0; i<8; i++) dapz_rad(cr2, .125*M_PI*(double)i, -545.0, 545.0);
        cairo_stroke(cr2); cairo_destroy(cr2); return p;
}

static void dapz_draw (ww_t * ww, cairo_t * cr2) {
	DA_XY01; int i, j, xi=0, yi=0, n2 = 2*DAPZ_N(ww), jp = 0, jz = n2, zf=0;
	const short * p = (short*) ww->etc; // 768-yi<y1+9  -yi<y1-759 yi>759-y1
	double x[n2], y[n2];              // 768-yi>y0-9  -yi>y0-777 yi<777-y9
	for (i=0; i<n2; i++) if ((i&1) ? (yi=768-yi) : (xi=p[i], zf=xi&2048, xi&=1023, yi=p[i+1]&1023),
				 (xi<x1+9 && xi>x0-9) && (yi<y1+9 && yi>y0-9))
					j = zf ? --jz : jp++, x[j]=(double)xi+.5, y[j]=(double)yi+.5;
	LOG("pz/draw: (%d+%d)/%d", jp, n2-jz, n2);
	cairo_set_line_width (cr2, 1.0); cairo_set_source_rgb (cr2, .0, .7, 1.0); for (i=0; i<jp; i++) 
		cairo_move_to(cr2, x[i]-5.0, y[i]-5.0), cairo_line_to(cr2, x[i]+5.0, y[i]+5.0),
		cairo_move_to(cr2, x[i]-5.0, y[i]+5.0), cairo_line_to(cr2, x[i]+5.0, y[i]-5.0);
	cairo_stroke(cr2); cairo_set_source_rgb (cr2, 1.0, .2, .2); for (i=jz; i<n2; i++)
		cairo_new_sub_path(cr2), cairo_arc(cr2, x[i], y[i], 6.0, 0.0, 2.0*M_PI);
	cairo_stroke(cr2);
}

static void dapz_skel(struct _ww_t * ww, const char **pp) {
	da_skel(ww, 767, 767);
	DAPZ_N(ww) = DAPZ_A(ww) = 0;
}

static void dapz_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return dapz_draw (ww, (cairo_t*)ww->arg[0].p);
	int i, k, n = 0;
	LOG("pz_cmd: arg = \"%s\"", arg);
	while (arg[0] && arg[1] && arg[2] && arg[3]) arg +=4, n++;
	if ((k = DAPZ_A(ww)) < n) {
		if (!k) k = 4; else free(ww->etc);
		do k += k; while (k<n);
		ww->etc = malloc(4*(DAPZ_A(ww)=k));
	}
	LOG("pz_cmd: n=%d, k=%d", n, k);
	DAPZ_N(ww) = n; arg -= 4*n;
	short * p = (short*) ww->etc;
	for (i=0; i<n; i++) p[2*i  ] = 64 * arg[4*i  ] + arg[4*i+1] - 3120,
			    p[2*i+1] = 64 * arg[4*i+2] + arg[4*i+3] - 3120;
	da_fullre(ww);
}

static void pz_skel (struct _topwin * tw, char * arg) {
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, "{PP}");
}

///////////////// errmsg /////////////////////////////////////////////////////

GtkWidget * err_vbl (struct _ww_t * ww, int ix) { return parse_w_s(ww->top, "{E0}"); }

static int err_nl = 10, err_red = 31, err_annoy = 1, err_annoy_t = 0;
static sthg err_tab[32];

static void local_error(int ec) {
	int buf[4]; buf[0]='>'<<24; buf[1]=qh4((ec>>16)&255); buf[2] = qh4(ec&65535); buf[3] = 0;
	topwin_mk(0x47, 'E', (char*)buf+3); }

static const char * daerr_txt(int i, int j) {
	static char buf[128]; if (!buf[2]) buf[2] = buf[5] = ':', buf[8] = ' ';
	time_t t = j; struct tm loct; localtime_r(&t, &loct);
	d59(buf, loct.tm_hour); d59(buf+3, loct.tm_min); d59(buf+6, loct.tm_sec);
	const char * s = (i&0x800000) ? err_str(i|0xfff00000) : strerror(i&0x7fffff);
	int k = 9, l = min_i(100, strlen(s)), c = ((i>>24)&127);
	if (c>1) { if (c>9) buf[k++]=48+(c/10), c%=10;  buf[k++] = 48+c; buf[k++]='x'; buf[k++]=' '; }
	memcpy(buf+k, s, l); buf[k+l] = 0; return buf;
}

static void err_upd_flg(struct _topwin * tw) {
	ww_t * ww = widg_lookup_pci(tw, 'A', -1); DABOOL(ww) = err_annoy; da_fullre(ww); }

static void err_upd_nl(struct _topwin * tw) {
	dacnt_set_x (widg_lookup_pci(tw, 'N', -1), err_nl, 256);
	vbox_show_bv(widg_lookup_pci(tw, 'E', -1), (1<<err_nl)-1); }

static void err_cmd (struct _topwin * tw, char * arg) {
	ww_t *p, *ww = widg_lookup_pci(tw, 'E', -1);
	int i, k, n, af = 0, wi = VB_WBASE(ww), nl = err_nl, red = err_red, t = time(NULL);
	switch(*arg) {
		case 0: return;
		case 'N': if (arg[1]=='+') nl += (4*arg[2] - 191);
			  else if (arg[1]=='-') nl -= (4*arg[2] - 191);
			  else nl = arg[1] - 48;
			  if (nl<3) nl = 3; else if (nl>30) nl = 30;
			  if (nl==err_nl) return;
			  if (nl>err_nl) for (i=err_nl; i<nl; i++) err_tab[i].i[0] = 0;
			  else if (nl<=red && red<31) memmove(err_tab, err_tab+(red+1-nl), 8*nl), red = nl-1;
			  break;
		case 'X': red = 31; for (i=0; i<nl; i++) err_tab[i].i[0] = 0; break;
		case 'A': err_annoy ^= 1; err_upd_flg(tw); break;
		case '>': ++arg; for (n=0; arg[n]; n++);
			  for (i=0; i<n; i+=8) if ((k = 65536*qh4rs(arg+i)+qh4rs(arg+4+i))!=0x01FFFFEA)
				  ++red, red&=-(red<nl), err_tab[red].i[0]=k,err_tab[red].i[1]=t; else af = 1;
			  break;
		default: LOG("err_cmd: undef: 0x%x '%c'", *arg, *arg); break;
	}
	if (red!=err_red) err_tab[err_red].i[0] &= ~0x80000000, err_tab[err_red=red].i[0] |= 0x80000000;
	if (nl!=err_nl) err_nl = nl, err_upd_nl(tw);
	for (i=0; i<nl; i++) if (p = widg_p(tw, wi+i), memcmp(p->arg+2, err_tab+i, 8))
						       memcpy(p->arg+2, err_tab+i, 8), da_fullre(p);
	if (af && !err_annoy) gtk_window_present(GTK_WINDOW(tw->w));
}

static void err_skel (struct _topwin * tw, char * arg) {
	if (!tw->state) {
		tw->arg[0].p = parse_w_s(tw,
		"[({B_clear$$>X}{8N#L$2$>N}{YAannoy$$>A}3()0{B_console$_c-1}{B_?$$?}){:EeN10}]");
		memcpy(tw->title, "Errors\0", 8); 
		err_upd_nl(tw); err_upd_flg(tw);
	} else if (err_annoy) {
		int t = time(NULL); 
		if (t!=err_annoy_t) err_annoy_t = t, gtk_window_present(GTK_WINDOW(tw->w));
	}
	if (arg) err_cmd(tw, arg);
}

///////////////// main config ////////////////////////////////////////////////

static void mcfg_skel (struct _topwin * tw, char * arg) {
	static const char *d[7], *v[6] = {"HOME", "LF_USERDIR", "LF_TMPROOT", "LF_DIR", "LF_TMPDIR", "LF_TLOG"};
	if (!d[0]) {
		int i,l; char *s; for (i=0; i<6; i++) if (!(d[i] = getenv(v[i]))) d[i] = "<<BUG!!!>>";
		l = strlen(d[1]); s = malloc(l+10); memcpy(s, d[1], l); memcpy(s+l, "/__asv.lf", 10); d[6]=s;
	}
	const char *ws="[(3[{C7,?}{C6,?}{C5,?}]0[{Yssv.exec$cs}{8amin/asv$2ca}{Ytau.tlog$ct}]"
	"[{8Ssv.bkup$2cS}{8Aas.bkup$2cA}{8Ttl.bkup$2cT}][3{B_?$$?}0{Yddev$cd}])"
	"([{B_wav-dir$$!fwW}{B_atmpdir$$!fkW}{L_homedir}{L_userdir}{L_tmp.dir}{L_instdir}{L_workdir}]"
	"3[(3[{Cw,/1234/6789/1234/6789/1234/6789/1234/6789/1234/6789}{Ck,?}]0[{MK$cK|k0}{8Llim$3cL}])"
	"{C0,?}{C1,?}{C2,?}{C3,?}{C4,?}]){B>saveCfg$c>}]";
	if (!tw->state) { 
		tw->arg[0].p = parse_w_s(tw, ws);
		memcpy(tw->title, "config", 6);
		int i; const char *s;
		for (i=0; i<7; i++) s = d[i], daclb_set(widg_lookup_pci(tw, 48+i, 0), &s, 3);
	} else {  gtk_window_present(GTK_WINDOW    (tw->w)); }
}

///////////////// audio config ///////////////////////////////////////////////

static void acfg_skel (struct _topwin * tw, char * arg) {
	const char *ws="[{CC%%%XXX(no audio output)}({!sspd$163A0s}{!rrsv$1faA0r}{!ttry$114A0t}{!wt/w$1faA0w}"
	"[{B_?$$?win.audio}{M_name:$A0n|_01}{en10$A0N}{M_chan.c:$A0o|S2}{eo10$A0O}{L##out: 0}{L_clock:}"
	"{Mc$A0c|S1}(3{Y0kill PA$A00}{Y1-9$A01}){B_restart$A0R}{B_saveCfg$_K}])]";
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, ws), memcpy(tw->title, "audio", 6);
	else gtk_window_present(GTK_WINDOW(tw->w)); 
}

///////////////// input box //////////////////////////////////////////////////

static void in01_skel (struct _topwin * tw, char * arg) {
	if (tw->state) return LOG("BUG: create called again for input win 0x%x", tw->id);
	int i, l = arg ? strlen(arg-1) : 0, n = arg ? ivlim(*arg-48, 1, 16) : 1;
	sprintf(tw->cmdpref, "%04x", (tw->id>>4)&65535); tw->cmdp_len = 4;
	char ws[256], *p = arg + 1, *q = ws+1;
	memcpy(tw->title,"input box", 12);
	ws[0] = '('; 
	if (l>15) hxdoub_str(tw->title, p, 15), p+=16, l-=16;
	for (i=0; i<n; i++, q+=14) {
		q[0] = '{'; q[1] = '!';
		memcpy(q+3, (l>15) ? (p+=16, l-=16, hxdoub_lbl(p-16)) : "---", 3);
		memcpy(q+6, "$564L0", 6); q[12] = q[2] = hexc1(i);
		q[13] = '}';
	}
	q[0] = ')'; q[1] = 0;
	tw->arg[0].p = parse_w_s(tw, ws);
	int ni = min_i(n, l>>1);
	for (i=0; i<ni; i++) LOG("setarg %d %d", i, hex2(p+2*i)), dacnt_set_x(widg_lookup_pci(tw, hexc1(i), 0), hex2(p+2*i), 512);
}

///////////////// audio file dialog //////////////////////////////////////////

static void a20_skel (struct _topwin * tw, char * arg) {
	const char * ws = "({B_>>wav$w3}{B_>>flac$w7}3[({L_skip:}3{eF6}0{L_len:}3{eL6})" 
	"3(3{B_keep$wf}{B_wav$w0$ %F$ %L}{B_flac$w4$ %F$ %L}{B_del$wd}{B_config$~cW})]0{B_?$$?})";
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, ws); 
	else   gtk_window_present(GTK_WINDOW(tw->w)); 
	int l=0; if (arg && (l=strlen(arg)) == 8) {
		memcpy(tw->cmdpref, arg, 8); tw->cmdpref[8] = '$'; tw->cmdp_len = 9;
		int i; for (i=0; i<8; i++) tw->title[i] = arg[i]|32; tw->title[8] = 0;
	} else { LOG("auconv BUG: l=%d", l); }
}
///////////////// calculator /////////////////////////////////////////////////

GtkWidget * calc_vbl (struct _ww_t * ww, int ix) {
	topwin * tw = ww->top;
	GtkWidget * rw = parse_w_s(tw, "({L0w88:}3{e130$Xx}0{L2W}");
	ww_t * lbl = widg_p(tw, VB_WBASE(ww) + 3*ix);
	char * s = DALBL_TXT(lbl);
	s[0] = 'z'-VB_ARG(ww);
	if (ix<10)      s[1] = 48+ix, s[2] = ':', s[3] = 0;
	else if (ix<20) s[1] = 49, s[2] = 38+ix;
	else		s[1] = 50, s[2] = 28+ix;
	ww_t * ent = widg_p(tw, VB_WBASE(ww) + 3*ix + 1);
	char * p = ent->cmd; 
	p[1] = s[0]; p[2] = s[1]; 
	if (s[2]==':') p[3] = 36, p[4] = 0; else p[3] = s[2], p[4] = 36, p[5] = 0;
        return rw;
}

static void calc_skel (struct _topwin * tw, char * arg) {
	const char * str = "[" TW_TOPH
		"({8x#in$2XX}{8z#tmp$2XZ}{8y#out$2XY}){:ZcN30}{:YcN31}]";
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, str);
}

///////////////// iter. box //////////////////////////////////////////////////

static void itb_skel (struct _topwin * tw, char * arg) {
	const char * str = "[" TW_TOPH
		"({M0[$Xb|_03}{M1B$Xb|_013}{M2]$Xb|_023}3{C3280$1$eeeeee333333(...bx1...)})"
		"({L_scale:}{Ms$Xs|i0}{Yrz1=z0zr$Xr})]";
	if (!tw->state) tw->arg[0].p = parse_w_s(tw, str);
}

///////////////// clipboard //////////////////////////////////////////////////

#define DACLIP_SEL(x) ((x)->arg[2].c[0])

static void daclip_inv(ww_t * ww, unsigned int bv) {
	if (!bv || !ww->w) return; if (bitcnt(bv)>2) return da_fullre(ww);
	GdkWindow * w = gtk_widget_get_window(ww->w);
	struct _GdkRectangle rct; rct.width = 29; rct.height = 57;
	BVFOR_JMC(bv) rct.x = 28*(j&7), rct.y = 7*(j&24), gdk_window_invalidate_rect(w, &rct, 0);
}

static unsigned int daclip_sel(ww_t * ww, int i) {
	int i0 = DACLIP_SEL(ww); return (i0==i) ? 0u : (DACLIP_SEL(ww)=i, (1u<<i)|(1u<<i0)); }

static void daclip_lb32(cairo_t * cr2, int i) {
	int x0 = 28*(i&7)+1, y0 = 7*(i&24)+1, ch = i + (i<10?48:55);
	cairo_rectangle(cr2, (double)x0, (double)y0, 14.0, 18.0); cairo_fill(cr2);
	cairo_set_source_rgb(cr2, .0, .0, .0); tx_box(cr2, x0, y0, 12, 16, conf_lbfs, (const char*)&ch); }

static cairo_surface_t * daclip_lzs(ww_t * ww) {
	LAZYSURF_HEAD(225, 225); int i;
	cairo_set_antialias(cr2, CAIRO_ANTIALIAS_NONE);
	cairo_set_source_rgb(cr2, .2, .2, .2); cairo_paint(cr2);
	cairo_set_line_width (cr2, 1.0); cairo_set_source_rgb(cr2, .8, .8, .8);
	for (i=0; i<9; i++) cairo_move_to(cr2, .5+(double)(28*i), .0), cairo_rel_line_to(cr2, .0, 225.0);
	for (i=0; i<5; i++) cairo_move_to(cr2, .0, .5+(double)(56*i)), cairo_rel_line_to(cr2, 225.0, .0);
	cairo_stroke(cr2);
	for (i=0; i<32; i++) cairo_set_source_rgb(cr2, .8, .8, .8), daclip_lb32(cr2, i);
	cairo_destroy(cr2); return p;
}

static void daclip_draw (ww_t * ww, cairo_t * cr2) {
	DA_XY01; cairo_set_line_width (cr2, 1.0);
	int i,f,sel=DACLIP_SEL(ww); for (i=0; i<32; i++) {
		char * s = ((char*)ww->etc) + 14*i;
		if (!(f = (s[2]!=48) + 2*(i==sel))) continue;
		int wx0 = 28*(i&7)+10; if (wx0>x1 || wx0<x0-19) continue;
		int wy0 = 7*(i&24)+12; if (wy0>y1 || wy0<y0-45) continue;
		if (f&2) cairo_set_source_rgb(cr2, 1.0, .0, .0), 
			 cairo_rectangle(cr2,(double)wx0-9.5,(double)wy0-11.5, 28.0,56.0), cairo_stroke(cr2),
			 daclip_lb32(cr2, i);
		if (f&1) da_wr18(cr2, wx0, wy0, s, 18);
	}}

static int daclip_get(void * to, struct _ww_t * ww, int ty) {
	return *(char*)to = i_to_b32(DACLIP_SEL(ww)), 1; }

static void daclip_cmd(struct _ww_t * ww, const char * arg) {
	if (!arg) return daclip_draw(ww, (cairo_t*)ww->arg[0].p);
	int k; unsigned int re = 0;
	while (*arg) {
		if ((k=b32_to_i(*arg))>=0) {
			if (arg[3]=='0') arg += 4, ((char*)ww->etc)[14*k+2] = '0'; 
			else memcpy((char*)ww->etc+14*k, arg+1, 14), arg += 15;
			re |= (1u<<k);
		} else if (*arg=='+') {
			if ((k=b32_to_i(arg[1]))>=0) { re |= daclip_sel(ww, k), arg += 2; }
			else { LOG("daclip/+: invalid ix 0x%x", arg[1]); break; }
		} else { LOG("daclip/+: invalid c0 0x%x", *arg); break; } }
	daclip_inv(ww, re);
}

static void daclip_clk(struct _ww_t * ww, int b9, int cx, int cy, GdkEventButton * ev) {
	char buf[16]; buf[0] = '~'; buf[1] = 'K';
	int i = 2 + tw_defcmd(ww->top, buf+2);
	if (b9<0) buf[i] = 'K'+32*(b9&1), buf[i+1] = 48+(cx>>4), buf[i+2] = hexc1(cx&15), i += 3;
	else buf[i] = 48+b9, buf[i+1] = 48+8*ivlim((cy-2)/56,0,3)+ivlim((cx-2)/28,0,7), i += 2;
	buf[i++] = 10; write(1, buf, i); 
}

static void daclip_skel(struct _ww_t * ww, const char **pp) {
	da_skel(ww, 226, 226);
	char * s = ww -> etc = malloc(448);
	DACLIP_SEL(ww) = 33;
	int i; for (i=0; i<32; i++) s[14*i+2] = '0';
}

static void clip_setflg(struct _topwin * tw, int flg) {
	const char * s;
	for (s="DACPX"; *s; s++,flg>>=1) {
		char nm[2]; nm[0] = *s; nm[1] = 0;
		ww_t * ww = widg_lookup_ps(tw, nm);
		if (!ww) { LOG("clip_setflg: %c not found"); return; }
		ww->arg[3].c[0] = (flg&1);
		da_fullre(ww);
	}}

GtkWidget * clip_vbl (struct _ww_t * ww, int ix) { return parse_w_s(ww->top, ix?"({L0sorry}{L1not yet}{C2,implemented})":"([])"); }

static void clip_skel (struct _topwin * tw, char * arg) {
	const char * str = "[{C.228$16$kkk000...}"
		"(3{YCcp$KC}{YPps$KP}{YD2x$KD}{YAau$KA}{YXxc$KX}{Brre$KW}{Ma+$|K0})"
		"{KK}{:WK230}]";
	ww_t * cl;
	if (tw->state) {
		cl = widg_lookup_ps(tw, "K"); if (!cl) {
			LOG("clip_skel: da not found"); return; }
		int i; for (i=0; i<32; i++) ((char*)cl->etc)[14*i+2] = 48;
	} else {
		tw->arg[0].p = parse_w(tw, &str);
		cl = widg_lookup_ps(tw, "K");
	}
	const char *s = arg;
	ww_t *ww = widg_lookup_pci(tw, 'W', -1); vbox_show_bv(ww, 1);
	if (!s || !*s || (daclb_set(widg_lookup_ps(tw, "."), &s, 1), !*s)) DACLIP_SEL(cl) = 0;
	else if ((DACLIP_SEL(cl)=b32_to_i(*s), *++s) && (clip_setflg(tw, *s - 48), *++s)) daclip_cmd(cl, s);
}

///////////////// main window ////////////////////////////////////////////////

static void mwin_skel (struct _topwin * tw, char * arg) {
	const char * str = "[{CN100$20$%%%uuu/1234/67890123/56789}3({!vvol$163vG}3[3{B_kussb+$r}"
	"(3{B_unplot$_g}{B_?$$?})({22}3[3{B_LR$$/}{Mm++$|.0}]){Bssave$s}{YWrec$_w}0{LV v??.??}])]";
        tw->arg[0].p = parse_w(tw, &str);
	memcpy(tw->title, "lf\0", 4);
}

///////////////// tree view //////////////////////////////////////////////////

static GtkTreeStore * tstore[2];
static GtkTreeView * tview[2];
static GtkWidget * tview_sc[2];
static GtkTreeIter troot[2];
static GtkTreeIter tsel[2];
#define FORD(i) (GTK_TREE_MODEL(tstore[i]))

static const char l_panel[] = "[{__}{CG/0$0$eeeeeeff0000hee}{__}{eE24}{__}"
"(3{BR|rnm$N@l$m%E}{BLlcp$N@l$c%E}{MY(W)$##|?1}{BC>+>$N@l$c@R.%E}{BM>>>$N@l$m@R.%E}{BDX$N@l$d}{MTnew$C@L$?%E$|/0})]";
static char r_panel[sizeof(l_panel)];

#define TVIEW_IX(i, nm) int i=(tw==tview[1]); if (!i && tw!=tview[0]) { LOG("%s: invalid treeview", #nm); }
#define GVZ(nm) GValue nm ; memset(& nm , 0, sizeof(GValue))

static ob_dir * tree_lookup (int id, int force);

static int tree_nd_id(int lr, GtkTreeIter *it) {
	GVZ(val); gtk_tree_model_get_value(FORD(lr), it, 1, &val);
	return g_value_get_int(&val);
}

static const char * tree_nd_name(int lr, GtkTreeIter *it) {
	GVZ(val); gtk_tree_model_get_value(FORD(lr), it, 0, &val);
	return g_value_get_string(&val);
}

static int tree_recdel(int lr, ob_dir *p, GtkTreeIter * it, int flg, char ** pp) { // 1:keep0th("...")
	int r = 0;
	if (dflg & DF_NCOLL) LOG("recdel: %x, %x", p?p->id:0xdeadbeef, tree_nd_id(lr, it));
	if (p) **pp = 36, ++*pp, *pp += hx5(*pp, p->id>>4);
	if (p && !(p->flg & (lr+1))) p = NULL;
	GtkTreeIter it2;
	GtkTreeModel * mdl = FORD(lr);
	if (!gtk_tree_model_iter_nth_child(mdl, &it2, it, 0)) return 0;
	dirtab_del(p, lr+1);
	while (1) {
		int id2 = tree_nd_id(lr, &it2);
		ob_dir * q = (id2>=0) ? tree_lookup(16 * id2 + lr + 1, 0) : NULL;
		r += tree_recdel(lr, q, &it2, 0, pp) + 1;
		if (flg&1) {
			gtk_tree_store_set(tstore[lr], &it2, 0, "...", 1, -1, -1);
			if (!gtk_tree_model_iter_next(mdl, &it2)) return r;
			flg = 0;
		} else {
			if (!gtk_tree_store_remove(tstore[lr], &it2)) return r;
		}
	}
}

#define NDC_REC(S) if (dflg&DF_REC) CMD("QRN`^%05x`%c$" S "`%05x`", id, 49+lr, id);
static void node_expand(GtkTreeView *tw, GtkTreeIter *it, GtkTreePath *path, gpointer _) {
	TVIEW_IX(lr, node_expand); int id = tree_nd_id(lr, it); NDC_REC("<");
	if (dflg & DF_NEXP) LOG("expand: %c %x, %c", "LR"[lr], id, 48+expand_cmd_flg);
	if (expand_cmd_flg) CMD("D#%x$%c", id, 32*lr+'E');
}

static void node_coll2(int lr, ob_dir * p, GtkTreeIter * it, int keepflg) {
	if (dflg & DF_NCOLL) LOG("collapse: %c %x p%p", "LR"[lr], p->id, p);
	char s0[1544], *s = s0 + 3;
	int r = -999;
	if (p->flg&(lr+1)) {
		s0[0] = '~'; s0[1] = 'x'; s0[2] = 49+lr;
		r = tree_recdel(lr, p, it, keepflg, &s);
		*(s++) = '\n'; write(1, s0, s-s0);
	}
	dirtab_del(p, lr+1);
	if (dflg & DF_NCOLL) LOG("collapse: %c 0x%x --> %d nodes deleted", "LR"[lr], p->id, r);
}

static void node_collapse(GtkTreeView *tw, GtkTreeIter *it, GtkTreePath *path, gpointer _) {
	TVIEW_IX(lr, node_collapse);
	int id = tree_nd_id(lr, it); NDC_REC(">");
	ob_dir * p = tree_lookup(16*id+lr+1, 0);
	if (p) node_coll2(lr, p, it, 1);
}

static void node_sel(GtkTreeView *tw, gpointer user_data)  {
	TVIEW_IX(lr, node_sel);
	GtkTreePath * path; gtk_tree_view_get_cursor(tw, &path, NULL);
	if (!gtk_tree_model_get_iter(FORD(lr), tsel+lr, path)) return LOG("node_sel: lookup failed");
	int id = 0xfffff&tree_nd_id(lr, tsel+lr);  NDC_REC("*");
	CMD("N#%x$%c", 0xfffff&tree_nd_id(lr, tsel+lr), "LR"[lr]);
}

static void node_act(GtkTreeView *tw, GtkTreePath *path, GtkTreeViewColumn * col, gpointer _) {
	TVIEW_IX(lr, node_act); int id = 0xfffff&tree_nd_id(lr, tsel+lr); 
	if (!gtk_tree_model_get_iter(FORD(lr), tsel+lr, path)) return LOG("node_act: lookup failed");
	NDC_REC("W"); CMD("D#%x$%c", id, 32*lr+'W');
}

static void mk_subnd(GtkTreeIter *to, int lr, GtkTreeIter *nd, int ix, int id, const char * nm) {
	static GtkTreeIter tmpit; if (!to) to = &tmpit;
	gtk_tree_store_insert(tstore[lr], to, nd, ix);
	gtk_tree_store_set(tstore[lr], to, 0, nm, 1, id, -1);
}

static void mk_tree_2(topwin * tw) { int i; for (i=0; i<2; i++) {
	tview[i] = GTK_TREE_VIEW(gtk_tree_view_new());
	gtk_tree_view_set_enable_tree_lines(tview[i], TRUE);
	tstore[i] = gtk_tree_store_new(2, G_TYPE_STRING, G_TYPE_INT);
	gtk_tree_store_insert(tstore[i], troot+i, NULL, 0);
	gtk_tree_store_set(tstore[i], troot+i, 0, "/", 1, 0, -1);
	mk_subnd(NULL, i, troot+i, 0, 1<<20, "...");
	memcpy(tsel+i, troot+i, sizeof(GtkTreeIter));
	gtk_tree_view_set_model(tview[i], GTK_TREE_MODEL(tstore[i])); 
	GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
	GtkTreeViewColumn *col = gtk_tree_view_column_new_with_attributes(
			8*i+"tree(L)\0tree(R)", renderer, "text", 0, NULL);
	gtk_tree_view_column_set_widget(col, parse_w_s(tw, "()"));
	gtk_tree_view_append_column(GTK_TREE_VIEW(tview[i]), col);
	g_signal_connect(tview[i], "row-expanded", G_CALLBACK(node_expand), NULL);
	g_signal_connect(tview[i], "row-collapsed", G_CALLBACK(node_collapse), NULL);
	g_signal_connect(tview[i], "row-activated", G_CALLBACK(node_act), NULL);
	g_signal_connect(tview[i], "cursor-changed", G_CALLBACK(node_sel), NULL);
	tview_sc[i] = gtk_scrolled_window_new(NULL, NULL);
	gtk_container_add(GTK_CONTAINER(tview_sc[i]), GTK_WIDGET(tview[i]));
	gtk_widget_show(GTK_WIDGET(tview[i])); gtk_widget_show(tview_sc[i]);
}}

static void lr_swap(char * p) {
	switch (*p) {
		case 'L': *p = 'R'; return;
		case 'R': *p = 'L'; return;
		case 'l': *p = 'r'; return;
		case 'r': *p = 'l'; return;
		default: return;
	}}

GtkWidget * t2_widg(int t) {
	const char * s; switch(t) {
		case 'L': return tview_sc[0];
		case 'R': return tview_sc[1];
		case 'l': return s=l_panel, parse_w(ot_etc+1,&s);
		case 'r': return s=r_panel, parse_w(ot_etc+1,&s);
		default: return OOPS;
	}
}

static const char * lastname(const char * s, int l) {
	if (l<2) return s;
	for (--l; l && s[l]&&(s[l]-46); l--);
	return s+l+1;
}

static void t2sel_upd(int lr, int nid, int ty, const char * rgb, const char * s) {
	int l = strlen(s)+1;
	ww_t * ww = widg_lookup_ps(ot_etc+1, "Y\0y\0"+2*lr);
	DALBL_TXT(ww)[1] = ty; da_fullre(ww);
	ww = widg_lookup_ps(ot_etc+1, "G\0g\0"+2*lr);
	if (ww->etc) free(ww->etc);
	ww->etc = malloc(l); memcpy(ww->etc, s, l);
	if (rgb && *rgb)  memcpy(ww->arg[2].c, rgb, 6), ww->arg[2].c[6]=',';
	else if (ty=='d') memcpy(ww->arg[2].c, "uuu000,", 8);
	else if (ty=='k') memcpy(ww->arg[2].c, "uu%000,", 8);
	else memcpy(ww->arg[2].c, "uuu%%z,", 8);
	da_fullre(ww);
	const char *s1 = lastname(s, l-1), *s2 = s1;
	char *q = lr_dirname + 24*lr;
	entry_set(widg_lookup_ps(ot_etc+1, "E\0e"+2*lr), s1);
	if (ty!='d'&&ty!='k') s2 = lastname(s, s2-s-1);	
	while(*s2&&(*s2-'.')) *(q++) = *(s2++);   *q = 0;
	return;
}

static void t2win_skel (struct _topwin * tw, char * arg) {
	static int iflg = 0; if (!iflg) { mk_tree_2(tw); iflg = 1; }
	const char * str = "[(3{C1,(here comes the filename)}0{M2+$|/1}{B_?$W.!b.?.win.tree})"
	"{__2}3(3[3TL0Tl]0{__4}3[3TR0Tr])]";
	memcpy(r_panel, l_panel, sizeof(l_panel));
	int i; for (i=0; i<sizeof(l_panel)-1; i++) {
		switch(r_panel[i]) {
			case '@': ++i; lr_swap(r_panel+i); break;
			case '{': ++i;
			case '%': ++i; if (r_panel[i]!='_') r_panel[i] |= 32; break;
			case '>': r_panel[i] = '<'; break;
		}
	}
	tw->arg[0].p = parse_w(tw, &str);
	gtk_window_set_default_size(GTK_WINDOW(tw->w), 320, 480);
	memcpy(tw->title, "LR\0", 4);
	t2sel_upd(0, 0, 'd', NULL, ".");
	t2sel_upd(1, 0, 'd', NULL, ".");
}

static ob_dir * tree_lookup (int id, int force) {
	int lr = (id&15) - 1; if (lr&~1) return NULL;
	int k = obidtab_lookup(&oi_dir, id, force&1); if (k<0) return NULL;
	ob_dir * p = ot_dir + (k & OI_ID);
	if (k&OI_NEW) ob_dir_ini(p, id&0xfffff0);
	if (force&1) p->flg |= (lr+1);
	return p;
}

static GtkTreeIter * tree_lookup_it (int id, int f) {
	ob_dir * p = tree_lookup(id, f);
	return (p && (p->flg&id&3)) ? p->nd + (id&3) - 1 : NULL;
}

// -1:not found -2:empty 64: only one
static int tree_find_child(GtkTreeIter * trg, int lr, GtkTreeIter * par, int id) {
	GtkTreeModel * mdl = FORD(lr); id>>=4;
	if (!gtk_tree_model_iter_nth_child(mdl, trg, par, 0)) return -1;
	if (tree_nd_id(lr, trg) & (1<<20)) return -2;
	int i; for (i=0; 1; i++) {
		if (tree_nd_id(lr, trg) == id)
			return (i||gtk_tree_model_iter_n_children(mdl, par)>1) ? i : 64;
		if ( !gtk_tree_model_iter_next(mdl, trg)) return -1;
	}
}

static void tree_rmch(int lr, GtkTreeIter * nd, int i0) {
	GtkTreeIter it2;
	if (gtk_tree_model_iter_nth_child(FORD(lr), &it2, nd, i0))
		do {} while (gtk_tree_store_remove(tstore[lr], &it2));
}

static void tree_fill_node(GtkTreeIter * nd, int lr, const char * s, int xp) {
	GtkTreeModel * mdl = FORD(lr);
	int i = 0, j, id2;
	char nmbuf[64];
	GtkTreeIter it2;
	int cflg = gtk_tree_model_iter_nth_child(mdl, &it2, nd, 0);
	int dirflg = 1;
	while (1) {
		if (*s=='x') dirflg = 0, ++s;
		if (!*s) {
			if (i) break;
			memcpy(nmbuf, "(...empty)", 11);
			id2 = 1048576|tree_nd_id(lr, nd); dirflg = 0;
		} else {
			id2 = 0; while (*s & 80) id2<<=4, id2 += hxd2i(*s), s++;
			if (!id2) break; if (*s == '$') ++s;
			j=0; while (*s && *s!=36 && j<63) nmbuf[j++] = *(s++);
			nmbuf[j] = 0; if (*s == '$') ++s;
		}
		if (cflg) tree_rmch(lr, &it2, 0);
		else gtk_tree_store_insert_before(tstore[lr], &it2, nd, NULL);
		gtk_tree_store_set(tstore[lr], &it2, 0, nmbuf, 1, id2, -1);
		if (dirflg) mk_subnd(NULL, lr, &it2, 0, id2|1048576, "...");
		if (cflg) cflg = gtk_tree_model_iter_next(mdl, &it2);
		i++;
	}
	if (cflg) tree_rmch(lr, nd, i);
	if (xp) {
		GtkTreePath * path = gtk_tree_model_get_path(mdl, nd);
		if (!gtk_tree_view_row_expanded(tview[lr], path)) {
			expand_cmd_flg = 0;
			gtk_tree_view_expand_row(tview[lr], path, FALSE);
			expand_cmd_flg = 1;
		}
		gtk_tree_path_free(path);
	}
}

static void tree_ins1(int lr, ob_dir * oe, int id, int dflg, const char *s) {
	GtkTreeIter it1, it2;
	GtkTreeModel * mdl = FORD(lr);
	int itf = gtk_tree_model_iter_nth_child(mdl, &it1, oe->nd+lr, 0);
	for (; itf; itf=gtk_tree_model_iter_next(mdl, &it1)) {
		int df2 = gtk_tree_model_iter_nth_child(mdl, &it2, &it1, 0);
		if (dflg && !df2) break;
		if (!dflg && df2) continue;
		int k = strcmp(tree_nd_name(lr, &it1), s);
		if (!k) LOG("ins1(%x): name \"%s\" already exists", id, s);
		if (k>0) break;
	}
	gtk_tree_store_insert_before(tstore[lr], &it2, oe->nd+lr, itf ? &it1 : NULL);
	gtk_tree_store_set(tstore[lr], &it2, 0, s, 1, id, -1);
	if (dflg) mk_subnd(NULL, lr, &it2, 0, id|1048576, "...");
}

static void dirtab_debug() {
	obidtab_debug(&oi_dir);
	int i,j; for (j=0; j<oi_dir.n; j++) {
		i = oi_dir.otab[j] & 255;
		int flg = ot_dir[i].flg;
		fprintf(stderr, " (id:%x flg:%x", ot_dir[i].id, flg);
		if (flg&1) fprintf(stderr, " L:%x", tree_nd_id(0, ot_dir[i].nd));
		if (flg&2) fprintf(stderr, " R:%x", tree_nd_id(1, ot_dir[i].nd+1));
		fprintf(stderr, ")");
	}
	fprintf(stderr, "\n"); fflush(stderr);
}


///////////////// command ////////////////////////////////////////////////////

static void ww_gencmd(ww_t * ww, const char *arg) { int c,k, cl=ww->cl->ch; switch(*arg) {
	case '?': return (void) ww_debug(ww, 1);
	case 'c': { ww_clk_fun cf = ww->cl->clk; int xy = (arg[1] && arg[2]==',') ? atoi_h(arg+3) : 0;
		    return cf ? (*cf)(ww, arg[1]-48, xy&4095, xy>>12, NULL) 
			      : LOG("ww_gencmd/c/'%c': no clkfun", cl); }
	case 'e': return (cl=='e') ? (void) gtk_entry_set_text(GTK_ENTRY(ww->w), arg+1)
		  	           : LOG("ww_gencmd/e/'%c': entry exp.", cl);
	case 'k': return ((cl|8)==43) ? (k=atoi_h(arg+1), dagrid_kcmd(ww, (k>>16)|32, k&255, (k>>8)&255))
		  	              : LOG("ww_gencmd/k/'%c': grid exp.", cl);
	case 'g': return (cl=='g') ? (c=arg[1], k=c?atoi_h(arg+2):0, graph_sendcmd(ww, c, k>>8, k&255))
		  	           : LOG("ww_gencmd/e/'%c': entry exp.", cl);
	default:  LOG("ww_gencmd: invalid arg \"%s\" for '%c'", arg, cl);
}}

static void cmd_tree(int id, char * s) {
	GtkTreeIter * it = tree_lookup_it(id, 0);
	if (!it) { LOG("T: lookup(0x%x) failed", id); return; }
	int lr = (id & 1) ^ 1, id2 = 0;
	while (*s & 80) id2<<=4, id2 += hxd2i(*s), s++;
	if (*s == 36) ++s;
	if ((id^id2)&15) { LOG("T: lr mismatch(0x%x,0x%x)", id, id2); return; }
	GtkTreeIter * it2 = tree_lookup_it (id2, 1);
	if (!it2) { LOG("TR: t_lu error(0x%x,0x%x)", id, id2); return; }
	if (!((id|id2)&0xfffff0)) {
		memcpy(it2, troot + lr, sizeof(GtkTreeIter));
	} else if (tree_find_child(it2, lr, it, id2) < 0) {
		LOG("T: t_find error(0x%x,0x%x)", id, id2); return; 
	}
	tree_fill_node(it2, lr, s, 1);
}

static void emp2nd(int lr, GtkTreeIter * it, const char * nm, int id, int dflg) {
	gtk_tree_store_set(tstore[lr], it, 0, nm, 1, id, -1);
	if (dflg) mk_subnd(NULL, lr, it, 0, id|1048576, "...");
}

typedef void (*path_fun) (GtkTreeView *v, GtkTreePath *p);
void cmd_pf_sel(GtkTreeView *v, GtkTreePath *p) { gtk_tree_view_set_cursor  (v, p, NULL, 0); }
void cmd_pf_exp(GtkTreeView *v, GtkTreePath *p) { gtk_tree_view_expand_row  (v, p, FALSE);   }
void cmd_pf_cls(GtkTreeView *v, GtkTreePath *p) { gtk_tree_view_collapse_row(v, p);          }
void cmd_pf_act(GtkTreeView *v, GtkTreePath *p) { gtk_tree_view_row_activated(v, p, 
							  gtk_tree_view_get_column(v, 0));   }

static void cmd_path(int lr, int id, int snid, GtkTreeIter * it, path_fun f) {
	if (!( (id>>4)|snid )) it = troot + lr; 
	else if (!it) return LOG("N*: %x.%x not found", id, snid);
	GtkTreePath * path = gtk_tree_model_get_path(FORD(lr), it);
	(*f)(tview[lr], path); gtk_tree_path_free(path); 
}

static void cmd_node(int id, char * s) {
	int k = obidtab_lookup(&oi_dir, id, 0);
	if (k<0) { LOG("N: node %x not found", id); return; }
	int lr = (id-1)&1, snid = 0;
	ob_dir * oe = ot_dir + k; int flg = oe->flg;
	/*LOG("oe->flg:%d", flg);*/ if (!(flg & (lr+1))) return;
	int c0 = *(s++);
	GtkTreeIter sn, *it = &sn;
	while (*s & 80) snid<<=4, snid += hxd2i(*s), s++;
	if (*s==36) ++s;
	int found = snid ? tree_find_child(it, lr, oe->nd+lr, 16*snid+1+lr) : -1, itflg = (found>=0);
	// LOG("cmd_node: itflg:%c s:\"%s\"", "ny"[itflg], s);
	switch(c0) {
		case '*': return cmd_path(lr, id, snid, itflg?it:0, cmd_pf_sel);
		case '<': return cmd_path(lr, id, snid, itflg?it:0, cmd_pf_exp);
		case '>': return cmd_path(lr, id, snid, itflg?it:0, cmd_pf_cls);
		case 'W': return cmd_path(lr, id, snid, itflg?it:0, cmd_pf_act);
		case '+':
		case ',':
			if (!snid) LOG("N%c: cannot rename root node", c0);
			else if (itflg) gtk_tree_store_set(tstore[lr], it, 0, s, -1);
			else if (found!=-2) tree_ins1(lr, oe, snid, c0-'+', s);
			else emp2nd(lr, it, s, snid, c0-'+');
			return;
		case '-':
			if (!snid) { LOG("N-: cannot delete root node"); return; }
			if (!itflg) { LOG("N-: subnode %x not found", snid); return; }
			k = obidtab_lookup(&oi_dir, 16*snid, 0);
			if (k>=0) node_coll2(lr, ot_dir+k, it, 0);
			if (found!=64) gtk_tree_store_remove(tstore[lr], it); 
			else gtk_tree_store_set(tstore[lr], it, 0, "(...empty)", 1, id|(1<<20), -1);
			return;
		default: LOG("N: unknown cmd 0x%x", c0); return;
	}
}

static void cmd_ob(int c0, int id, char * arg) {
	switch (c0) {
		case 'T': cmd_tree(id, arg); return; 
		case 'N': cmd_node(id, arg); return; 
		case 'C': topwin_mk(id, arg[0], arg+1); return;
		default: break;
	}
	topwin * tw; ww_t * ww;
	if (!(tw = tw_lookup(id))) { if ((dflg&DF_CWLU)|((id-0x57)&~16)) LOG("%c: lookup(0x%x) failed",c0,id);
				     return;}
	if (!*arg) { LOG("%c: missing twclass", c0); return; }
	if (*arg!='?' && *arg!=tw->cl->ch) {
		LOG("%c: exp:%c got:%c", c0, *arg, tw->cl->ch); return; }
	const char * s = arg + 1;
	switch(c0) {
		case 'W': case 'V':
			if (!(ww = wlu_any_pp(tw, &s))) LOG("%c: widget not found \"%s\"", c0, arg);
			else if (!ww->cl) LOG("%c: widget: zero class \"%s\"", c0, arg);
			else (*((c0&1)?ww->cl->cmd:&ww_gencmd)) (ww, s);
			return;
		case 'U': return (*tw->cl->cmd)(tw, arg+1);
		case 'P': gtk_window_present(GTK_WINDOW(tw->w)); return;
		case 'Z': return (id<32) ? LOG("cowardly refusing to destroy window 0x%x",id) : tw_close(tw);
		default: LOG("unknown cmd 0x%x \"%s\"", c0); return;
	}
}

static int last_obid = 0;

static void debug_wid(int fs) {
	if (fs<6 || fs>25) { LOG("debug_wid: invalid fs %d", fs); return; }
	char buf[96], *p = FONT_WTAB(fs);
	int i; for (i=0; i<95; i++) buf[i] = p[i+32] + 48;
	buf[95] = 0; LOG("debug_wid(%d):\"%s\"", fs, buf);
}

static void cmd1(char * str) {
	int i, t, k = *str;
	char * s = str + 1;
	if ((unsigned int)(k-65) < 26u) {
		int id = 0; while (*s & 80) id<<=4, id += hxd2i(*s), s++;
		if (*s==36) ++s;
		if (id) cmd_ob(k, last_obid = id, s);
		else if (last_obid) cmd_ob(k, last_obid, s);
		else LOG("cmd_ob(%c): missing obid", k);
		return;
	}
	switch(k) {
		case 0: return;
		case 'q': bye(*s);
		case '#': return;
		case '^': return CMD("%s", s);
		case '+': 
			  i = *(s++) - 48; t = *(s++);
			  for (k=0; *s&80; s++) k=16*k + hxd2i(*s);
			  if (*s==36) ++s;
			  t2sel_upd(i&1, k, t, *s=='.'?"":s, s+6*(*s!='.'));
			  return;
		case '=':
			   if ((unsigned int)(i = *(s++)-60) > 2u) LOG("=: exp. <,=,>");
			   else lsr_upd(i, s);
			   return;
		case 'v': memcpy(help_vpos, s, 4); return;
		case 't': switch(*s) { case 'w': case 0: write_tlog(); return;
				       case 'c': tlog_c_onq = (s[1]>63), tlog_c_bk = s[1]&15; return;
				       default: LOG("t: invalid subcmd"); return; }
		case 'f': choo_cmd(choo_tab+chtab_get(&choo_ch, *s), s+1); return; 
		case 'd': switch(*s) {
				case 'S': sleep(1); return;
				case 'B': obidtab_debug(&oi_box); return;
				case 'D': dirtab_debug(); return;
				case 'E': obidtab_debug(&oi_etc); return;
				case 'T': debug_wid(atoi(s+1)); return;
				case 't':
					 for (++s, k=0; 47<*s&&*s<58; ++s) k = 10*k + *s - 48;
					 if (*s==36) ++s;
					 LOG("wid: %d", tx_len(k, s));
					 return;
				case '=': lsr_debug();
				case 'f': if (s[1]=='?') LOG("dlfg: %d %s", dflg, dflg_s); 
					  else dflg = atoi_h(s+1);	
					  return;
				default: LOG("unknown debug-cmd 0x%x(%c)",*s, *s); return;
			}
		case 'm': switch(*s) {
				  case ':': return menu_act_s(s+1);
				  default:  return LOG("unknown menu-cmd 0x%x(%c)",*s,*s);
			  }
		default:LOG("unknown cmd 0x%x",*str); return;
	}
}

static ww_t * dot_mkwin(char ** pp) {
	if (**pp!='C') { LOG("dot_mkwin: exp:'C' got:0x%x'%c'", **pp, **pp); return NULL; }
	topwin * tw; ww_t * ww;
	int k = 0; ++*pp; while (**pp && **pp!='_') k=16*k+hxd2i(**pp), ++*pp;
	if (**pp!='_') goto ex_; else ++*pp;
	char wnm[4]; wnm[0] = 0;
	int tt = 16*hxd2i(**pp) + hxd2i((*pp)[1]); *pp += 2;
	int wt = 16*hxd2i(**pp) + hxd2i((*pp)[1]); *pp += 2;
	if (!(tw = tw_lookup(k)) && !(tw = topwin_mk(k, tt, wnm))) {
		LOG("dot_mkwin: create window(0x%x) failed", k); return NULL; }
	wnm[0]=wt; wnm[1] = ((**pp|32)=='z')?'.':**pp; wnm[2]=(*pp)[1]; wnm[3]=0;
	*pp += 2; if (**pp!='_') goto ex_; else ++*pp;
	ww = wlu_any_s(tw, wnm);
	if (ww && ww->cl->ch=='g') return ww;
	LOG("dot_mkwin: no%s", " graph widget found"+6*!ww); return NULL; 
ex_:	LOG("dot_mkwin: exp:'_' got:0x%x'%c'", **pp, **pp); return NULL;
}

static void dot_line(char * s) {
	static int state='G'; // G N S g n J x
	static ww_t * ww = NULL;
	if (!s) { ww = NULL; state = 'x'; return; }
	while (*s==' ' || *s=='\t' || *s=='\n') ++s;
	if (!*s) return;
	switch (state) {
		case 'G':
			if (memcmp(s, "grap", 4)) break;
			state = 'N'; return;
		case 'N': // node C<id24>_TTWWRij##<arg>
			if (memcmp(s, "node", 4)) break;
			s += 4; while (*s==' ' || *s=='\t') ++s;
			if (!(ww = dot_mkwin(&s))) break;
			dagraph_ini(ww, 1, s); state = 'S'; return;
		case 'S': case 'J':
			if (!memcmp(s, "stop", 4)) { state = 'g'; return; }
			if (state=='J') return;
			LOG("junk after def-node: \"%s\"", s); state = 'J'; return;
		case 'g':
			if (memcmp(s, "grap", 4)) break;
			dagraph_d_bbox(ww, s+6);
			state = 'n'; return;
		case 'n':
			if (!memcmp(s, "node", 4)) return dagraph_cd_node(ww, 1, s+5);
			else if (!memcmp(s, "edge", 4)) return dagraph_d_edge(ww, s+5);
			else if (memcmp(s, "stop", 4)) break;
			dagraph_done(ww, 1); ww = NULL; state = 'G';
			return;
		case 'x':
			if (!memcmp(s, "grap", 4)) state = 'N';
			else if (!memcmp(s, "stop", 4)) state = 'G'; 
			return;
		default:
			break;
	}
	LOG("dot_line: ?? st:%c in:\"%s\"", state, s);
	state = 'x';
}

static gboolean cmd_in (GIOChannel *src, GIOCondition condition, gpointer data) {
	static int zr_c = 0;
	char *str;
	gsize siz = 0;
	GError * err = 0;
	g_io_channel_read_line (src, &str, &siz, NULL, &err);
	if (!str) { 
		LOG("zero cmdstr #%d", ++zr_c);;
		if (zr_c>9) bye(0); else return TRUE;
	} else { zr_c = 0; }
	if (!siz) { LOG("cmdstr length 0 -- this should not happen!"); return TRUE; }
	last_obid = 0;
	if (str[siz-1]=='\n') str[--siz]=0;
	else LOG("warning: string ends with char 0x%x(%c) instead of nl", str[siz-1], str[siz-1]);
	if (*str == 9) {
		char *s, *toks;
		for (s=strtok_r(str+1, "\t", &toks); s; s=strtok_r(0, "\t", &toks))
			cmd1(s);
	} else {
		cmd1(str);
	}
        g_free(str);
	return TRUE;
}

static gboolean from_dot (GIOChannel *src, GIOCondition condition, gpointer data) {
	char *str;
	static int cont = 0;
	gsize siz = 0;
	GError * err = 0;
	g_io_channel_read_line (src, &str, &siz, NULL, &err);
	if (!str) { LOG("zero dot_str"); return TRUE; }
	if (!siz) { LOG("empty dot_str"); return TRUE; }
	if (str[siz-1]=='\n') str[--siz]=0;
	int c2 = str[siz-1]=='\\';
	// LOG("from_dot %02o l=%d \"%s\"", 8*cont+c2, siz, str);
	if (!cont) dot_line(str);
	cont = c2; g_free(str);
	return TRUE;
}
									                                                          
///////////////// guess what /////////////////////////////////////////////////

static void add_in(int fd, gboolean (*fun)(GIOChannel*,GIOCondition,gpointer), int flg) {
	if (fd<0) { const char * s = tpipe_name(-fd);
		     if ((fd = open(s, O_RDWR))<0) return LOG("%s: %s", s, strerror(errno)); }
	GError * err = 0;
	GIOChannel * chan = g_io_channel_unix_new(fd); 
	g_io_channel_set_encoding(chan, NULL, &err);
	if (flg&1) g_io_channel_set_buffered (chan, FALSE);
	if (flg&2) g_io_add_watch(chan, G_IO_HUP, &fifo_bye, pf_buf+(2*fd)),
		   g_io_add_watch(chan, G_IO_ERR, &fifo_bye, pf_buf+(2*fd+1));
	g_io_add_watch(chan, G_IO_IN, fun, 0);
}	

int main(int ac, char **av) {
	if (getenv("LF_TMPDIR")) write(2, "&CMD1\n", 6), signal(SIGHUP, SIG_IGN), signal(SIGINT, SIG_IGN);
	else LOG("(gui test mode)");
	int tpipe_fd = ac<2 ? -1 : qh4r(*(int*)av[1]);
	obidtab_ini(&oi_box, 0);
	obidtab_ini(&oi_dir, 1); ot_dir[0].flg = 3; ot_dir[0].id = 0; ot_dir[0].clip = NULL;
	obidtab_ini(&oi_etc, 0); // ot_etc[0].state = ot_etc[1].state = 0;
	cltab_init(); txtm_init(); choo_init(); pf_buf[0] = '~';
	gtk_init (&ac, &av);
	const char * s = getenv("LF_NOMWIN"); conf_nomwin = s && *s && (*s|32)-'n'; 
	const char * rcfn = getenv("LF_GTKRC"); if (!rcfn) rcfn = "lf.gtk.rc";
	char rch[20];
	int r, fd = open(rcfn, O_RDONLY);
	if (fd<0) LOG("%s: error %d", rcfn, errno);
	else if ((r=read(fd, rch, 16))!=16) LOG("%s: ret %d, errno %d", rcfn, r, errno);
	else    conf_lbh    = 10*rch[1] + rch[2] - 528,
		conf_lbfh   = conf_lbh  - rch[3] + 48,
		conf_lbfh_s = conf_lbfh - rch[4] + 48,
		close(fd), gtk_rc_parse (rcfn);
	conf_lbfs = get_fontsize(conf_lbfh, 0, NULL);
	conf_lbfs_s = get_fontsize(conf_lbfh_s, 0, NULL);
	conf_portwid = 6*FONT_NWID(conf_lbfs_s) + 6;
	menutab_init();
	help_readbuf((s=getenv("LF_HLP"))?s:"help.txt");
	if (!conf_nomwin) { 
		topwin_mk(0x07, '.', NULL); 
		topwin_mk(0x17, '/', NULL); 
		memcpy(&(ot_dir[0].nd), troot, 2*sizeof(GtkTreeIter));
	}
	add_in(0,        &cmd_in,   2);
	add_in(-'T',     &from_dot, 0);
	add_in(tpipe_fd, &tpipe_in, 3);
	gtk_main (); write(2, "whaaat???\n&CMD0\n", 16); return 0;
}

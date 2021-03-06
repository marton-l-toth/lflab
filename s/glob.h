#ifndef __qwe_glob_h__
#define __qwe_glob_h__

#define GLF_SAVED     0x20000
#define GLF_HITHERE   0x40000
#define GLF_AVOLSHR   0x80000
#define GLF_LIBMODE   0x100000
#define GLF_EMPTY     0x200000
#define GLF_RECOVER   0x400000
#define GLF_FSTATE    (GLF_EMPTY|GLF_RECOVER|GLF_SAVED)
#define GLF_RECORD    0x1000000
#define GLF_FINAL_ASV 0x2000000
#define GLF_SILENCE   0x4000000
#define GLF_INI0      0x8000000
#define GLF_INI1      0x10000000
#define GLF_GUIOK     0x20000000

#define DFLG_RFINST 1
#define DFLG_MODELDEL 2
#define DFLG_CMDECHO 4
#define DFLG_VOLTAB 8
#define DFLG_AUDIO 16
#define DFLG_AOET  32
#define DFLG_EXPR 64
#define DFLG_TRK  128
#define DFLG_GUI2 256
#define DFLG_SAVE 512
#define DFLG_MX 0x400
#define DFLG_JQ 0x800 
#define DFLG_GRTMP  0x1000
#define DFLG_RWMIDI 0x2000
#define DFLG_WRAP   0x4000
#define DFLG_TRK_C  0x8000
#define DFLG_PT     0x10000
#define DFLG_FE     0x20000
#define DFLG_MIDIEV 0x40000
#define DFLGZ       0x7ffff
#define DFLG_HELP "1:recfilt-state 2:model-del 4:guicmd 8:voltab 10:audio 20-sel_loop 40-expr 80-trk "\
		  "100-gui2 200:save 400:mixer 800:jobq 1000:graph(>workdir) 2000:midi-rw(ini only) "\
		  "4000-wrap 8000-trkcalc 10000:pt 20000:fe 40000:midiev"
#define MAX31 0x7fffffff
#define DEBUG_UTXT(J) (debug_utxt_buf + 64*(J))
#define IFDBGX(F) if (debug_flags&DFLG_##F)
extern int glob_flg;
extern int debug_flags, zero_sec;
extern int sample_rate;
extern double sample_length, natural_bpm, natural_bpm_r;
extern int killer_fd;
extern int qstat_size;

#define zeroblkD ((double*)(zeroblkC))
#define junkbuf  ((double*)(junkbufC))
extern char zeroblkC[32768], junkbufC[32768];
extern char save_file_name[1024], debug_utxt_buf[1024];

extern int v_major, v_minor;
extern const char * patch_string;

void ini_err(int ec);
void log(const char * fmt, ...);
void log_n(const char * fmt, ...);

void gui_errq_add(int ec, const char * s = 0), gui_errq_add2(int e1, int e2, const char*s = 0),
     gui_closewin(int oid), gui_sliderwin(int oid, int n, const double * lbl, const unsigned char * v0);
int  gui_acv_op(int j, int opw = -1);
#endif // __qwe_glob_h__

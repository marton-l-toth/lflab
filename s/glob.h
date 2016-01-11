#ifndef __qwe_glob_h__
#define __qwe_glob_h__

#define GLF_LIBMODE   0x100000
#define GLF_FINAL_ASV 0x200000

#define DFLG_RFINST 1
#define DFLG_MODELDEL 2
#define DFLG_GUICMD 4
#define DFLG_VOLTAB 8
#define DFLG_AUDIO 16
#define DFLG_PLOT2 32
#define DFLG_EXPR 64
#define DFLG_TRK 128
#define DFLG_GUI2 256
#define DFLG_SAVE 512
#define DFLG_MX 0x400
#define DFLG_JQ 0x800
#define DFLG_HELP "1:recfilt-state 2:model-del 4:guicmd 8:voltab 10:audio 20-plot2 40-expr 80-trk\n"\
		  "100-gui2        200:save   400:mixer 800:jobq"

extern int glob_flg;
extern int debug_flags;
extern int sample_rate;
extern double sample_length;
extern int killer_fd;

extern char mostly_zero[0x8080];
extern double junkbuf[4096];
extern char save_file_name[1024];

extern int v_major, v_minor;

#define zeroblkC (mostly_zero+128)
#define zeroblkD ((double*)(mostly_zero+128))
#define imp4097  ((double*)(mostly_zero+120))

void ini_err(int ec);
void log(const char * fmt, ...);
void log_n(const char * fmt, ...);

void gui_errq_add(int ec), gui_closewin(int oid),
     gui_sliderwin(int oid, int n, const double * lbl, const unsigned char * v0);
#endif // __qwe_glob_h__

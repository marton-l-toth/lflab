// This is an example source file. illustrating how to implement new primitive boxes
// The name of the contrib folder will be the name of the source file (without the ".cc")
// This name has to begin with a letter and be a valid C identifier
// It will be placed under ".!c.<first_letter>" in the object tree
// Please create only one source file for your collection of contrib boxes, in s/c
// It is recommended to use ../../t/prod.sh for compiling (change SDIR and possibly PDIR)

// first, some includes

#include <math.h>

#include "../box0.h" // needed for all primitive boxes
#include "../pzrf.h" // only if you want to make pole/zero filters
#include "../util.h" // error codes, NaN-lists and Scale01

// help texts: lines beginning with "//?" are collected from source files
// help texts are identified by {{{<package>.<cl_id>}}}
// where package (for contrib boxes) is the main folder name
// please include a short description of what the box does, what the input
// ports are for, and (unless it is small and constant) some hints for CPU usage
// (see all the builtin boxes for examples)
// because this is free software :) and also because it is internally used as a 
// command/field separator character, please avoid the dollar sign ($)

// since these example boxes are not intended for actual use (they are copies of ones in
// the "!b" tree), all boxes will have the same help text, referring to this file:

//?{{{zB._zB}}}
//?This is an example contrib box, please do not use it.
//?(It is a copy of a box in the "!b" tree, and may change/
//?disappear in future versions) The purpose of the !c.zB subtree
//?is to illustrate how to implement primitive boxes in C/C++.
//?If you are interested, please see the source file, c/zB.cc

// for any new box you have to
// (1) define a class (subcl. of BoxInst) with a calc() function
//     (either explicitly or using a convenience macro)
// (2) insert a node referring to the class into the object tree
//     (do this in the ci_<foldername>() init function, see at the end of this file)

// the simplest case: stateless box (1 in, 1 out) computing a single-variable function
// just use the following macro FUN1_BOX(<classname>, <expression>) to define the class
// in the expression "x" is a variable of type "double"
// double-precision is used throughout the program, for reducing worries about
// numerical stability

FUN1_BOX(ZBAbs, fabs(x))

// For anything more complicated, you need a calc() funtion:
// virtual int calc(int inflg, double** inb, double** outb, int n)
// inflg: for all k=0...29, iflg&(1<<k) is nonzero if inb[k] is non-constant
// inb: pointers to input buffers, inb[k][0] for const, inb[k][0]...inb[k][n-1] otherwise
// outb: same for output buffers
// n: number of samples (0...4095)
// return value: OK: a bit-vector for output, analogous to "inflg"
//               error: negative value (see errtab.txt, esp. RTE_*)

// It is possible that more input buffers point to the same address
// Output buffers are always distinct (except when all point to the junk buffer)
// Whether or not an input buffer may have the same address as an output buffer
// depends on the configuration (see qmk_box below)
// Any two non-junk buffers (input/output) are either idential or disctinct -- they
// never overlap.

// For stateless boxes, there is a macro for class definition and calc() head,
// you only have to write the function body (this is the two input "+" box, it
// outputs the sum of the two inputs (which is constant if both inputs are constant)

STATELESS_BOX_0(ZBAdd) { 
	double *o = outb[0]; switch(inflg&3) {
		case 0: return **outb = *inb[0] + *inb[1], 0;
		case 1: { double *p=inb[0],y=*inb[1]; for (int i=0;i<n;i++) o[i] = p[i]+y; return 1; }
		case 2: { double x=*inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = x+q[i]; return 1; }
		case 3: { double *p=inb[0],*q=inb[1]; for (int i=0;i<n;i++) o[i] = p[i]+q[i]; return 1; }}}

// Primitive boxes have no configurable static properties, they can only be configured via inputs.
// If you want to implement a box with a variable number of inputs/outputs, the only possibility
// is to make a small _collection_ of boxes, one box for each #in/#out configuration (defined by
// an integer parameter). In stateless parametrized boxes defined using the following macro,
// this parameter can be accessed as member variable "m_arg". 

// This example simply copies one of the inputs to output, using the last input as index;
// here "m_arg" is the number of selectable inputs
// (the macro BVFOR_JM iterates through the nonzero bits of a 32-bit integer)

inline int ixround(double x, int n) { int r = (int)lround(x);
	        return ((unsigned int)r<(unsigned int)n) ? r : (r<0?0:n-1); }

STATELESS_BOX_1(ZBSel) {
	double *pv, *q = *outb, *pi = inb[m_arg];
	int ixm = 1<<m_arg;
	if (!(inflg&ixm)) {
		int j = ixround(*pi, m_arg), r = (inflg>>j)&1;
		if (!r) *q = *inb[j]; else if ((pv=inb[j])!=q) for (int i=0;i<n;i++) q[i]=pv[i];
		return r;
	}
	double con[m_arg]; BVFOR_JM((ixm-1)&~inflg) con[j] = inb[j][0];
	for (int j,i=0; i<n; i++) j    = ixround(pi[i], m_arg),
		q[i] = (inflg&(1<<j)) ? inb[j][i] : con[j];
	return 1;
}

// If you want to implement a box with a state (other than pole/zero filters, see below),
// you have to explicitly define a class. Beside the calc() function you also have to
// define a contructor to initialize the state (default (no args) constructor for
// singular boxes, and one with an int argument for parametrized boxes. )
// (for the stateless boxes above, the contructors are defined by the macros.)
// Unless you explicitly allocate heap memory, you don't have to write a destructor.
// The following is a copy of the simple accumulator filter "=acc"
// sample_rate(int) and sample_length(double) are global variables
// sample_length is ( 1.0/(double)sample_rate ), sample_rate is currently fixed at 44100

class ZBAcc : public BoxInst {
        public: 
                ZBAcc() : m_v(0.0) {}
                virtual int calc(int inflg, double** inb, double** outb, int n);
        protected:
                double m_v;
};

int ZBAcc::calc(int inflg, double** inb, double** outb, int n) {
        double x, y, *p = *outb, *q0, *q1; switch (inflg & 3) {
        case 0 : x = sample_length * **inb, y = exp(-sample_length * inb[1][0]);
                 for (int i=0; i<n; i++) p[i] = ((m_v+=x)*=y); return 1;
        case 1 : q0 = *inb, y = exp(-sample_length * inb[1][0]);
                 for (int i=0; i<n; i++) p[i] = ((m_v+=sample_length*q0[i])*=y); return 1;
        case 2 : x = sample_length * **inb, q1 = inb[1];
                 for (int i=0; i<n; i++) p[i] = ((m_v+=x)*=exp(-sample_length*q1[i])); return 1;
        case 3 : q0 = *inb, q1 = inb[1];
                 for (int i=0;i<n;i++) p[i]= ((m_v+=sample_length*q0[i])*=exp(-sample_length*q1[i])); return 1;
        }} // no it doesn't

// For filters composed of pole/zero pairs (either sequentially or parallelly combined),
// you only have to write the fuction that generates the pole/zero pairs, the rest is taken
// care of by the following macros and the base classes.
// The function head is "int NM::mk_filter(double **inb)", hidden in the macro
// The pole-zero pairs have to be filled out from an array of RECF_PZItem
// You should allocate this array with "new" and not free it (it could be used for the
// command "filter display" from the main menu (which can be BTW quite useful for checking
// if your , and only freed when replaced with a new one)
// You can fill them manually (pr, pi, zr, zi and c (constant multipler)), or you can use
// some convenience functions (see decl. in pzrf.h)
// You can find more examples in pzrf.cc (e.g. QUICK_PZFILT_1 for parametrized filter boxes,
// or QUICK_PPZFILT for parallel (added) filters)

// In the following example (harmonic equalizer filter), one of the arguments in a NaN list
// the macro NAN_UNPK_32(name, ptr, defval) unpacks at most 16 32-bit-integers from the NaN
// pointed by "ptr", to name_x (length to name_l); the rest of the array (all in the case
// of empty list or number) is filled with defval (NAN_UNPK_8 does the same with 8-bit chars)
// ipow() and ipows() calculate power with unsigned/signed integer exponents.
// fq_warp() does "pre-warping" for frequency -- it converts value in Hz to the "frequency"
// required by filters.

QUICK_PZFILT(ZBEqLF) {
        double fq = inb[0][0], step = fq*inb[1][0], wid1 = inb[2][0],
               amp1 = exp(M_LN2 / inb[3][0]), wid2 = wid1 * fq;
        NAN_UNPK_32(amp, inb[4], 0);
        int i, j, k, np = amp_n; 
	set_n(np); // sets number poles, allocates m_ab
        RECF_PZItem * pzn = new RECF_PZItem[np]; // pole-zero pairs
        for (i=j=0; i<np; i++, fq+=step) if ((k=amp_x[i])) pzn[j].eql(fq_warp(fq),wid2/fq,ipows(amp1,k)), j++;
        rfpz_transform(m_ab, pzn, j); return 0; // pzn not freed now, available for main menu/filter disp.
}

// And finally, the init function (called before files are read or GUI started)
// You start with an empty folder node (the parameter), and you can add more folders
// and boxes with the qmk_dir() and qmk_box() functions.

// folders are easy: you pass a node pointer and a name, the function:
// ANode * qmk_dir(ANode * up, const char * nm);
// returns the address of the newly created folder node, which you can use in further
// qmk_dir() or qmk_box() calls.

// boxes are a little more complicated, the function is:
// ANode * qmk_box(ANode * up, const char * nm, qmb_arg_t qa, int k, int ni, int no, 
//				               const char * cl, const char * fmt, ...);
// up: parent node (must be a folder)
// nm: name
// qa: QMB_ARG0(simple_box_class) or QMB_ARG1(param_box_class), type is ptr-to-function
// k:  parameter (for single boxes it is ignored)
// ni: number of input ports
// no: number of output ports (add 32 for io-alias, see below)
// cl: class-id (used for help text lookup)
// fmt: GUI config format string -- if NULL or "", the box will have the default GUI config,
//      and no further arguments are expected, if non-empty, it can contain the following chars:
//      i I o O r R  -- target select (no arg)  iI: input oO: output rR: rgb/all
//      1..7  -- load(small trg) / store(cap trg) to slot 1..7 (here r/R means all)  (no arg)
//      * -- set all rgb (string, see rgb below) or all input/output (string, see iolbl below)
//      + -- rgb fg (string)
//      - -- rgb bg (string) or i/o interval (int:256*start+len, string)
//      arguments: rgb is 3x or 6x char ( '%'(37)...'z'(122), 8-bit val is 3*(chr-37) )
//	           iolbl string: $-separated labels (max. 4 char/label)
// although the function returns the address of the new box node, this value is not really useful
// for contrib init functions

// io-alias: when allowed (no&32) input buffers(s) may have the same address as output buffer(s)
// that can cause a problem if you use the output buffer as a temporary store while still
// needing the input buffers -- if unsure, leave it disabled (its only disadvantage is that
// graph boxes containing this box might allocate one more temporary buffer)

// In the early initialization there is no fancy error handling: when an operation fails,
// the program quits with a core dump. (since you start with an empty folder, nothing
// surprising should happen.) In case something impossible happens, you too can abort the
// program with the bug() function; if you only want to write a debug message to log/console
// use the log() function -- both log() and bug() take printf-style arguments.

void ci_zB(ANode * r) {
	log("lflab is easy to extend, see %s for details", __FILE__);
	qmk_box(r, "abs", QMB_ARG0(ZBAbs), 0, 1, 33, "zB", "i*o*R*1", "X", "absX", "Pz%%0%"); 
	qmk_box(r, "+", QMB_ARG0(ZBAdd), 0, 2, 1, "zB", "1i*o*", "x$y", "x+y");

	ANode * dsel = qmk_dir(r, "sel"); qmb_arg_t qa = QMB_ARG1(ZBSel);
	char nm[8]; memcpy(nm, "sel02", 6); 
	qmk_box(dsel, nm, qa, 2, 3, 33, "zB", "R*1i-", "uu%99%", 513, "sel"); // store color only
	for (int i=3; i<30; i++) nm[3] = '0'+i/10, nm[4] = '0'+i%10,   // keep default labels (exc. "sel")
		qmk_box(dsel, nm, qa, i, i+1, 33, "zB", "1i-", 1+256*i, "sel");

	qmk_box(r, "=acc", QMB_ARG0(ZBAcc), 0, 2, 1, "zB", "i*", "in$loss");
	qmk_box(r, "=eql-ls", QMB_ARG0(ZBEqLF),0,6,33,"zB", "i*R*", "in$fq1$step$wid$2div$[a]", "uuu3%F");
}

// you can find more examples (although with less comments) in the main source files.
// coding style:
// - a tab is 8 characters
// - lines shall be no longer than 111 characters
// - names of member variables (but not member functions) begin with "m_"
// - otherwise I don't expect you to imitate my coding style :)
// Have fun!

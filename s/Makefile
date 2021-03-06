CC=g++

SRC1 = main.cc node.cc util.cc cmd.cc box.cc util2.cc asnd.cc job.cc mx.cc pt.cc midi.cc
INC1 = glob.h  node.h  util.h  cmd.h  box.h  util2.h  asnd.h  job.h  mx.h  pt.h  midi.h
OBJ1 = main.o  node.o  util.o  cmd.o  box.o  util2.o  asnd.o  job.o  mx.o  pt.o  midi.o

SRC2 = expr.cc wrap.cc guistub.cc combo.cc pzrf.cc trk.cc
INC2 = expr.h  wrap.h  guistub.h  combo.h  pzrf.h  trk.h
OBJ2 = expr.o  wrap.o  guistub.o  combo.o  pzrf.o  trk.o

SRC3 = nz.cc rf2.cc fe.cc arbox.cc grbox.cc it.cc in.cc b0.cc v.cc
OBJ3 = nz.o  rf2.o  fe.o  arbox.o  grbox.o  it.o  in.o  b0.o  v.o 

INCX = aobuf.h lwarr.h uc0.h uc1.h tab7.h box0.h contrib.h distr.h

SRCC = $(shell c/upd_contrib.sh cc c/)
OBJC = $(shell c/upd_contrib.sh o  co/)

SRCS = $(SRC1) $(SRC2) $(SRC3) $(SRCC)
INCS = $(INC1) $(INC2) $(INCX)
OBJS = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJC) contrib.o hlpn.o
export LF_BINS = lf.gui lf.bb lf.pump lf.bin

GTKINC = -I/usr/include/gtk-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include/ -I/usr/include/glib-2.0 -I/usr/include/cairo -I/usr/include/pango-1.0/ -I/usr/include/gdk-pixbuf-2.0 -I/usr/lib/x86_64-linux-gnu/gtk-2.0/include -I/usr/include/atk-1.0
GTKLIB = -lgtk-x11-2.0 -lgdk-x11-2.0 -lcairo -lgobject-2.0 -lglib-2.0 -lm

all: $(LF_BINS) lf help.txt lf.con lf.acv lf.lic lf.ed

lf.ed:
	ln -sf lf.bb lf.ed
lf.con:
	ln -sf lf.bb lf.con
lf.acv:
	ln -sf lf.bb lf.acv
lf.lic:
	ln -sf lf.bb lf.lic
lf:
	ln -sf lf.bin lf

contrib.cc: c c/upd_contrib.sh
	c/upd_contrib.sh

cfgtab.inc: cfgtab.txt cfg2inc.sh
	./cfg2inc.sh cfgtab.txt > cfgtab.inc

errtab.inc: errtab.txt errtab.sh
	./errtab.sh errtab.txt > errtab.inc

co/%.o: c/%.cc $(INCS) errtab.inc cfgtab.inc
	$(CC) -g -c $(CFLAGS) -o$@ $<
%.o: %.cc $(INCS) errtab.inc cfgtab.inc
	$(CC) -g -c $(CFLAGS) -o$@ $<
%.s: %.cc $(INCS) errtab.inc cfgtab.inc
	$(CC) -S $(CFLAGS) -o$@ $<

lf.bin: $(OBJS)  
	$(CC) -g $(OBJS) -lasound -olf.bin

lf.gui: gui.c uc0.h uc1.h errtab.inc cfgtab.inc
	gcc -g $(CFLAGS) gui.c $(GTKINC) $(GTKLIB) -olf.gui

lf.bb:  bb.c uc0.h
	gcc -g $(CFLAGS) bb.c -olf.bb -lm

lf.pump: asnd.cc
	ln -sf asnd.cc pump.c && gcc -g -Wall $(CFLAGS) -olf.pump pump.c -lm -lasound

hlpt.txt: hlpt.sh helpm.txt
	./hlpt.sh helpm.txt > hlpt.txt

hlpn.cc:  hlpc.sh helpm.txt
	./hlpc.sh helpm.txt > hlpn.cc

help.txt: hlpt.txt $(SRCS)
	tail -n +2 hlpt.txt > help.txt && sed -n '/^\/\/?/{s/^\/\/? *//;p}' $(SRCS) >> help.txt

deb:	$(LF_BINS) help.txt ex.lf lf.dsc.txt COPYING deb.sh pkg.sh
	./deb.sh -v

rpm:	$(LF_BINS) help.txt ex.lf lf.dsc.txt COPYING rpm.sh lf.spec pkg.sh
	./rpm.sh 

clean:
	rm -f $(OBJS) lf.gui lf.bin lf.bb lf help.txt cfgtab.inc errtab.inc deb lflab-*.deb contrib.cc

stat:
	wc $(SRCS) $(INCS) gui.c bb.c | sort -n


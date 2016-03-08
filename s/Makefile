CC=g++
CFLAGS ?= -O2 -fpredictive-commoning -fgcse-after-reload

SRC1 = main.cc node.cc util.cc cmd.cc box.cc util2.cc asnd.cc job.cc mx.cc pt.cc midi.cc
INC1 = glob.h  node.h  util.h  cmd.h  box.h  util2.h  asnd.h  job.h  mx.h  pt.h  midi.h
OBJ1 = main.o  node.o  util.o  cmd.o  box.o  util2.o  asnd.o  job.o  mx.o  pt.o  midi.o

SRC2 = expr.cc wrap.cc nz.cc guistub.cc gp.cc combo.cc pzrf.cc trk.cc
INC2 = expr.h  wrap.h  nz.h  guistub.h  gp.h  combo.h  pzrf.h  trk.h
OBJ2 = expr.o  wrap.o  nz.o  guistub.o  gp.o  combo.o  pzrf.o  trk.o

SRC3 = rf2.cc fe.cc arbox.cc grbox.cc in.cc b0.cc v.cc
OBJ3 = rf2.o  fe.o  arbox.o  grbox.o  in.o  b0.o  v.o 

INCX = aobuf.h lwarr.h uc0.h uc1.h tab7.h box0.h contrib.h

SRCC = $(shell c/upd_contrib.sh cc c/)
OBJC = $(shell c/upd_contrib.sh o)

SRCS = $(SRC1) $(SRC2) $(SRC3) $(SRCC)
INCS = $(INC1) $(INC2) $(INCX)
OBJS = $(OBJ1) $(OBJ2) $(OBJ3) $(OBJC) contrib.o

GTKINC = -I/usr/include/gtk-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include/ -I/usr/include/glib-2.0 -I/usr/include/cairo -I/usr/include/pango-1.0/ -I/usr/include/gdk-pixbuf-2.0 -I/usr/lib/x86_64-linux-gnu/gtk-2.0/include -I/usr/include/atk-1.0
GTKLIB = -lgtk-x11-2.0 -lgdk-x11-2.0 -lcairo -lgobject-2.0 -lglib-2.0 -lm

all: lf.gui lf.bb lf.bin lf lf.ini help.txt lf.io lf.con lf.acv lf.lic

lf.io:
	ln -sf lf.bb lf.io
lf.con:
	ln -sf lf.bb lf.con
lf.acv:
	ln -sf lf.bb lf.acv
lf.lic:
	ln -sf lf.bb lf.lic

contrib.cc: c c/upd_contrib.sh
	c/upd_contrib.sh

lf: cfgtab.txt lf.src cfg2sh.sh
	./cfg2sh.sh cfgtab.txt lf.src > lf; chmod 755 lf

cfgtab.inc: cfgtab.txt cfg2inc.sh
	./cfg2inc.sh cfgtab.txt > cfgtab.inc

lf.ini: cfgtab.txt cfg2ini.awk
	awk '-F|' -f cfg2ini.awk < cfgtab.txt > lf.ini

errtab.inc: errtab.txt errtab.sh
	./errtab.sh errtab.txt > errtab.inc

%.o: c/%.cc $(INCS) errtab.inc cfgtab.inc
	$(CC) -g -Wall -Wno-unused-function -c $(CFLAGS) -o$@ $<
%.o: %.cc $(INCS) errtab.inc cfgtab.inc
	$(CC) -g -Wall -Wno-unused-function -c $(CFLAGS) -o$@ $<

lf.bin: $(OBJS)  
	$(CC) -g $(OBJS) -lasound -olf.bin

lf.gui: gui.c uc0.h uc1.h errtab.inc
	gcc -g -Wall -Wno-unused-function $(CFLAGS) gui.c $(GTKINC) $(GTKLIB) -olf.gui

lf.bb:  bb.c uc0.h
	gcc -g -Wall $(CFLAGS) bb.c -olf.bb -lm

help.txt: helpm.txt $(SRCS)
	cp helpm.txt help.txt && sed -n '/^\/\/?/{s/^\/\/? *//;p}' $(SRCS) >> help.txt

deb:	lf.gui lf.bb lf.bin lf help.txt lf.ini ex.lf lf.dsc.txt COPYING deb.sh pkg.sh
	./deb.sh -v

clean:
	rm -f $(OBJS) lf.gui lf.bin lf.bb lf help.txt cfgtab.inc errtab.inc deb lflab-*.deb contrib.cc

stat:
	wc $(SRCS) $(INCS) gui.c bb.c


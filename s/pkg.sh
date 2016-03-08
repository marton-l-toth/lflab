nm=lflab
tmpr=/run/shm
([[ -d $tmpr ]] && [[ -w $tmpr ]]) || tmpr=/tmp
vers2=$(sed 's/^.*jor=//;s/,v_mi.*;//' v.cc).$(sed 's/^.*minor=//;s/;//' v.cc)
vers=$vers2-0
trgdir=/usr/share/$nm
tmptree=$tmpr/$$
tmplf=$tmptree$trgdir
bin="lf.bb lf.gui lf.bin"
misc="lf lf.ini help.txt ex.lf lf.gtk.rc.def COPYING"
who="Marton Laszlo Toth (initial dot initial dot fullsurname at gmail dot com)"
dsc0="audio lab with emphasis on linear filters"

mkdir -p $tmplf
for a in $bin; do strip -o $tmplf/$a $a || exit 1; done
cp $misc $tmplf

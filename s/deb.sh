#!/bin/bash
set -e
rm -f deb
nm=lflab
tmpr=/run/shm
([[ -d $tmpr ]] && [[ -w $tmpr ]]) || tmpr=/tmp
vers=$(sed 's/^.*jor=//;s/,v_mi.*;//' v.cc).$(sed 's/^.*minor=//;s/;//' v.cc)-0
arch=$(dpkg-architecture -qDEB_HOST_ARCH)
pkgnm=$nm-$vers-$arch
pkgdir=$tmpr/$pkgnm
debdir=$pkgdir/DEBIAN
trgdir=/usr/share/$nm
pktrg=$pkgdir/$trgdir
bin="lf.bb lf.gui lf.bin"
misc="lf lf.ini help.txt ex.lf lf.gtk.rc.def lf.cleanup COPYING"
ctl=$debdir/control
dep0="graphviz, gnuplot-x11"
who="Marton Laszlo Toth (initial dot initial dot fullsurname at gmail dot com)"
dsc0="audio lab with emphasis on linear filters"
bblnk="/usr/bin/lf.acv $trgdir/lf.con $trgdir/lf.io $trgdir/lf.lic"

[[ -d $pkgdir ]] && rm -r $pkgdir
mkdir -p $debdir
mkdir -p $pktrg
for a in $bin; do strip -o $pktrg/$a $a || exit 1; done
cp $misc $pktrg

echo -en "Package: $nm\nVersion: $vers\nArchitecture: $arch\nDepends: $dep0" > $ctl
dpkg -S $(ldd $bin | awk '/=>.*\//{print $3}' | sort -u) |  cut -d: -f1 | sort -u |
	(while read a; do echo -n ", $a" >> $ctl; done)
echo -e "\nMaintainer: $who\nDescription: $dsc0" >> $ctl

(echo -en "#!/bin/sh\nln -sf $trgdir/lf /usr/bin/$nm"
for t in $bblnk; do echo -n " && ln -s $trgdir/lf.bb $t"; done; echo) > $debdir/postinst
echo -e "#!/bin/sh\nrm /usr/bin/$nm $bblnk" >> $debdir/prerm
chmod 755 $debdir/postinst $debdir/prerm

dpkg-deb -b $pkgdir && cp $tmpr/$pkgnm.deb ./
dpkg -I ./$pkgnm.deb > deb || rm deb
[[ "$1" != "-v" ]] || cat deb

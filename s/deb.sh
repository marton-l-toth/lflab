#!/bin/bash
# does not work when dirnames contain spaces, sorry
set -e
rm -f deb
. pkg.sh
libdep=$(ldd $bin | awk '/=>.*\//{print $3}' | sort -u | grep -Ev 'libffi|libpango|libgraphite|libharfbuzz')
arch=$(dpkg-architecture -qDEB_HOST_ARCH)
pkgnm=$nm-$vers-$arch
pkgdir=$tmpr/$pkgnm
debdir=$pkgdir/DEBIAN
ctl=$debdir/control
dep0="graphviz, gnuplot-x11"
bblnk="/usr/bin/lf.acv $trgdir/lf.con $trgdir/lf.ed $trgdir/lf.lic"

[[ -d $pkgdir ]] && rm -r $pkgdir
mv $tmptree $pkgdir
mkdir $debdir

echo -en "Package: $nm\nVersion: $vers\nArchitecture: $arch\nDepends: $dep0" > $ctl
dpkg -S $libdep |  cut -d: -f1 | sort -u | (while read a; do echo -n ", $a" >> $ctl; done)
echo -en "\nMaintainer: $who\nDescription:" >> $ctl
grep -v '^#' lf.dsc.txt | sed 's/^/ /' >> $ctl

(echo -en "#!/bin/sh\nln -sf $trgdir/lf.bin /usr/bin/$nm"
for t in $bblnk; do echo -n " && ln -s $trgdir/lf.bb $t"; done; echo) > $debdir/postinst
echo -e "#!/bin/sh\nrm /usr/bin/$nm $bblnk" >> $debdir/prerm
chmod 755 $debdir/postinst $debdir/prerm

fakeroot dpkg-deb -b $pkgdir && cp $tmpr/$pkgnm.deb ./
dpkg -I ./$pkgnm.deb > deb || rm deb
[[ "$1" != "-v" ]] || cat deb

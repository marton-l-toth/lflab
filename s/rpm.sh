#!/bin/bash
# does not work when dirnames contain spaces, sorry
set -e
rm -f rpm
. pkg.sh

nmv=$nm-$vers2
mv $tmptree $tmpr/$nmv
rpmdir=$tmpr/$nm-$vers2-rpm
rm -rf $rpmdir
mkdir -p $rpmdir/SOURCES
tar -C$tmpr -cf $rpmdir/SOURCES/lf_pk.tar $nmv
dscf=$tmpr/$nm-$vers2-dsc.txt
sum=$(echo -n $(head -2 $dscf | tail -1))
tail -n +4 lf.dsc.txt > $dscf

rpmbuild -bb -D "_topdir $rpmdir" -D "LF_DIR $trgdir" -D "LF_DSC $dscf" -D "LF_NM $nm" \
-D "LF_VER $vers2" -D "LF_SUM $sum" -D "LF_URL https://github.com/marton-l-toth/lflab" lf.spec 2>&1 | tee -a rpm
cp -v $rpmdir/RPMS/*/*.rpm ./ 2>&1 | tee -a rpm
echo $0 done at $(date) 2>&1 | tee -a rpm

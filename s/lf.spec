Name: %{LF_NM}	
Version: %{LF_VER}	
Release: 0
Summary: %{LF_SUM}
Group: Applications/Multimedia
License: GPL 2.0
URL: %{LF_URL}
Source0: lf_pk.tar

Requires: gnuplot graphviz
%define debug_package %{nil}

%description
%include %{LF_DSC}

%prep
%setup -q

%build
%install
cp -a usr %{buildroot}
mkdir %{buildroot}/usr/bin

%post
ln -sf %{LF_DIR}/lf /usr/bin/%{LF_NM}
ln -sf %{LF_DIR}/lf.bb /usr/bin/lf.acv
ln -sf %{LF_DIR}/lf.bb %{LF_DIR}/lf.con
ln -sf %{LF_DIR}/lf.bb %{LF_DIR}/lf.io
ln -sf %{LF_DIR}/lf.bb %{LF_DIR}/lf.lic

%postun
rm /usr/bin/%{LF_NM} /usr/bin/lf.acv %{LF_DIR}/lf.con %{LF_DIR}/lf.io %{LF_DIR}/lf.lic
rmdir %{LF_DIR}

%files
%{LF_DIR}

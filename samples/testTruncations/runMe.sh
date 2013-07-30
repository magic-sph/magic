#!/bin/bash
export OMP_NUM_THREADS=8
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100m
export KMP_AFFINITY=noverbose,granularity=core,scatter

# nphi=96 (default)
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test96m4",/g' input.nml
./magic.exe input.nml
# nphi=128
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =128,/g' input.nml
sed -i 's/tag.*/tag         ="test128",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test128m4",/g' input.nml
./magic.exe input.nml
# nphi=192
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =192,/g' input.nml
sed -i 's/tag.*/tag         ="test192",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test192m4",/g' input.nml
./magic.exe input.nml
# nphi=256
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =256,/g' input.nml
sed -i 's/tag.*/tag         ="test256",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test256m4",/g' input.nml
./magic.exe input.nml
# nphi=288
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =288,/g' input.nml
sed -i 's/tag.*/tag         ="test288",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test288m4",/g' input.nml
./magic.exe input.nml
# nphi=320
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =320,/g' input.nml
sed -i 's/tag.*/tag         ="test320",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test320m4",/g' input.nml
./magic.exe input.nml
# nphi=384
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =384,/g' input.nml
sed -i 's/tag.*/tag         ="test384",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test384m4",/g' input.nml
./magic.exe input.nml
# nphi=400
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =400,/g' input.nml
sed -i 's/tag.*/tag         ="test400",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test400m4",/g' input.nml
./magic.exe input.nml
# nphi=512
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =512,/g' input.nml
sed -i 's/tag.*/tag         ="test512",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test512m4",/g' input.nml
./magic.exe input.nml
# nphi=640
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =640,/g' input.nml
sed -i 's/tag.*/tag         ="test640",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test640m4",/g' input.nml
./magic.exe input.nml
# nphi=768
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =768,/g' input.nml
sed -i 's/tag.*/tag         ="test768",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test768m4",/g' input.nml
./magic.exe input.nml
# nphi=800
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =800,/g' input.nml
sed -i 's/tag.*/tag         ="test800",/g' input.nml
./magic.exe input.nml
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         ="test800m4",/g' input.nml
./magic.exe input.nml
# nphi=864
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =864,/g' input.nml
sed -i 's/tag.*/tag         ="test864m4",/g' input.nml
./magic.exe input.nml
# nphi=1024
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =1024,/g' input.nml
sed -i 's/tag.*/tag         ="test1024m4",/g' input.nml
./magic.exe input.nml

# Restore initial values
sed -i 's/n_phi_tot.*/n_phi_tot   =96,/g' input.nml
sed -i 's/tag.*/tag         ="test96",/g' input.nml
sed -i 's/minc.*/minc        =1,/g' input.nml

# Build the equivalent output
cat e_kin.test96* e_kin.test128* e_kin.test192* e_kin.test256* e_kin.test288* e_kin.test320* e_kin.test384* e_kin.test400* e_kin.test512* e_kin.test640* e_kin.test768* e_kin.test800* e_kin.test864m4 e_kin.test1024m4 > e_kin.test

# Cleaning
rm scond.dat
rm *.test96* *.test128* *.test192* *.test256* *.test288* *.test320* *.test384* *.test400* *.test512* *.test640* *.test768* *.test800* *.test864m4 *.test1024m4

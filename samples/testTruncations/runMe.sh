#!/bin/bash
export OMP_NUM_THREADS=16
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100m
export KMP_AFFINITY=noverbose,granularity=core,scatter

sed -i "s/nThreadsRun *= *[0-9]*/nThreadsRun = 16/" input.nml

# nphi=96 (default)
tag='test96'
sed -i 's/n_phi_tot.*/n_phi_tot   =96,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
sed -i 's/minc.*/minc        =1,/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag="test96m4"
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=128
tag='test128'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =128,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test128m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=192
tag='test192'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =192,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test192m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=256
tag='test256'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =256,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test256m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=288
tag='test288'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =288,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test288m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=320
tag='test320'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =320,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test320m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=384
tag='test384'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =384,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test384m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=400
tag='test400'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =400,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test400m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=512
tag='test512'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =512,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test512m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=640
tag='test640'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =640,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test640m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=768
tag='test768'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =768,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test768m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=800
tag='test800'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =800,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
tag='test800m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=864
tag='test864m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =864,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag
# nphi=1024
tag='test1024m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =1024,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
./magic.exe input.nml >output.$tag

# Restore initial values
tag='test96'
sed -i 's/n_phi_tot.*/n_phi_tot   =96,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
sed -i 's/minc.*/minc        =1,/g' input.nml

# Build the equivalent output
cat e_kin.test96* e_kin.test128* e_kin.test192* e_kin.test256* e_kin.test288* e_kin.test320* e_kin.test384* e_kin.test400* e_kin.test512* e_kin.test640* e_kin.test768* e_kin.test800* e_kin.test864m4 e_kin.test1024m4 > e_kin.test

# Cleaning
rm scond.dat
rm *.test96* *.test128* *.test192* *.test256* *.test288* *.test320* *.test384* *.test400* *.test512* *.test640* *.test768* *.test800* *.test864m4 *.test1024m4

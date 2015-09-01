#!/bin/bash
if [ $# -ge 1 ]; then
    hybrid=$1
    nomp=$2
else
    hybrid="no"
fi

export F_UFMTENDIAN=big
if [ "$hybrid" == "hybrid" ]; then
    nmpi=2
    #nomp=4
    echo "Running hybrid mode with $nomp OpenMP threads per $nmpi MPI processes."
    export I_MPI_PIN_DOMAIN=socket
    export MP_BINDPROC=no
    export MP_TASK_AFFINITY=core:$nomp
    export KMP_STACKSIZE=1g
    export KMP_AFFINITY=noverbose,granularity=core,compact
else
    nomp=1
    nmpi=8
    echo "Running pure MPI code with $nmpi MPI processes."
    export I_MPI_PIN_PROCESSOR_LIST=allcores
fi
export OMP_NUM_THREADS=$nomp

keep=no

exitcode=0
# nphi=96 (default)
rm e_kin.test
tag='test96'
sed -i 's/n_phi_tot.*/n_phi_tot   =96,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
sed -i 's/minc.*/minc        =1,/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
    exitcode=42
    echo "Run failed for tag $tag"
fi

tag="test96m4"
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=128
tag='test128'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =128,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test128m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=192
tag='test192'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =192,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test192m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=256
tag='test256'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =256,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test256m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=288
tag='test288'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =288,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test288m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=320
tag='test320'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =320,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test320m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=384
tag='test384'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =384,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test384m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=400
tag='test400'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =400,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test400m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=512
tag='test512'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =512,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test512m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=640
tag='test640'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =640,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test640m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=768
tag='test768'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =768,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test768m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=800
tag='test800'
sed -i 's/minc.*/minc        =1,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =800,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
tag='test800m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=864
tag='test864m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =864,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi
# nphi=1024
tag='test1024m4'
sed -i 's/minc.*/minc        =4,/g' input.nml
sed -i 's/n_phi_tot.*/n_phi_tot   =1024,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
echo "Running $tag"
rm -f *.$tag
mpiexec -n $nmpi ./magic.exe input.nml >output.$tag
if [ $? -eq 0 ]; then
    # successful run
    cat e_kin.$tag >>e_kin.test
    if [ $keep != "yes" ]; then
	rm *.$tag
    fi
else
exitcode=42
    echo "Run failed for tag $tag"
fi

# Restore initial values
tag='test96'
sed -i 's/n_phi_tot.*/n_phi_tot   =96,/g' input.nml
sed -i 's/tag.*/tag         =\"'"$tag"'\",/g' input.nml
sed -i 's/minc.*/minc        =1,/g' input.nml

exit $exitcode
# Build the equivalent output
#cat e_kin.test96* e_kin.test128* e_kin.test192* e_kin.test256* e_kin.test288* e_kin.test320* e_kin.test384* e_kin.test400* e_kin.test512* e_kin.test640* e_kin.test768* e_kin.test800* e_kin.test864m4 e_kin.test1024m4 > e_kin.test

# Cleaning
#rm *.test96* *.test128* *.test192* *.test256* *.test288* *.test320* *.test384* *.test400* *.test512* *.test640* *.test768* *.test800* *.test864m4 *.test1024m4

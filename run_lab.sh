#!/bin/bash
#YBATCH -r epyc-7502_8
#SBATCH -N 1
#SBATCH -J ELPA
#SBATCH --time=24:00:00

source ~/.bashrc

source /etc/profile.d/modules.sh
module purge
module load gcc/12.2 cuda intel/2022/mkl cmake openmpi/4.0.5

# git clone https://gitlab.mpcdf.mpg.de/elpa/elpa.git
# cd elpa
# mkdir build

# ./autogen.sh
# CC=mpicc CXX=mpicxx LDFLAGS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl" ./configure --prefix=$PWD/build --enable-generic-kernels --disable-sse --disable-avx --disable-avx512-kernels --enable-avx2-kernels
# make -j
# make install

# cd ..

ELSES_ROOT=/home/sameer.deshmukh/ELSES_mat_calc-master
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PWD/elpa/build/lib/pkgconfig
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/elpa/build/lib

rm a.out

mpicxx -g $(pkg-config --cflags elpa) \
       -I$ELSES_ROOT -m64  -I"${MKLROOT}/include" \
        main.cpp -c -o main.o

mpicxx main.o $ELSES_ROOT/src/src.a $ELSES_ROOT/xmlf90-1.2g-elses/macros/lib/libflib.a \
       -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -Wl,--no-as-needed \
       -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 \
       -lgomp -lpthread -lm -ldl $(pkg-config --libs elpa) -lgfortran -o a.out

# ELSES matrix generation.
exec_supercell=$ELSES_ROOT/make_supercell_C60_FCCs_w_noise/a.out
exec_elses_xml_generate=$ELSES_ROOT/bin/elses-xml-generate

# Generate the points for the ELSES matrix.
ny=1
nz=1

mol_folder=$ELSES_ROOT/sample/sample_non_geno/C60_fcc2x2x2_disorder_expand_1x1x1
source_file=$mol_folder/C60_fcc2x2x2_20220727.xyz
fcc_xml_file=C60_fcc.xml

# Run the ELPA example with various matrix sizes.

export ELPA_DEFAULT_solver=ELPA_SOLVER_2STAGE
export OMP_NUM_THREADS=64
export OMP_PROC_BIND=close

for nx in 3 4 5 6 7 8 9 10 11 12 13; do
    # Generate the xml file from the source geometry depenending on the number of repetitions specified.
    $exec_supercell $nx $ny $nz $source_file

    cp C60_fcc.xyz $ELSES_ROOT/make_supercell_C60_FCCs_w_noise

    # generate config.xml.
    $exec_elses_xml_generate $ELSES_ROOT/make_supercell_C60_FCCs_w_noise/generate.xml $fcc_xml_file

    N=$(($nx * $ny * $nz * 1 * 1 * 1 * 32 * 60 * 4))

    mpirun --mca opal_warn_on_missing_libcuda 0  -n 1 ./a.out $N
done

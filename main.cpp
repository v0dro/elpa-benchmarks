#include <iostream>

#include "mpi.h"

extern "C" {
#include "elses.h"
}

extern "C" {
  /* Cblacs declarations: https://netlib.org/blacs/BLACS/QRef.html */
  void Cblacs_pinfo(int*, int*);
  void Cblacs_get(int CONTEXT, int WHAT, int*VALUE);
  void Cblacs_gridinit(int*, const char*, int, int);
  // returns the co-ordinates of the process number PNUM in PROW and PCOL.
  void Cblacs_pcoord(int CONTEXT, int PNUM, int* PROW, int* PCOL);
  void Cblacs_gridexit(int CONTEXT);
  void Cblacs_barrier(int, const char*);
  void Cblacs_exit(int CONTINUE);
  void Cblacs_gridmap(int* CONTEXT, int* USERMAP, const int LDUMAP,
                      const int NPROW, const int NPCOL);
  void Cblacs_gridinfo(int CONTEXT, int *NPROW, int *NPCOL,
                       int *MYROW, int *MYCOL);

  // calculate the number of rows and cols owned by process IPROC.
  // IPROC :: (local input) INTEGER
  //          The coordinate of the process whose local array row or
  //          column is to be determined.
  // ISRCPROC :: The coordinate of the process that possesses the first
  //             row or column of the distributed matrix. Global input.
  int numroc_(const int* N, const int* NB, const int* IPROC, const int* ISRCPROC,
              const int* NPROCS);
}

int MPIGRID[2];
int MPIRANK, MPISIZE, MYROW, MYCOL;
int BLACS_CONTEXT;

int main(int argc, char* argv[]) {
  long int N = atoll(argv[1]);

  std::cout << "N: " << N << std::endl;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIRANK);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISIZE);
  MPI_Dims_create(MPISIZE, 2, MPIGRID);

  Cblacs_get(-1, 0, &BLACS_CONTEXT );
  Cblacs_gridinit(&BLACS_CONTEXT, "Row", MPIGRID[0], MPIGRID[1]);
  Cblacs_pcoord(BLACS_CONTEXT, MPIRANK, &MYROW, &MYCOL);

  Cblacs_gridexit(BLACS_CONTEXT);
  Cblacs_exit(1);
  MPI_Finalize();

  return 0;
}

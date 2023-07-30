#include <iostream>
#include <vector>
#include <chrono>

#include <elpa/elpa.h>

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

  // init descriptor for scalapack matrices.
  void descinit_(int *desc,
                 const int *m,  const int *n, const int *mb, const int *nb,
                 const int *irsrc, const int *icsrc,
                 const int *BLACS_CONTEXT,
                 const int *lld, int *info);
}

int MPIGRID[2];
int MPIRANK, MPISIZE, MYROW, MYCOL;
int BLACS_CONTEXT;
const int ZERO = 0;

void assert_elpa_ok(const int error, const std::string& error_str,
                    bool& success, std::string& error_msg) {
  if (error != ELPA_OK) {
    success = false;
    error_msg = error_str;
  }
}

int
indxl2g(int indxloc, int nb, int iproc, int nprocs) {
  return nprocs * nb * ((indxloc - 1) / nb) +
    (indxloc-1) % nb + ((nprocs + iproc) % nprocs) * nb + 1;
}

int main(int argc, char* argv[]) {
  int N = atoi(argv[1]);
  int NB = 240;                 // default scalapack block size.

  std::cout << "N: " << N << std::endl;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIRANK);
  MPI_Comm_size(MPI_COMM_WORLD, &MPISIZE);
  MPI_Dims_create(MPISIZE, 2, MPIGRID);

  Cblacs_get(-1, 0, &BLACS_CONTEXT );
  Cblacs_gridinit(&BLACS_CONTEXT, "Row", MPIGRID[0], MPIGRID[1]);
  Cblacs_pcoord(BLACS_CONTEXT, MPIRANK, &MYROW, &MYCOL);

  const int local_nrows = numroc_(&N, &NB, &MYROW, &ZERO, &MPIGRID[0]);
  const int local_ncols = numroc_(&N, &NB, &MYCOL, &ZERO, &MPIGRID[1]);

  // Initialize local matrix
  init_elses_state();
  int info, desc[9];
  descinit_(desc, &N, &N, &NB, &NB, &ZERO, &ZERO, &BLACS_CONTEXT, &local_nrows, &info);
  std::vector<double> A((int64_t)local_nrows * (int64_t)local_ncols);

  for (int i = 0; i < local_nrows; i++) {
    for (int j = 0; j < local_ncols; j++) {
      // std::cout << "i: " << i << " j: " << j << std::endl;
      long int g_row = indxl2g(i + 1, NB, MYROW, MPIGRID[0]);
      long int g_col = indxl2g(j + 1, NB, MYCOL, MPIGRID[1]);
      double val;
      get_elses_matrix_value(&g_row, &g_col, &val);

      A[i + j * local_nrows] = val;
    }
  }


  std::cout << "----- DONE GEN ------ \n";

  std::vector<double> ev(N);
  bool elpa_success = true;
  std::string elpa_error_str = "";
  elpa_t handle;
  int error;

  error = elpa_init(ELPA_API_VERSION);
  assert_elpa_ok(error, "ELPA API version not supported",
                 elpa_success, elpa_error_str);
  handle = elpa_allocate(&error);
  assert_elpa_ok(error, "ELPA instance allocation failed",
                 elpa_success, elpa_error_str);

  /* Set parameters the matrix and it's MPI distribution */
  elpa_set(handle, "na", N, &error);
  elpa_set(handle, "local_nrows", local_nrows, &error);
  elpa_set(handle, "local_ncols", local_ncols, &error);
  elpa_set(handle, "nblk", NB, &error);
  elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);
  elpa_set(handle, "process_row", MYROW, &error);
  elpa_set(handle, "process_col", MYCOL, &error);
  assert_elpa_ok(error, "ELPA matrix initialization failed",
                 elpa_success, elpa_error_str);
  error = elpa_setup(handle);
  assert_elpa_ok(error, "ELPA setup failed",
                 elpa_success, elpa_error_str);

  double elpa_time;
  if (elpa_success) {
    std::cout << "Now calculating eigenvalues...\n";
    auto start_time = std::chrono::system_clock::now();
    elpa_eigenvalues(handle, A.data(), ev.data(), &error);
    auto stop_time = std::chrono::system_clock::now();
    elpa_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count();
  }

  elpa_deallocate(handle, &error);
  elpa_uninit(&error);
  Cblacs_gridexit(BLACS_CONTEXT);
  Cblacs_exit(1);
  MPI_Finalize();

  std::cout << "ELPA_TIME: " << N << "," << elpa_time << std::endl;

  return 0;
}

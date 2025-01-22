
#ifndef SPECLE_GLOBALS_HPP
#define SPECLE_GLOBALS_HPP

#include <cstddef>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <hdf5.h>

#include <list>
#include <functional>

static const double clight = 299792458.0;

extern const std::complex<double> ii;

extern ptrdiff_t N0, N1;
extern int NOm;
extern double lDomain;
extern double laserLambda;
extern double laserK;
extern double omMin, omMax;
extern double om;

extern double intensity;
extern std::complex<double> phase;

extern fftw_complex *data;

extern fftw_plan plan, planbw;
extern ptrdiff_t local_n0, local_0_start;

// iteration counter for loop commands
extern int local_om_pos;
extern int local_x_pos, local_y_pos;
extern int local_xy_index;

extern int mpi_rank, mpi_size;

extern hid_t hdf5FileId;
extern hid_t hdf5DxplId;
extern hid_t hdfDatasetRealId;
extern hid_t hdfDatasetImagId;
extern hid_t hdfDataSpaceId;

extern hid_t hdf5ReadFileId;
extern hid_t hdf5ReadDxplId;
extern hid_t hdfReadDatasetRealId;
extern hid_t hdfReadDatasetImagId;
extern hid_t hdfReadDataSpaceId;

extern std::list<std::function<void()>> stepOmCommands;

extern bool debugFlag;

/**
 * A simple 2d array of doubles
 */
class Array2D {
    public:
        Array2D(int N0, int N1);
        ~Array2D();

        double **data;
        int N0, N1;
};

#endif
#include "globals.hpp"


const std::complex<double> ii(0, 1);

fftw_plan plan, planbw;
ptrdiff_t local_n0, local_0_start;
int local_om_pos = 0;
int local_x_pos = 0, local_y_pos = 0;
int local_xy_index = 0;

ptrdiff_t N0 = 16000, N1 = 16000;
int NOm = 3000;
double lDomain = 500.0;
double laserLambda = 1.06e-3;
double laserK = 2*M_PI/laserLambda;
double omMin = 2*M_PI*clight/850e-6;
double omMax = 2*M_PI*clight/750e-6;

double om = omMin;

double intensity = 0;
std::complex<double> phase{1.0, 0.0};

fftw_complex *data = nullptr;

hid_t hdf5FileId = -1;
hid_t hdf5DxplId = -1;
hid_t hdfDatasetRealId = -1;
hid_t hdfDatasetImagId = -1;
hid_t hdfDataSpaceId = -1;

hid_t hdf5ReadDxplId = -1;
hid_t hdf5ReadFileId = -1;
hid_t hdfReadDatasetRealId = -1;
hid_t hdfReadDatasetImagId = -1;
hid_t hdfReadDataSpaceId = -1;

int mpi_rank, mpi_size;

int counter = 0;

std::list<std::function<void()>> stepOmCommands;

bool debugFlag = false;

Array2D::Array2D(int N0, int N1) : N0(N0), N1(N1) {
    data = new double*[N0];
    for (int i = 0; i < N0; ++i) {
        data[i] = new double[N1];
    }
}

Array2D::~Array2D() {
    for (int i = 0; i < N0; ++i) {
        delete[] data[i];
    }
    delete[] data;
}
#include "../globals.hpp"
#include "../range.hpp"
#include "../transfer.hpp"
#include "../parse.hpp"

#include <list>
#include <string>
#include <iostream>

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

XYRange range{N0, N1, -lDomain, lDomain, -lDomain, lDomain};

void allocateData() {
    if (data != nullptr) {
        fftw_destroy_plan(plan);
        fftw_destroy_plan(planbw);
        fftw_free(data);
    }

    /* get local data size and allocate */
    ptrdiff_t alloc_local = fftw_mpi_local_size_2d(
        N0, N1, MPI_COMM_WORLD, &local_n0, &local_0_start
    );
    data = fftw_alloc_complex(alloc_local);

    /* create plan for in-place forward DFT */
    plan = fftw_mpi_plan_dft_2d(
        N0, N1, data, data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE
    );    

    planbw = fftw_mpi_plan_dft_2d(
        N0, N1, data, data, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE
    );  

    range = XYRange(N0, N1, -lDomain, lDomain, -lDomain, lDomain);
}

/**
 * Fill the field with a function
 * 
 * The function is called with the coordinates of the point in the field.
 */
template<typename Func>
void fillField(fftw_complex *data, Func func, XYRange range) {
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            std::complex<double> z = func(range.getCoord(i, j));
            double *c = data[i*N1 + j];
            c[0] = z.real();
            c[1] = z.imag();
        }
    }
}

void setupSGBeamCommand(TokenStream& tokens) {
    allocateData();

    double cx = tokens.nextDouble();
    double cy = tokens.nextDouble();
    double radius = tokens.nextDouble();
    int order = tokens.nextInt();

    if (mpi_rank == 0) {
        std::cout << "creating super-Gaussian beam: " << cx << ", " << cy << ", " << radius << ", " << order << std::endl;
    }
    SuperGaussian beam{{cx, cy}, radius, order};
    auto envelopedBeam = [&beam](Coord c) {
        return beam(c) * intensity * phase;
    };
    fillField(data, envelopedBeam, range);
}
#include "../globals.hpp"
#include "../range.hpp"
#include "../transfer.hpp"
#include "../io.hpp"
#include "../parse.hpp"

#include "../data/data-util.hpp"

#include <list>
#include <string>
#include <iostream>

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>


template<typename Func>
void multiplyField(fftw_complex *data, Func func, XYRange range) {
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            std::complex<double> z = func(range.getCoord(i, j));
            double *c = data[i*N1 + j];
            z = z*std::complex<double>{c[0], c[1]};
            c[0] = z.real();
            c[1] = z.imag();
        }
    }
}

template<typename Func>
void multiplyFieldK(fftw_complex *data, Func func, XYRange range) {
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            std::complex<double> z = func(range.getFFTCoord(i, j));
            double *c = data[i*N1 + j];
            z = z*std::complex<double>{c[0], c[1]};
            c[0] = z.real();
            c[1] = z.imag();
        }
    }
}

template<typename Func>
void propagate(fftw_complex *data, Func func, XYRange range) {
    fftw_execute(plan);
    multiplyFieldK(data, func, range);
    fftw_execute(planbw);
    renorm(data);
}

void renormSingle(fftw_complex *data) {
    double norm = 1.0/sqrt(N0*N1);
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double *c = data[i*N1 + j];
            c[0] *= norm;
            c[1] *= norm;
        }
    }
}

void renormOm(fftw_complex *data) {
    double norm = 1.0/sqrt(NOm);
    for (ptrdiff_t i = 0; i < NOm; ++i) {
        double *c = data[i];
        c[0] *= norm;
        c[1] *= norm;
    }
}

void preshift(fftw_complex *data) {
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = (i+local_0_start)%2; j < N1; j+=2) {
            double *c = data[i*N1 + j];
            c[0] = -c[0];
            c[1] = -c[1];
        }
    }
}

void fraunhofer(fftw_complex *data, double distance) {
    preshift(data);
    fftw_execute(planbw);
    renormSingle(data);
    double lSpotDomain = N0*laserLambda*distance / (4*lDomain);

    if (mpi_rank == 0) {
        std::cout << "Fraunhofer spot domain size: " << lSpotDomain << std::endl;
    }
}

void inverseFFTOm(fftw_complex *data) {
    fftw_execute(planbw);
    renormOm(data);
}

void mutateParabolaCommand(TokenStream& tokens) {
    double focalLength = tokens.nextDouble();

    if (mpi_rank == 0) {
        std::cout << "parabola: " << focalLength << std::endl;
    }
    TransMirrorParabola transMirror{focalLength};
    multiplyField(data, transMirror, range);
}

void mutateGridNearestNeighbourCommand(TokenStream& tokens) {
    std::string key = tokens.next();
    double size = tokens.nextDouble();
    double multiplier = tokens.nextDouble();
    std::string mode = tokens.next();
    if (mode == "path") {
        multiplier = multiplier*om/clight;
    } else if (mode != "simple") {
        throw std::runtime_error("Unknown mode: " + mode);
    }
    if (mpi_rank == 0) {
        std::cout << "Nearest Neighbour: " << key << ", " << mode << " -- " << multiplier << std::endl;
    }
    if (dataStore.find(key) == dataStore.end()) {
        throw std::runtime_error("Data store - key not found: " + key);
    }
    Array2D *array = dataStore[key].get();
    XYRange arange{array->N0, array->N1, -size, size, -size, size};
    TransDataNearestNeighbour transDataNN{array, arange, multiplier};
    multiplyField(data, transDataNN, range);
}

void mutatePropagateCommand(TokenStream& tokens) {
    double distance = tokens.nextDouble();

    if (mpi_rank == 0) {
        std::cout << "propagating: " << distance << std::endl;
    }

    TransFree transFree{distance};
    propagate(data, transFree, range);
}

void mutateFocalFarField(TokenStream& tokens) {
    double distance = tokens.nextDouble();

    if (mpi_rank == 0) {
        std::cout << "Fraunhofer focal spot: " << distance << std::endl;
    }

    fraunhofer(data, distance);
}

void mutateOmToTime(TokenStream& tokens) {
    if (mpi_rank == 0) {
        std::cout << "Inverse Fourier transform: om -> time" << std::endl;
    }

    inverseFFTOm(data);
}
#ifndef SPECLE_RANGE_HPP
#define SPECLE_RANGE_HPP

#include "globals.hpp"

#include <cstddef>
#include <cmath>
// #include <algorithm>

struct Coord {
    double x, y;
};

struct XYRange {
    ptrdiff_t Nx, Ny;
    ptrdiff_t Nxh, Nyh;
    double Xmin, Xmax;
    double Ymin, Ymax;
    double Dx, Dy;
    double Kx, Ky;
    double KxMin, KyMin;

    XYRange(
        ptrdiff_t Nx, 
        ptrdiff_t Ny,
        double Xmin, // -lDomain
        double Xmax, //  lDomain
        double Ymin, 
        double Ymax
    ) : Nx(Nx), Ny(Ny), Nxh(Nx/2), Nyh(Ny/2), Xmin(Xmin), Xmax(Xmax), Ymin(Ymin), Ymax(Ymax)
    {
        Dx = (Xmax - Xmin)/double(Nx); // (2*lDomain)/N
        Dy = (Ymax - Ymin)/double(Ny);
        Kx = 2.0*M_PI/(Xmax - Xmin); // pi/(lDomain)
        Ky = 2.0*M_PI/(Ymax - Ymin);
        KxMin = -M_PI/Dx; // -pi*N/(2*lDomain)
        // KxMax = KxMin + Kx*N = -pi*N/(2*lDomain) + pi*N/(lDomain) 
        //   = pi*N/(2*lDomain) = -KxMin
        KyMin = -M_PI/Dy;
    }

    Coord getCoord(ptrdiff_t i, ptrdiff_t j) {
        return Coord{Dx*(i + local_0_start) + Xmin, Dy*j + Ymin};
    }

    /**
     * Get the nearest neighbor to a coordinate.
     * 
     * The coordinate is assumed to be for a cell-centered grid.
     * This function is NOT the inverse of getCoord.
     * 
     * For all coordinates in the range Xmin <= c.x <= Xmax and Ymin <= c.y <= Ymax
     * this function returns a pair of indices i, j such that 
     * 0 <= i < Nx and 0 <= j < Ny
     */
    std::pair<ptrdiff_t, ptrdiff_t> getNearestNeighbour(Coord c) {
        ptrdiff_t i = std::max(0l, std::min(Nx-1, ptrdiff_t((c.x - Xmin)/Dx)));
        ptrdiff_t j = std::max(0l, std::min(Ny-1, ptrdiff_t((c.y - Ymin)/Dy)));
        return std::pair<ptrdiff_t, ptrdiff_t>{i, j};
    }

    bool isInside(Coord c) {
        return Xmin <= c.x && c.x <= Xmax && Ymin <= c.y && c.y <= Ymax;
    }

    Coord getFFTCoord(ptrdiff_t i, ptrdiff_t j) {
        ptrdiff_t is = (i + local_0_start + Nxh) % (2*Nxh);
        ptrdiff_t js = (j + Nyh) % (2*Nyh);
        return Coord{Kx*is + KxMin, Ky*js + KyMin};
        // return Coord{double(is), double(js)};
    }
};

// defined in setup.cpp
extern XYRange range;

#endif
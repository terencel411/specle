#ifndef SPECLE_TRANSFER_HPP
#define SPECLE_TRANSFER_HPP

#include "globals.hpp"
#include "range.hpp"

#include <complex>
#include <iostream>

// =======================
// data provided by Robbie
//
// units are mm


// ======================

extern const double cavityLength;

extern const double mirrorRadCurvature;
extern const double mirrorRadius;
extern const double mirrorRadCurvature_2;

struct Distance {
    double cx, cy;
    double operator()(double x, double y) {
        return (x-cx)*(x-cx) + (y-cy)*(y-cy);
    }
};

struct SuperGaussian {
    Distance dist;
    double radius;
    int order;
    std::complex<double> operator()(Coord c) {
        return exp(
            - 2.0*pow(dist(c.x, c.y)/(radius*radius), order)
        );
    }
};

struct MaskMirror {
    Distance dist;
    double radius;
    std::complex<double> operator()(Coord c) {
        return dist(c.x, c.y) > radius*radius ? 0 : 1;
    }
};

struct TransFree {
    double distance;
    std::complex<double> operator()(Coord c) {
        return exp(
            ii*sqrt(laserK*laserK - c.x*c.x - c.y*c.y)*distance
        );
    }
};

struct TransMirror {
    std::complex<double> operator()(Coord c) {
        double R2 = c.x*c.x + c.y*c.y;
        return exp(-ii*laserK*(
            2*R2/(mirrorRadCurvature*(1+sqrt(1 - R2/mirrorRadCurvature_2)))
        ));
    }
};

struct TransMirrorParabola {
    double focalLength;
    std::complex<double> operator()(Coord c) {
        double R2 = c.x*c.x + c.y*c.y;
        return exp(-ii*laserK*(
            R2/(2*focalLength)
        ));
    }
};

struct TransMirrorDefocus {
    double maxHeight;
    std::complex<double> operator()(Coord c) {
        double R2 = c.x*c.x + c.y*c.y;
        return exp(-ii*laserK*(
            2*R2/(mirrorRadCurvature*(1+sqrt(1 - R2/mirrorRadCurvature_2)))
            + maxHeight*R2/(mirrorRadius*mirrorRadius)
        ));
    }
};

struct TransMirrorPerturb {
    double lambda;
    double maxHeight;
    double phx = 0.0;
    double phy = 0.0;
    std::complex<double> operator()(Coord c) {
        double R2 = c.x*c.x + c.y*c.y;
        return exp(-ii*laserK*(
            2*R2/(mirrorRadCurvature*(1+sqrt(1 - R2/mirrorRadCurvature_2)))
            + maxHeight*sin(2*M_PI*c.x/lambda + phx)*sin(2*M_PI*c.y/lambda + phy)
        ));
    }
};

/**
 * Uses the data stored in an array to provide the phase shift at a given point.
 * The phase shift is calculated by finding the nearest neighbour to the point
 * and using the phase shift of that point.
 * 
 */
struct TransDataNearestNeighbour
{
    Array2D *data;
    XYRange range;
    double phaseMultiplier;

    std::complex<double> operator()(Coord c) {
        if (!range.isInside(c)) {
            return 0;
        }
        auto p = range.getNearestNeighbour(c);
        return exp(ii*phaseMultiplier*data->data[p.first][p.second]);
    }
};


#endif // SPECLE_TRANSFER_HPP
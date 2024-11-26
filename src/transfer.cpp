#include "transfer.hpp"

// =======================
// data provided by Robbie
//
// units are mm

double mrr_l_rad = 6.35;
double mrr_r_rad = 150;
int sg_exp = 20;
double sg_radius_l = 0.80*mrr_l_rad;

// =======================

const double sourceCentreX = 106.076017177982;
const double sourceCentreY = 0;
const double sourceRadius = sg_radius_l;
const double cavityLength = 300e3;
const int superGaussianOrder = sg_exp/2;

const double mirrorRadCurvature = 300e3;
const double mirrorRadius = mrr_r_rad;

const double mirrorRadCurvature_2 = mirrorRadCurvature*mirrorRadCurvature;


#include "../globals.hpp"
#include "../parse.hpp"

#include <list>
#include <string>
#include <iostream>

void updateLocalXYIndex() {
    local_x_pos = local_xy_index % N0;
    local_y_pos = local_xy_index / N0;
}

void setSizeCommand(TokenStream& tokens) {
    N0 = tokens.nextInt();
    N1 = tokens.nextInt();
    updateLocalXYIndex();

    if (mpi_rank == 0) {
        std::cout << "setting size: " << N0 << ", " << N1 << std::endl;
    }
}

void setOmSizeCommand(TokenStream& tokens) {
    NOm = tokens.nextInt();
    local_om_pos = 0;
    if (mpi_rank == 0) {
        std::cout << "setting size in omega: " << NOm << std::endl;
    }
}

void setDomainCommand(TokenStream& tokens) {
    lDomain = tokens.nextDouble();
    if (mpi_rank == 0) {
        std::cout << "setting domain: " << lDomain << std::endl;
    }
}

void setLambdaCommand(TokenStream& tokens) {
    laserLambda = tokens.nextDouble();
    laserK = 2*M_PI/laserLambda;
    if (mpi_rank == 0) {
        std::cout << "setting lambda: " << laserLambda << " " << laserK << std::endl;
    }
}

void setLambdaDomainCommand(TokenStream& tokens) {
    double lambdaMin = tokens.nextDouble();
    double lambdaMax = tokens.nextDouble();

    omMin = 2*M_PI*clight/lambdaMax;
    omMax = 2*M_PI*clight/lambdaMin;
    om = omMin;
    if (mpi_rank == 0) {
        std::cout << "setting lambda domain: " << lambdaMin << " " << lambdaMax << std::endl;
        std::cout << "     => omega domain: " << omMin << " " << omMax << std::endl;
    }
}

double sgSpectrumOm0 = 2*M_PI*clight/800e-6;
double sgSpectrumDOm = 2*M_PI*clight/40e-6;
int sgSpectrumOrder = 8;
double spectrumCohBW = 2*sgSpectrumDOm;

void stepOmSGSpectrum() {
    // std::cout << "stepOmSGSpectrum" << std::endl;
    // std::cout << "  om " << om << std::endl;
    // std::cout << "  sgSpectrumOm0 " << sgSpectrumOm0 << std::endl;
    // std::cout << "  sgSpectrumDOm " << sgSpectrumDOm << std::endl;
    // std::cout << "  sgSpectrumOrder " << sgSpectrumOrder << std::endl;
    // std::cout << "  delta " << om - sgSpectrumOm0 << std::endl;
    // std::cout << "  deltaNorm " << (om - sgSpectrumOm0) / sgSpectrumDOm << std::endl;
    intensity =  exp(-pow(double(om - sgSpectrumOm0)/sgSpectrumDOm, 2*sgSpectrumOrder));
    if (mpi_rank == 0) {
       std::cout << "intensity[" << om << "] = " << intensity << std::endl;
    }
}

void stepOmRWSpectrum() {
    double dOm = (omMax - omMin) / double(NOm);
    double dPhase = 2*M_PI * rand() / double(RAND_MAX) * sqrt(dOm / spectrumCohBW);
    phase *= exp(ii * dPhase);
    if (mpi_rank == 0) {
        std::cout << "phase[" << om << "] = " << phase << std::endl;
    }
}

void setSGSpectrumCommand(TokenStream& tokens) {
    double lambdaMin = tokens.nextDouble();
    double lambdaMax = tokens.nextDouble();
    sgSpectrumOrder = tokens.nextInt();

    double sgOmMin = 2*M_PI*clight/lambdaMax;
    double sgOmMax = 2*M_PI*clight/lambdaMin;
    sgSpectrumOm0 = 0.5*(sgOmMin + sgOmMax);
    sgSpectrumDOm  = 0.5*(sgOmMax - sgOmMin);
    if (mpi_rank == 0) {
        std::cout << "setting SG spectrum: " << sgSpectrumOm0 << " " << sgSpectrumDOm << " " << sgSpectrumOrder << std::endl;
    }
    stepOmSGSpectrum();
    stepOmCommands.push_back(stepOmSGSpectrum);
}

void setRWSpectrumPhaseCommand(TokenStream& tokens) {
    double spectrumCohBWLambda = tokens.nextDouble();
    spectrumCohBW = 2*M_PI*clight/spectrumCohBWLambda;
    if (mpi_rank == 0) {
        std::cout << "setting random walk spectrum (thermal): " << spectrumCohBWLambda << " " << spectrumCohBW << std::endl;
    }
    stepOmRWSpectrum();
    stepOmCommands.push_back(stepOmRWSpectrum);
}

void initLoopOmCommand() {
    local_om_pos = 0;
    for (auto &cmd : stepOmCommands) {
        cmd();
    }
}

bool loopOmCommand(TokenStream& tokens) {
    if (++local_om_pos < NOm) {
        om = omMin + (omMax - omMin) * local_om_pos / double(NOm);
        if (mpi_rank == 0) {
            std::cout << "omega[" << local_om_pos << "] = " << om << std::endl;
        }
        for (auto &cmd : stepOmCommands) {
            cmd();
        }
        return true;
    } else {
        local_om_pos = 0;
        om = omMin;
        return false;
    }
}

bool loopXYCommand(TokenStream& tokens) {
    if (++local_xy_index < N0*N1) {
        updateLocalXYIndex();
        // if (mpi_rank == 0) {
        //     std::cout << "(x,y) = (" << local_xy_index << ", " << local_xy_index << ")" << std::endl;
        // }
        for (auto &cmd : stepOmCommands) {
            cmd();
        }
        return true;
    } else {
        local_xy_index = mpi_rank;
        updateLocalXYIndex();
        return false;
    }
}
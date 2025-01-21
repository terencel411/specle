#include "range.hpp"
#include "globals.hpp"
#include "parse.hpp"
#include "io.hpp"
#include "transfer.hpp"
#include "commands/basicio.hpp"
#include "commands/hdfio.hpp"
#include "commands/mutate.hpp"
#include "commands/setters.hpp"
#include "commands/setup.hpp"

#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <cmath>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <complex>
#include <cstdlib>
#include <functional>

// process with
// ./bin/dataToAscii beam-out 7700 6000 600 600


double randomDouble() {
    return double(rand())/double(RAND_MAX);
}

double randomPhase() {
    return 2.0*M_PI*randomDouble();
}

double spectrum(Coord coord) {
    double sum = 0.0;
    for (int kx=0; kx<N0/2; ++kx) {
        double kxn = (double)kx;
        double kyn = kxn*kxn/(double)N0;
        sum += 1/(kxn*kxn + 1) * sin((coord.x*kxn + coord.y*kyn)*2.0*M_PI);
    }
    return sum;
}


double real(fftw_complex c) {
    return c[0];
}

double imag(fftw_complex c) {
    return c[1];
}

double absolute(fftw_complex c) {
    return c[0]*c[0] + c[1]*c[1];
}

double logAbsolute(fftw_complex c) {
    return log(c[0]*c[0] + c[1]*c[1]);
}

void parseInstructions(std::string filename) {
    Parser parser;
    
    parser.addCommand("size", setSizeCommand);
    parser.addCommand("omSize", setOmSizeCommand);
    parser.addCommand("domain", setDomainCommand);
    parser.addCommand("lambda", setLambdaCommand);
    parser.addCommand("lambdaDomain", setLambdaDomainCommand);
    parser.addCommand("sgspectrum", setSGSpectrumCommand);
    parser.addCommand("rwphase", setRWSpectrumPhaseCommand);
    parser.addLoopCommand("loopOm", loopOmCommand);
    parser.addLoopCommand("loopXY", loopXYCommand);

    parser.addCommand("sgbeam", setupSGBeamCommand);

    parser.addCommand("parabola", mutateParabolaCommand);
    parser.addCommand("gridNN", mutateGridNearestNeighbourCommand);
    parser.addCommand("propagate", mutatePropagateCommand);
    parser.addCommand("focal", mutateFocalFarField);
    parser.addCommand("omToTime", mutateOmToTime);

    parser.addCommand("readAscii", readAsciiCommand);
    parser.addCommand("write", writeBinaryCommand);
    parser.addCommand("writeFFT", writeFFTCommand);

    parser.addCommand("createHDF", createHdfCommand);
    parser.addCommand("openHDF", openHdfCommand);
    parser.addCommand("writeHDF", writeHdfCommand);
    parser.addCommand("writeTimeHDF", writeHdfTimelineCommand);
    parser.addCommand("readTimeHDF", readHdfTimelineCommand);
    parser.addCommand("closeHDF", closeHdfCommand);
    parser.addCommand("closeHDFRead", closeHdfReadCommand);

    parser.parse(filename);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    local_xy_index = mpi_rank;

    fftw_mpi_init();

    parseInstructions(argv[1]);
    fftw_destroy_plan(plan);

    MPI_Finalize();
}

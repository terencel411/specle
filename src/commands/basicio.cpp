
#include "../globals.hpp"
#include "../io.hpp"
#include "../parse.hpp"

#include "../data/data-util.hpp"

#include <list>
#include <string>
#include <iostream>

std::string mpiFilename(std::string base) {
    std::ostringstream name;
    name << "output/" << base << "_" << mpi_rank << ".out";
    return name.str();
}

std::string mpiMetaFilename(std::string base) {
    std::ostringstream name;
    name << "output/" << base << "_descr.out";
    return name.str();
}

void writeFFT(fftw_complex *data, std::string filename) {
    fftw_execute(plan);
    
    writeBinary(data, mpiFilename(filename));
    writeMetaInfo(mpiMetaFilename(filename));
    
    fftw_execute(planbw);
    renorm(data);
}

void readAsciiCommand(TokenStream& tokens) {
    std::string filename = tokens.next();
    std::string key = tokens.next();

    if (mpi_rank == 0) {
        std::cout << "reading: " << key << " --> " << filename << std::endl;
    }

    readAscii(filename, key);
}

void writeBinaryCommand(TokenStream& tokens) {
    std::string filename = tokens.next();

    if (mpi_rank == 0) {
        std::cout << "writing: " << filename << std::endl;
    }

    writeBinary(data, mpiFilename(filename));
    writeMetaInfo(mpiMetaFilename(filename));
}

void writeFFTCommand(TokenStream& tokens) {
    std::string filename = tokens.next();

    if (mpi_rank == 0) {
        std::cout << "writing FFT: " << filename << std::endl;
    }
    writeFFT(data, filename);
}
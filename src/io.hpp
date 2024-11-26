#ifndef SPECLE_IO_HPP
#define SPECLE_IO_HPP

#include "range.hpp"
#include "globals.hpp"

#include <fftw3.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>

extern std::map<std::string, std::shared_ptr<Array2D>> dataStore;

void writeAscii(fftw_complex *data, XYRange range, std::string filename);
void writeBinary(fftw_complex *data, std::string filename);
void writeMetaInfo(std::string filename);
void readAscii(std::string filename, std::string key);
void createHdf5File(std::string fname);
void openHdf5File(std::string fname);
void writeHdf5Dataset(fftw_complex *data, XYRange range);
void writeHdf5TimelineDataset(fftw_complex *data);
void readHdf5TimelineDataset(fftw_complex *data);
void closeHdf5File();
void closeHdf5ReadFile();

template<typename Func>
void writePGM(fftw_complex *data, Func func, std::string filename) 
{
    double maxVal = 0.0;
    double minVal = 1e100;
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double amp = func(data[i*N1 + j]);
            maxVal = fmax(maxVal, amp);
            minVal = fmin(minVal, amp);
        }
    }

    std::cerr << "PNM: min = " << minVal << ";  max = " << maxVal << std::endl;

    std::ofstream out(filename);

    out << "P2" << std::endl
        << local_n0 << " " << N1 << std::endl
        << 255 << std::endl;
    
    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double amp = func(data[i*N1 + j]);
            amp = 255*(amp - minVal)/(maxVal - minVal);
            out << int(amp) << " ";
        }
    }
}


#endif // SPECLE_IO_HPP

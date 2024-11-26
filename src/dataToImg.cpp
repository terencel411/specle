
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

const size_t N0 = 16000, N1 = 16000;
const size_t local_n0 = 667;
const size_t local_n0f = N0 - 23*local_n0;

double *data;

std::string mpiFilename(std::string base, size_t count) {
    std::ostringstream name;
    name << "output/" << base << "_" << count << ".out";
    return name.str();
}

std::string pgmFilename(std::string base) {
    std::ostringstream name;
    name << base << ".pgm";
    return name.str();
}

template<typename Func>
void writePGM(Func func, std::string filename) 
{
    double maxVal = 0.0;
    double minVal = 1e100;
    for (ptrdiff_t i = 0; i < N0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double amp = func(&data[2*(i*N1 + j)]);
            maxVal = fmax(maxVal, amp);
            minVal = fmin(minVal, amp);
        }
    }

    std::cerr << "PNM: min = " << minVal << ";  max = " << maxVal << std::endl;

    std::ofstream out(filename);

    out << "P2" << std::endl
        << N0 << " " << N1 << std::endl
        << 255 << std::endl;
    
    for (ptrdiff_t i = 0; i < N0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double amp = func(&data[2*(i*N1 + j)]);
            amp = 255*(amp - minVal)/(maxVal - minVal);
            out << int(amp) << " ";
        }
    }
}

void readSlice(size_t count, size_t length, std::string filename) 
{
    std::cout << "reading slice " << count << std::endl;
    std::ifstream in(filename, std::ios::out | std::ios::binary);

    in.read((char *) &(data[2*count*local_n0*N1]), 2*length*N1*sizeof(double));

    in.close();
}

double absolute(double *c) {
    return c[0]*c[0] + c[1]*c[1];
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cout << "Must specify filename" << std::endl;
    }

    data = new double[2*N0*N1];

    for (size_t i=0; i<23; ++i) {
        readSlice(i, local_n0, mpiFilename(argv[1], i));
    }
    readSlice(23, local_n0f, mpiFilename(argv[1], 23));

    writePGM(absolute, pgmFilename(argv[1]));
}

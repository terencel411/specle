
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

long N0 = 16000, N1 = 16000;
long *local_n0;
long *local_n0_start;

int mpi_size;

size_t xmin, xmax;
size_t ymin, ymax;

double *data;

std::string mpiFilename(std::string base, size_t count) {
    std::ostringstream name;
    name << "output/" << base << "_" << count << ".out";
    return name.str();
}

std::string mpiMetaFilename(std::string base) {
    std::ostringstream name;
    name << "output/" << base << "_descr.out";
    return name.str();
}

std::string outFilename(std::string base) {
    std::ostringstream name;
    name << base << ".dat";
    return name.str();
}

template<typename Func>
void writeAscii(size_t sliceXMin, size_t sliceXMax, Func func, std::ofstream &out) 
{
    size_t rowLo = std::max(xmin, sliceXMin);
    size_t rowHi = std::min(xmax, sliceXMax);
    std::cout << "  writing rows " << rowLo << " -- " << rowHi << std::endl;
    for (ptrdiff_t i = rowLo; i < rowHi; ++i) {
        for (ptrdiff_t j = ymin; j < ymax; ++j) {
            double amp = func(&data[2*((i - sliceXMin)*N1 + j)]);
            out << i << " " << j << " " << amp << std::endl;
        }
        out << std::endl;
    }
}

void readSlice(size_t count, size_t length, std::string filename) 
{
    std::cout << "  reading slice " << count << std::endl;
    std::ifstream in(filename, std::ios::out | std::ios::binary);

    in.read((char *)data, 2*length*N1*sizeof(double));

    in.close();
}

/**
 * Reads the data about the domain decomposition from all processes.
 * 
 * The data is stored in a file with the following format:
 * 
 * # Comment
 * N0 N1
 * mpi_size
 * local_n0_start[0] ... local_n0_start[mpi_size-1]
 * local_n0[0] ... local_n0[mpi_size-1]
 * 
 * @param filename 
 */
void readDataDescription(std::string filename) 
{
    std::cout << "Reading data description from " << filename << std::endl;
    std::ifstream in(filename);
    std::string line;
    std::getline(in, line);
    while (line[0] == '#') { // skip comment
        std::getline(in, line);
    }
    // read dimensions
    std::istringstream dim(line);
    dim >> N0 >> N1;
    std::cout << "  N0 = " << N0 << ", N1 = " << N1 << std::endl;
    std::getline(in, line);
    std::istringstream mpiSize(line);
    mpiSize >> mpi_size;
    std::cout << "  mpi_size = " << mpi_size << std::endl;
    local_n0 = new long[mpi_size];
    local_n0_start = new long[mpi_size];

    std::getline(in, line);
    std::istringstream local(line);
    for (int i=0; i<mpi_size; ++i) {
        local >> local_n0_start[i];
    }

    std::getline(in, line);
    std::istringstream localN0(line);
    for (int i=0; i<mpi_size; ++i) {
        localN0 >> local_n0[i];
    }
}

double absolute(double *c) {
    return c[0]*c[0] + c[1]*c[1];
}

double phase(double *c) {
    return atan2(c[1], c[0]);
}

int main(int argc, char **argv)
{
    if (argc < 6) {
        std::cout << "Format dataToAscii <filename-base> <xmin> <ymin> <sizex> <sizey> [abs|phase]" << std::endl;
    }

    std::string filename{argv[1]};
    xmin = std::stoi(argv[3]);
    ymin = std::stoi(argv[2]);
    xmax = xmin + std::stoi(argv[5]);
    ymax = ymin + std::stoi(argv[4]);
    bool writeAbs = true;
    if (argc > 6) {
        writeAbs = std::string(argv[6]) == "abs";
        if (!writeAbs && std::string(argv[6]) != "phase") {
            std::cout << "Unknown option " << argv[6] << std::endl;
            return 1;
        }
    }

    readDataDescription(mpiMetaFilename(filename));
    std::ofstream outAscii(outFilename(filename));

    std::cout << "Reading from " << mpi_size << " slices" << std::endl;

    std::vector<size_t> sliceOrder(mpi_size);
    for (size_t i=0; i<mpi_size; ++i) {
        sliceOrder[i] = i;
    }

    // sort slices by local_n0_start
    std::sort(sliceOrder.begin(), sliceOrder.end(), [&](size_t a, size_t b) {
        return local_n0_start[a] < local_n0_start[b];
    });

    // get maximum slice length and allocate memory
    long maxSliceLength = 0;
    for (size_t i=0; i<mpi_size; ++i) {
        maxSliceLength = std::max(maxSliceLength, local_n0[sliceOrder[i]]);
    }
    data = new double[2*maxSliceLength*N1];

    for (size_t i=0; i<mpi_size; ++i) {
        long sliceXMin = local_n0_start[sliceOrder[i]];
        long sliceXMax = local_n0_start[sliceOrder[i]] + local_n0[sliceOrder[i]];
        long sliceXLength = local_n0[sliceOrder[i]];

        std::cout << i << ": slice " << sliceOrder[i] << std::endl;

        if (xmin < sliceXMax && xmax >= sliceXMin) {
            readSlice(sliceOrder[i], sliceXLength, mpiFilename(filename, sliceOrder[i]));
            if (writeAbs)
                writeAscii(sliceXMin, sliceXMax, absolute, outAscii);
            else
                writeAscii(sliceXMin, sliceXMax, phase, outAscii);
        }
    }
}

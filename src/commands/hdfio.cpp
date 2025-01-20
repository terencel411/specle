#include "../globals.hpp"
#include "../io.hpp"
#include "../parse.hpp"

#include <list>
#include <string>
#include <iostream>

void allocateTimelineData() {
    if (data != nullptr) {
        fftw_destroy_plan(plan);
        fftw_destroy_plan(planbw);
        fftw_free(data);
    }

    data = fftw_alloc_complex(NOm);

    /* create plan for in-place forward DFT */
    plan = fftw_plan_dft_1d(
        NOm, data, data, FFTW_FORWARD, FFTW_ESTIMATE
    );    

    planbw = fftw_plan_dft_1d(
        NOm, data, data, FFTW_BACKWARD, FFTW_ESTIMATE
    );  
}

void createHdfCommand(TokenStream& tokens) {
    std::string filename = tokens.next();

    if (mpi_rank == 0) {
        std::cout << "create HDF5 file: " << filename << std::endl;
    }
    createHdf5File(filename);
}


void openHdfCommand(TokenStream& tokens) {
    std::string filename = tokens.next();

    if (mpi_rank == 0) {
        std::cout << "opening HDF5 file: " << filename << std::endl;
    }
    openHdf5File(filename);
}

void writeHdfCommand(TokenStream& tokens) {
    if (mpi_rank == 0) {
        std::cout << "writing data to HDF" << std::endl;
    }

    writeHdf5Dataset(data, range);
}

void writeHdfTimelineCommand(TokenStream& tokens) {
    // if (mpi_rank == 0) {
    //     std::cout << "writing timeline to HDF5" << std::endl;
    // }

    hsize_t dims[3] = {1, 1, (hsize_t)NOm};
    hsize_t offset[3] = {(hsize_t)local_x_pos, (hsize_t)local_y_pos, 0};

    std::cout << "Rank "<<mpi_rank<<": writing timeline to HDF5" << std::endl;

    // MPI_Barrier(MPI_COMM_WORLD);
    // std::cout << "Rank "<<mpi_rank<<": All Ranks have reached writeHdfTimelineCommand "<< offset[0] << " " << offset[1] << " " << offset[2] << std::endl;

    writeHdf5TimelineDataset(data);

    // MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Rank "<<mpi_rank<<": All Ranks have executed writeHdfTimelineCommand "<< offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
}

void readHdfTimelineCommand(TokenStream& tokens) {
    // if (mpi_rank == 0) {
    //     std::cout << "reading timeline from HDF5" << std::endl;
    // }
    std::cout << "Rank "<<mpi_rank<<": reading timeline from HDF5" << std::endl;
    allocateTimelineData();
    readHdf5TimelineDataset(data);
}

void closeHdfCommand(TokenStream& tokens) {
    // if (mpi_rank == 0) {
    //     std::cout << "closing HDF file" << std::endl;
    // }
    std::cout <<"Rank "<<mpi_rank<<": closing HDF file" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Rank "<<mpi_rank<<": All Ranks have reached closeHdfCommand" << std::endl;

    closeHdf5File();
    std::cout << "Rank "<<mpi_rank<<": All Ranks have executed closeHdfCommand" << std::endl;
}

void closeHdfReadCommand(TokenStream& tokens) {
    // if (mpi_rank == 0) {
    //     std::cout << "closing HDF file for reading" << std::endl;
    // }
    std::cout <<"Rank "<<mpi_rank<<": closing HDF file for reading" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout <<"Rank "<<mpi_rank<<": closing HDF file for reading (after barrier)" << std::endl;
    closeHdf5ReadFile();
    std::cout <<"Rank "<<mpi_rank<<": HDF file for reading has been closed" << std::endl;
}
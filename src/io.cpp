#include "io.hpp"

#include <iostream>
#include <fstream>

#include <mpi.h>
#include <hdf5.h>

std::map<std::string, std::shared_ptr<Array2D>> dataStore;

void wrapErr(herr_t err) {
    if (err < 0) {
        throw std::runtime_error("HDF5 error");
    }
}

hid_t wrapInvalid(hid_t id) {
    if (id == H5I_INVALID_HID) {
        throw std::runtime_error("HDF5 error: invalid ID");
    }
    return id;
}

std::ofstream createOfstream(std::string filename, std::ios_base::openmode mode = std::ios_base::out) {
    std::ofstream out(filename, mode);
    if (!out.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    return out;
}

void writeAscii(fftw_complex *data, XYRange range, std::string filename) 
{
    std::ofstream out = createOfstream(filename);

    for (ptrdiff_t i = 0; i < local_n0; ++i) {
        for (ptrdiff_t j = 0; j < N1; ++j) {
            double *c = data[i*N1 + j];
            Coord pos = range.getCoord(i, j);
            out << pos.x << " " << pos.y << " " << c[0] << " " << c[1] << std::endl;
        }
        out << std::endl;
    }
    out.close();
}

void writeBinary(fftw_complex *data, std::string filename) 
{
    std::ofstream out = createOfstream(filename, std::ios::out | std::ios::binary);

    out.write((char *) data, local_n0*N1*sizeof(fftw_complex));

    out.close();
}

/**
 * Gathers the data about the domain decomposition from all processes.
 * The master process writes this data to an ASCII file.
 */
void writeMetaInfo(std::string filename) {
    std::ofstream out;

    int mpi_size;
    int root = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == root) {
        out = createOfstream(filename);
        out << "# Cavity Data Description" << std::endl;
        out << N0 << " " << N1 << std::endl;
        out << mpi_size << std::endl;
    }
    long *gatherBuffer = new long[mpi_size];
    long local_n0_start_long = local_0_start;
    long local_n0_long = local_n0;

    MPI_Gather(&local_n0_start_long, 1, MPI_LONG, gatherBuffer, 1, MPI_LONG, root, MPI_COMM_WORLD);

    if (mpi_rank == root) {
        for (int i=0; i<mpi_size; ++i) {
            out << gatherBuffer[i] << " ";
        }
        out << std::endl;
    }

    MPI_Gather(&local_n0_long, 1, MPI_LONG, gatherBuffer, 1, MPI_LONG, root, MPI_COMM_WORLD);

    if (mpi_rank == root) {
        for (int i=0; i<mpi_size; ++i) {
            out << gatherBuffer[i] << " ";
        }
        out << std::endl;
    }

    delete[] gatherBuffer;
    if (mpi_rank == root) {
       out.close();
    }
}

/**
 * Read a 2D array from an ascii file.
 * The first line contains the dimensions of the array.
 * Each subsequent line of the file should contain X Y Value.
 * 
 * The data is stored in the dataStore map.
 * 
 * @param filename The name of the file to read from.
 * @param key The key to store the data under in the dataStore map.
 */
void readAscii(std::string filename, std::string key) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    int N0, N1;
    in >> N0 >> N1;

    std::shared_ptr<Array2D> arrayPtr(new Array2D(N0, N1));
    Array2D *array = arrayPtr.get();
    int x, y;
    double value;
    for (int i = 0; i < N0; ++i) {
        for (int j = 0; j < N1; ++j) {
            in >> x >> y >> value;
            array->data[x-1][y-1] = value;
        }
    }
    dataStore[key] = arrayPtr;
}

/**
 * Create an HDF5 file with a 3D dataset.
 * The dataset is created with the dimensions of the data N0 x N1 x NOm.
 *
 * @param range The range of the data.
 * @param filename The name of the file to write to.
 * @param datasetName The name of the dataset to write.
 */
void createHdf5File(std::string fname) {
    const int sieve_buf_size = 262144;
    const int align_threshold = 1; //524288;
    const int alignment = 1; // 262144;

    MPI_Info mpi_info;

    int rank, size;
    MPI_Group worldGroup;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

    /* setup file access template */
    hid_t plist_id = wrapInvalid(H5Pcreate (H5P_FILE_ACCESS));

    wrapErr(H5Pset_sieve_buf_size(plist_id, sieve_buf_size));
    wrapErr(H5Pset_alignment(plist_id, align_threshold, alignment));

    MPI_Info_create(&mpi_info);

    // MPI
    const char access_style[]         = "access_style";
    const char write_once[]           = "write_once";
    const char collective_buffering[] = "collective_buffering";
    const char strue[]                = "true";
    const char cb_block_size[]        = "cb_block_size";
    const char n1048576[]             = "1048576";
    const char cb_buffer_size[]       = "cb_buffer_size";
    const char n4194304[]             = "4194304";

    MPI_Info_set(mpi_info, access_style, write_once);
    MPI_Info_set(mpi_info, collective_buffering, strue);
    MPI_Info_set(mpi_info, cb_block_size, n1048576);
    MPI_Info_set(mpi_info, cb_buffer_size, n4194304);

    /* set Parallel access with communicator */
    wrapErr(H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, mpi_info));

    /* open the file collectively */
    hdf5FileId = wrapInvalid(H5Fcreate (fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id));

    if (hdf5FileId < 0) {
        throw std::runtime_error("Could not open HDF5 file: " + fname);
    }

    /* Release file-access template */
    MPI_Info_free(&mpi_info);
    wrapErr(H5Pclose(plist_id));

    hdf5DxplId = wrapInvalid(H5Pcreate(H5P_DATASET_XFER));
    wrapErr(H5Pset_dxpl_mpio(hdf5DxplId, H5FD_MPIO_COLLECTIVE));

    // Create the dataset
    if (N0 < 0 || N1 < 0 || NOm < 0) {
        throw std::runtime_error("Data dimensions are invalid");
    }

    hsize_t dims[3] = {(hsize_t)N0, (hsize_t)N1, (hsize_t)NOm};
    hdfDataSpaceId = wrapInvalid(H5Screate_simple(3, dims, NULL));
    hdfDatasetRealId = wrapInvalid(H5Dcreate2(hdf5FileId, "real", H5T_NATIVE_DOUBLE, hdfDataSpaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    hdfDatasetImagId = wrapInvalid(H5Dcreate2(hdf5FileId, "imag", H5T_NATIVE_DOUBLE, hdfDataSpaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // Label 1.1
    // std::cout<<"File Created : " <<fname<<endl;
}

/**
 * Open an existing HDF5 file with a 3D dataset.
 */
void openHdf5File(std::string fname) {
    const int sieve_buf_size = 262144;
    const int align_threshold = 1; //524288;
    const int alignment = 1; // 262144;

    MPI_Info mpi_info;

    int rank, size;
    MPI_Group worldGroup;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

    /* setup file access template */
    hid_t plist_id = wrapInvalid(H5Pcreate (H5P_FILE_ACCESS));

    wrapErr(H5Pset_sieve_buf_size(plist_id, sieve_buf_size));
    wrapErr(H5Pset_alignment(plist_id, align_threshold, alignment));

    MPI_Info_create(&mpi_info);

    // MPI
    const char access_style[]         = "access_style";
    const char write_once[]           = "write_once";
    const char collective_buffering[] = "collective_buffering";
    const char strue[]                = "true";
    const char cb_block_size[]        = "cb_block_size";
    const char n1048576[]             = "1048576";
    const char cb_buffer_size[]       = "cb_buffer_size";
    const char n4194304[]             = "4194304";

    MPI_Info_set(mpi_info, access_style, write_once);
    MPI_Info_set(mpi_info, collective_buffering, strue);
    MPI_Info_set(mpi_info, cb_block_size, n1048576);
    MPI_Info_set(mpi_info, cb_buffer_size, n4194304);

    /* set Parallel access with communicator */
    wrapErr(H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, mpi_info));

    /* open the file collectively */
    hdf5ReadFileId = wrapInvalid(H5Fopen (fname.c_str(), H5F_ACC_RDONLY, plist_id));

    if (hdf5ReadFileId < 0) {
        throw std::runtime_error("Could not open HDF5 file: " + fname);
    }

    /* Release file-access template */
    MPI_Info_free(&mpi_info);
    wrapErr(H5Pclose(plist_id));

    hdf5ReadDxplId = wrapInvalid(H5Pcreate(H5P_DATASET_XFER));
    wrapErr(H5Pset_dxpl_mpio(hdf5ReadDxplId, H5FD_MPIO_COLLECTIVE));

    // Create the dataset
    if (N0 < 0 || N1 < 0 || NOm < 0) {
        throw std::runtime_error("Data dimensions are invalid");
    }

    hdfReadDatasetRealId = wrapInvalid(H5Dopen(hdf5ReadFileId, "real", H5P_DEFAULT));
    hdfReadDatasetImagId = wrapInvalid(H5Dopen(hdf5ReadFileId, "imag", H5P_DEFAULT));

    hsize_t dims[3] = {(hsize_t)N0, (hsize_t)N1, (hsize_t)NOm};
    hdfReadDataSpaceId = wrapInvalid(H5Screate_simple(3, dims, NULL));

    // Label 1.1
    // std::cout<<"File Opened : " <<fname<<endl;
}


/**
 * Write the local section of the 2D FFTW data array to the current HDF5 file.
 * The index in the third dimension is given by local_om_pos.
 */
void writeHdf5Dataset(fftw_complex *data, XYRange range) {
    if (hdf5FileId < 0 || hdf5DxplId < 0) {
        throw std::runtime_error("HDF5 file has not been initialised");
    }

    hsize_t dims[3] = {(hsize_t)local_n0, (hsize_t)N1, 1};
    hsize_t offset[3] = {(hsize_t)local_0_start, 0, (hsize_t)local_om_pos};
    // std::cout << "Writing dataset at " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
    // std::cout << "Dims: " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    // std::cout << "Selecting hyperslab" << std::endl;
    wrapErr(H5Sselect_hyperslab(hdfDataSpaceId, H5S_SELECT_SET, offset, NULL, dims, NULL));

    // std::cout << "Creating memspace" << std::endl;
    hsize_t mem_dims[2] = {(hsize_t)local_n0, (hsize_t)(2*N1)};
    hsize_t mem_offset[2] = {0, 0};
    hsize_t mem_stride[2] = {1, 2};
    hid_t memspace_id = wrapInvalid(H5Screate_simple(2, mem_dims, NULL));
    wrapErr(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offset, mem_stride, dims, NULL));

    // std::cout << "Writing real data" << std::endl;
    wrapErr(H5Dwrite(hdfDatasetRealId, H5T_NATIVE_DOUBLE, memspace_id, hdfDataSpaceId, hdf5DxplId, data));

    mem_offset[1] = 1;
    wrapErr(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offset, mem_stride, dims, NULL));
    // std::cout << "Writing imag data" << std::endl;
    wrapErr(H5Dwrite(hdfDatasetImagId, H5T_NATIVE_DOUBLE, memspace_id, hdfDataSpaceId, hdf5DxplId, data));

    // std::cout << "Data written" << std::endl;
    wrapErr(H5Sclose(memspace_id));
    // std::cout << "Memspace closed" << std::endl;
}

/**
 * Write the timeline data to the HDF5 file.
 * The indices in the x and y dimensions are given by local_x_pos and local_y_pos.
 */
void writeHdf5TimelineDataset(fftw_complex *data) {
    if (hdf5FileId < 0 || hdf5DxplId < 0) {
        throw std::runtime_error("HDF5 file has not been initialised");
    }

    hsize_t dims[3] = {1, 1, (hsize_t)NOm};
    hsize_t offset[3] = {(hsize_t)local_x_pos, (hsize_t)local_y_pos, 0};
    // std::cout << "Writing dataset at " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
    // std::cout << "Dims: " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    // std::cout << "Selecting hyperslab" << std::endl;
    wrapErr(H5Sselect_hyperslab(hdfDataSpaceId, H5S_SELECT_SET, offset, NULL, dims, NULL));

    // std::cout << "Creating memspace" << std::endl;
    hsize_t mem_dims[1] = {(hsize_t)(2*NOm)};
    hsize_t mem_select_dims[1] = {(hsize_t)NOm};
    hsize_t mem_offset[1] = {0};
    hsize_t mem_stride[1] = {2};
    hid_t memspace_id = wrapInvalid(H5Screate_simple(1, mem_dims, NULL));
    wrapErr(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offset, mem_stride, mem_select_dims, NULL));

    // std::cout << "Writing real data" << std::endl;
    wrapErr(H5Dwrite(hdfDatasetRealId, H5T_NATIVE_DOUBLE, memspace_id, hdfDataSpaceId, hdf5DxplId, data));

    mem_offset[0] = 1;
    wrapErr(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offset, mem_stride, mem_select_dims, NULL));
    // std::cout << "Writing imag data" << std::endl;
    wrapErr(H5Dwrite(hdfDatasetImagId, H5T_NATIVE_DOUBLE, memspace_id, hdfDataSpaceId, hdf5DxplId, data));

    // std::cout << "Data written" << std::endl;
    wrapErr(H5Sclose(memspace_id));
    // std::cout << "Memspace closed" << std::endl;
}


/**
 * Read the timeline data from the HDF5 file.
 * The indices in the x and y dimensions are given by local_x_pos and local_y_pos.
 */
void readHdf5TimelineDataset(fftw_complex *data) {
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (hdf5ReadFileId < 0 || hdf5ReadDxplId < 0) {
        throw std::runtime_error("HDF5 file has not been initialised");
    }

    hsize_t dims[3] = {1, 1, (hsize_t)NOm};
    hsize_t offset[3] = {(hsize_t)local_x_pos, (hsize_t)local_y_pos, 0};
    std::cout <<"Rank "<<rank<<": Reading dataset at " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
    std::cout <<"Rank "<<rank<<": Dims: " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    std::cout <<"Rank "<<rank<<": Selecting hyperslab" << std::endl;
    wrapErr(H5Sselect_hyperslab(hdfReadDataSpaceId, H5S_SELECT_SET, offset, NULL, dims, NULL));

    // std::cout << "Creating memspace" << std::endl;
    hsize_t mem_dims[1] = {(hsize_t)(2*NOm)};
    hsize_t mem_select_dims[1] = {(hsize_t)NOm};
    hsize_t mem_offset[1] = {0};
    hsize_t mem_stride[1] = {2};
    std::cout <<"Rank "<<rank<<": Creating memspace" << std::endl;
    std::cout <<"Rank "<<rank<<": Mem dims: " << mem_dims[0] << std::endl;
    std::cout <<"Rank "<<rank<<": Mem offset: " << mem_offset[0] << std::endl;
    std::cout <<"Rank "<<rank<<": Mem stride: " << mem_stride[0] << std::endl;

    hid_t memspace_id = wrapInvalid(H5Screate_simple(1, mem_dims, NULL));
    wrapErr(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offset, mem_stride, mem_select_dims, NULL));

    std::cout <<"Rank "<<rank<<": Reading real data" << std::endl;
    wrapErr(H5Dread(hdfReadDatasetRealId, H5T_NATIVE_DOUBLE, memspace_id, hdfReadDataSpaceId, hdf5ReadDxplId, data));

    mem_offset[0] = 1;
    wrapErr(H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, mem_offset, mem_stride, mem_select_dims, NULL));
    // std::cout << "Reading imag data" << std::endl;
    wrapErr(H5Dread(hdfReadDatasetImagId, H5T_NATIVE_DOUBLE, memspace_id, hdfReadDataSpaceId, hdf5ReadDxplId, data));

    // std::cout << "Data written" << std::endl;
    wrapErr(H5Sclose(memspace_id));
    // std::cout << "Memspace closed" << std::endl;
}


void closeHdf5File() {
    wrapErr(H5Sclose(hdfDataSpaceId));
    wrapErr(H5Dclose(hdfDatasetRealId));
    wrapErr(H5Dclose(hdfDatasetImagId));
    wrapErr(H5Fclose(hdf5FileId));
    wrapErr(H5Pclose(hdf5DxplId));
}

void closeHdf5ReadFile() {
    wrapErr(H5Sclose(hdfReadDataSpaceId));
    wrapErr(H5Dclose(hdfReadDatasetRealId));
    wrapErr(H5Dclose(hdfReadDatasetImagId));
    wrapErr(H5Fclose(hdf5ReadFileId));
    wrapErr(H5Pclose(hdf5ReadDxplId));
}
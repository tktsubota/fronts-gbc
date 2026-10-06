#pragma once
#include <H5Cpp.h>
#include <armadillo>
#include <string>
#include <stdexcept>

class H5DatasetBatchLoader {
private:
    std::string dataset_name;
    size_t batch_size;
    size_t nrows, ncols;

    H5::DataSet dataset;
    arma::mat cache;
    size_t batch_start, batch_cols;

    void load_batch(size_t k) {
        batch_start = k;
        batch_cols = std::min(batch_size, ncols - k);
        cache.set_size(nrows, batch_cols);

        // Define hyperslab in the file
        std::vector<hsize_t> offset = {batch_start, 0};
        std::vector<hsize_t> count  = {batch_cols, nrows};

        H5::DataSpace filespace = dataset.getSpace();
        filespace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data());

        // Define memory space for the buffer
        H5::DataSpace memspace(2, count.data());

        // Flattened buffer
        std::vector<double> buffer(nrows * batch_cols);

        // Read hyperslab into buffer
        dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);

        // Wrap into Armadillo (copy)
        cache = arma::mat(buffer.data(), nrows, batch_cols, /*copy_aux_mem=*/true);
    }

public:
    H5DatasetBatchLoader(H5::H5File& file,
                         const std::string& dataset_name,
                         size_t batch_size)
        : dataset_name(dataset_name),
          batch_size(batch_size),
          batch_start(0),
          batch_cols(0) {
        dataset = file.openDataSet(dataset_name);
        H5::DataSpace dataspace = dataset.getSpace();

        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims);

        nrows = dims[0];
        ncols = dims[1];

        cache.set_size(nrows, 0);
    }

    arma::vec get_column(size_t k) {
        if (k >= ncols) throw std::out_of_range("Column index out of range");
        if (k < batch_start || k >= batch_start + batch_cols) {
            load_batch(k);
        }
        return cache.col(k - batch_start);
    }

};

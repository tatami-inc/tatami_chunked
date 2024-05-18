#ifndef TATAMI_TEST_MOCK_BLOB_HPP
#define TATAMI_TEST_MOCK_BLOB_HPP

#include <vector>
#include <algorithm>

template<bool row_major_>
struct MockDenseBlob {
    typedef double value_type;

    bool is_row_major() const {
        return row_major_;
    }

private:
    int nrows, ncols;
    std::vector<double> contents;

public:
    MockDenseBlob() = default;

    MockDenseBlob(int nr, int nc, std::vector<double> c) : nrows(nr), ncols(nc), contents(std::move(c)) {}

    int nrow() const {
        return nrows;
    }

    int ncol() const {
        return ncols;
    }

    void inflate(std::vector<double>& buffer) const {
        buffer.resize(contents.size());
        std::copy(contents.begin(), contents.end(), buffer.begin());
    }
};

template<bool csr_>
struct MockSparseBlob {
    typedef int index_type;
    typedef double value_type;

    bool is_csr() const {
        return csr_;
    }

private:
    int nrows, ncols;
    std::vector<double> vcontents;
    std::vector<int> icontents;
    std::vector<size_t> pcontents;

public:
    MockSparseBlob() = default;

    MockSparseBlob(int nr, int nc, std::vector<double> v, std::vector<int> i, std::vector<size_t> p) : 
        nrows(nr), ncols(nc), vcontents(std::move(v)), icontents(std::move(i)), pcontents(std::move(p)) {}

    int nrow() const {
        return nrows;
    }

    int ncol() const {
        return ncols;
    }

    void inflate(std::vector<double>& vbuffer, std::vector<int>& ibuffer, std::vector<size_t>& pbuffer) const {
        vbuffer.resize(vcontents.size());
        std::copy(vcontents.begin(), vcontents.end(), vbuffer.begin());
        ibuffer.resize(icontents.size());
        std::copy(icontents.begin(), icontents.end(), ibuffer.begin());
        pbuffer.resize(pcontents.size());
        std::copy(pcontents.begin(), pcontents.end(), pbuffer.begin());
    }
};

#endif

#pragma once

#include "../Define.h"
#include "../SparseStorage.h"

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cusolverSp.h>

namespace EQlib {

class LinearSolver
{
public:
    virtual bool analyze(const SparseStorage<double, int, true>& structure) = 0;

    virtual bool factorize(Ref<const Vector> a) = 0;

    virtual bool solve(Ref<const Vector> b, Ref<Vector> x) = 0;
};

class CudaCholeskyLinearSolver : public LinearSolver
{
private:    // variables
    cusolverSpHandle_t m_handle;
    SparseStorage<double, int, true> m_structure;
    double* m_device_a;
    int* m_device_ia;
    int* m_device_ja;
    double* m_device_x;
    double* m_device_b;
    cusparseMatDescr_t m_descr_a;

public:     // constructors

public:     // methods
    ~CudaCholeskyLinearSolver()
    {
        free();
    }

    void free()
    {
        if (m_device_a) cudaFree(m_device_a);
        if (m_device_ia) cudaFree(m_device_ia);
        if (m_device_ja) cudaFree(m_device_ja);
        if (m_device_x) cudaFree(m_device_x);
        if (m_device_b) cudaFree(m_device_b);
    }

    bool analyze(const SparseStorage<double, int, true>& structure) override
    {
        free();

        m_structure = structure;

        cusolverSpCreate(&m_handle);

        cudaMalloc((void**)&m_device_a, sizeof(double) * m_structure.nnz());
        cudaMalloc((void**)&m_device_ia, sizeof(int) * m_structure.ia().size());
        cudaMalloc((void**)&m_device_ja, sizeof(int) * m_structure.ja().size());
        cudaMalloc((void**)&m_device_x, sizeof(double) * m_structure.cols());
        cudaMalloc((void**)&m_device_b, sizeof(double) * m_structure.rows());

        cudaMemcpy(m_device_ia, m_structure.ia().data(), sizeof(int) * m_structure.ia().size(), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_ja, m_structure.ja().data(), sizeof(int) * m_structure.ja().size(), cudaMemcpyHostToDevice);
        
        cusparseCreateMatDescr(&m_descr_a);
        cusparseSetMatType(m_descr_a, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatFillMode(m_descr_a, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatIndexBase(m_descr_a, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatDiagType(m_descr_a, CUSPARSE_DIAG_TYPE_NON_UNIT);
    }

    bool factorize(Ref<const Vector> a) override
    {
        cudaMemcpy(m_device_a, a.data(), sizeof(double) * m_structure.nnz(), cudaMemcpyHostToDevice);
    }

    bool solve(Ref<const Vector> b, Ref<Vector> x) override
    {
        cudaMemcpy(m_device_x, x.data(), sizeof(double) * m_structure.cols(), cudaMemcpyHostToDevice);
        cudaMemcpy(m_device_b, b.data(), sizeof(double) * m_structure.rows(), cudaMemcpyHostToDevice);

        double tol = 1.e-6;
        int reorder = 1;
        int singularity = 0;
        
        cusolverStatus_t error = cusolverSpDcsrlsvchol(
            m_handle,
            m_structure.rows(),
            m_structure.nnz(),
            m_descr_a,
            m_device_a,
            m_device_ia,
            m_device_ja,
            m_device_b,
            tol,
            reorder,
            m_device_x,
            &singularity
        );

        cudaMemcpy(x.data(), m_device_x, sizeof(double) * h_x.size(), cudaMemcpyDeviceToHost);
    }
};

} // namespace EQlib
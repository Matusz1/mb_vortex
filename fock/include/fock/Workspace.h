#ifndef _WORKSPACE_H_
#define _WORKSPACE_H_

#define NDEBUG
#include <Eigen/Sparse>
#include "FockCollection.h"

namespace Mbs {

class Workspace {
private:
    FockCollection m_basis;

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> m_energy;
    Eigen::SparseMatrix<double> m_potential;
    Eigen::SparseMatrix<double> m_hamiltonian;

    Eigen::MatrixXd m_eigenvectors;
    Eigen::VectorXd m_eigenvalues;
public:
    explicit Workspace(FockCollection basis);
    
    void solve(double g);
    [[nodiscard]] const auto& eigenvalues() const { return m_eigenvalues; };
    [[nodiscard]] const auto& eigenvectors() const { return m_eigenvectors; };
    FockLinear<double> eigenvector_combination(uint i) const;
};

} // namespace Mbs

#endif // _WORKSPACE_H_

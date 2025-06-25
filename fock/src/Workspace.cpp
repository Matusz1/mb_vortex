#include "fock/FockState.h"
#include "fock/Workspace.h"
#include "fock/FockCollection.h"

#include <Eigen/src/SparseCore/SparseMatrix.h>
/*#include <Spectra/SymEigsShiftSolver.h>*/
#include <Spectra/SymEigsSolver.h>
/*#include <Spectra/MatOp/SparseSymShiftSolve.h>*/
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <chrono>
#include <iostream>

namespace Mbs {

Workspace::Workspace(FockCollection basis) :
    m_basis(std::move(basis)),
    m_energy(Eigen::DiagonalMatrix<double, Eigen::Dynamic>(m_basis.size())),
    m_potential(Eigen::SparseMatrix<double>(m_basis.size(), m_basis.size())),
    m_hamiltonian(Eigen::SparseMatrix<double>(m_basis.size(), m_basis.size()))
{

    // Time
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Eigen::Triplet<double>> triplets;
    // Temporary vector for speedup access
    std::vector<Fock> temp;
    for (const auto& state : m_basis.states())
        temp.push_back(state);

    for (uint i = 0; i < temp.size(); ++i)
        m_energy.diagonal()[i] = temp[i].energy();

    for (uint i = 0; i < temp.size(); ++i) {
        for (uint j = i; j < temp.size(); ++j) {
            const double val = potential_element(temp[i], temp[j]);
            if (val != 0.0)
                triplets.emplace_back(j, i, val);
        }
    }
    m_potential.setFromTriplets(triplets.begin(), triplets.end());
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    std::cout << "  Fill factor:                " << m_potential.nonZeros() / double(m_basis.size() * m_basis.size()) << std::endl;
    std::cout << "  Matrix filling time:        " << elapsed.count() << " s" << std::endl;
}

void Workspace::solve(double g)
{
    m_hamiltonian = g * m_potential;
    for (uint i = 0; i < m_basis.size(); ++i)
        m_hamiltonian.coeffRef(i, i) += m_energy.diagonal()[i];
    m_hamiltonian.makeCompressed();

    /*auto start = std::chrono::high_resolution_clock::now();*/
    /*Spectra::SparseSymShiftSolve<double> op(m_hamiltonian);*/
    /*Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> solver(op, 1, 4, 0.0);*/
    /*solver.init();*/
    /*auto stop = std::chrono::high_resolution_clock::now();*/
    /*std::cout << "  Solver initialization time: " << std::chrono::duration<double>(stop - start).count() << " s" << std::endl;*/
    /**/
    /*start = std::chrono::high_resolution_clock::now();*/
    /*solver.compute(Spectra::SortRule::LargestMagn, 1000, 1e-10, Spectra::SortRule::SmallestAlge);*/
    /*stop = std::chrono::high_resolution_clock::now();*/
    /*std::chrono::duration<double> elapsed = stop - start;*/
    /*std::cout << "  Solving time:               " << elapsed.count() << " s" << std::endl;*/

    auto start = std::chrono::high_resolution_clock::now();

    Spectra::SparseSymMatProd<double> op(m_hamiltonian);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> solver(op, 1, 4);
    solver.init();
    start = std::chrono::high_resolution_clock::now();
    solver.compute(Spectra::SortRule::SmallestMagn, 1000, 1e-10, Spectra::SortRule::SmallestAlge);

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    std::cout << "  Solving time:               " << elapsed.count() << " s" << std::endl;

    if (solver.info() != Spectra::CompInfo::Successful)
        throw std::runtime_error("Failed to solve the eigenvalue problem");
    m_eigenvalues = solver.eigenvalues();
    m_eigenvectors = solver.eigenvectors();
}

FockLinear<double> Workspace::eigenvector_combination(uint i) const
{
    FockLinear<double> linear;

    const auto& vec = m_eigenvectors.col(i);
    auto it = m_basis.states().begin();
    for (uint j = 0; j < m_basis.size(); ++j, ++it)
        linear.m_states[*it] = vec[j];
    return linear;
}

} // namespace Mbs

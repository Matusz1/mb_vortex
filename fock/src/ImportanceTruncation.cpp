#include "fock/AnnichilationAdvancer.h"
#include "fock/ImportanceTruncation.h" 
#include "fock/DiskState.h"
#include "fock/FockCollection.h"
#include "fock/FockState.h"
#include "fock/Integrator.h"
#include "fock/Workspace.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <unordered_map>


namespace Mbs {

using FockMap = std::unordered_map<Fock, double, Fock::Hash>;

static inline void _construct_reachable(const u16* snpart, const u16* rest, u16* dst, uint npart)
{
    std::merge(snpart, snpart+npart, rest, rest+2, dst);
}

// Some function that will help find all the overlap elements
static FockMap _find_overlap_elements(
    const Fock& fock,
    double cutoff
) {
    const double ener = fock.energy();
    Fock sf(fock);
    Fock rf(fock);
    sf.m_size = fock.size() - 2;

    // <ij|V|kl>
    std::array<u16,4> ijkl = {0, 0, 0, 0};

    FockMap reachables;

    // To find possible matrix elemnts, check what can be annichillated
    for (auto s0 = AnnichilationAdvancer<false>(sf.m_state.data(), fock.size()); s0.valid(); s0.advance()) {
        const double e0 = ener - energy(s0.index());
        const int m0 = angular_momentum(s0.index());
        const uint v1 = s0.count();
        ijkl[3] = s0.index();

        // Annichilating the second particle
        for (auto s1 = AnnichilationAdvancer<true>(sf.m_state.data(), fock.size()-1); s1.valid(); s1.advance()) {
            const double e1 = e0 - energy(s1.index());
            const int m1 = m0 + angular_momentum(s1.index());
            const uint v2 = s1.count();
            ijkl[2] = s1.index();

            // All the possibilities states
            ijkl[0] = 0;
            for (; e1 + energy(ijkl[0]) < cutoff; ++ijkl[0]) {
                const double e2 = e1 + energy(ijkl[0]);
                const int m2 = m1 - angular_momentum(ijkl[0]);

                for (ijkl[1] = ijkl[0]; e2 + energy(ijkl[1]) < cutoff; ++ijkl[1]) {
                    const int m3 = m2 - angular_momentum(ijkl[1]);
                    if (m3 != 0) // Angular momentum conservation
                        continue;
                    _construct_reachable(sf.m_state.data(), ijkl.data(), rf.m_state.data(), sf.size());
                    if (ijkl[0] == ijkl[1]) {
                        const uint v3 = sf.count(ijkl[0]) + 1;
                        double current = potential_integrate_cached_unchecked(ijkl);
                        reachables[rf] += 0.5 * std::sqrt(v1*v2*v3*(v3+1)) * current;
                    } else {
                        const uint v3 = sf.count(ijkl[0]) + 1;
                        const uint v4 = sf.count(ijkl[1]) + 1;
                        double current = potential_integrate_cached_unchecked(ijkl);
                        std::swap(ijkl[0], ijkl[1]);
                        current += potential_integrate_cached_unchecked(ijkl);
                        std::swap(ijkl[0], ijkl[1]);
                        reachables[rf] += 0.5 * std::sqrt(v1*v2*v3*v4) * current;
                    }
                }
            }
        }
    }
    return reachables;
}

Workspace importance_truncation_scheme(double g, uint N, int mom, TruncationParameters params)
{
    std::cout << "### Starting importance truncation scheme ###" << std::endl;
    std::cout << "=============================================" << std::endl;

    // Construct beggining basis
    // It is desirable to start from slightly larger basis,
    // so that the first iteration is not too expensive
    // and gives us a decent approximation for the ground state.
    FockCollection basis = FockCollection::from_size(N, 10, mom);
    const auto& basis_set = basis.states();
    const double ener_low = basis.states().begin()->energy();
    double ret = 0.0;
    double prev = 0.0;

    for (uint iteration = 0;; ++iteration) {
        std::cout << "\n### Iteration " << iteration << " ###" << std::endl;

        // Find a solution for the current basis
        std::cout << "\nComputational basis size: " << basis.size() << std::endl;
        std::cout << "Computing the solution..." << std::endl;
        Workspace workspace(basis);
        workspace.solve(g);
        std::cout << "Eigenvalue: " << std::setprecision(10) << workspace.eigenvalues()[0] << std::endl;
        Eigen::VectorXd v0 = workspace.eigenvectors().col(0);
        ret = workspace.eigenvalues()[0];
        if (std::abs((ret - prev)/ret) < 1e-5)
            return workspace;
        prev = ret;

        // Keep only the coefficients above the threshold
        for (auto& val : v0)
            if (std::abs(val) < params.C_min)
                val = 0.0;
        v0.normalize();

        // Determine the reachable elements
        std::cout << "Finding the reachable elements..." << std::endl;
        FockMap reachables;
        uint i = 0;
        for (auto it = basis_set.begin(); it != basis_set.end(); ++i, ++it) {
            if (v0[i] == 0.0)
                continue;
            auto next = _find_overlap_elements(*it, params.ener_cutoff);
            const double val = v0[i];
            for (auto jt = next.begin(); jt != next.end(); ++jt)
                reachables[jt->first] += val * jt->second;
        }
        std::cout << "Found " << reachables.size() << " reachable elements." << std::endl;

        // Delete elements present in basis
        std::cout << "Constructing the new basis..." << std::endl;
        for (auto it = basis_set.begin(); it != basis.states().end(); ++it) {
            auto found = reachables.find(*it);
            if (found != reachables.end())
                reachables.erase(found);
        }

        // Determine new basis
        for (auto it = reachables.begin(); it != reachables.end(); ++it)
            if (std::abs(g * it->second / (ener_low - it->first.energy())) > params.kappa_min)
                basis.insert(it->first);
    }
}

}

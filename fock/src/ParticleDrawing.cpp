#include <numeric>
#include "fock/ParticleDrawing.h"
#include "fock/Numeric.h"

namespace Mbs {

static constexpr uint _n_r_bins = 40;
static constexpr uint _n_phi_bins = 40;
static constexpr double dr = 1.0 / _n_r_bins;
static constexpr double dphi = 2.0 * M_PI / _n_phi_bins;

static RPhi _draw_pos(const std::vector<double> r_dist, const std::vector<double> phi_dist)
{
    // Draw r
    const double r_rand = r_dist.back() * (double)rand() / RAND_MAX;
    auto it = std::lower_bound(r_dist.begin(), r_dist.end(), r_rand);
    uint ir = it - r_dist.begin();
    double r = ir * dr + dr*(double)rand() / RAND_MAX;
    assert(ir < _n_r_bins);

    // Draw phi
    const double phi_rand = phi_dist[(ir+1)*_n_phi_bins - 1] * (double)rand() / RAND_MAX;
    it = std::lower_bound(phi_dist.begin() + ir*_n_phi_bins, phi_dist.begin() + (ir+1)*_n_phi_bins, phi_rand);
    uint iphi = it - (phi_dist.begin() + ir*_n_phi_bins);
    assert(iphi < _n_phi_bins);
    double phi = iphi * dphi + dphi*(double)rand() / RAND_MAX;

    return {r, phi};
}

static std::pair<std::vector<double>, std::vector<double>> _comp_dens_dist(const FockLinear<Complex>& state)
{
    auto rho = state.density_matrix();
    std::vector<double> phi_dist(_n_r_bins*_n_phi_bins);
    std::vector<double> r_dist(_n_r_bins);
    for (uint i = 0; i < _n_r_bins; ++i) {
        const double r = (i + 0.5) * dr;
        for (uint j = 0; j < _n_phi_bins; ++j) {
            const double phi = (j + 0.5) * dphi;
            phi_dist[i*_n_phi_bins + j] = r * density_matrix_compute(r, phi, rho);
            assert(phi_dist[i*_n_phi_bins + j] >= 0.0);
        }
        std::partial_sum(phi_dist.begin() + i*_n_phi_bins, phi_dist.begin() + (i+1)*_n_phi_bins, phi_dist.begin() + i*_n_phi_bins);
        r_dist[i] = phi_dist[(i+1)*_n_phi_bins - 1];
    }
    std::partial_sum(r_dist.begin(), r_dist.end(), r_dist.begin());

    return {r_dist, phi_dist};
}

static void _rec_draw_pos(const FockLinear<Complex>& state, std::vector<RPhi>& ret, double c_cutoff)
{
    auto [r_dist, phi_dist] = _comp_dens_dist(state);
    RPhi rphi = _draw_pos(r_dist, phi_dist);
    ret.push_back(rphi);
    if (state.npart() > 2) {
        auto next = apply_psi(rphi.r, rphi.phi, state);
        next.normalize();
        next.reduce(c_cutoff);
        next.normalize();
        _rec_draw_pos(next, ret, c_cutoff);
    }
}

std::vector<RPhi> draw_positions(const FockLinear<double>& state, double c_cutoff)
{
    FockLinear<Complex> psi = state.complexify();
    psi.reduce(c_cutoff);
    psi.normalize();
    std::vector<RPhi> ret;
    _rec_draw_pos(psi, ret, c_cutoff);
    return ret;
}

std::vector<Complex> fill_wave_function(const std::vector<RPhi>& positions, const FockLinear<double>& fs, uint npts, double cutoff)
{
    FockLinear<double> state = fs;
    state.reduce(cutoff);
    state.normalize();

    const uint npart = state.npart();
    const uint state_len = state.size();

    Eigen::MatrixXcd sub_matrix(npart-1, npart);
    Eigen::MatrixXcd permanents(npart, state_len);

    auto it = state.begin();
    for (uint i = 0; i < state_len; ++i, ++it) {
        // Fill the submatrix
        for (uint j = 0; j < npart-1; ++j)
            for (uint k = 0; k < npart; ++k)
                sub_matrix(j, k) = wave_function(it->first[k], positions[j].r, positions[j].phi);

        // Compute the norm for the permanent given sum
        const double snorm = it->first.permanent_norm();
        for (uint j = 0; j < npart; ++j) {
            const Complex perm = permanent(&sub_matrix(0, 1), npart-1);
            sub_matrix.col(0).swap(sub_matrix.col(j));
            permanents(j, i) = perm * snorm;
        }
    }

    // Fill the wave function
    std::vector<Complex> ret(npts*npts);
    for (uint i = 0; i < npts; ++i) {
        const double y = -1.0 + 2.0 * i / (npts-1);
        for (uint j = 0; j < npts; ++j) {
            const double x = -1.0 + 2.0 * j / (npts-1);
            const RPhi rphi = cartesian_to_polar(x, y);

            // Compute at point (x, y)
            Complex psi = 0.0;
            it = state.begin();
            for (uint k = 0; k < state_len; ++k, ++it) {
                Complex psi_k = 0.0;
                for (uint n = 0; n < npart; ++n) {
                    const Complex wf = wave_function(it->first[n], rphi.r, rphi.phi);
                    psi_k += permanents(n, k) * wf;
                }
                psi += it->second * psi_k;
            }
            ret[i*npts + j] = psi;
        }
    }

    return ret;
}

}

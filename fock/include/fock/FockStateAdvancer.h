#ifndef _FOCK_STATE_ADVANCER_H_
#define _FOCK_STATE_ADVANCER_H_

#include "DiskState.h"
#include "FockState.h"
#include <array>
#include <vector>

namespace Mbs {

class FockStateAdvancer {
    /*std::vector<uint> m_state;*/
    /*uint npart;*/
    /**/
    /*std::vector<double> ener_accum;*/
    /*std::vector<double> ener_sp;*/
    /*double ener_cutoff;*/
    /**/
    /*std::vector<int>    amom_accum;*/
    /*std::vector<int>    amom_sp;*/
    /*int tot_amom;*/
    /**/
    /*uint i;*/
    Fock m_state;

    std::array<double, Mbs::CONFING_MAX_FOCK_SIZE+1> ener_accum;
    std::vector<double> ener_sp{};
    double ener_cutoff;

    std::array<int, Mbs::CONFING_MAX_FOCK_SIZE+1> amom_accum;
    std::vector<int> amom_sp{};
    int tot_amom;

    uint i;
public:
    FockStateAdvancer(uint npart, double ener_cutoff, int tot_amom);
    void advance() noexcept;
    bool valid() const noexcept { return i < m_state.m_size; }
    uint size() const noexcept { return m_state.m_size; }
    const Fock& state() const noexcept { return m_state; }
};

inline FockStateAdvancer::FockStateAdvancer(uint npart, double ener_cutoff, int tot_amom) :
    m_state(npart),
    ener_cutoff(ener_cutoff),
    tot_amom(tot_amom),
    i(npart - 1)
{
    ener_accum.fill(0.0);
    amom_accum.fill(0);

    for (u16 i = 0; energy(i) < ener_cutoff; ++i) {
        ener_sp.push_back(energy(i));
        amom_sp.push_back(angular_momentum(i));
    }

    m_state[npart] = ener_sp.size();
    if (tot_amom != 0 || npart * ener_sp[0] >= ener_cutoff)
        advance();
}

inline void FockStateAdvancer::advance() noexcept
{
    while (i < m_state.size()) {
        if (m_state[i] == m_state[i+1]) {
            m_state[i] = 0;
            ++i;
            continue;
        }
        ++m_state[i];

        // This ensures that the energy is below the cutoff after this point
        ener_accum[i] = ener_sp[m_state[i]] + ener_accum[i+1];
        if (ener_accum[i] + i*ener_sp[0] >= ener_cutoff) {
            m_state[i] = 0;
            ++i;
            continue;
        }

        const int amom_curr = amom_sp[m_state[i]] + amom_accum[i+1];
        if (i > 0) {
            amom_accum[i] = amom_curr;
            --i;
        }

        if (amom_curr == tot_amom) // At this point, we know energy is below cutoff
            break;
    };
}

}

#endif // _FOCK_STATE_ADVANCER_H_

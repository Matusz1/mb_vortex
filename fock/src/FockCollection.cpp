#include "fock/FockState.h"
#include "fock/FockCollection.h"
#include "fock/FockStateAdvancer.h"

#include <iostream>
#include <iomanip>

namespace Mbs {

FockCollection::FockCollection(uint npart, double ener_cutoff, int ang_mom)
{
    for (FockStateAdvancer it(npart, ener_cutoff, ang_mom); it.valid(); it.advance())
        m_states.insert(it.state());
}

uint FockCollection::calculate_size(uint npart, double ener_cutoff, int ang_mom)
{
    uint count = 0;
    for (FockStateAdvancer it(npart, ener_cutoff, ang_mom); it.valid(); it.advance())
        ++count;
    return count;
}

[[nodiscard]] FockCollection FockCollection::from_size(uint npart, uint size, int ang_mom)
{
    constexpr double delta_cutoff = 10.0;
    double init_cutoff = 0.0;
    double curr_size = 0;
    while (curr_size < size) {
        init_cutoff += delta_cutoff;
        curr_size = FockCollection::calculate_size(npart, init_cutoff, ang_mom);
    }
    return FockCollection(npart, init_cutoff-delta_cutoff, ang_mom);
}

// TODO: Use an operator<< for Fock class.
std::ostream& operator<<(std::ostream& os, const FockCollection& basis)
{
    const uint size = basis.size();
    const uint npart = basis.npart();
    const uint format_size = std::log10(size) + 1;
    auto it = basis.m_states.begin();
    for (uint i = 0; i < size; ++i, ++it) {
        os << std::setw(format_size) << i << ": {";
        for (uint j = 0; j < npart; ++j)
            os << (*it)[j] << (j == npart - 1 ? "" : ", ");
        os << "}" << '\n';
    }
    return os;
}

} // namespace Mbs

#ifndef _FOCK_BASIS_H_
#define _FOCK_BASIS_H_

#include "Core.h"
#include "FockState.h"
#include <set>
#include <Eigen/Dense>

namespace Mbs {

class FockCollection {
private:
    std::set<Fock> m_states;

public:
    FockCollection() = default;
    FockCollection(uint npart, double ener_cutoff, int ang_mom = 0);

    [[nodiscard]] uint size() const { return m_states.size(); }
    [[nodiscard]] uint npart() const { return m_states.empty() ? 0 : m_states.begin()->size(); }
    
    const std::set<Fock>& states() const { return m_states; }

    /* Inserts into the basis.
     * If the state is already in the basis, it is not added. */
    bool contains(const Fock& fock) const { return m_states.find(fock) != m_states.end(); }

    void insert(const Fock& fock) { m_states.insert(fock); }

    void merge(std::set<Fock>& states) { m_states.merge(states); }
    void merge(std::set<Fock>&& states) { m_states.merge(std::move(states)); }

    // TODO: Function name is misleading. Something better?
    // Maybe also free the function.
    [[nodiscard]] static uint calculate_size(uint npart, double ener_cutoff, int ang_mom = 0);
    [[nodiscard]] static FockCollection from_size(uint npart, uint size, int ang_mom = 0);

    friend std::ostream& operator<<(std::ostream& os, const FockCollection& basis);
};

std::ostream& operator<<(std::ostream& os, const FockCollection& basis);

} // namespace Mbs

#endif // _FOCK_BASIS_H_

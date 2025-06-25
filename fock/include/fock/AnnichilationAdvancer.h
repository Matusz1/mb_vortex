#ifndef _ANNIHILATION_ADVANCER_H_
#define _ANNIHILATION_ADVANCER_H_

#include <algorithm>
#include "Core.h"

namespace Mbs {

/*
 * Initialized with some pointer to MB state with npart particles,
 * this class allows to iterate over all possible states with one less particle.
 * That is, annihilates one particle and after using the advance() method, another.
 *
 * Because _AnnichilationAdvancer modifies underlying state, we introduce a Reset
 * template parameter to allow for the state to be reset to its original state,
 * after all the possible particles have been annihilated.
 *
 * TODO: Use standart iterator interface.
 */
template<bool Reset = true>
class AnnichilationAdvancer {
    u16* m_state;
    uint m_size;
    int i, j;
public:
    AnnichilationAdvancer(u16* state, uint size);
      
    u16 index() const noexcept { return m_state[m_size]; }
    uint count() const noexcept { return j-i; }
    bool valid() const noexcept { return count() != 0; }
    void advance() noexcept;
};

template<bool Reset>
inline AnnichilationAdvancer<Reset>::AnnichilationAdvancer(u16* state, uint size) :
    m_state(state),
    m_size(size-1),
    i(size-1),
    j(size-1)
{
    if (size == 0)
        return;
    while (--i >= 0 && m_state[i] == m_state[m_size])
        ; // Do nothing
}

template<bool Reset>
inline void AnnichilationAdvancer<Reset>::advance() noexcept
{
    if (i < 0) { // i < 0 impies last possible particle was alredy annihilated
        if (j != i) { // But j != i means no advance was made after last annihilation
            if (Reset)
                std::rotate(m_state, m_state + m_size, m_state + m_size + 1);
            j = i;
        }
        return;
    }
    std::swap(m_state[i], m_state[m_size]);
    j = i;
    while (--i >= 0 && m_state[i] == m_state[m_size])
        ; // Do nothing
}

} // namespace Mbs


#endif // _ANNIHILATION_ADVANCER_H_

#pragma once

namespace kt84 {
template <class Container>
struct ReversedRange {
protected:
    Container& container;
public:
    ReversedRange(Container& container) : container(container) {}
    auto begin() -> decltype(container.rbegin()) { return container.rbegin(); }
    auto end  () -> decltype(container.rend  ()) { return container.rend  (); }
};

template <class Container>
inline ReversedRange<Container> make_ReversedRange(Container& container) { return ReversedRange<Container>(container); }

}

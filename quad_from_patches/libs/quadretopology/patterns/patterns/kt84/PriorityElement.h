#pragma once

namespace kt84 {
    template <typename PriorityType, typename ElementType>
    struct PriorityElement {
        PriorityType priority;
        ElementType element;
        bool operator<(const PriorityElement& rhs) const { return priority < rhs.priority; }
    };
    template <typename PriorityType, typename ElementType>
    inline PriorityElement<PriorityType, ElementType> make_PriorityElement(PriorityType priority, const ElementType& element) {
        PriorityElement<PriorityType, ElementType> result = { priority, element };
        return result;
    }
}

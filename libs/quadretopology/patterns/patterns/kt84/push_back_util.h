#pragma once
namespace kt84 {
    namespace internal {
        template <class Container>
        struct PushBackUtil {
            PushBackUtil(Container& container) : container(container) {}
            template <class Element>
            PushBackUtil& operator<<(const Element& element) {
                container.push_back(element);
                return *this;
            }
        protected:
            Container& container;
        };
    }
    template <class Container>
    inline internal::PushBackUtil<Container> push_back_util(Container& container) { return internal::PushBackUtil<Container>(container); }
}

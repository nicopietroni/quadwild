#pragma once
#include <boost/lexical_cast.hpp>

namespace kt84 {
    namespace internal {
        template <class IStream>
        struct SafeIStream {
            SafeIStream(IStream& is) : is(is) {}
            template <typename T>
            SafeIStream& operator>>(T& value) {
                std::string value_str;
                is >> value_str;
                try {
                    value = boost::lexical_cast<T>(value_str);
                    fail_flag = false;
                } catch (const boost::bad_lexical_cast&){
                    fail_flag = true;
                }
                return *this;
            }
            bool fail() const { return fail_flag; }
        private:
            IStream& is;
            bool fail_flag;
        };
    }
    template <class IStream>
    inline internal::SafeIStream<IStream> make_safe_istream(IStream& is) { return internal::SafeIStream<IStream>(is); }
}

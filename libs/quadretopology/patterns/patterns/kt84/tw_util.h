#pragma once
#include <cassert>
#include <functional>
#include <AntTweakBar.h>

namespace kt84 {
    namespace tw_util {
        // AddButton |
        //-----------+
        namespace internal {
            inline void TW_CALL AddButton_sub(void *clientData) {
                auto& func = *static_cast<std::function<void(void)>*>(clientData);
                func();
            }
        }
        inline int AddButton(TwBar *bar, const char *name, std::function<void(void)> func, const char *def) {
            static const int buffer_size = 1000;
            static std::function<void(void)> buffer[buffer_size];
            static int counter = 0;
            assert(counter < buffer_size);
            buffer[counter] = func;
            return TwAddButton(bar, name, internal::AddButton_sub, &buffer[counter++], def);
        }
        // AddVarCB |
        //----------+
        namespace internal {
            template <typename T>
            struct FuncPair {
                std::function<void(const T&)> set_func;
                std::function<void(      T&)> get_func;
            };
            template <typename T>
            inline FuncPair<T> make_FuncPair(std::function<void(const T&)> set_func, std::function<void(T&)> get_func) {
                FuncPair<T> result = { set_func, get_func };
                return result;
            }
            template <typename T>
            inline void TW_CALL AddVar_sub_set(const void *value, void *clientData) {
                static_cast<FuncPair<T>*>(clientData)->set_func(*static_cast<const T*>(value));
            }
            template <typename T>
            inline void TW_CALL AddVar_sub_get(void *value, void *clientData) {
                static_cast<FuncPair<T>*>(clientData)->get_func(*static_cast<      T*>(value));
            }
        }
        template <typename T>
        inline int AddVarCB(TwBar *bar, const char *name, TwType type, std::function<void(const T&)> set_func, std::function<void(T&)> get_func, const char *def) {
            static const int buffer_size = 1000;
            static internal::FuncPair<T> buffer[buffer_size];
            static int counter = 0;
            assert(counter < buffer_size);
            buffer[counter] = internal::make_FuncPair(set_func, get_func);
            return TwAddVarCB(bar, name, type, internal::AddVar_sub_set<T>, internal::AddVar_sub_get<T>, &buffer[counter++], def);
        }
        // AddVarCB_default: default get func and custom func after default set func |
        //---------------------------------------------------------------------------+
        //template <typename T> // This doesn't work in Visual Studio 2012 probably because of a compiler bug... (search keyword: _VARIADIC_EXPAND_P1_0)
        //inline int AddVarCB_default(TwBar *bar, const char *name, TwType type, T& client_value, std::function<void(void)> func_after_set, const char *def) {
        //    return AddVarCB<T>(bar, name, type, [&](const T& value) { client_value = value; func_after_set(); },
        //                                        [&](      T& value) { value = client_value; }, def);
        //}
        namespace internal {
            template <typename T>
            struct ValueFunc {
                T* value;
                std::function<void(void)> func;
            };
            template <typename T>
            inline ValueFunc<T> make_ValueFunc(T& value, std::function<void(void)> func) {
                ValueFunc<T> result = { &value, func };
                return result;
            }
            template <typename T>
            inline void TW_CALL AddVar_default_sub_set(const void *value, void *clientData) {
                ValueFunc<T>& value_func = *static_cast<ValueFunc<T>*>(clientData);
                *value_func.value = *static_cast<const T*>(value);
                value_func.func();
            }
            template <typename T>
            inline void TW_CALL AddVar_default_sub_get(void *value, void *clientData) {
                *static_cast<T*>(value) = *static_cast<ValueFunc<T>*>(clientData)->value;
            }
        }
        template <typename T>
        inline int AddVarCB_default(TwBar *bar, const char *name, TwType type, T& client_value, std::function<void(void)> func_after_set, const char *def) {
            static const int buffer_size = 1000;
            static internal::ValueFunc<T> buffer[buffer_size];
            static int counter = 0;
            assert(counter < buffer_size);
            buffer[counter] = internal::make_ValueFunc(client_value, func_after_set);
            return TwAddVarCB(bar, name, type, internal::AddVar_default_sub_set<T>, internal::AddVar_default_sub_get<T>, &buffer[counter++], def);
        }
        // API with TW_TYPE_*** omitted for certain types
        namespace internal {
            template <typename T> struct ETwTypeT;
            template <> struct ETwTypeT<bool  > { static const ETwType Value = TW_TYPE_BOOLCPP; };
            template <> struct ETwTypeT<int   > { static const ETwType Value = TW_TYPE_INT32  ; };
            template <> struct ETwTypeT<float > { static const ETwType Value = TW_TYPE_FLOAT  ; };
            template <> struct ETwTypeT<double> { static const ETwType Value = TW_TYPE_DOUBLE ; };
        }
        template <typename T>
        inline int AddVarCB(TwBar *bar, const char *name, std::function<void(const T&)> set_func, std::function<void(T&)> get_func, const char *def) {
            return AddVarCB<T>(bar, name, internal::ETwTypeT<T>::Value, set_func, get_func, def);
        }
        template <typename T>
        inline int AddVarCB_default(TwBar *bar, const char *name, T& client_value, std::function<void(void)> func_after_set, const char *def) {
            return AddVarCB_default(bar, name, internal::ETwTypeT<T>::Value, client_value, func_after_set, def);
        }
        // little shortcut for for custom enums
        template <typename T>
        inline int AddVarCB(TwBar *bar, const char *name, const char *enum_type_name, std::function<void(const T&)> set_func, std::function<void(T&)> get_func, const char *def) {
            return AddVarCB<T>(bar, name, TwDefineEnum(enum_type_name, 0, 0), set_func, get_func, def);
        }
        template <typename T>
        inline int AddVarCB_default(TwBar *bar, const char *name, const char *enum_type_name, T& client_value, std::function<void(void)> func_after_set, const char *def) {
            return AddVarCB_default(bar, name, TwDefineEnum(enum_type_name, 0, 0), client_value, func_after_set, def);
        }
    }
}

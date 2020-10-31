#pragma once

namespace kt84 {

namespace loop_util {
    // find_if for two synchronized sequences
    template<class InputIterator1, class InputIterator2, class Predicate>
    inline InputIterator1 find_if(
        InputIterator1 first1,
        InputIterator1 last1,
        InputIterator2 first2,
        Predicate pred)
    {
        auto i1(first1);
        auto i2(first2);
        
        for ( ; i1 != last1; ++i1, ++i2)
            if (pred(*i1, *i2)) break;
        
        return i1;
    }
    
    // accumulate for two synchronized sequences
    template<class InputIterator1, class InputIterator2, class Type, class TernaryOperation> 
    inline Type accumulate(
        InputIterator1 first1,
        InputIterator1 last1,
        InputIterator2 first2,
        Type val,
        TernaryOperation ternary_op)
    {
        Type partial_result(val);
        
        auto i1(first1);
        auto i2(first2);
        
        for ( ; i1 != last1; ++i1, ++i2)
            partial_result = ternary_op(partial_result, *i1, *i2);
        
        return partial_result;
    }
    
    // for_each for two synchronized sequences
    template <class InputIterator1, class InputIterator2, class BinaryFunction>
    inline BinaryFunction for_each(
        InputIterator1 first1,
        InputIterator1 last1,
        InputIterator2 first2,
        BinaryFunction func)
    {
        auto i1(first1);
        auto i2(first2);
        
        for ( ; i1 != last1; ++i1, ++i2)
            func(*i1, *i2);
        
        return func;
    }
    
    // for_each for three synchronized sequences
    template <class InputIterator1, class InputIterator2, class InputIterator3, class TernaryFunction>
    inline TernaryFunction for_each(
        InputIterator1 first1,
        InputIterator1 last1,
        InputIterator2 first2,
        InputIterator3 first3,
        TernaryFunction func)
    {
        auto i1(first1);
        auto i2(first2);
        auto i3(first3);
        
        for ( ; i1 != last1; ++i1, ++i2, ++i3)
            func(*i1, *i2, *i3);
        
        return func;
    }
    
    // tansform for three synchronized sequences
    template<class InputIterator1, class InputIterator2, class InputIterator3, class OutputIterator, class TernaryFunction> 
    inline OutputIterator transform(
        InputIterator1 first1,
        InputIterator1 last1,
        InputIterator2 first2,
        InputIterator3 first3,
        OutputIterator result,
        TernaryFunction func)
    {
        auto i1(first1);
        auto i2(first2);
        auto i3(first3);
        auto o (result);
        
        for ( ; i1 != last1; ++i1, ++i2, ++i3, ++o)
            *o = func(*i1, *i2, *i3);
        
        return o;
    }
    
    // *_with_index versions, just for convenience... |
    //------------------------------------------------+
    
    // accumulate with index itself
    template<class IndexType, class Type, class BinaryOperation> 
    inline Type accumulate_with_index( 
        IndexType first_index,
        IndexType last_index,
        Type val,
        BinaryOperation binary_op)
    {
        Type partial_result(val);
        
        auto index(first_index);
        
        for ( ; index != last_index; ++index)
            partial_result = binary_op(partial_result, index);
        
        return partial_result;
    }
    
    // accumulate over a sequence with index
    template<class IndexType, class InputIterator, class Type, class TernaryOperation> 
    inline Type accumulate_with_index( 
        IndexType first_index,
        IndexType last_index,
        InputIterator first,
        Type val,
        TernaryOperation ternary_op)
    {
        Type partial_result(val);
        
        auto index(first_index);
        auto i(first);
        
        for ( ; index != last_index; ++index, ++i)
            partial_result = ternary_op(partial_result, index, *i);
        
        return partial_result;
    }
    
    // for_each with index itself
    template<class IndexType, class Function> 
    inline Function for_each_with_index( 
        IndexType first_index,
        IndexType last_index,
        Function func)
    {
        auto index(first_index);
        
        for ( ; index != last_index; ++index)
            func(index);
        
        return func;
    }

    // for_each over a sequence with index
    template<class IndexType, class InputIterator, class Function> 
    inline Function for_each_with_index( 
        IndexType first_index,
        IndexType last_index,
        InputIterator first,
        Function func)
    {
        auto index(first_index);
        auto i(first);
        
        for ( ; index != last_index; ++index, ++i)
            func(index, *i);
        
        return func;
    }
    
    // generate
    template<class IndexType, class ForwardIterator, class Generator> 
    inline void generate_with_index(
        IndexType first_index,
        IndexType last_index,
        ForwardIterator first,
        Generator gen)
    {
        auto index(first_index);
        auto i(first);
        
        for ( ; index != last_index; ++index, ++i)
            *i = gen(index);
    }

    // tansform over a sequence with index
    template<class IndexType, class InputIterator, class OutputIterator, class TernaryFunction> 
    inline OutputIterator transform_with_index(
        IndexType first_index,
        IndexType last_index,
        InputIterator first,
        OutputIterator result,
        TernaryFunction func)
    {
        auto index(first_index);
        auto i(first);
        auto o(result);
        
        for ( ; index != last_index; ++index, ++i, ++o)
            *o = func(index, *i);
        
        return o;
    }
    
    // tansform_with_index for two synchronized sequences
    template<class IndexType, class InputIterator1, class InputIterator2, class OutputIterator, class TernaryFunction> 
    inline OutputIterator transform_with_index(
        IndexType first_index,
        IndexType last_index,
        InputIterator1 first1,
        InputIterator2 first2,
        OutputIterator result,
        TernaryFunction func)
    {
        auto index(first_index);
        auto i1(first1);
        auto i2(first2);
        auto o (result);
        
        for ( ; index != last_index; ++index, ++i1, ++i2, ++o)
            *o = func(index, *i1, *i2);
        
        return o;
    }
    
    // simple for loop with index
    template<class IndexType, class Function> 
    inline Function loop_with_index( 
        IndexType first_index,
        IndexType last_index,
        Function func)
    {
        auto index(first_index);
        
        for ( ; index != last_index; ++index)
            func(index);
        
        return func;
    }
    
    // simple for loop without index
    template<class Function> 
    inline Function loop_without_index( 
        int count,
        Function func)
    {
        for (int i = 0; i < count; ++i)
            func();
        
        return func;
    }
}

}

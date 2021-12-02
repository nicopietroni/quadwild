#pragma once
#include <Eigen/Core>

namespace Eigen {

//-------------------------------------+
// mostly copied from src/Core/Array.h |
//                                     |
#define EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)   \
/** \ingroup arraytypedefs */                                    \
typedef Array<Type, Size, Size> Array##SizeSuffix##SizeSuffix##TypeSuffix;  \
/** \ingroup arraytypedefs */                                    \
typedef Array<Type, Size, 1>    Array##SizeSuffix##TypeSuffix;

#define EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, Size)         \
/** \ingroup arraytypedefs */                                    \
typedef Array<Type, Size, Dynamic> Array##Size##X##TypeSuffix;  \
/** \ingroup arraytypedefs */                                    \
typedef Array<Type, Dynamic, Size> Array##X##Size##TypeSuffix;

#define EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 2, 2) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 3, 3) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 4, 4) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 5, 5) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 6, 6) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 7, 7) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 8, 8) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 9, 9) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 10, 10) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, Dynamic, X) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 2) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 3) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 4) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 5) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 6) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 7) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 8) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 9) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 10)

EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(bool, b)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(char, c)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(unsigned int, ui)
EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES(unsigned char, uc)

#define EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10(Type, TypeSuffix) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 5, 5) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 6, 6) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 7, 7) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 8, 8) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 9, 9) \
EIGEN_MAKE_ARRAY_TYPEDEFS(Type, TypeSuffix, 10, 10) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 5) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 6) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 7) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 8) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 9) \
EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS(Type, TypeSuffix, 10)

EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10(int,                  i)
EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10(float,                f)
EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10(double,               d)
EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10(std::complex<float>,  cf)
EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10(std::complex<double>, cd)

#undef EIGEN_MAKE_ARRAY_TYPEDEFS_ALL_SIZES
#undef EIGEN_MAKE_ARRAY_TYPEDEFS_SIZES_5_10
#undef EIGEN_MAKE_ARRAY_FIXED_TYPEDEFS
#undef EIGEN_MAKE_ARRAY_TYPEDEFS
//                                     |
// mostly copied from src/Core/Array.h |
//-------------------------------------+

//--------------------------------------+
// mostly copied from src/Core/Matrix.h |
//                                      |
#define EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Size, SizeSuffix)   \
/** \ingroup matrixtypedefs */                                    \
typedef Matrix<Type, Size, Size> Matrix##SizeSuffix##TypeSuffix;  \
/** \ingroup matrixtypedefs */                                    \
typedef Matrix<Type, Size, 1>    Vector##SizeSuffix##TypeSuffix;  \
/** \ingroup matrixtypedefs */                                    \
typedef Matrix<Type, 1, Size>    RowVector##SizeSuffix##TypeSuffix;

#define EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, Size)         \
/** \ingroup matrixtypedefs */                                    \
typedef Matrix<Type, Size, Dynamic> Matrix##Size##X##TypeSuffix;  \
/** \ingroup matrixtypedefs */                                    \
typedef Matrix<Type, Dynamic, Size> Matrix##X##Size##TypeSuffix;

#define EIGEN_MAKE_TYPEDEFS_ALL_SIZES(Type, TypeSuffix) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 2, 2) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 3, 3) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 4, 4) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 5, 5) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 6, 6) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 7, 7) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 8, 8) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 9, 9) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 10, 10) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, Dynamic, X) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 2) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 3) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 4) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 5) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 6) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 7) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 8) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 9) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 10)

EIGEN_MAKE_TYPEDEFS_ALL_SIZES(bool, b)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(char, c)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(unsigned char, uc)
EIGEN_MAKE_TYPEDEFS_ALL_SIZES(unsigned int, ui)

#define EIGEN_MAKE_TYPEDEFS_SIZES_5_10(Type, TypeSuffix) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 5, 5) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 6, 6) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 7, 7) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 8, 8) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 9, 9) \
EIGEN_MAKE_TYPEDEFS(Type, TypeSuffix, 10, 10) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 5) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 6) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 7) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 8) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 9) \
EIGEN_MAKE_FIXED_TYPEDEFS(Type, TypeSuffix, 10)

EIGEN_MAKE_TYPEDEFS_SIZES_5_10(int,                  i)
EIGEN_MAKE_TYPEDEFS_SIZES_5_10(float,                f)
EIGEN_MAKE_TYPEDEFS_SIZES_5_10(double,               d)
EIGEN_MAKE_TYPEDEFS_SIZES_5_10(std::complex<float>,  cf)
EIGEN_MAKE_TYPEDEFS_SIZES_5_10(std::complex<double>, cd)

#undef EIGEN_MAKE_TYPEDEFS_ALL_SIZES
#undef EIGEN_MAKE_TYPEDEFS_SIZES_5_10
#undef EIGEN_MAKE_TYPEDEFS
#undef EIGEN_MAKE_FIXED_TYPEDEFS
//                                      |
// mostly copied from src/Core/Matrix.h |
//--------------------------------------+

}


namespace kt84 {

#define EIGEN_MAKE_VECTOR_5(Type, TypeSuffix) \
    inline Eigen::Vector5##TypeSuffix make_Vector5##TypeSuffix (const Type& v0, const Type& v1, const Type& v2, const Type& v3, const Type& v4) { \
        Eigen::Vector5##TypeSuffix vec; \
        vec(0) = v0; \
        vec(1) = v1; \
        vec(2) = v2; \
        vec(3) = v3; \
        vec(4) = v4; \
        return vec; \
    }
#define EIGEN_MAKE_VECTOR_6(Type, TypeSuffix) \
    inline Eigen::Vector6##TypeSuffix make_Vector6##TypeSuffix (const Type& v0, const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5) { \
        Eigen::Vector6##TypeSuffix vec; \
        vec(0) = v0; \
        vec(1) = v1; \
        vec(2) = v2; \
        vec(3) = v3; \
        vec(4) = v4; \
        vec(5) = v5; \
        return vec; \
    }
#define EIGEN_MAKE_VECTOR_7(Type, TypeSuffix) \
    inline Eigen::Vector7##TypeSuffix make_Vector7##TypeSuffix (const Type& v0, const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5, const Type& v6) { \
        Eigen::Vector7##TypeSuffix vec; \
        vec(0) = v0; \
        vec(1) = v1; \
        vec(2) = v2; \
        vec(3) = v3; \
        vec(4) = v4; \
        vec(5) = v5; \
        vec(6) = v6; \
        return vec; \
    }
#define EIGEN_MAKE_VECTOR_8(Type, TypeSuffix) \
    inline Eigen::Vector8##TypeSuffix make_Vector8##TypeSuffix (const Type& v0, const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5, const Type& v6, const Type& v7) { \
        Eigen::Vector8##TypeSuffix vec; \
        vec(0) = v0; \
        vec(1) = v1; \
        vec(2) = v2; \
        vec(3) = v3; \
        vec(4) = v4; \
        vec(5) = v5; \
        vec(6) = v6; \
        vec(7) = v7; \
        return vec; \
    }
#define EIGEN_MAKE_VECTOR_9(Type, TypeSuffix) \
    inline Eigen::Vector9##TypeSuffix make_Vector9##TypeSuffix (const Type& v0, const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5, const Type& v6, const Type& v7, const Type& v8) { \
        Eigen::Vector9##TypeSuffix vec; \
        vec(0) = v0; \
        vec(1) = v1; \
        vec(2) = v2; \
        vec(3) = v3; \
        vec(4) = v4; \
        vec(5) = v5; \
        vec(6) = v6; \
        vec(7) = v7; \
        vec(8) = v8; \
        return vec; \
    }
#define EIGEN_MAKE_VECTOR_10(Type, TypeSuffix) \
    inline Eigen::Vector10##TypeSuffix make_Vector10##TypeSuffix (const Type& v0, const Type& v1, const Type& v2, const Type& v3, const Type& v4, const Type& v5, const Type& v6, const Type& v7, const Type& v8, const Type& v9) { \
        Eigen::Vector10##TypeSuffix vec; \
        vec(0) = v0; \
        vec(1) = v1; \
        vec(2) = v2; \
        vec(3) = v3; \
        vec(4) = v4; \
        vec(5) = v5; \
        vec(6) = v6; \
        vec(7) = v7; \
        vec(8) = v8; \
        vec(9) = v9; \
        return vec; \
    }

#define EIGEN_MAKE_VECTOR_5_10(Type, TypeSuffix) \
EIGEN_MAKE_VECTOR_5(Type, TypeSuffix) \
EIGEN_MAKE_VECTOR_6(Type, TypeSuffix) \
EIGEN_MAKE_VECTOR_7(Type, TypeSuffix) \
EIGEN_MAKE_VECTOR_8(Type, TypeSuffix) \
EIGEN_MAKE_VECTOR_9(Type, TypeSuffix) \
EIGEN_MAKE_VECTOR_10(Type, TypeSuffix)

EIGEN_MAKE_VECTOR_5_10(int                 , i )
EIGEN_MAKE_VECTOR_5_10(float               , f )
EIGEN_MAKE_VECTOR_5_10(double              , d )
EIGEN_MAKE_VECTOR_5_10(std::complex<float> , cf)
EIGEN_MAKE_VECTOR_5_10(std::complex<double>, cd)
EIGEN_MAKE_VECTOR_5_10(bool                , b )
EIGEN_MAKE_VECTOR_5_10(char                , c )
EIGEN_MAKE_VECTOR_5_10(unsigned char       , uc)
EIGEN_MAKE_VECTOR_5_10(unsigned int        , ui)

#undef EIGEN_MAKE_VECTOR_5_10
#undef EIGEN_MAKE_VECTOR_5
#undef EIGEN_MAKE_VECTOR_6
#undef EIGEN_MAKE_VECTOR_7
#undef EIGEN_MAKE_VECTOR_8
#undef EIGEN_MAKE_VECTOR_9
#undef EIGEN_MAKE_VECTOR_10

#define EIGEN_MAKE_MATRIX_2(Type, TypeSuffix) \
    inline Eigen::Matrix2##TypeSuffix make_Matrix2##TypeSuffix ( \
        const Type& v00, const Type& v01, \
        const Type& v10, const Type& v11) \
    { \
        Eigen::Matrix2##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; \
        mat(1,0)=v10; mat(1,1)=v11; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_3(Type, TypeSuffix) \
    inline Eigen::Matrix3##TypeSuffix make_Matrix3##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, \
        const Type& v10, const Type& v11, const Type& v12, \
        const Type& v20, const Type& v21, const Type& v22) \
    { \
        Eigen::Matrix3##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_4(Type, TypeSuffix) \
    inline Eigen::Matrix4##TypeSuffix make_Matrix4##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33) \
    { \
        Eigen::Matrix4##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_5(Type, TypeSuffix) \
    inline Eigen::Matrix5##TypeSuffix make_Matrix5##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, const Type& v04, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, const Type& v14, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, const Type& v24, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33, const Type& v34, \
        const Type& v40, const Type& v41, const Type& v42, const Type& v43, const Type& v44) \
    { \
        Eigen::Matrix5##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; mat(0,4)=v04; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; mat(1,4)=v14; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; mat(2,4)=v24; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; mat(3,4)=v34; \
        mat(4,0)=v40; mat(4,1)=v41; mat(4,2)=v42; mat(4,3)=v43; mat(4,4)=v44; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_6(Type, TypeSuffix) \
    inline Eigen::Matrix6##TypeSuffix make_Matrix6##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, const Type& v04, const Type& v05, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, const Type& v14, const Type& v15, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, const Type& v24, const Type& v25, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33, const Type& v34, const Type& v35, \
        const Type& v40, const Type& v41, const Type& v42, const Type& v43, const Type& v44, const Type& v45, \
        const Type& v50, const Type& v51, const Type& v52, const Type& v53, const Type& v54, const Type& v55) \
    { \
        Eigen::Matrix6##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; mat(0,4)=v04; mat(0,5)=v05; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; mat(1,4)=v14; mat(1,5)=v15; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; mat(2,4)=v24; mat(2,5)=v25; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; mat(3,4)=v34; mat(3,5)=v35; \
        mat(4,0)=v40; mat(4,1)=v41; mat(4,2)=v42; mat(4,3)=v43; mat(4,4)=v44; mat(4,5)=v45; \
        mat(5,0)=v50; mat(5,1)=v51; mat(5,2)=v52; mat(5,3)=v53; mat(5,4)=v54; mat(5,5)=v55; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_7(Type, TypeSuffix) \
    inline Eigen::Matrix7##TypeSuffix make_Matrix7##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, const Type& v04, const Type& v05, const Type& v06, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, const Type& v14, const Type& v15, const Type& v16, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, const Type& v24, const Type& v25, const Type& v26, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33, const Type& v34, const Type& v35, const Type& v36, \
        const Type& v40, const Type& v41, const Type& v42, const Type& v43, const Type& v44, const Type& v45, const Type& v46, \
        const Type& v50, const Type& v51, const Type& v52, const Type& v53, const Type& v54, const Type& v55, const Type& v56, \
        const Type& v60, const Type& v61, const Type& v62, const Type& v63, const Type& v64, const Type& v65, const Type& v66) \
    { \
        Eigen::Matrix7##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; mat(0,4)=v04; mat(0,5)=v05; mat(0,6)=v06; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; mat(1,4)=v14; mat(1,5)=v15; mat(1,6)=v16; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; mat(2,4)=v24; mat(2,5)=v25; mat(2,6)=v26; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; mat(3,4)=v34; mat(3,5)=v35; mat(3,6)=v36; \
        mat(4,0)=v40; mat(4,1)=v41; mat(4,2)=v42; mat(4,3)=v43; mat(4,4)=v44; mat(4,5)=v45; mat(4,6)=v46; \
        mat(5,0)=v50; mat(5,1)=v51; mat(5,2)=v52; mat(5,3)=v53; mat(5,4)=v54; mat(5,5)=v55; mat(5,6)=v56; \
        mat(6,0)=v60; mat(6,1)=v61; mat(6,2)=v62; mat(6,3)=v63; mat(6,4)=v64; mat(6,5)=v65; mat(6,6)=v66; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_8(Type, TypeSuffix) \
    inline Eigen::Matrix8##TypeSuffix make_Matrix8##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, const Type& v04, const Type& v05, const Type& v06, const Type& v07, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, const Type& v14, const Type& v15, const Type& v16, const Type& v17, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, const Type& v24, const Type& v25, const Type& v26, const Type& v27, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33, const Type& v34, const Type& v35, const Type& v36, const Type& v37, \
        const Type& v40, const Type& v41, const Type& v42, const Type& v43, const Type& v44, const Type& v45, const Type& v46, const Type& v47, \
        const Type& v50, const Type& v51, const Type& v52, const Type& v53, const Type& v54, const Type& v55, const Type& v56, const Type& v57, \
        const Type& v60, const Type& v61, const Type& v62, const Type& v63, const Type& v64, const Type& v65, const Type& v66, const Type& v67, \
        const Type& v70, const Type& v71, const Type& v72, const Type& v73, const Type& v74, const Type& v75, const Type& v76, const Type& v77) \
    { \
        Eigen::Matrix8##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; mat(0,4)=v04; mat(0,5)=v05; mat(0,6)=v06; mat(0,7)=v07; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; mat(1,4)=v14; mat(1,5)=v15; mat(1,6)=v16; mat(1,7)=v17; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; mat(2,4)=v24; mat(2,5)=v25; mat(2,6)=v26; mat(2,7)=v27; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; mat(3,4)=v34; mat(3,5)=v35; mat(3,6)=v36; mat(3,7)=v37; \
        mat(4,0)=v40; mat(4,1)=v41; mat(4,2)=v42; mat(4,3)=v43; mat(4,4)=v44; mat(4,5)=v45; mat(4,6)=v46; mat(4,7)=v47; \
        mat(5,0)=v50; mat(5,1)=v51; mat(5,2)=v52; mat(5,3)=v53; mat(5,4)=v54; mat(5,5)=v55; mat(5,6)=v56; mat(5,7)=v57; \
        mat(6,0)=v60; mat(6,1)=v61; mat(6,2)=v62; mat(6,3)=v63; mat(6,4)=v64; mat(6,5)=v65; mat(6,6)=v66; mat(6,7)=v67; \
        mat(7,0)=v70; mat(7,1)=v71; mat(7,2)=v72; mat(7,3)=v73; mat(7,4)=v74; mat(7,5)=v75; mat(7,6)=v76; mat(7,7)=v77; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_9(Type, TypeSuffix) \
    inline Eigen::Matrix9##TypeSuffix make_Matrix9##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, const Type& v04, const Type& v05, const Type& v06, const Type& v07, const Type& v08, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, const Type& v14, const Type& v15, const Type& v16, const Type& v17, const Type& v18, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, const Type& v24, const Type& v25, const Type& v26, const Type& v27, const Type& v28, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33, const Type& v34, const Type& v35, const Type& v36, const Type& v37, const Type& v38, \
        const Type& v40, const Type& v41, const Type& v42, const Type& v43, const Type& v44, const Type& v45, const Type& v46, const Type& v47, const Type& v48, \
        const Type& v50, const Type& v51, const Type& v52, const Type& v53, const Type& v54, const Type& v55, const Type& v56, const Type& v57, const Type& v58, \
        const Type& v60, const Type& v61, const Type& v62, const Type& v63, const Type& v64, const Type& v65, const Type& v66, const Type& v67, const Type& v68, \
        const Type& v70, const Type& v71, const Type& v72, const Type& v73, const Type& v74, const Type& v75, const Type& v76, const Type& v77, const Type& v78, \
        const Type& v80, const Type& v81, const Type& v82, const Type& v83, const Type& v84, const Type& v85, const Type& v86, const Type& v87, const Type& v88) \
    { \
        Eigen::Matrix9##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; mat(0,4)=v04; mat(0,5)=v05; mat(0,6)=v06; mat(0,7)=v07; mat(0,8)=v08; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; mat(1,4)=v14; mat(1,5)=v15; mat(1,6)=v16; mat(1,7)=v17; mat(1,8)=v18; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; mat(2,4)=v24; mat(2,5)=v25; mat(2,6)=v26; mat(2,7)=v27; mat(2,8)=v28; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; mat(3,4)=v34; mat(3,5)=v35; mat(3,6)=v36; mat(3,7)=v37; mat(3,8)=v38; \
        mat(4,0)=v40; mat(4,1)=v41; mat(4,2)=v42; mat(4,3)=v43; mat(4,4)=v44; mat(4,5)=v45; mat(4,6)=v46; mat(4,7)=v47; mat(4,8)=v48; \
        mat(5,0)=v50; mat(5,1)=v51; mat(5,2)=v52; mat(5,3)=v53; mat(5,4)=v54; mat(5,5)=v55; mat(5,6)=v56; mat(5,7)=v57; mat(5,8)=v58; \
        mat(6,0)=v60; mat(6,1)=v61; mat(6,2)=v62; mat(6,3)=v63; mat(6,4)=v64; mat(6,5)=v65; mat(6,6)=v66; mat(6,7)=v67; mat(6,8)=v68; \
        mat(7,0)=v70; mat(7,1)=v71; mat(7,2)=v72; mat(7,3)=v73; mat(7,4)=v74; mat(7,5)=v75; mat(7,6)=v76; mat(7,7)=v77; mat(7,8)=v78; \
        mat(8,0)=v80; mat(8,1)=v81; mat(8,2)=v82; mat(8,3)=v83; mat(8,4)=v84; mat(8,5)=v85; mat(8,6)=v86; mat(8,7)=v87; mat(8,8)=v88; \
        return mat; \
    }
#define EIGEN_MAKE_MATRIX_10(Type, TypeSuffix) \
    inline Eigen::Matrix10##TypeSuffix make_Matrix10##TypeSuffix ( \
        const Type& v00, const Type& v01, const Type& v02, const Type& v03, const Type& v04, const Type& v05, const Type& v06, const Type& v07, const Type& v08, const Type& v09, \
        const Type& v10, const Type& v11, const Type& v12, const Type& v13, const Type& v14, const Type& v15, const Type& v16, const Type& v17, const Type& v18, const Type& v19, \
        const Type& v20, const Type& v21, const Type& v22, const Type& v23, const Type& v24, const Type& v25, const Type& v26, const Type& v27, const Type& v28, const Type& v29, \
        const Type& v30, const Type& v31, const Type& v32, const Type& v33, const Type& v34, const Type& v35, const Type& v36, const Type& v37, const Type& v38, const Type& v39, \
        const Type& v40, const Type& v41, const Type& v42, const Type& v43, const Type& v44, const Type& v45, const Type& v46, const Type& v47, const Type& v48, const Type& v49, \
        const Type& v50, const Type& v51, const Type& v52, const Type& v53, const Type& v54, const Type& v55, const Type& v56, const Type& v57, const Type& v58, const Type& v59, \
        const Type& v60, const Type& v61, const Type& v62, const Type& v63, const Type& v64, const Type& v65, const Type& v66, const Type& v67, const Type& v68, const Type& v69, \
        const Type& v70, const Type& v71, const Type& v72, const Type& v73, const Type& v74, const Type& v75, const Type& v76, const Type& v77, const Type& v78, const Type& v79, \
        const Type& v80, const Type& v81, const Type& v82, const Type& v83, const Type& v84, const Type& v85, const Type& v86, const Type& v87, const Type& v88, const Type& v89, \
        const Type& v90, const Type& v91, const Type& v92, const Type& v93, const Type& v94, const Type& v95, const Type& v96, const Type& v97, const Type& v98, const Type& v99) \
    { \
        Eigen::Matrix10##TypeSuffix mat; \
        mat(0,0)=v00; mat(0,1)=v01; mat(0,2)=v02; mat(0,3)=v03; mat(0,4)=v04; mat(0,5)=v05; mat(0,6)=v06; mat(0,7)=v07; mat(0,8)=v08; mat(0,9)=v09; \
        mat(1,0)=v10; mat(1,1)=v11; mat(1,2)=v12; mat(1,3)=v13; mat(1,4)=v14; mat(1,5)=v15; mat(1,6)=v16; mat(1,7)=v17; mat(1,8)=v18; mat(1,9)=v19; \
        mat(2,0)=v20; mat(2,1)=v21; mat(2,2)=v22; mat(2,3)=v23; mat(2,4)=v24; mat(2,5)=v25; mat(2,6)=v26; mat(2,7)=v27; mat(2,8)=v28; mat(2,9)=v29; \
        mat(3,0)=v30; mat(3,1)=v31; mat(3,2)=v32; mat(3,3)=v33; mat(3,4)=v34; mat(3,5)=v35; mat(3,6)=v36; mat(3,7)=v37; mat(3,8)=v38; mat(3,9)=v39; \
        mat(4,0)=v40; mat(4,1)=v41; mat(4,2)=v42; mat(4,3)=v43; mat(4,4)=v44; mat(4,5)=v45; mat(4,6)=v46; mat(4,7)=v47; mat(4,8)=v48; mat(4,9)=v49; \
        mat(5,0)=v50; mat(5,1)=v51; mat(5,2)=v52; mat(5,3)=v53; mat(5,4)=v54; mat(5,5)=v55; mat(5,6)=v56; mat(5,7)=v57; mat(5,8)=v58; mat(5,9)=v59; \
        mat(6,0)=v60; mat(6,1)=v61; mat(6,2)=v62; mat(6,3)=v63; mat(6,4)=v64; mat(6,5)=v65; mat(6,6)=v66; mat(6,7)=v67; mat(6,8)=v68; mat(6,9)=v69; \
        mat(7,0)=v70; mat(7,1)=v71; mat(7,2)=v72; mat(7,3)=v73; mat(7,4)=v74; mat(7,5)=v75; mat(7,6)=v76; mat(7,7)=v77; mat(7,8)=v78; mat(7,9)=v79; \
        mat(8,0)=v80; mat(8,1)=v81; mat(8,2)=v82; mat(8,3)=v83; mat(8,4)=v84; mat(8,5)=v85; mat(8,6)=v86; mat(8,7)=v87; mat(8,8)=v88; mat(8,9)=v89; \
        mat(9,0)=v90; mat(9,1)=v91; mat(9,2)=v92; mat(9,3)=v93; mat(9,4)=v94; mat(9,5)=v95; mat(9,6)=v96; mat(9,7)=v97; mat(9,8)=v98; mat(9,9)=v99; \
        return mat; \
    }

#define EIGEN_MAKE_MATRIX_5_10(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_2(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_3(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_4(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_5(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_6(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_7(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_8(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_9(Type, TypeSuffix) \
EIGEN_MAKE_MATRIX_10(Type, TypeSuffix)

EIGEN_MAKE_MATRIX_5_10(int                 , i )
EIGEN_MAKE_MATRIX_5_10(float               , f )
EIGEN_MAKE_MATRIX_5_10(double              , d )
EIGEN_MAKE_MATRIX_5_10(std::complex<float> , cf)
EIGEN_MAKE_MATRIX_5_10(std::complex<double>, cd)
EIGEN_MAKE_MATRIX_5_10(bool                , b )
EIGEN_MAKE_MATRIX_5_10(char                , c )
EIGEN_MAKE_MATRIX_5_10(unsigned char       , uc)
EIGEN_MAKE_MATRIX_5_10(unsigned int        , ui)

#undef EIGEN_MAKE_MATRIX_5_10
#undef EIGEN_MAKE_MATRIX_2
#undef EIGEN_MAKE_MATRIX_3
#undef EIGEN_MAKE_MATRIX_4
#undef EIGEN_MAKE_MATRIX_5
#undef EIGEN_MAKE_MATRIX_6
#undef EIGEN_MAKE_MATRIX_7
#undef EIGEN_MAKE_MATRIX_8
#undef EIGEN_MAKE_MATRIX_9
#undef EIGEN_MAKE_MATRIX_10

}

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "RBF.h"
#include "MLS.h"
#include "HermiteRBF.h"
#include "GeneralizedMLS.h"
#include <boost/mpl/list.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/int.hpp>
using namespace std;
using namespace kt84;
namespace mpl = boost::mpl;

// utility for interpolation algorithms |
//--------------------------------------+
template <template <int, int, class, int> class FuncT>
struct FuncUtil;
template <> struct FuncUtil<RBF> {
    static const char* getName() { return "RBF"; }
    template <class Func> static void add_constraint(Func& f, const typename Func::Point& p, const typename Func::Value& v, const typename Func::Gradient& g) { f.add_constraint(p, v); }
    template <class Func> static void preprocess(Func& f) { f.factorize_and_solve(); }
    template <class Func> static typename Func::Gradient gradient(const Func& f, const typename Func::Point& p) { return f.gradient(p); }
    template <class Func>
    static pair<double, double> getMaxError(const Func& f, double epsilon) {
        pair<double, double> result(0, 0);
        for (size_t i = 0; i < f.constraints.size(); ++i) {
            double error_v = (f(f.constraints[i].first) - f.constraints[i].second).norm();
            result.first = max<double>(result.first, error_v);
        }
        return result;
    }
};
template <> struct FuncUtil<MLS> {
    static const char* getName() { return "MLS"; }
    template <class Func> static void add_constraint(Func& f, const typename Func::Point& p, const typename Func::Value& v, const typename Func::Gradient& g) { f.add_constraint(p, v); }
    template <class Func> static void preprocess(Func& f) { }
    template <class Func> static typename Func::Gradient gradient(const Func& f, const typename Func::Point& p) { return typename Func::Gradient::Zero(); }
    template <class Func>
    static pair<double, double> getMaxError(const Func& f, double epsilon) {
        pair<double, double> result(0, 0);
        for (size_t i = 0; i < f.constraints.size(); ++i) {
            double error_v = (f(f.constraints[i].first) - f.constraints[i].second).norm();
            result.first = max<double>(result.first, error_v);
        }
        return result;
    }
};
template <> struct FuncUtil<HermiteRBF> {
    static const char* getName() { return "HRBF"; }
    template <class Func> static void add_constraint(Func& f, const typename Func::Point& p, const typename Func::Value& v, const typename Func::Gradient& g) { f.add_constraint(p, v, g); }
    template <class Func> static void preprocess(Func& f) { f.factorize_and_solve(); }
    template <class Func> static typename Func::Gradient gradient(const Func& f, const typename Func::Point& p) { return f.gradient(p); }
    template <class Func>
    static pair<double, double> getMaxError(const Func& f, double epsilon) {
        pair<double, double> result(0, 0);
        for (size_t i = 0; i < f.constraints.size(); ++i) {
            double error_v = (f(f.constraints[i].get<0>()) - f.constraints[i].get<1>()).norm();
            double error_g = (f.gradient_fd(f.constraints[i].get<0>(), epsilon) - f.constraints[i].get<2>()).norm();
            result.first  = max<double>(result.first , error_v);
            result.second = max<double>(result.second, error_g);
        }
        return result;
    }
};
template <> struct FuncUtil<GeneralizedMLS> {
    static const char* getName() { return "GMLS"; }
    template <class Func> static void add_constraint(Func& f, const typename Func::Point& p, const typename Func::Value& v, const typename Func::Gradient& g) { f.add_constraint(p, v, g); }
    template <class Func> static void preprocess(Func& f) { }
    template <class Func> static typename Func::Gradient gradient(const Func& f, const typename Func::Point& p) { return typename Func::Gradient::Zero(); }
    template <class Func>
    static pair<double, double> getMaxError(const Func& f, double epsilon) {
        pair<double, double> result(0, 0);
        for (size_t i = 0; i < f.constraints.size(); ++i) {
            double error_v = (f(f.constraints[i].get<0>()) - f.constraints[i].get<1>()).norm();
            double error_g = (f.gradient_fd(f.constraints[i].get<0>(), epsilon) - f.constraints[i].get<2>()).norm();
            result.first  = max<double>(result.first , error_v);
            result.second = max<double>(result.second, error_g);
        }
        return result;
    }
};

// utility for RBF kernels |
//-------------------------+
template <class RBFKernel_Core>
struct KernelUtil;
template <> struct KernelUtil<RBFKernel_Gaussian      > { static const char* getName() { return "Gauss"; } static RBFKernel_Gaussian      ::Param getParam() { return RBFKernel_Gaussian      ::Param(10); } };
template <> struct KernelUtil<RBFKernel_SquaredInverse> { static const char* getName() { return "SqInv"; } static RBFKernel_SquaredInverse::Param getParam() { return RBFKernel_SquaredInverse::Param(20); } };
template <> struct KernelUtil<RBFKernel_Wendland      > { static const char* getName() { return "Wendl"; } static RBFKernel_Wendland      ::Param getParam() { return RBFKernel_Wendland      ::Param(1.5); } };
template <> struct KernelUtil<RBFKernel_Cubed         > { static const char* getName() { return "Cubed"; } static RBFKernel_Cubed         ::Param getParam() { return RBFKernel_Cubed         ::Param()   ; } };
template <> struct KernelUtil<RBFKernel_Identity      > { static const char* getName() { return "Ident"; } static RBFKernel_Identity      ::Param getParam() { return RBFKernel_Identity      ::Param()   ; } };
template <> struct KernelUtil<RBFKernel_SquaredLog    > { static const char* getName() { return "SqLog"; } static RBFKernel_SquaredLog    ::Param getParam() { return RBFKernel_SquaredLog    ::Param()   ; } };

// test functions |
//----------------+
template <int Dim, template <int, int, class, int> class FuncT>
struct Tester;

// 1D test
template <template <int, int, class, int> class FuncT>
struct Tester<1, FuncT> {
    template <class RBFKernel_Core, int DegreePolynomial>
    static void go(const typename RBFKernel_Core::Param& param) {
        typedef FuncT<1, 1, RBFKernel_Core, DegreePolynomial> Func;
        typedef Func::Point Point;
        typedef Func::Value Value;
        typedef Func::Gradient Gradient;
        
        Func f;
        f.kernel.param() = param;
        
        ifstream fin("input1d.txt");
        if (!fin)
            throw exception("input1d.txt not found!");
        while (!fin.eof()) {
            Point p;
            Value v;
            Gradient g;
            
            fin >> p[0] >> v[0] >> g[0];
            if (!fin.eof())
                FuncUtil<FuncT>::add_constraint(f, p, v, g);
        }
        
        FuncUtil<FuncT>::preprocess(f);
        
        stringstream fname;
        fname << "test1d_" << FuncUtil<FuncT>::getName() << "_" << KernelUtil<RBFKernel_Core>::getName() << "_p" << DegreePolynomial << ".txt";
        ofstream fout(fname.str().c_str());
        // report error first
        pair<double, double> maxError = FuncUtil<FuncT>::getMaxError(f, 0.00001);
        fout << "# max error (value, gradient): (" << maxError.first << ", " << maxError.second << ")" << endl << endl;
        // plot data
        const int N = 100;
        for (double i = 0; i <= N; ++i) {
            Func::Point point = Func::Point::Constant(i / N);
            Func::Value value = f(point);
            Func::Gradient gradient = FuncUtil<FuncT>::gradient(f, point);
            Func::Gradient gradient_fd = f.gradient_fd(point, 0.00001);
            fout
                << setw(12) << point[0] << " "
                << setw(12) << value[0] << " "
                << setw(12) << gradient[0] << " "
                << setw(12) << gradient_fd[0] << " "
                << endl;
        }
    }
    // generate gnuplot commands |
    //---------------------------+
    template <class RBFKernel_Core>
    static void genCmd(int DegreePolynomial) {
        stringstream fname;
        const char* funcName = FuncUtil<FuncT>::getName();
        const char* kernelName = KernelUtil<RBFKernel_Core>::getName();
        fname << "cmd1d_" << funcName << "_" << kernelName << ".txt";
        ofstream fout(fname.str().c_str());
        fout << "plot\\" << endl;
        for (int i = 0; i <= DegreePolynomial; ++i) {
            fout << "    \"test1d_" << funcName << "_" << kernelName << "_p" << i << ".txt\" w linespoints";
            if (i < DegreePolynomial)
                fout << ",\\";
            fout << endl;
        }
        fout << "pause -1" << endl;
        fout.close();
        fname.str("");
        fname << "cmd1d_" << funcName << "_" << kernelName << "_g.txt";
        fout.open(fname.str().c_str());
        fout << "plot\\" << endl;
        for (int i = 0; i <= DegreePolynomial; ++i) {
            fout << "    \"test1d_" << funcName << "_" << kernelName << "_p" << i << ".txt\" u 1:4 w linespoints";
            if (i < DegreePolynomial)
                fout << ",\\";
            fout << endl;
        }
        fout << "pause -1" << endl;
    }
};

// 2D test
template <template <int, int, class, int> class FuncT>
struct Tester<2, FuncT> {
    template <class RBFKernel_Core, int DegreePolynomial>
    static void go(const typename RBFKernel_Core::Param& param) {
        typedef FuncT<2, 1, RBFKernel_Core, DegreePolynomial> Func;
        typedef Func::Point Point;
        typedef Func::Value Value;
        typedef Func::Gradient Gradient;
        
        Func f;
        f.kernel.param() = param;
        
        ifstream fin("input2d.txt");
        if (!fin)
            throw exception("input2d.txt not found!");
        while (!fin.eof()) {
            Point p;
            Value v;
            Gradient g;
            
            fin >> p[0] >> p[1] >> v[0] >> g[0] >> g[1];
            if (!fin.eof())
                FuncUtil<FuncT>::add_constraint(f, p, v, g);
        }
        
        FuncUtil<FuncT>::preprocess(f);
        
        stringstream fname;
        fname << "test2d_" << FuncUtil<FuncT>::getName() << "_" << KernelUtil<RBFKernel_Core>::getName() << "_p" << DegreePolynomial << ".txt";
        ofstream fout(fname.str().c_str());
        // report error first
        pair<double, double> maxError = FuncUtil<FuncT>::getMaxError(f, 0.00001);
        fout << "# max error (value, gradient): (" << maxError.first << ", " << maxError.second << ")" << endl << endl;
        // plot data
        const int N = 32;
        for (double j = 0; j <= N; ++j) {
            for (double i = 0; i <= N; ++i) {
                Func::Point point = Func::Point(i / N, j / N);
                Func::Value value = f(point);
                Func::Gradient gradient = FuncUtil<FuncT>::gradient(f, point);
                Func::Gradient gradient_fd = f.gradient_fd(point, 0.00001);
                fout
                    << setw(12) << point[0] << " "
                    << setw(12) << point[1] << " "
                    << setw(12) << value[0] << " "
                    << setw(12) << gradient[0] << " "
                    << setw(12) << gradient[1] << " "
                    << setw(12) << gradient_fd[0] << " "
                    << setw(12) << gradient_fd[1] << " "
                    << endl;
            }
            fout << endl;
        }
    }
    // generate gnuplot commands |
    //---------------------------+
    template <class RBFKernel_Core>
    static void genCmd(int DegreePolynomial) {
        const char* funcName = FuncUtil<FuncT>::getName();
        const char* kernelName = KernelUtil<RBFKernel_Core>::getName();
        for (int i = 0; i <= DegreePolynomial; ++i) {
            // value
            stringstream fname;
            fname << "cmd2d_" << funcName << "_" << kernelName << "_p" << i << ".txt";
            ofstream fout(fname.str().c_str());
            fout << "set pm3d; unset surface; set pm3d hidden3d 100;set view 60, 320\n";
            fout << "splot \"test2d_" << funcName << "_" << kernelName << "_p" << i << ".txt\" u 1:2:3\n";
            fout << "pause -1\n";
            fout.close();
            // gradient
            fname.str("");
            fname << "cmd2d_" << funcName << "_" << kernelName << "_p" << i << "_g.txt";
            fout.open(fname.str().c_str());
            fout << "set pm3d; unset surface; set pm3d hidden3d 100;set view 60, 320;\n";
            fout << "splot \"test2d_" << funcName << "_" << kernelName << "_p" << i << "_g.txt\" u 1:2:6\n";
            fout << "pause -1\n";
        }
    }
};

// list of valid kernels for each algorithm |
//------------------------------------------+
template <template <int, int, class, int> class FuncT>
struct KernelList;
template <> struct KernelList<RBF> {
    typedef mpl::list<
        RBFKernel_Gaussian      ,
        RBFKernel_SquaredInverse,
        RBFKernel_Wendland      ,
        RBFKernel_Cubed         ,
        RBFKernel_Identity      ,
        RBFKernel_SquaredLog    
    > Value;
};
template <> struct KernelList<HermiteRBF> {
    typedef mpl::list<
        RBFKernel_Gaussian      ,
        RBFKernel_SquaredInverse,
        RBFKernel_Cubed         
    > Value;
};
template <> struct KernelList<MLS> {
    typedef mpl::list<
        RBFKernel_Gaussian      ,
        RBFKernel_SquaredInverse,
        RBFKernel_Wendland      
    > Value;
};
template <> struct KernelList<GeneralizedMLS> {
    typedef mpl::list<
        RBFKernel_Gaussian      ,
        RBFKernel_SquaredInverse,
        RBFKernel_Wendland      
    > Value;
};

// utility for looping over parameters |
//-------------------------------------+
template <template <int, int, class, int> class FuncT>
struct TestLooper {
    static void go() {
        typedef mpl::list<
            mpl::int_<1>,
            mpl::int_<2>
        > DimList;
        mpl::for_each<DimList>(LoopDim());
    }
    struct LoopDim {
        template <class IntDim>
        void operator()(const IntDim&) const {
            mpl::for_each<KernelList<FuncT>::Value>(LoopKernel<IntDim::value>());
        }
        template <int Dim>
        struct LoopKernel {
            template <class RBFKernel>
            void operator()(const RBFKernel&) const {
                typedef mpl::list<
                    mpl::int_<0>,
                    mpl::int_<1>,
                    mpl::int_<2>
                > DegreeList;
                mpl::for_each<DegreeList>(LoopDegree<RBFKernel>());
                Tester<Dim, FuncT>::genCmd<RBFKernel>(2);
            }
            template <class RBFKernel>
            struct LoopDegree {
                template <class IntDegree>
                void operator()(const IntDegree&) const {
                    Tester<Dim, FuncT>::go<RBFKernel, IntDegree::value>(KernelUtil<RBFKernel>::getParam());
                }
            };
        };
    };
};

int main() {
    TestLooper<RBF           >::go();       // list of templates cannot be handled by Boost.MPL
    TestLooper<MLS           >::go();
    TestLooper<HermiteRBF    >::go();
    TestLooper<GeneralizedMLS>::go();
}

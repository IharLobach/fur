#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <stdio.h>
#include <stdlib.h> /* atoi */
#include <random>
#define BOOST_PYTHON_MAX_ARITY 17
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/numpy.hpp>
#include <boost/math/special_functions/sinc.hpp>
namespace p = boost::python;
namespace np = boost::python::numpy;

double N_periods = 10;
double lambda_wiggler_m = 0.055;
double alpha = 1.0 / 137.0;
int bessel_cutoff = 10;
std::vector<int> harmonics = {1,2,3};

template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U> &l, const std::pair<T, U> &r)
{
    return {l.first + r.first, l.second + r.second};
}

double besselj(int n, double x)
{
    double aux1 = 1.0;
    if (x < 0)
    {
        aux1 = pow(-1, n);
    }
    double aux2 = 1.0;
    if (n < 0)
    {
        aux2 = pow(-1, n);
    }
    return aux1 * aux2 * std::cyl_bessel_j(int(abs(n)), abs(x));
}

bool isinEllipse(double x, double y, double xmax, double ymax)
{
    return (pow(x/xmax,2)+pow(y/ymax,2)) < 1.0;
}

double myrand(double LO, double HI)
{
    return LO + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (HI - LO)));
}

std::pair<double, double> CalcAmpOneHarmonic(int harmonic, double x, double y, double lam, double gam, double K_peak)
{
    double aux_const = 1.0 + pow(K_peak, 2) / 2.0;
    double lambda1_um = 1e6 * lambda_wiggler_m / 2 / pow(gam, 2) * aux_const;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double r2 = gam * gam * (x * x + y * y);
    double A = aux_const + r2;
    double Y = harmonic * pow(K_peak, 2) / 4.0 / A;
    double X = 2 * harmonic * gam * K_peak * x / A;
    double jv2pm1 = besselj(harmonic + 2 * (-bessel_cutoff) - 1, X);
    for (int p = -bessel_cutoff; p < bessel_cutoff + 1; p++)
    {
        double jvpY = besselj(p, Y);
        sum1 += besselj(harmonic + 2 * p, X) * jvpY;
        sum2 += jv2pm1 * jvpY;
        double jv2pp1 = besselj(harmonic + 2 * p + 1, X);
        sum3 += jv2pp1 * jvpY;
        jv2pm1 = jv2pp1;
    }
    //std::cout << "sum1 = " << sum1 << ", sum2 = " << sum2
    //<< ", sum3 = " << sum3 << std::endl;
    double aux_factor = sqrt(alpha) * harmonic * gam * N_periods / A;
    double bessel_part_x = aux_factor * (2 * gam * x * sum1 - K_peak * (sum2 + sum3));
    //std::cout << "bessel_part_x = " << bessel_part_x << std::endl;
    double bessel_part_y = aux_factor * 2 * gam * y * sum1;
    double dw = lambda1_um / lam - harmonic;
    double L = boost::math::sinc_pi(M_PI * N_periods * (harmonic * r2 + dw * A) / aux_const) / sqrt(lam);

    return std::make_pair(bessel_part_x * L, bessel_part_y * L);
}

std::pair<double, double> CalcAmp(std::vector<int> harmonics, double x, double y, double lam, double gam, double K_peak)
{
    std::pair<double, double> res = std::make_pair(0.0, 0.0);
    for (int harmonic : harmonics)
    {
        res = res + CalcAmpOneHarmonic(harmonic, x, y, lam, gam, K_peak);
    }
    return res;
}

np::ndarray CalcSum(np::ndarray params, np::ndarray spectral_transmission, double lam_step)
{
    const double *st = reinterpret_cast<double *>(spectral_transmission.get_data());
    const double *prms = reinterpret_cast<double *>(params.get_data());
    double xmin = prms[0];
    double xmax = prms[1];
    double ymin = prms[2];
    double ymax = prms[3];
    double lmin = prms[4];
    double lmax = prms[5];
    int64_t m0 = (int64_t)prms[6];
    int mfold = (int)prms[7];
    int seed = (int)prms[8];
    double gam = prms[9];
    double K_peak = prms[10];

    double V = (lmax - lmin) * (xmax - xmin) * (ymax - ymin);
    std::cout << "Starting the loop" << std::endl;
    double s = 0;
    int64_t imax = m0 * pow(2, mfold) + 1;
    std::vector<double> res(mfold + 1);
    int res_idx = 0;
    for (int64_t i = 1; i < imax; i++)
    {
        double lam = myrand(lmin, lmax);
        double tr = st[int((lam-lmin)/lam_step)];
        double x = myrand(xmin, xmax);
        double y = myrand(ymin, ymax);
        double amp = CalcAmp(harmonics, x, y, lam, gam, K_peak).first;
        s = s * (double(i - 1) / i) + (1.0 / i) * pow(amp, 2)*tr;
        if (i == m0)
        {
            res[res_idx] = V * s;
            std::time_t t = std::time(nullptr);
            std::cout << std::put_time(std::localtime(&t), "%c %Z") << std::endl;
            std::cout << "n points = " << i << ", photons_per_electron = " << res[res_idx] << std::endl;
            m0 = 2 * m0;
            res_idx += 1;
        }
    }
    Py_intptr_t sh[1] = {res.size()};
    np::ndarray result = np::zeros(1, sh, np::dtype::get_builtin<double>());
    std::copy(res.begin(), res.end(), reinterpret_cast<double*>(result.get_data()));

    return result;
}

np::ndarray CalcM(np::ndarray params, np::ndarray spectral_transmission, double lam_step)
{
    const double *st = reinterpret_cast<double *>(spectral_transmission.get_data());
    std::complex<double> s, ar;
    const double *prms = reinterpret_cast<double*>(params.get_data());

    double Sx = prms[0];
    double Sy = prms[1];
    double dx = prms[2];
    double dy = prms[3];
    double sxp = prms[4];
    double syp = prms[5];
    double xmin = prms[6];
    double xmax = prms[7];
    double ymin = prms[8];
    double ymax = prms[9];
    double lmin = prms[10];
    double lmax = prms[11];
    double sz = prms[12];
    double tot = prms[13];
    int64_t m0 = (int64_t)prms[14];
    int mfold = (int)prms[15];
    int seed =  (int)prms[16];
    double gam = prms[17];
    double K_peak = prms[18];

    double k0, re, im, V, xp, yp, fm1, fn1, fm2, fn2, lam;
    double fm1Mixp, fn1Miyp, fm2Mixp, fn2Miyp;
    std::srand(seed);
    std::default_random_engine generator;
    std::normal_distribution<double> distx(0, sqrt(2.0)*sxp);
    std::normal_distribution<double> disty(0, sqrt(2.0)*syp);
    std::normal_distribution<double> distfm(0, 1.0 / sqrt(2.0)  / Sx);
    std::normal_distribution<double> distfn(0, 1.0 / sqrt(2.0) / Sy);


    std::cout << "tot = " << tot << std::endl;
    V = (lmax-lmin) * (xmax-xmin) * (ymax-ymin);
    std::cout << "Starting the loop" << std::endl;
    s = 0;
    int64_t imax = m0*pow(2, mfold)+1;
    ar.real(0.0);
    std::vector<std::complex<double>> res(mfold+1);
    int res_idx = 0;
    for (int64_t i=1;i<imax;i++)
    {
        lam = myrand(lmin, lmax);
        double tr = st[int((lam-lmin)/lam_step)];
        k0 = 2 * M_PI / lam;
        double fm1mfm2 = distfm(generator)/k0;
        fm2 = myrand(xmin, xmax);
        fm1 = fm2+fm1mfm2;
        double fn1mfn2 = distfn(generator)/k0;
        fn2 = myrand(ymin, ymax);
        fn1 = fn2+fn1mfn2;
        
        xp = distx(generator);
        yp = disty(generator);
        fm1Mixp = fm1 - xp;
        fn1Miyp = fn1 - yp;
        fm2Mixp = fm2 - xp;
        fn2Miyp = fn2 - yp;
        
        if ( isinEllipse(fm1,fn1, xmax, ymax)
                & isinEllipse(fm2, fn2, xmax, ymax)
                & isinEllipse(fm1Mixp, fn1Miyp, xmax, ymax)
                & isinEllipse(fm2Mixp, fn2Miyp, xmax, ymax)
        )
        {
            ar.imag(k0 * dx * (fm1 - fm2) * xp + k0 * dy * (fn1 - fn2) * yp);
            s = s * (double(i - 1) / i) + (1.0 / i) * (2 * M_PI /2 / k0 / k0 / Sx / Sy) * pow(lam, 2)
            * CalcAmp(harmonics, fm1Mixp, fn1Miyp, lam, gam, K_peak).first
            * CalcAmp(harmonics, fm1, fn1, lam, gam, K_peak).first
            * CalcAmp(harmonics, fm2Mixp, fn2Miyp, lam, gam, K_peak).first
            * CalcAmp(harmonics, fm2, fn2, lam, gam, K_peak).first
            * exp(-ar)
            * pow(tr, 2);
        }
        if (i == m0)
        {
            std::complex<double> M = pow(tot, 2) / (1.0 / 2.0 / sqrt(M_PI) / sz * V * s);
            res[res_idx] = M;
            std::time_t t = std::time(nullptr);
            std::cout << std::put_time(std::localtime(&t), "%c %Z") << std::endl;
            std::cout << "n points = " << i << ", M = " << M << std::endl;
            m0 = 2 * m0;
            res_idx+=1;
        }
    }
    Py_intptr_t sh[1] = {res.size()};
    np::ndarray result = np::zeros(1, sh, np::dtype::get_builtin<std::complex<double>>());
    std::copy(res.begin(), res.end(), reinterpret_cast<std::complex<double>*>(result.get_data()));

    return result;
}

int nx, ny, nz;

double ex3d(const double *Ex3d, int z, int y, int x)
{
    return Ex3d[z * (nx * ny) + y * (nx) + x];
};

bool isinRange(int val, int r)
{
    if ((val < 0) || (val > r - 1))
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool isin(int a, int b, int c, int d, int l)
{
    if (not isinRange(a, nx))
    {
        return false;
    }
    if (not isinRange(b, ny))
    {
        return false;
    }
    if (not isinRange(c, nx))
    {
        return false;
    }
    if (not isinRange(d, ny))
    {
        return false;
    }
    if (not isinRange(l, nz))
    {
        return false;
    }
    return true;
}

double CalcSumTensor(const double *Ex3d, int nx, int ny, int nz)
{
    double res = 0;
#pragma simd
    for (int i = 0; i < nx * ny * nz; i++)
    {
        res += pow(abs(Ex3d[i]), 2);
    }
    return res;
}

np::ndarray CalcMFromPrecalculatedFieldAmps(np::ndarray input_npy, np::ndarray spectral_transmission, np::ndarray params)
{
    std::complex<double> s, ar;
    const double *st = reinterpret_cast<double *>(spectral_transmission.get_data());
    const double *prms = reinterpret_cast<double *>(params.get_data());

    double Sx = prms[0];
    double Sy = prms[1];
    double dx = prms[2];
    double dy = prms[3];
    double sxp = prms[4];
    double syp = prms[5];
    double xmin = prms[6];
    double xmax = prms[7];
    double ymin = prms[8];
    double ymax = prms[9];
    double lmin = prms[10];
    double lmax = prms[11];
    double sz = prms[12];
    double sum = prms[13];
    int64_t m0 = (int64_t)prms[14];
    int mfold = (int)prms[15];
    int seed = (int)prms[16];

    double dax, day, dl, k0, re, im, V, xp, yp, fm1, fn1, fm2, fn2, lam, xplim, yplim;
    int ixp, iyp, m1, n1, m2, n2, il;
    int m1Mixp, n1Miyp, m2Mixp, n2Miyp;
    std::srand(seed);
    std::default_random_engine generator;
    std::normal_distribution<double> distx(0, sqrt(2.0) * sxp);
    std::normal_distribution<double> disty(0, sqrt(2.0) * syp);
    std::normal_distribution<double> distfm(0, 1.0 / sqrt(2.0) / Sx);
    std::normal_distribution<double> distfn(0, 1.0 / sqrt(2.0) / Sy);

    auto shape = input_npy.get_shape();
    nx = (int)(shape[2]);
    std::cout << "nx = " << nx << std::endl;
    ny = (int)(shape[1]);
    std::cout << "ny = " << ny << std::endl;
    nz = (int)(shape[0]);
    dax = (xmax - xmin) / (nx - 1);
    day = (ymax - ymin) / (ny - 1);
    dl = (lmax - lmin) / (nz - 1);
    std::cout << "nz = " << nz << std::endl;
    const double *Ex3d = reinterpret_cast<double *>(input_npy.get_data());
    double tot;
    if (sum == 0.0)
    {
        tot = dax*day*dl*CalcSumTensor(Ex3d, nx, ny, nz);
    }
    else
    {
        tot = sum;
    }
    std::cout << "tot = " << tot << std::endl;
    V = (lmax - lmin) * pow((xmax - xmin) * (ymax - ymin), 1);
    std::cout << "Starting the loop" << std::endl;
    s = 0;
    int64_t imax = m0 * pow(2, mfold) + 1;
    ar.real(0.0);
    std::vector<std::complex<double>> res(mfold + 1);
    int res_idx = 0;
    for (int64_t i = 1; i < imax; i++)
    {
        lam = myrand(lmin, lmax);
        k0 = 2 * M_PI / lam;
        double fm1mfm2 = distfm(generator) / k0;
        fm2 = myrand(xmin, xmax);
        fm1 = fm2 + fm1mfm2;
        double fn1mfn2 = distfn(generator) / k0;
        fn2 = myrand(ymin, ymax);
        fn1 = fn2 + fn1mfn2;
        m1 = int((fm1 - xmin) / dax);
        m2 = int((fm2 - xmin) / dax);
        if (m2 == nx)
        {
            m2 = nx - 1;
        }
        n1 = int((fn1 - ymin) / day);
        n2 = int((fn2 - ymin) / day);
        if (n2 == ny)
        {
            n2 = ny - 1;
        }
        il = int((lam - lmin) / dl);
        if (il == nz)
        {
            il = nz - 1;
        }
        xp = distx(generator);
        yp = disty(generator);
        m1Mixp = int((fm1 - xmin - xp) / dax + 0.5);
        n1Miyp = int((fn1 - ymin - yp) / day + 0.5);
        m2Mixp = int((fm2 - xmin - xp) / dax + 0.5);
        n2Miyp = int((fn2 - ymin - yp) / day + 0.5);

        s = s * (double(i - 1) / i);
        if (
            ((pow(fm1/xmax,2)+pow(fn1/ymax,2)) < 1)
            & ((pow(fm2/xmax,2)+pow(fn2/ymax,2)) < 1)
            & isin(m1, n1, m2, n2, il)
            & isin(m1Mixp, n1Miyp, m2Mixp, n2Miyp, il)
            )
        {
            ar.imag(k0 * dx * (fm1 - fm2) * xp + k0 * dy * (fn1 - fn2) * yp);
            s = s + (1.0 / i) * (2 * M_PI / 2 / k0 / k0 / Sx / Sy) * pow(lam, 2) * ex3d(Ex3d, il, n1Miyp, m1Mixp) * ex3d(Ex3d, il, n1, m1) * ex3d(Ex3d, il, n2Miyp, m2Mixp) * ex3d(Ex3d, il, n2, m2) * exp(-ar) * pow(st[il],2);
        }
        if (i == m0)
        {
            std::complex<double> M = pow(tot, 2) / (1.0 / 2.0 / sqrt(M_PI) / sz * V * s);
            res[res_idx] = M;
            std::time_t t = std::time(nullptr);
            std::cout << std::put_time(std::localtime(&t), "%c %Z") << std::endl;
            std::cout << "n points = " << i << ", M = " << M << std::endl;
            m0 = 2 * m0;
            res_idx += 1;
        }
    }
    Py_intptr_t sh[1] = {res.size()};
    np::ndarray result = np::zeros(1, sh, np::dtype::get_builtin<std::complex<double>>());
    std::copy(res.begin(), res.end(), reinterpret_cast<std::complex<double> *>(result.get_data()));

    return result;
}

BOOST_PYTHON_MODULE(coherent_modes_cpp)
{
    using namespace boost::python;
    Py_Initialize();
    np::initialize();
    def("CalcM", CalcM);
    def("CalcSum", CalcSum);
    def("CalcMFromPrecalculatedFieldAmps", CalcMFromPrecalculatedFieldAmps);
}
#include <cmath>
#include <iostream>
#include <boost/math/special_functions/sinc.hpp>

template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U> &l, const std::pair<T, U> &r)
{
    return {l.first + r.first, l.second + r.second};
}

double besselj(int n, double x){
    double aux1 = 1.0;
    if (x<0){
        aux1 = pow(-1, n);
    }
    double aux2 = 1.0;
    if (n<0){
        aux2 = pow(-1, n);
    }
    return aux1*aux2*std::cyl_bessel_j(int(abs(n)), abs(x));
}

double gam = 96.4 / 0.511;
double K_peak = 1.0;

double N_periods = 10;
double lambda_wiggler_m = 0.055;
double alpha = 1.0 / 137.0;
int bessel_cutoff = 10;




std::pair<double, double> CalcAmpOneHarmonic(int harmonic, double x, double y, double lam, double gam, double K_peak){
    double aux_const = 1.0 + pow(K_peak, 2) / 2.0;
    double lambda1_um = 1e6 * lambda_wiggler_m / 2 / pow(gam, 2) * aux_const;
    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double r2 = gam*gam*(x*x+y*y);
    double A = aux_const+r2;
    double Y = harmonic*pow(K_peak,2)/4.0/A;
    double X = 2*harmonic*gam*K_peak*x/A;
    double jv2pm1 = besselj(harmonic + 2 * (-bessel_cutoff) - 1, X);
    for (int p = -bessel_cutoff;p<bessel_cutoff+1;p++){
        double jvpY = besselj(p, Y);
        sum1 += besselj(harmonic+2*p, X)*jvpY;
        sum2 += jv2pm1*jvpY;
        double jv2pp1 = besselj(harmonic+2*p+1, X);
        sum3 += jv2pp1*jvpY;
        jv2pm1 = jv2pp1;
    }
    //std::cout << "sum1 = " << sum1 << ", sum2 = " << sum2 
    //<< ", sum3 = " << sum3 << std::endl;
    double aux_factor = sqrt(alpha)*harmonic*gam*N_periods/A;
    double bessel_part_x = aux_factor*(2*gam*x*sum1-K_peak*(sum2+sum3));
    //std::cout << "bessel_part_x = " << bessel_part_x << std::endl;
    double bessel_part_y = aux_factor*2*gam*y*sum1;
    double dw = lambda1_um/lam-harmonic;
    double L = boost::math::sinc_pi(M_PI*N_periods*(harmonic*r2+dw*A)/aux_const)/sqrt(lam);
    
    return std::make_pair(bessel_part_x*L, bessel_part_y*L);
}

std::pair<double, double> CalcAmp(std::vector<int> harmonics, double x, double y, double lam, double gam, double K_peak)
{
    std::pair<double,double> res = std::make_pair(0.0,0.0);
    for (int harmonic : harmonics){
        res = res + CalcAmpOneHarmonic(harmonic, x, y, lam, gam, K_peak);
    }
    return res;
}

    int main()
{
    std::vector<int> harmonics = {1};
    std::cout << "Result = "
              << CalcAmp(
                     harmonics,
                     0.0,
                     0.0,
                     1.075752508361204,
                     gam,
                     K_peak)
                     .first
              << std::endl;

    std::cout << besselj(3,-3) << std::endl;
    // std::cout << boost::math::sinc_pi(0.5) << std::endl;
    // std::cout << "lambda1_um = " << lambda1_um << std::endl;
}
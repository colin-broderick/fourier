#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "fourier.hpp"

Fourier::FourierFit::FourierFit(const Points &points_, const int order_)
: points(points_), order(order_)
{
}

/** \brief Gets the needed mathematical data for the lines that make the shape. */
void Fourier::FourierFit::create_line_segments()
{
    for (unsigned int i = 0; i < points.size() - 1; i++)
    {
        lines.push_back({points[i], points[i + 1]});
    }
    lines.push_back({points.back(), points.front()});
}

void Fourier::FourierFit::process()
{
    create_line_segments();
    
    p = lines.size();

    // Finds the coefficients for both coordinates
    double ax_0 = get_a0(0);
    double ay_0 = get_a0(1);
    std::vector<double> ax_coefficients;
    std::vector<double> bx_coefficients;
    std::vector<double> ay_coefficients;
    std::vector<double> by_coefficients;

    for (int i = 0; i < order; i++)
    {
        int n = i + 1;
        ax_coefficients.push_back(get_an(0, n));
        bx_coefficients.push_back(get_bn(0, n));
        ay_coefficients.push_back(get_an(1, n));
        by_coefficients.push_back(get_bn(1, n));
    }

    std::vector<double> timestamps = arange(0, p, p / 1000.0);

    // Get plot data
    for (double &timestamp : timestamps)
    {
        fx.push_back(end_fun(ax_0, ax_coefficients, bx_coefficients, p, timestamp));
        fy.push_back(end_fun(ay_0, ay_coefficients, by_coefficients, p, timestamp));
    }

    processed = true;
}

/** \brief Given initial coordinate I, final coordinate F, returns contribution to a0 coefficient. */
double Fourier::FourierFit::a0_segment(const double I, const double F)
{
    return (I + F) / 2.0;
}

/** \brief Given initial coordinate I, final coordinate F, line number L,and n, returns contribution to an coefficient. */
double Fourier::FourierFit::an_segment(const double I, const double F, const double L, const double n)
{
    double N = 2 * pi * n / p;
    double D = F - I;
    double term1 = ((I - L * D) / N) * (std::sin(N * (L + 1)) - std::sin(N * L));
    double term2 = (D / N / N) * (N * (L + 1) * std::sin(N * (L + 1)) - N * L * std::sin(N * L) + std::cos(N * (L + 1)) - std::cos(N * L));
    return term1 + term2;
}

/** \brief Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to bn coefficient. */
double Fourier::FourierFit::bn_segment(const double I, const double F, const double L, const double n)
{
    double N = 2 * pi * n / p;
    double D = F - I;
    double term1 = (-(I - L * D) / N) * (std::cos(N * (L + 1)) - std::cos(N * L));
    double term2 = (D / N / N) * (-N * (L + 1) * std::cos(N * (L + 1)) + N * L * std::cos(N * L) + std::sin(N * (L + 1)) - std::sin(N * L));
    return term1 + term2;
}

double Fourier::FourierFit::get_a0(const double c)
{
    double a = 0;
    for (int i = 0; i < p; i++)
    {
        a += a0_segment(lines[i][0][c], lines[i][1][c]);
    }
    double a0 = a / static_cast<double>(p);

    return a0;
}

double Fourier::FourierFit::get_an(const double c, const double n)
{
    double a = 0;
    for (int i = 0; i < p; i++)
    {
        a += an_segment(lines[i][0][c], lines[i][1][c], i, n);
    }
    double an = 2 * a / static_cast<double>(p);
    return an;
}

double Fourier::FourierFit::get_bn(const double c, const double n)
{
    double b = 0;
    for (int i = 0; i < p; i++)
    {
        b += bn_segment(lines[i][0][c], lines[i][1][c], i, n);
    }
    double bn = 2 * b / static_cast<double>(p);
    return bn;
}

double Fourier::FourierFit::sine_part(const double n, const double T, const double x)
{
    return std::sin(2 * pi * n * x / T);
}

double Fourier::FourierFit::cosine_part(const double n, const double T, const double x)
{
    return std::cos(2 * pi * n * x / T);
}

/** \brief Puts together the final Fourier series and returns y for given x. */
double Fourier::FourierFit::end_fun(const double a0, const std::vector<double> &a, const std::vector<double> &b, const double T, const double x)
{
    double y = a0;
    for (int i = 0; i < order; i++)
    {
        y += a[i] * cosine_part(i + 1, T, x);
        y += b[i] * sine_part(i + 1, T, x);
    }
    return y;
}

std::vector<double> Fourier::FourierFit::arange(const double start, const double end, const double step)
{
    std::vector<double> t;

    double s = start;
    while (s < end)
    {
        t.push_back(s);
        s += step;
    }

    return t;
}

/** \brief Saves the result of .process() to a CSV-formatted file.
 * 
 * If .process() has not been run, it will be run before saving.
 * 
 * \param filename The location to save the output.
 */
void Fourier::FourierFit::save(const std::string &filename)
{
    if (!processed)
    {
        process();
    }

    std::ofstream output{filename};

    output << std::setprecision(15);
    for (unsigned int i = 0; i < fx.size(); i++)
    {
        output << fx[i] << "," << fy[i] << "\n";
    }

    output.close();
}

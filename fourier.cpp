#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

class FourierFit
{
    using Points = std::vector<std::array<double, 2>>;
    using Lines = std::vector<std::array<std::array<double, 2>, 2>>;
public:
    FourierFit(const Points& pts) : points(pts)
    {
    }
    void process()
    {
        std::cout << "Processing started" << std::endl;

        std::array<int, 5> orders{2, 4, 6, 10, 15};

        // Gets highest order in 'orders'
        int high_ord = -1;
        for (auto &order : orders)
        {
            if (order > high_ord)
            {
                high_ord = order;
            }
        }

        // Gets the needed mathematical data for the lines that make the shape
        for (unsigned int i = 0; i < points.size() - 1; i++)
        {
            lines.push_back({points[i], points[i + 1]});
        }
        lines.push_back({points.back(), points.front()});

        p = lines.size();

        // Finds the coefficients for both coordinates
        double ax0 = get_a0(0);
        double ay0 = get_a0(1);
        std::vector<double> ax;
        std::vector<double> bx;
        std::vector<double> ay;
        std::vector<double> by;

        for (int i = 0; i < high_ord; i++)
        {
            int n = i + 1;
            ax.push_back(get_an(0, n));
            bx.push_back(get_bn(0, n));
            ay.push_back(get_an(1, n));
            by.push_back(get_bn(1, n));
        }

        // Prepare storage for plot data
        std::vector<double> fx0;
        std::vector<double> fy0;
        std::vector<double> fx1;
        std::vector<double> fy1;
        std::vector<double> fx2;
        std::vector<double> fy2;
        std::vector<double> fx3;
        std::vector<double> fy3;
        std::vector<double> fx4;
        std::vector<double> fy4;

        std::vector<double> timestamps = arange(0, p, p / 10000000.0);

        // Get plot data
        for (auto &timestamp : timestamps)
        {
            fx0.push_back(end_fun(ax0, ax, bx, p, timestamp, orders[0]));
            fy0.push_back(end_fun(ay0, ay, by, p, timestamp, orders[0]));
            fx1.push_back(end_fun(ax0, ax, bx, p, timestamp, orders[1]));
            fy1.push_back(end_fun(ay0, ay, by, p, timestamp, orders[1]));
            fx2.push_back(end_fun(ax0, ax, bx, p, timestamp, orders[2]));
            fy2.push_back(end_fun(ay0, ay, by, p, timestamp, orders[2]));
            fx3.push_back(end_fun(ax0, ax, bx, p, timestamp, orders[3]));
            fy3.push_back(end_fun(ay0, ay, by, p, timestamp, orders[3]));
            fx4.push_back(end_fun(ax0, ax, bx, p, timestamp, orders[4]));
            fy4.push_back(end_fun(ay0, ay, by, p, timestamp, orders[4]));
        }

        // Output data to log file.
        std::cout << "Writing output to file" << std::endl;
        std::ofstream output{"data.txt"};
        output << std::setprecision(15);
        for (unsigned int i = 0; i < fx4.size(); i++)
        {
            output << fx4[i] << " " << fy4[i] << "\n";
        }
        output.close();

        std::cout << "Processing complete, results in data.txt" << std::endl;
    }

private:
    // Given initial coordinate I, final coordinate F and line number L, returns contribution to a0 coefficient
    double a0_segment(const double line_stard_coord, const double line_end_coord)
    {
        return (line_stard_coord + line_end_coord) / 2.0;
    }

    // Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to an coefficient
    double an_segment(const double line_start_coord, const double line_end_coord, const double L, const double n)
    {
        double N = 2 * pi * n / p;
        double D = line_end_coord - line_start_coord;
        double term1 = ((line_start_coord - L * D) / N) * (std::sin(N * (L + 1)) - std::sin(N * L));
        double term2 = (D / N / N) * (N * (L + 1) * std::sin(N * (L + 1)) - N * L * std::sin(N * L) + std::cos(N * (L + 1)) - std::cos(N * L));
        return term1 + term2;
    }

    // Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to bn coefficient
    double bn_segment(const double line_start_coord, const double line_end_coord, const double L, const double n)
    {
        double N = 2 * pi * n / p;
        double D = line_end_coord - line_start_coord;
        double term1 = (-(line_start_coord - L * D) / N) * (std::cos(N * (L + 1)) - std::cos(N * L));
        double term2 = (D / N / N) * (-N * (L + 1) * std::cos(N * (L + 1)) + N * L * std::cos(N * L) + std::sin(N * (L + 1)) - std::sin(N * L));
        return term1 + term2;
    }

    // c defines which coordinate to get coefficient for: 0=x, 1=y
    double get_a0(const double c)
    {
        double a = 0;
        for (int i = 0; i < p; i++)
        {
            a += a0_segment(lines[i][0][c], lines[i][1][c]);
        }
        double a0 = a / static_cast<double>(p);

        return a0;
    }

    double get_an(const double c, const double n)
    {
        double a = 0;
        for (int i = 0; i < p; i++)
        {
            a += an_segment(lines[i][0][c], lines[i][1][c], i, n);
        }
        double an = 2 * a / static_cast<double>(p);
        return an;
    }

    double get_bn(const double c, const double n)
    {
        double b = 0;
        for (int i = 0; i < p; i++)
        {
            b += bn_segment(lines[i][0][c], lines[i][1][c], i, n);
        }
        double bn = 2 * b / static_cast<double>(p);
        return bn;
    }

    double sine_part(const double n, const double T, const double x)
    {
        return std::sin(2 * pi * n * x / T);
    }

    double cosine_part(const double n, const double T, const double x)
    {
        return std::cos(2 * pi * n * x / T);
    }

    // Puts together the final Fourier series and returns y for given x
    double end_fun(const double a0, const std::vector<double> &a, const std::vector<double> &b, const double T, const double x, const double o)
    {
        double y = a0;
        for (int i = 0; i < o; i++)
        {
            y += a[i] * cosine_part(i + 1, T, x);
            y += b[i] * sine_part(i + 1, T, x);
        }
        return y;
    }

    std::vector<double> arange(const double start, const double end, const double step)
    {
        std::vector<double> t;
        double s = start;
        double e = end;
        while (s < e)
        {
            t.push_back(s);
            s += step;
        }
        return t;
    }

    Points points;
    Lines lines;
    double pi = 3.14159;
    int p;
};

int main()
{
    FourierFit f{{
        {1, 1},
        {0, 0},
        {1, -1},
        {-1, -1},
        {0, 0},
        {-1, 1}
    }};
    f.process();
}

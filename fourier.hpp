#ifndef _CB_MP_FOURIER_HPP
#define _CB_MP_FOURIER_HPP

#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <cmath>

namespace Fourier
{
    using Points = std::vector<std::array<double, 2>>;
    using Lines = std::vector<std::array<std::array<double, 2>, 2>>;

    class FourierFit
    {
        public:
            FourierFit(const Points& pts, const int order);
            void process();
            void save(const std::string& filename);

        private:
            /** \brief Given initial coordinate I, final coordinate F and line number L, returns contribution to a0 coefficient. */
            double a0_segment(const double line_stard_coord, const double line_end_coord);

            /** \brief Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to an coefficient. */
            double an_segment(const double line_start_coord, const double line_end_coord, const double L, const double n);

            /** \brief Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to bn coefficient. */
            double bn_segment(const double line_start_coord, const double line_end_coord, const double L, const double n);

            /** \brief c defines which coordinate to get coefficient for: 0=x, 1=y */
            double get_a0(const double c);
            std::tuple<double, double> get_a0();


            double get_an(const double c, const double n);

            double get_bn(const double c, const double n);

            double sine_part(const double n, const double T, const double x);

            double cosine_part(const double n, const double T, const double x);

            void create_line_segments();

            /** \brief Puts together the final Fourier series and returns y for given x. */
            double end_fun(const double a0, const std::vector<double> &a, const std::vector<double> &b, const double T, const double x);

            std::vector<double> arange(const double start, const double end, const double step);

            Points points;
            std::vector<double> fx, fy;
            Lines lines;
            int order;
            const double pi = 3.14159;
            int p;
            bool processed = false;
    };
}

#endif

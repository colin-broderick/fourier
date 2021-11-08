#include "fourier.hpp"

int main()
{
    Fourier::Points pts = {{1, 1}, {0, 0}, {1, -1}, {-1, -1}, {0, 0}, {-1, 1}};
    Fourier::FourierFit fourier{pts, 100};

    fourier.process();
    fourier.save("data.txt");
}

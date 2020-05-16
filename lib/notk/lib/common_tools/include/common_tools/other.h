template <class T>
T round_mantissa(T a, unsigned n)
{
    if (a == 0)
    {
        return a;
    }
    int    x = std::floor(std::log10(std::abs(a)));
    double t = std::pow(10.0, (double)n - (double)x);
    return std::round(a * t) / t;
}
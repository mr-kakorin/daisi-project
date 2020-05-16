#include <array>
#include <cmath>

#include "GeneralTools.h"
#include "ellipse.h"

std::array<float, 5> est_ellipse(const std::vector<float>& x, const std::vector<float>& y)
{
    std::array<float, 5> result;

    int   n  = x.size();   // ���������� ����� � ���������
    float mx = sum(x) / n; // ���������� ������ ���� �� x
    float my = sum(y) / n; // ���������� ������ ���� �� y

    float Dx  = disp(x);
    float Dy  = disp(y);
    float Kxy = K(x, y);

    result[0] = 0.5 * atan(2 * Kxy / (Dx - Dy));
    result[1] = 4 * sqrt(std::abs(Dx * std::cos(result[0]) * std::cos(result[0]) + Kxy * std::sin(2 * result[0]) +
                         Dy * std::sin(result[0]) * std::sin(result[0])));
    result[2] = 4 * sqrt(std::abs(Dx * std::sin(result[0]) * std::sin(result[0]) - Kxy * std::sin(2 * result[0]) +
                         Dy * std::cos(result[0]) * std::cos(result[0])));

    result[3] = mx;
    result[4] = my;

    return result;
}

std::array<float, 5> calc_ellipse(const std::array<float, 5>& input, const float target, const std::vector<float>& x,
                                  const std::vector<float>& y)
{
    std::array<float, 5> result;

    float k1 = 0.01;
    float k2 = 1.0;
    // k[i] = (k1 + k2) / 2;
    // f[i] = inEllipse(k[i] * ra, k[i] * rb, alpha, x, y, mx, my);
    float k;
    float f;

    while (1)
    {
        k = (k1 + k2) / 2;
        f = inEllipse(k * input[1], k * input[2], input[0], x, y, input[3], input[4]);

        if (f < target)
            k1 = k;
        else
            k2 = k;

        if (std::abs(f - target) < 0.01 || std::abs(k1 - k2) < 1e-3)
            break;
    };

    result[3] = k;
    result[4] = f;

    result[0] = input[0];
    result[1] = k * input[1];
    result[2] = k * input[2];

    return result;
}

std::pair<std::vector<float>, std::vector<float>>
calc_alpha_beta(const std::vector<std::vector<float>>& x, const std::vector<std::vector<float>>& y, const float target)
{
    std::pair<std::vector<float>, std::vector<float>> result;

    for (size_t i = 0; i < x[0].size(); i++)
    {
        std::vector<float> x_cur(x.size());
        std::vector<float> y_cur(x.size());
        for (size_t j = 0; j < x.size(); j++)
        {
            x_cur[j] = x[j][i];
            y_cur[j] = y[j][i];
        }
        auto al_a_b  = est_ellipse(x_cur, y_cur);
        auto ell_par = calc_ellipse(al_a_b, target, x_cur, y_cur);

        auto phi     = ell_par[0];
        auto ra      = ell_par[1];
        auto rb      = ell_par[2];
        auto rb2     = rb * rb;
        auto ra2     = ra * ra;
        auto sinphi2 = std::sin(phi) * std::sin(phi);

        auto alpha = std::tan(2.0 * phi) * (-2 * rb2 * sinphi2 + 2 * ra2 * sinphi2 + rb2 - ra2) / (2 * ra * rb);
        auto beta  = (rb2 * sinphi2 - ra2 * sinphi2 + ra2) / (ra * rb);

        result.first.push_back(alpha);
        result.second.push_back(beta);
    }

    return result;
}

std::pair<std::vector<float>, std::vector<float>> calc_cm(const std::vector<std::vector<float>>& x,
                                                          const std::vector<std::vector<float>>& y)
{
    std::pair<std::vector<float>, std::vector<float>> result;

    for (size_t i = 0; i < x[0].size(); i++)
    {
        std::vector<float> x_cur(x.size());
        std::vector<float> y_cur(x.size());
        for (size_t j = 0; j < x.size(); j++)
        {
            x_cur[j] = x[j][i];
            y_cur[j] = y[j][i];
        }
        auto al_a_b = est_ellipse(x_cur, y_cur);

        result.first.push_back(al_a_b[3]);
        result.second.push_back(al_a_b[4]);
    }

    return result;
}

void calculateEllipsesParameters(std::vector<std::string>& out, std::vector<std::vector<float>>& x1El,
                                 std::vector<std::vector<float>>& x2El, std::vector<std::vector<float>>& yEl1,
                                 std::vector<std::vector<float>>& yEl2, const std::vector<float>& x,
                                 const std::vector<float>& y, const std::vector<float>& targets)
{
    int nEll = targets.size();
    out.resize(nEll);
    x1El.resize(nEll);
    x2El.resize(nEll);
    yEl1.resize(nEll);
    yEl2.resize(nEll);

    auto al_a_b = est_ellipse(x, y);

    float alpha = al_a_b[0];
    float ra    = al_a_b[1];
    float rb    = al_a_b[2];
    float mx    = al_a_b[3];
    float my    = al_a_b[4];

    std::vector<float> k(nEll);
    std::vector<float> f(nEll);

    for (int i = 0; i < nEll; i++)
    {
        auto ell_par = calc_ellipse(al_a_b, targets[i], x, y);

        f[i]            = roundf(ell_par[4] * 1000) / 1000;
        k[i]            = ell_par[3];
        std::string ss1 = std::to_string(f[i] * 100);
        std::string s1  = std::string(ss1.begin(), ss1.begin() + 4);

        float em = k[i] * k[i] * ra * rb / 1e-6;
        em       = roundf(em * 10) / 10;

        std::string ss2 = std::to_string(em);
        std::string s2  = std::string(ss2.begin(), ss2.begin() + 4);
        plotEllipse(k[i] * ra, k[i] * rb, alpha, x1El[i], x2El[i], yEl1[i], yEl2[i], mx, my);

        float xmf = vectorMax(x1El[i]);
        xmf       = roundf(xmf * 1000) / 1000;
        ;

        std::string xm   = std::to_string(xmf);
        std::string xms2 = std::string(xm.begin(), xm.begin() + 5);

        float ymf = vectorMax(yEl1[i]);
        ymf       = roundf(ymf * 1000) / 1000;

        std::string ym   = std::to_string(ymf);
        std::string yms2 = std::string(ym.begin(), ym.begin() + 6);

        if (nEll == 1)
            out[i] = s2 + " pi*mm*mrad, Xmax= " + xms2 + ", Ymax= " + yms2 + " Inside " + s1 + "% ";
        else
            out[i] = s2 + " pi*mm*mrad, Inside " + s1 + "% ";
    }
};
void plotEllipse(float ra, float rb, float alpha, std::vector<float>& x1, std::vector<float>& x2,
                 std::vector<float>& y1, std::vector<float>& y2, float mx, float my)
{
    int   n  = 500;
    float dx = 2 * ra / n;
    x1.reserve(n);
    x2.reserve(n);

    y1.reserve(n);
    y2.reserve(n);
    float xCur = -ra;
    float dxC;
    while (1)
    {
        if (std::abs(ra - std::abs(xCur)) < 0.1 * dx)
            dxC = 0.04 * dx;
        else
            dxC = dx;

        float x1tmp = rb * sqrt(1 - xCur * xCur / (ra * ra));
        float x2tmp = -rb * sqrt(1 - xCur * xCur / (ra * ra));

        x1.push_back(xCur * std::cos(alpha) - x1tmp * std::sin(alpha) + mx);
        x2.push_back(xCur * std::cos(alpha) - x2tmp * std::sin(alpha) + mx);

        y1.push_back(xCur * std::sin(alpha) + x1tmp * std::cos(alpha) + my);
        y2.push_back(xCur * std::sin(alpha) + x2tmp * std::cos(alpha) + my);

        if (xCur == ra)
            break;

        xCur = xCur + dxC;

        if (xCur > ra)
            xCur = ra;
    };
};
float inEllipse(float ra, float rb, float alpha, const std::vector<float>& x, const std::vector<float>& y, float mx,
                float my)
{
    float inEl = 0;
    for (int i = 0; i < x.size(); i++)
    {
        float xtmp1 = x[i] - mx;
        float ytmp1 = y[i] - my;
        float xtmp  = xtmp1 * std::cos(alpha) + ytmp1 * std::sin(alpha);
        float ytmp  = -xtmp1 * std::sin(alpha) + ytmp1 * std::cos(alpha);
        if ((xtmp * xtmp / (ra * ra) + ytmp * ytmp / (rb * rb)) < 1)
            inEl++;
    };
    inEl = inEl / float(x.size());
    return inEl;
};
float sum(const std::vector<float>& v)
{
    float res = 0;
    for (int i = 0; i < v.size(); i++)
        res = res + v[i];

    return res;
};
float disp(const std::vector<float>& v)
{
    int   n   = v.size();
    float m   = sum(v) / n;
    float res = 0;

    for (int i = 0; i < v.size(); i++)
        res = res + (v[i] - m) * (v[i] - m);

    res = res / n;
    return res;
};
float K(const std::vector<float>& vx, const std::vector<float>& vy)
{
    float Kxy = 0;
    int   n   = vx.size();
    float mx  = sum(vx) / n;
    float my  = sum(vy) / n;

    for (int i = 0; i < vx.size(); i++)
        Kxy = Kxy + (vx[i] - mx) * (vy[i] - my);

    Kxy = Kxy / n;
    return Kxy;
}
#include "../base/AccelFlow.h"
#include "RFQSolver.h"

#include <common_tools/constants.h>

#include <boost/math/special_functions/bessel.hpp>

double f2(double m, double k, double r, double t)
{
    return t -
           (commtools::PI() / 4) * (m * m - 1) /
               (m * m * boost::math::cyl_bessel_i(0, 2 * k * r / (m + 1)) +
                boost::math::cyl_bessel_i(0, 2 * m * k * r / (m + 1)));
};
int GenerateRFQForFlowPr(int& succes, const std::vector<double>& RFQParameters,
                         std::shared_ptr<AccelFlow>& flow, const std::vector<double>& Modulations,
                         const std::vector<double>& SyncPhases,
                         const std::vector<double>& MatcherRadii, std::vector<double>& L,
                         std::vector<double>&       MinumalRadii,
                         std::vector<double>&       MinumalRadiiRegular,
                         const std::vector<double>& OutRadii, std::vector<double>& AccEff,
                         std::vector<double>& AvEnergies)
{

    /*0 "Frequency (Hz)"
    1 "Channel Radius",
    2 "Number of buncher cells",
    3 "Buncher phase",
    4 "Buncher modulation",
    5 "Number of former cells",
    6 "Former phase",
    7 "Former modulation",
    8 "Number of accel cells",
    9 "Matcher max radius",
    10 "Number of matcher cells"*/
    int nCellMatcher = MatcherRadii.size();
    int nCellRegular = Modulations.size() - nCellMatcher;
    int nCellOut     = OutRadii.size();

    succes        = 1;
    double lambda = commtools::LIGHT_VELOCITY() / RFQParameters[0];

    double Bo = flow->GetStartVelocity();
    double Lo = Bo * lambda;

    double eUTmax = flow->GetVoltage();
    double Wo     = flow->GetRestEnergy();

    std::vector<double> velocities;
    std::vector<double> dz(200);
    L.clear();
    AccEff.clear();
    MinumalRadii.clear();
    MinumalRadiiRegular.clear();
    L.resize(nCellRegular + nCellMatcher);
    MinumalRadii.resize(nCellRegular + nCellMatcher);

    MinumalRadiiRegular.resize(nCellRegular);
    AccEff.resize(nCellRegular + nCellMatcher);

    AvEnergies.clear();
    double T = 0;

    double periodMomentum = flow->GetStartMomentum();
    double periodVelocity = flow->GetStartVelocity();
    double z              = 0;
    double wfreq          = 2 * commtools::PI() / lambda;
    double alpha          = flow->GetAlpha();
    double startMomentum  = periodMomentum;
    double startZ         = 0;

    for (int period = 0; period < nCellRegular + nCellMatcher; period++)
    {
        velocities.push_back(periodVelocity);
        AvEnergies.push_back(flow->GetMass() * commtools::LIGHT_VELOCITY() * commtools::LIGHT_VELOCITY() *
                             periodVelocity * periodVelocity / (2 * std::abs(commtools::ELECTRON_CHARGE())));

        double phi1 = SyncPhases[period];
        double phi2 = SyncPhases[period + 1];

        double T;

        T = lambda * (phi2 - phi1 + commtools::PI()) / (2 * commtools::PI());

        L[period] = periodVelocity * T;

        double t = 0;
        double kk;
        if (period % 2 == 0)
            kk = commtools::PI();
        else
            kk = 0;

        double dt = T / 400;
        double gamma;

        double dzC;
        double Teff;
        for (int iter = 0; iter < 2; iter++)
        {
            double tt      = 0;
            z              = startZ;
            int s          = 0;
            periodMomentum = startMomentum;

            Teff = std::abs(f2(Modulations[period], commtools::PI() / L[period], RFQParameters[1], 0));

            while (1)
            {
                double Ez = (2 * Teff * eUTmax / L[period]) * std::cos(commtools::PI() * z / L[period] + kk) *
                            std::cos(wfreq * tt + phi1 + kk);
                periodMomentum = periodMomentum + dt * alpha * Ez;
                gamma          = sqrt(1 + periodMomentum * periodMomentum);
                //	dz[s] = dt * periodMomentum / gamma;
                dzC = dt * periodMomentum / gamma;
                z   = z + dzC;
                tt  = tt + dt;
                s++;
                if (z > L[period])
                    break;
            }
            tt        = tt - dt * (z - L[period]) / dzC;
            double dL = (periodMomentum / gamma) * (T - tt);
            L[period] = L[period] + dL;
        }
        AccEff[period] = Teff;
        startZ         = z - L[period];
        startZ         = 0;
        startMomentum  = periodMomentum;
        periodVelocity = periodMomentum / gamma;
    }

    /*std::vector<double> modulations;
    //double = lambda*lambda / (commtools::PI()*commtools::PI())*.5*ke / km / (Wo*NA);

    double kMatcher = (RFQParameters[9] - RFQParameters[1]) / (nCellMatcher *nCellMatcher);

    for (int i = 0; i < nCellMatcher; i++)
    {
    MatcherRadii[i] = RFQParameters[1] + (nCellMatcher - i)*(nCellMatcher - i)*kMatcher;
    AllRadii[i] = MatcherRadii[i];
    }*/

    /*for (int period = RFQParameters[5]; period < RFQParameters[2] + RFQParameters[5]; period++)
    {
    double k = commtools::PI() / L[period];
    double st = u1c[period] * RFQParameters[6];

    if (st == 0)
    {
    modulations.push_back(1);
    continue;
    }

    double qm0;
    if (period == RFQParameters[5])
    qm0 = 1;
    else
    qm0 = modulations[period - RFQParameters[5] - 1];

    modulations.push_back(Dmath::zerosSearchRFQ(qm0, 1e-6, 1e-2, f2, k, RFQParameters[1], st));
    if (modulations.back() == -1)
    {
    succes = 0;
    return period - RFQParameters[5];
    }
    }*/

    for (int period          = 0; period < nCellMatcher; period++)
        MinumalRadii[period] = MatcherRadii[period];

    for (int period = nCellMatcher; period < nCellRegular + nCellMatcher - nCellOut; period++)
    {

        MinumalRadii[period] = (2 / (Modulations[period] + 1)) * RFQParameters[1];
        MinumalRadiiRegular[period - nCellMatcher] = MinumalRadii[period];
        //	AllRadii[period + nCellMatcher] = MinumalRadii[period];
    };

    for (int period = nCellRegular + nCellMatcher - nCellOut; period < nCellRegular + nCellMatcher;
         period++)
    {
        MinumalRadii[period] = OutRadii[period - (nCellRegular + nCellMatcher - nCellOut)];
    }

    FILE* fp = fopen("u1.txt", "w");

    fprintf(fp, "Acceleration efficiency \n");
    FILE* fp1 = fopen("u2.txt", "w");
    fprintf(fp1, "Sync phases \n");

    std::vector<double> LL(MinumalRadii.size() + 1);

    LL[0] = 0;
    std::vector<std::vector<float>> OutParams(7);

    for (int i = 0; i < MinumalRadii.size(); i++)
    {
        LL[i + 1] = LL[i] + L[i];
    };
    for (int i = 0; i < MinumalRadii.size(); i = i + 10)
    {
        fprintf(fp, "%lf\t%lf \n", LL[i], 2 * AccEff[i] / L[i]);
        fprintf(fp1, "%lf\t%lf \n", LL[i], SyncPhases[i]);
    };
    fclose(fp);
    fclose(fp1);

    return -1;
};

//
// Created by artoria on 5/3/20.
//
#include "EnergyDistribution.h"
#include <cmath>
#include <algorithm>
#include <random>

namespace
{
	constexpr double x[] = {-12.668877,-12.612855,-12.556834,-12.494587,-12.413667,-12.376319,-12.264276,-12.226928,-12.164682,-12.108660,-12.058863,-12.009066,-11.859675,-11.828552,-11.816103,-11.685386,-11.654263,-11.604465,-11.542219,-11.504871,-11.467524,-11.399053,-11.355480,-11.293234,-11.212314,-11.143843,-11.081597,-11.013126,-10.957104,-10.913532,-10.863735,-10.820162,-10.770365,-10.733018,-10.689445,-10.627199,-10.564953,-10.521380,-10.484032,-10.421786,-10.396888,-10.334641,-10.297294,-10.241272,-10.216373,-10.166576,-10.116779,-10.098106,-10.060758,-9.992287,-9.973613,-9.923816,-9.905142,-9.855345,-9.830447,-9.793099,-9.730852,-9.699729,-9.587686,-9.550338,-9.519215,-9.481867,-9.432070,-9.376049,-9.338701,-9.282679,-9.245332,-9.195535,-9.114614,-9.033694,-8.983897,-8.940325,-8.865629,-8.809608,-8.703789,-8.647767,-8.604195,-8.529499,-8.504601,-8.473478,-8.361434,-8.398782,-8.261840,-8.218268,-8.174696,-8.124899,-8.093775,-8.081326,-8.075102,-8.062652,-8.025305,-7.994181,-7.944384,-7.863464,-7.851015,-7.838566,-7.801218,-7.782544,-7.770095,-7.763870,-7.763870,-7.738972,-7.682950,-7.658051,-7.626928,-7.583356,-7.558457,-7.558457,-7.546008,-7.533559,-7.508660,-7.502436,-7.489986,-7.458863,-7.421516,-7.390392,-7.353045,-7.315697,-7.272125,-7.216103,-7.166306,-7.141407,-7.104059,-7.054263,-7.023139,-6.979567,-6.942219,-6.904871,-6.879973,-6.836400,-6.792828,-6.743031,-6.711908,-6.630988,-6.574966,-6.537618,-6.494046,-6.462923,-6.438024,-6.382003,-6.319756,-6.238836,-6.120568,-6.089445,-6.002300,-5.921380,-5.778214,-5.547903,-5.336265,-5.062382,-4.944114,-4.800947,-4.620433,-4.495940,-4.302977,-4.060217,-3.798782,-3.618268,-3.493775};
	constexpr double y[] = {0.30566038,0.21886792,0.23396226,0.24528302,0.26037736,0.27547170,0.27924528,0.28679245,0.30188679,0.31320755,0.32452830,0.34339623,0.32830189,0.32830189,0.32830189,0.34716981,0.35094340,0.35849057,0.37735849,0.38490566,0.38490566,0.41509434,0.42264151,0.42264151,0.41886792,0.41886792,0.41132075,0.41132075,0.40377358,0.39245283,0.36226415,0.35471698,0.35471698,0.34339623,0.33584906,0.32830189,0.32075472,0.32075472,0.30566038,0.29056604,0.29056604,0.29056604,0.28301887,0.28301887,0.28301887,0.27547170,0.27169811,0.26037736,0.20754717,0.19622642,0.18490566,0.14716981,0.13962264,0.13584906,0.13584906,0.12452830,0.12452830,0.10943396,0.09056604,0.08679245,0.07924528,0.06792453,0.06037736,0.06037736,0.06037736,0.04528302,0.04905660,0.05660377,0.06415094,0.07924528,0.09433962,0.09811321,0.11320755,0.12452830,0.14716981,0.18867925,0.21132075,0.24528302,0.25660377,0.27924528,0.30566038,0.29433962,0.36226415,0.38113208,0.40377358,0.43396226,0.46415094,0.49056604,0.50566038,0.54339623,0.56603774,0.56603774,0.58113208,0.58490566,0.58867925,0.59622642,0.61509434,0.62641509,0.64528302,0.67924528,0.66415094,0.69056604,0.72830189,0.76603774,0.78490566,0.79245283,0.81509434,0.80000000,0.82641509,0.84150943,0.86792453,0.88679245,0.90188679,0.94716981,0.95849057,0.97358491,0.98490566,0.98867925,0.98867925,0.98113208,0.95094340,0.92830189,0.91320755,0.88679245,0.87924528,0.84905660,0.82264151,0.80000000,0.77735849,0.73584906,0.66792453,0.63018868,0.59245283,0.54716981,0.49433962,0.47169811,0.42641509,0.39622642,0.35471698,0.30943396,0.28301887,0.21886792,0.17735849,0.15849057,0.14716981,0.13207547,0.11698113,0.12452830,0.12452830,0.12452830,0.12452830,0.12452830,0.12452830,0.12830189,0.12830189,0.12830189,0.12830189,0.12830189,0.12452830};
	constexpr double I[] = {0.005209883717460731, 0.009707544700980154, 0.01499658365796888, 0.02225121415053546, 0.02579943235307489, 0.03681880030086659, 0.04056691805796445, 0.04706358593186312, 0.05317301794281135, 0.05880348610606657, 0.06470048531290468, 0.08249143220719533, 0.08611457397401424, 0.08756380739801355, 0.1032183263434574, 0.1070705172783299, 0.1133341224944953, 0.1214549574472289, 0.1265024227393641, 0.1315997264517355, 0.1413114688391562, 0.1477832684164904, 0.1571118172547104, 0.1691848216112322, 0.179354664944719, 0.1885166324393835, 0.1985032353458855, 0.2065991696892976, 0.2127501532982532, 0.2194134292112478, 0.2249523567551635, 0.2312158361027839, 0.2358383909607567, 0.2410857960714826, 0.2484153702217183, 0.2555783631387586, 0.260534245739041, 0.2646821628039057, 0.2712621214047621, 0.2738274311387214, 0.2802409115391458, 0.2840389025729786, 0.2896610793079504, 0.2921598579756891, 0.2970906821936837, 0.3019215572261904, 0.3036830782494164, 0.3067815223052121, 0.3116832036663673, 0.3129450700059552, 0.3158769113970782, 0.3168264345626511, 0.3192585302912004, 0.3204578959017419, 0.3221820300978179, 0.3249306644920786, 0.3262216689941627, 0.3301946383984133, 0.3313690486373581, 0.3322852453647525, 0.3332597559829223, 0.3343925129181494, 0.3355918892333937, 0.3363914877164937, 0.3374409607255624, 0.3380656302812786, 0.3389984888455302, 0.3407309592185785, 0.3427882423898066, 0.3443207958021196, 0.3458075264455143, 0.3486061211363642, 0.3509673932073398, 0.3560648172924344, 0.3594006421644547, 0.3624907097792456, 0.3685376731090391, 0.3707531678511856, 0.373709984630043, 0.3853291497896145, 0.3813561449260364, 0.3972980615881599, 0.4030409231263615, 0.4091044520279417, 0.4165006882225011, 0.421456638976972, 0.4235638576680915, 0.4246631878085531, 0.426978815701191, 0.4343249299333276, 0.440571926780968, 0.4507001062605397, 0.4674290675924274, 0.4700193640989037, 0.4726346473999689, 0.4806556195180237, 0.4847660553543304, 0.4875729039628606, 0.4890347481792889, 0.4890347481792889, 0.4950149182217925, 0.5091078411884545, 0.5157046168553926, 0.524262727517674, 0.5364480884848355, 0.5435446198500052, 0.5435446198500052, 0.5471343876544563, 0.5508157736542321, 0.5583620851654684, 0.5602984052884513, 0.5642466341062715, 0.574449734363587, 0.5870680597598582, 0.5977296010141244, 0.6106977413801369, 0.6237661787903269, 0.639041607449531, 0.6586067824355832, 0.6756647688717327, 0.6839607140186065, 0.6961545906210143, 0.7120461845633992, 0.7217914995419766, 0.7351429237119022, 0.7462123649143424, 0.7569569692661313, 0.7639199527729794, 0.7756100051713682, 0.7864543934437629, 0.7979152280529606, 0.8046617678318325, 0.8210117558498048, 0.8313565609269356, 0.8377533485930851, 0.8446914248774258, 0.8492307633716392, 0.8525458097298541, 0.8591423792650511, 0.8656807975423226, 0.8728812889285578, 0.8811895697694159, 0.8830427859611243, 0.8877653963604949, 0.891771684777225, 0.8980934511545403, 0.9079550940064208, 0.9173003714644921, 0.9293941950030293, 0.9346165429481582, 0.9409383534823709, 0.948909291664863, 0.9544898076924898, 0.9632686568798799, 0.9743130199377246, 0.9862070018788498, 0.9944194839723731, 0.9999999999999999};
	constexpr double Q = 2.820132465156285;
	constexpr unsigned int spectre_length = 159;
	constexpr double a = 0.045283019;
	constexpr double b = 0.98867925;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(a,b);

	constexpr double trapz( unsigned int max_idx ) noexcept
	{
		double sum = 0;
		for(unsigned int i = 0; i < max_idx; ++i)
			sum += (x[i+1]-x[i])*(y[i+1]+y[i]);
		return sum/2;
	}

	constexpr bool internal_pdf( double alpha, double eps, int& idx ) noexcept
	{
		for ( unsigned int i = 0; i < spectre_length; ++i)
		{
			const bool cond = std::abs( I[i] - alpha ) > eps;
			if ( i + 1 == spectre_length && cond ) {
				idx = -1;
				return false;
			}
			if ( !cond ) {
				idx = i + 1;
				return true;
			}
		}
	}
}

double get_initial_energy_bimodal() noexcept
{
	double alpha = distribution(generator);
	double eps = 1e-3;
	int idx = -1;
	bool result = internal_pdf( alpha, eps, idx );
	while ( !result )
	{
		eps *= 10;
		result = internal_pdf( alpha, eps, idx );
	}
	return x[idx];
}
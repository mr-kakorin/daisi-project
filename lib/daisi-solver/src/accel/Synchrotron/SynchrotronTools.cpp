#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "../base/AccelFlow.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
#include "Tools.h"
#include <algorithm>

#include <common_tools/constants.h>

void save_errors(const std::vector<std::pair<std::string, std::array<double, 8>>>& errors,
                 const std::string& fileName)
{
    boost::property_tree::ptree pt;

    for (const auto& val : errors)
    {
        boost::property_tree::ptree pt_loc;
        pt_loc.put("z", val.second[0]);
        pt_loc.put("x", val.second[1]);
        pt_loc.put("y", val.second[2]);
        pt_loc.put("xy", val.second[3]);
        pt_loc.put("xz", val.second[4]);
        pt_loc.put("yz", val.second[5]);
		pt_loc.put("sigma", val.second[6]);
		pt_loc.put("L", val.second[7]);

        pt.put_child(val.first, pt_loc);
    }
    boost::property_tree::write_json(fileName, pt);
}

std::shared_ptr<std::vector<std::pair<std::string, std::array<double, 8>>>>
read_errors(const boost::property_tree::ptree& pt)
{
    auto result = std::make_shared<std::vector<std::pair<std::string, std::array<double, 8>>>>();
    auto read_errors = [&](const boost::property_tree::ptree& pt) -> bool {
        bool correct = true;

        for (auto& item : pt.get_child(""))
        {
            std::array<double, 8> tols;
            tols[0] = item.second.get<double>("z", 0);
            tols[1] = item.second.get<double>("x", 0);
            tols[2] = item.second.get<double>("y", 0);
            tols[3] = item.second.get<double>("xy", 0);
            tols[4] = item.second.get<double>("xz", 0);
            tols[5] = item.second.get<double>("yz", 0);
			tols[6] = item.second.get<double>("sigma", 0);
			tols[7] = item.second.get<double>("L", 0);

            result->emplace_back(item.first, tols);
        }

        return correct;
    };
    if (read_property_tree(pt, "Error read tolerances file", read_errors))
    {
        return result;
    }
    return nullptr;
}

OpticElement::OpticElement(){

};

std::vector<double> OpticElement::GetErrors() const
{
    return alignment_errors;
}

std::pair<std::vector<double>, std::vector<double>>
OpticElement::GetErrorsBorders(double RelError) const
{
    std::pair<std::vector<double>, std::vector<double>> result;

    if (type == "RBEND" || type == "QUADRUPOLE" || type == "SBEND")
    {
        for (size_t i = 0; i < 8; i++)
        {
            result.first.push_back(-RelError / 2.0);
            result.second.push_back(RelError / 2.0);
        }
    }

    return result;
}

std::pair<std::vector<double>, std::vector<double>>
OpticElement::GetErrorsBorders(const std::map<std::string, std::array<double, 8>>& errors) const
{

    std::pair<std::vector<double>, std::vector<double>> result;

    auto fillup = [&](const std::string& key) {
        for (size_t i = 0; i < 8; i++)
        {
            result.first.push_back(-errors.find(key)->second[i] / 2.0);
            result.second.push_back(errors.find(key)->second[i] / 2.0);
        }
    };

    if (type == "QUADRUPOLE")
    {
        fillup("quadrupole");
    }
    else if (type == "RBEND")
    {
        fillup("rbend");
    }
    else if (type == "SBEND")
    {
        fillup("sbend");
    }

    return result;
}

std::pair<size_t, size_t>
OpticElementsSequence::GetErrIndexes(const std::vector<int>& monitiors) const
{
    std::pair<size_t, size_t> result;
    result.first  = 0;
    result.second = 0;

    for (int i = 0; i < length(); i++)
    {
        if (*monitiors.begin() == i)
        {
            break;
        }
        if (elements[i].GetErrors().size() != 0)
        {
            result.first++;
        }
    }

    for (int i = 0; i < length(); i++)
    {
        if (monitiors.back() == i)
        {
            break;
        }
        if (elements[i].GetErrors().size() != 0)
        {
            result.second++;
        }
    }

    return result;
}
std::vector<double> OpticElementsSequence::GetErrors() const
{
    std::vector<double> result;

    for (int i = 0; i < length(); i++)
    {
        if (elements[i].GetErrors().size() != 0)
        {
            auto tmp = elements[i].GetErrors();
            result.insert(result.end(), tmp.begin(), tmp.end());
        }
    }

    return result;
}
void OpticElementsSequence::insert_shuffle(const std::string& type, 
	const std::vector<std::string>& shuffled, const std::map<std::string, double>& actual_vals)
{
	auto l_old = GetL();
	int j = 0;

	double LL = 0;

	for (auto v : actual_vals)
	{
		LL += v.second;
	}
	LL /= double(actual_vals.size());

	for (int i = 0; i < length(); i++)
	{
		if (GetType(i) == type)
		{
			SetLabel(i, shuffled[j]);

			double L_new = getParameter(i, "L");
				
			auto it = actual_vals.find(shuffled[j]);

			if (it != actual_vals.end())
			{
				L_new = it->second;
			}


			auto L_old = getParameter(i, "L");
			auto dL2 = (L_old - L_new) / 2.0;

			setParameter(i, "L", L_new);

			if(i>0)
				setParameter(i - 1, "L", getParameter(i - 1, "L") + dL2);

			if(i<length()-1)
				setParameter(i + 1, "L", getParameter(i + 1, "L") + dL2);

			auto cur_err = elements[i].GetErrors();

			cur_err[6] = L_new / LL - 1.0;

			elements[i].InsertErrors(cur_err);


			j++;
		}
	}
	auto l_new = GetL();

}

std::vector<std::pair<std::string, std::array<double, 8>>>
OpticElementsSequence::GetErrorsStruct() const
{
    std::vector<std::pair<std::string, std::array<double, 8>>> result;

    for (int i = 0; i < length(); i++)
    {
        if (elements[i].GetErrors().size() != 0)
        {
            auto tmp_v = elements[i].GetErrors();

            std::array<double, 8> tmp;

            for (size_t ind = 0; ind < tmp_v.size(); ind++)
            {
                tmp[ind] = tmp_v[ind];
            }

            result.emplace_back(elements[i].GetLabel(), tmp);
        }
    }

    return result;
}
std::vector<std::string>  OpticElementsSequence::get_types_array(const std::string& type)
{
	std::vector<std::string> result;
	for (int i = 0; i < length(); i++)
	{
		if (GetType(i) == type)
		{
			result.push_back(GetLabel(i));
		}
	}
	return result;
}

void OpticElementsSequence::mtrans()
{
    size_t len = length();
    T.resize(len);
    sT.resize(len);

    Tx.resize(len);
    Ty.resize(len);

    sTx.resize(len);
    sTy.resize(len);

    Bx.resize(len);
    By.resize(len);
    sBx.resize(len);
    sBy.resize(len);
    sB.resize(len);
    B.resize(len);

    for (int i = 0; i < length(); i++)
    {
        Tx[i] = arma::mat(3, 3, arma::fill::eye);
        Ty[i] = arma::mat(3, 3, arma::fill::eye);

        Bx[i] = arma::vec(3, arma::fill::zeros);
        By[i] = arma::vec(3, arma::fill::zeros);

        T[i] = arma::mat(5, 5, arma::fill::zeros);

        B[i] = arma::vec(5, arma::fill::zeros);

        std::array<arma::mat, 6> P;
        std::array<arma::vec, 6> p;

        for (auto& val : P)
        {
            val.zeros(5, 5);
        }
        for (auto& val : p)
        {
            val.zeros(5);
        }

        int j;
        for (j = 0; j < SynchrotronOptixsflags.size(); j++)
        {
            if (GetType(i) == SynchrotronOptixsflags[j])
                break;
        }
		double L = getParameter(i, "L");
		auto errors = elements[i].GetErrors();

		if (j != 0 && j != 6 && j != 7 && j != 3 && j != 5)
		{
			L = L * (1 + errors[7]);
			auto dL = L * errors[7];

			if (i - 1 >= 0)
			{
				double L_prev = getParameter(i - 1, "L");
				setParameter(i - 1, "L", L_prev - dL / 2);
			}
			if (i + 1 == length())
			{
				double L_next = getParameter(i + 1, "L");
				setParameter(i + 1, "L", L_next - dL / 2);
			}
		}

        switch (j)
        {
        case 0: //"DRIFT"
        {
            Tx[i](0, 0) = 1;
            Tx[i](0, 1) = L;
            Tx[i](0, 2) = 0;
            Tx[i](1, 0) = 0;
            Tx[i](1, 1) = 1;
            Tx[i](1, 2) = 0;
            Tx[i](2, 0) = 0;
            Tx[i](2, 1) = 0;
            Tx[i](2, 2) = 1;

            Ty[i](0, 0) = 1;
            Ty[i](0, 1) = L;
            Ty[i](0, 2) = 0;
            Ty[i](1, 0) = 0;
            Ty[i](1, 1) = 1;
            Ty[i](1, 2) = 0;
            Ty[i](2, 0) = 0;
            Ty[i](2, 1) = 0;
            Ty[i](2, 2) = 1;

            for (size_t d = 0; d < 5; d++)
                T[i](d, d) = 1;

            T[i](0, 1) = L;
            T[i](2, 3) = L;
        };
        break;
        case 1: //"RBEND"
        {
            arma::mat Msh(3, 3);
            arma::mat Msv(3, 3);
            arma::mat Me1h(3, 3);
            arma::mat Me2h(3, 3);
            arma::mat Me1v(3, 3);
            arma::mat Me2v(3, 3);

            arma::mat M;
            M.zeros(5, 5);
            arma::mat M1;
            M1.zeros(5, 5);
            arma::mat M2;
            M2.zeros(5, 5);
            /*	double eta = getParameter(i, "ANGLE");
            double eps1 = getParameter(i, "E1");
            double eps2 = getParameter(i, "E2");

            arma::mat Msh(3, 3); arma::mat Msv(3, 3); arma::mat Me1h(3, 3); arma::mat Me2h(3, 3);
            arma::mat Me1v(3, 3); arma::mat Me2v(3, 3);


            Msh(0, 0) = std::cos(eta); Msh(0, 1) = p*std::sin(eta); Msh(0, 2) = p*(1 - std::cos(eta));
            Msh(1, 0) = -std::sin(eta) / p; Msh(1, 1) = std::cos(eta); Msh(1, 2) = std::sin(eta);
            Msh(2, 0) = 0; Msh(2, 1) = 0; Msh(2, 2) = 1;

            Msv(0, 0) = 1; Msv(0, 1) = p*eta; Msv(0, 2) = 0;
            Msv(1, 0) = 0; Msv(1, 1) = 1; Msv(1, 2) = 0;
            Msv(2, 0) = 0; Msv(2, 1) = 0; Msv(2, 2) = 1;


            Me1h(0, 0) = 1; Me1h(0, 1) = 0; Me1h(0, 2) = 0;
            Me1h(1, 0) = tan(eps1) / p; Me1h(1, 1) = 1; Me1h(1, 2) = 0;
            Me1h(2, 0) = 0; Me1h(2, 1) = 0; Me1h(2, 2) = 1;

            Me2h(0, 0) = 1; Me2h(0, 1) = 0; Me2h(0, 2) = 0;
            Me2h(1, 0) = tan(eps2) / p; Me2h(1, 1) = 1; Me2h(1, 2) = 0;
            Me2h(2, 0) = 0; Me2h(2, 1) = 0; Me2h(2, 2) = 1;

            Me1v(0, 0) = 1; Me1v(0, 1) = 0; Me1v(0, 2) = 0;
            Me1v(1, 0) = -tan(eps1) / p; Me1v(1, 1) = 1; Me1v(1, 2) = 0;
            Me1v(2, 0) = 0; Me1v(2, 1) = 0; Me1v(2, 2) = 1;

            Me2v(0, 0) = 1; Me2v(0, 1) = 0; Me2v(0, 2) = 0;
            Me2v(1, 0) = -tan(eps2) / p; Me2v(1, 1) = 1; Me2v(1, 2) = 0;
            Me2v(2, 0) = 0; Me2v(2, 1) = 0; Me2v(2, 2) = 1;

            Tx[i] = Me2h*Msh*Me1h;
            Ty[i] = Me2v*Msv*Me1v;*/

            double eta  = getParameter(i, "ANGLE");
            double eps1 = getParameter(i, "E1");
            double eps2 = getParameter(i, "E2");

            double pp = L / eta;

            for (size_t d = 0; d < 5; d++)
                M(d, d) = 1;

            M(0, 0) = Msh(0, 0) = std::cos(eta);
            M(0, 1) = Msh(0, 1) = pp * std::sin(eta);
            M(0, 4) = Msh(0, 2) = pp * (1 - std::cos(eta));
            M(1, 0) = Msh(1, 0) = -std::sin(eta) / pp;
            M(1, 1) = Msh(1, 1) = std::cos(eta);
            M(1, 4) = Msh(1, 2) = std::sin(eta);
            Msh(2, 0) = 0;
            Msh(2, 1) = 0;
            Msh(2, 2) = 1;

            Msv(0, 0) = 1;
            M(2, 3) = Msv(0, 1) = pp * eta;
            Msv(0, 2) = 0;
            Msv(1, 0) = 0;
            Msv(1, 1) = 1;
            Msv(1, 2) = 0;
            Msv(2, 0) = 0;
            Msv(2, 1) = 0;
            Msv(2, 2) = 1;

            for (size_t d = 0; d < 5; d++)
                M1(d, d) = 1;

            for (size_t d = 0; d < 5; d++)
                M2(d, d) = 1;

            Me1h(0, 0) = 1;
            Me1h(0, 1) = 0;
            Me1h(0, 2) = 0;
            M1(1, 0) = Me1h(1, 0) = tan(eta / 2 + eps1) / pp;
            Me1h(1, 1) = 1;
            Me1h(1, 2) = 0;
            Me1h(2, 0) = 0;
            Me1h(2, 1) = 0;
            Me1h(2, 2) = 1;

            Me2h(0, 0) = 1;
            Me2h(0, 1) = 0;
            Me2h(0, 2) = 0;
            M2(1, 0) = Me2h(1, 0) = tan(eta / 2 + eps2) / pp;
            Me2h(1, 1) = 1;
            Me2h(1, 2) = 0;
            Me2h(2, 0) = 0;
            Me2h(2, 1) = 0;
            Me2h(2, 2) = 1;

            Me1v(0, 0) = 1;
            Me1v(0, 1) = 0;
            Me1v(0, 2) = 0;
            M1(3, 2) = Me1v(1, 0) = -tan(eta / 2 + eps1) / pp;
            Me1v(1, 1) = 1;
            Me1v(1, 2) = 0;
            Me1v(2, 0) = 0;
            Me1v(2, 1) = 0;
            Me1v(2, 2) = 1;

            Me2v(0, 0) = 1;
            Me2v(0, 1) = 0;
            Me2v(0, 2) = 0;
            M2(3, 2) = Me2v(1, 0) = -tan(eta / 2 + eps2) / pp;
            Me2v(1, 1) = 1;
            Me2v(1, 2) = 0;
            Me2v(2, 0) = 0;
            Me2v(2, 1) = 0;
            Me2v(2, 2) = 1;

            Tx[i] = Me2h * Msh * Me1h;
            Ty[i] = Me2v * Msv * Me1v;
            T[i]  = M2 * M * M1;
			B[i](0) = errors[6] * pp * (1 - std::cos(eta));
			B[i](1) = errors[6] * pp * std::sin(eta);
			B[i] = M2 * B[i];
        }
        break;
        case 2: //"QUADRUPOLE"
        {
            double Kx = std::abs(getParameter(i, "k1")*(1 + errors[6]));
            double Ky = std::abs(getParameter(i, "k1")*(1 + errors[6]));
            // FIXME hard code ANGLE
            double eta = 0.0659;
            // double eta = getParameter(i, "ANGLE");

            double pp = L / eta;

            double ksix = L * sqrt(Kx);
            double ksiy = L * sqrt(Ky);
            T[i].zeros(5, 5);
            if (getParameter(i, "k1") > 0) // ������������ �����
            {

                T[i](0, 0) = Tx[i](0, 0) = std::cos(ksix);
                T[i](0, 1) = Tx[i](0, 1) = 1.0 / sqrt(Kx) * std::sin(ksix);
                Tx[i](0, 2) = 0;
                T[i](1, 0) = Tx[i](1, 0) = -sqrt(Kx) * std::sin(ksix);
                T[i](1, 1) = Tx[i](1, 1) = std::cos(ksix);
                Tx[i](1, 2) = 0;
                Tx[i](2, 0) = 0;
                Tx[i](2, 1) = 0;
                Tx[i](2, 2) = 1;
                //   T[i](0, 4)               = (1.0 / (Kx * pp)) * (1 - std::cos(ksix));
                //   T[i](1, 4)               = (1.0 / (sqrt(Kx) * pp)) * (std::sin(ksix));

                T[i](2, 2) = Ty[i](0, 0) = cosh(ksiy);
                T[i](2, 3) = Ty[i](0, 1) = 1.0 / sqrt(Ky) * sinh(ksiy);
                Ty[i](0, 2) = 0;
                T[i](3, 2) = Ty[i](1, 0) = sqrt(Ky) * sinh(ksiy);
                T[i](3, 3) = Ty[i](1, 1) = cosh(ksiy);
                Ty[i](1, 2) = 0;
                Ty[i](2, 0) = 0;
                Ty[i](2, 1) = 0;
                Ty[i](2, 2) = 1;

                //  T[i](2, 4) = (1.0 / (Kx * pp)) * (cosh(ksiy) - 1);
                //   T[i](3, 4) = (1.0 / (sqrt(Kx) * pp)) * (sinh(ksix));
                T[i](4, 4) = 1;
            }
            else // �������������� �����
            {
                T[i](0, 0) = Tx[i](0, 0) = cosh(ksix);
                T[i](0, 1) = Tx[i](0, 1) = 1.0 / sqrt(Kx) * sinh(ksix);
                Tx[i](0, 2) = 0;
                T[i](1, 0) = Tx[i](1, 0) = sqrt(Kx) * sinh(ksix);
                T[i](1, 1) = Tx[i](1, 1) = cosh(ksix);
                Tx[i](1, 2) = 0;
                Tx[i](2, 0) = 0;
                Tx[i](2, 1) = 0;
                Tx[i](2, 2) = 1;
                //     T[i](0, 4)               = (1.0 / (Kx * pp)) * (cosh(ksix) - 1);
                //   T[i](1, 4)               = (1.0 / (sqrt(Kx) * pp)) * (sinh(ksix));

                T[i](2, 2) = Ty[i](0, 0) = std::cos(ksiy);
                T[i](2, 3) = Ty[i](0, 1) = 1.0 / sqrt(Ky) * std::sin(ksiy);
                Ty[i](0, 2) = 0;
                T[i](3, 2) = Ty[i](1, 0) = -sqrt(Ky) * std::sin(ksiy);
                T[i](3, 3) = Ty[i](1, 1) = std::cos(ksiy);
                Ty[i](1, 2) = 0;
                Ty[i](2, 0) = 0;
                Ty[i](2, 1) = 0;
                Ty[i](2, 2) = 1;

                //  T[i](2, 4) = (1.0 / (Kx * pp)) * (1 - std::cos(ksiy));
                // T[i](3, 4) = (1.0 / (sqrt(Kx) * pp)) * (std::sin(ksix));
                T[i](4, 4) = 1;
            }
        };
        break;

        case 3: //"KICKER"
        {
            double HKICK = getParameter(i, "HKICK");
            double VKICK = getParameter(i, "VKICK");

            Tx[i](0, 0) = 1;
            Tx[i](0, 1) = L;
            Tx[i](0, 2) = 0;
            Tx[i](1, 0) = 0;
            Tx[i](1, 1) = 1;
            Tx[i](1, 2) = 0;
            Tx[i](2, 0) = 0;
            Tx[i](2, 1) = 0;
            Tx[i](2, 2) = 1;

            Ty[i](0, 0) = 1;
            Ty[i](0, 1) = L;
            Ty[i](0, 2) = 0;
            Ty[i](1, 0) = 0;
            Ty[i](1, 1) = 1;
            Ty[i](1, 2) = 0;
            Ty[i](2, 0) = 0;
            Ty[i](2, 1) = 0;
            Ty[i](2, 2) = 1;

            for (size_t d = 0; d < 5; d++)
                T[i](d, d) = 1;

            T[i](0, 1) = L;
            T[i](2, 3) = L;

            B[i](0) = Bx[i](0) = 0.5 * L * HKICK;
            B[i](1) = Bx[i](1) = HKICK;

            B[i](2) = By[i](0) = 0.5 * L * VKICK;
            B[i](3) = By[i](1) = VKICK;

            /*double Ik = 1; // (�)���� ����
            double	alpha = lambdak*Ik*Lk / Brho;

            Tx[i](0, 0) = 1; Tx[i](0, 1) = alpha*Lk / 2; Tx[i](0, 2) = 0;
            Tx[i](1, 0) = 0; Tx[i](1, 1) = alpha; Tx[i](1, 2) = 0;
            Tx[i](2, 0) = 0; Tx[i](2, 1) = 0; Tx[i](2, 2) = 1;

            Ty[i](0, 0) = 1; Ty[i](0, 1) = alpha*Lk / 2; Ty[i](0, 2) = 0;
            Ty[i](1, 0) = 0; Ty[i](1, 1) = alpha; Ty[i](1, 2) = 0;
            Ty[i](2, 0) = 0; Ty[i](2, 1) = 0; Ty[i](2, 2) = 1;*/
        }
        break;

        case 4: //"SBEND"
        {
            arma::mat Msh(3, 3);
            arma::mat Msv(3, 3);
            arma::mat Me1h(3, 3);
            arma::mat Me2h(3, 3);
            arma::mat Me1v(3, 3);
            arma::mat Me2v(3, 3);

            arma::mat M;
            M.zeros(5, 5);

            double eta  = getParameter(i, "ANGLE");
            double eps1 = getParameter(i, "E1");
            double eps2 = getParameter(i, "E2");

            double pp = L / eta;
            // double pp = 14.09;

            for (size_t d = 0; d < 5; d++)
                M(d, d) = 1;

            M(0, 0) = Msh(0, 0) = std::cos(eta);
            M(0, 1) = Msh(0, 1) = pp * std::sin(eta);
            M(0, 4) = Msh(0, 2) = pp * (1 - std::cos(eta));
            M(1, 0) = Msh(1, 0) = -std::sin(eta) / pp;
            M(1, 1) = Msh(1, 1) = std::cos(eta);
            M(1, 4) = Msh(1, 2) = std::sin(eta);
            Msh(2, 0) = 0;
            Msh(2, 1) = 0;
            Msh(2, 2) = 1;

            Msv(0, 0) = 1;
            M(2, 3) = Msv(0, 1) = pp * eta;
            Msv(0, 2) = 0;
            Msv(1, 0) = 0;
            Msv(1, 1) = 1;
            Msv(1, 2) = 0;
            Msv(2, 0) = 0;
            Msv(2, 1) = 0;
            Msv(2, 2) = 1;

            Me1h(0, 0) = 1;
            Me1h(0, 1) = 0;
            Me1h(0, 2) = 0;
            Me1h(1, 0) = tan(eta / 2 + eps1) / pp;
            Me1h(1, 1) = 1;
            Me1h(1, 2) = 0;
            Me1h(2, 0) = 0;
            Me1h(2, 1) = 0;
            Me1h(2, 2) = 1;

            Me2h(0, 0) = 1;
            Me2h(0, 1) = 0;
            Me2h(0, 2) = 0;
            Me2h(1, 0) = tan(eta / 2 + eps2) / pp;
            Me2h(1, 1) = 1;
            Me2h(1, 2) = 0;
            Me2h(2, 0) = 0;
            Me2h(2, 1) = 0;
            Me2h(2, 2) = 1;

            Me1v(0, 0) = 1;
            Me1v(0, 1) = 0;
            Me1v(0, 2) = 0;
            Me1v(1, 0) = -tan(eta / 2 + eps1) / pp;
            Me1v(1, 1) = 1;
            Me1v(1, 2) = 0;
            Me1v(2, 0) = 0;
            Me1v(2, 1) = 0;
            Me1v(2, 2) = 1;

            Me2v(0, 0) = 1;
            Me2v(0, 1) = 0;
            Me2v(0, 2) = 0;
            Me2v(1, 0) = -tan(eta / 2 + eps2) / pp;
            Me2v(1, 1) = 1;
            Me2v(1, 2) = 0;
            Me2v(2, 0) = 0;
            Me2v(2, 1) = 0;
            Me2v(2, 2) = 1;

            Tx[i] = Me2h * Msh * Me1h;
            Ty[i] = Me2v * Msv * Me1v;

            Tx[i] = Msh;
            Ty[i] = Msv;
            T[i]  = M;
			B[i](0) = errors[6] * pp * (1 - std::cos(eta));
			B[i](1) = errors[6] * pp * std::sin(eta);

        };
        break;

        case 5: //"SEXTUPOLE"
        {
            Tx[i](0, 0) = 1;
            Tx[i](0, 1) = L;
            Tx[i](0, 2) = 0;
            Tx[i](1, 0) = 0;
            Tx[i](1, 1) = 1;
            Tx[i](1, 2) = 0;
            Tx[i](2, 0) = 0;
            Tx[i](2, 1) = 0;
            Tx[i](2, 2) = 1;

            Ty[i](0, 0) = 1;
            Ty[i](0, 1) = L;
            Ty[i](0, 2) = 0;
            Ty[i](1, 0) = 0;
            Ty[i](1, 1) = 1;
            Ty[i](1, 2) = 0;
            Ty[i](2, 0) = 0;
            Ty[i](2, 1) = 0;
            Ty[i](2, 2) = 1;
        };
        default:
            T[i] = arma::mat(5, 5, arma::fill::eye);
            break;
        };

        double m11 = T[i](0, 0);
        double m12 = T[i](0, 1);
        double m21 = T[i](1, 0);
        double m22 = T[i](1, 1);

        double m33 = T[i](2, 2);
        double m34 = T[i](2, 3);
        double m43 = T[i](3, 2);
        double m44 = T[i](3, 3);

        if (j == 2)
        {
            P[0](0, 0) = -m21;
            P[0](0, 1) = m11 - m22;
            P[0](1, 0) = 0;
            P[0](1, 1) = m21;
            P[0](2, 2) = -m43;
            P[0](2, 3) = m33 - m44;
            P[0](3, 3) = m43;

            p[1](0) = 1 - m11;
            p[1](1) = -m21;

            p[2](2) = 1 - m33;
            p[2](3) = -m43;

            P[3](0, 2) = m11 - m33;
            P[3](0, 3) = m12 - m34;
            P[3](1, 2) = m21 - m33;
            P[3](1, 3) = m22 - m44;
            P[3](2, 0) = m11 - m33;
            P[3](2, 1) = m12 - m34;
            P[3](3, 0) = m21 - m43;
            P[3](3, 1) = m22 - m44;

            p[4](0) = 0.5 * L * (1 + m11) - m12;
            p[4](1) = 0.5 * L * m21 + 1 - m22;

            p[5](2) = m34 - 0.5 * L * (1 + m33);
            p[5](3) = m44 - 1 - 0.5 * L * m43;
        }

        if (j == 1 || j == 4)
        {
            double phi = getParameter(i, "ANGLE");

			double Ro = L / phi;

            P[0](0, 0) = -m21 * std::cos(phi);
            P[0](0, 1) = m11 - m22 * std::cos(phi);
            P[0](1, 1) = m21;
            P[0](2, 2) = -m43 * std::cos(phi);
            P[0](3, 2) = m33 - m44 * std::cos(phi);
            P[0](3, 3) = m43;

            p[0](0) = -std::sin(phi);

            P[1](0, 0) = -m21 * std::sin(phi);
            P[1](0, 1) = -m22 * std::sin(phi);
            P[1](2, 2) = -m43 * std::sin(phi);
            P[1](3, 2) = -m44 * std::sin(phi);

            p[1](0) = std::cos(phi) - m11;
            p[1](1) = -m21;

            p[2](2) = 1 - m33;
            p[2](3) = -m43;

            P[3](2, 0) = m11 * std::cos(phi) - m33;
            P[3](3, 0) = m21 * std::cos(phi) - m43;
            P[3](2, 1) = m12 * std::cos(phi) - m34;
            P[3](3, 1) = m22 * std::cos(phi) - m44;
            P[3](0, 2) = m11 - m33 * std::cos(phi);
            P[3](0, 3) = m12 - m34 * std::cos(phi);
            P[3](1, 2) = m21 - m43 * std::cos(phi);
            P[3](1, 3) = m22 - m44 * std::cos(phi);

            p[3](2) = Ro * std::sin(phi) * std::tan(phi / 2.0);
            p[3](3) = std::sin(phi);

            p[4](0) = (1 + m11) * Ro * std::tan(phi / 2.0) - m12;
            p[4](1) = m21 * Ro * std::tan(phi / 2.0) + 1 - m22;

            P[5](2, 0) = m11 * std::sin(phi);
            P[5](3, 0) = m21 * std::sin(phi);
            P[5](2, 1) = m12 * std::sin(phi);
            P[5](3, 1) = m22 * std::sin(phi);
            P[5](0, 2) = -m33 * std::sin(phi);
            P[5](0, 3) = -m34 * std::sin(phi);
            P[5](1, 2) = -m43 * std::sin(phi);
            P[5](1, 3) = -m44 * std::sin(phi);

            p[5](2) = m34 - Ro * (std::cos(phi) + m33) * std::tan(phi / 2.0);
            p[5](3) = m44 - std::cos(phi) - Ro * m43 * std::tan(phi / 2.0);
        }


        /*if (GetType(i) == "QUADRUPOLE") //TODO remove
        {
                errors.resize(6);
                errors[0] = 0.1;
                errors[1] = 0.1;
                errors[2] = 0.1;
                errors[3] = 0.1;
                errors[4] = 0.1;
                errors[5] = 0.1;
        }*/

        if (errors.size() >= 6)
        {
            for (size_t sig = 0; sig < 6; sig++)
            {
                T[i] = T[i] + P[sig] * errors[sig];
            }
            for (size_t sig = 0; sig < 6; sig++)
            {
                B[i] = B[i] + p[sig] * errors[sig];
            }
        }

        if (i == 0)
        {
            sTx[i] = Tx[i];
            sTy[i] = Ty[i];
            sT[i]  = T[i];

            sBx[i] = Bx[i];
            sBy[i] = By[i];
            sB[i]  = B[i];
        }
        else
        {
            sTx[i] = Tx[i] * sTx[i - 1];
            sTy[i] = Ty[i] * sTy[i - 1];
            sT[i]  = T[i] * sT[i - 1];

            sBx[i] = Tx[i] * sBx[i - 1] + Bx[i];
            sBy[i] = Ty[i] * sBy[i - 1] + By[i];
            sB[i]  = T[i] * sB[i - 1] + B[i];
        }
    };
};

OpticElement::OpticElement(const std::vector<std::string>& input, std::string& error)
{
    label = input[0];
    type  = input[1];
    double tmp;
    int    j = 0;
    for (j = 0; j < SynchrotronOptixsflags.size(); j++)
    {
        if (type == SynchrotronOptixsflags[j])
            break;
    }
    if (j == SynchrotronOptixsflags.size() && input.size() > 3)
    {
        error = "Incorrect description of element " + type;
        return;
    };
    /*if (input.size() > 3)
    {
            for (int k = 0; k < CorrectParametersNames[j].size(); k++)
            {
                    int i;
                    for (i = 2; i < input.size() - 1; i = i + 2)
                    {
                            if (CorrectParametersNames[j][k] == input[i])
                                    break;
                    }
                    if (i == input.size() - 1)
                    {
                            error = "Incomplete description of element " + type;
                            return;
                    };
            };
    }*/
    for (int j = 3; j < input.size(); j = j + 2)
    {
        int er = sscanf(input[j].c_str(), "%lf", &tmp);
        if (!er)
        {
            error = "Incorrect input character " + input[j];
            return;
        };
        parameters.insert(std::make_pair(input[j - 1], tmp));
    }
    if (type == "RBEND" || type == "QUADRUPOLE" || type == "SBEND")
    {
        alignment_errors.resize(8);
    }
    std::fill(alignment_errors.begin(), alignment_errors.end(), 0);
};
void OpticElement::copy(const OpticElement& obj)
{
    type = obj.type;
    if (parameters.size() == 0)
        parameters = obj.parameters;
};
void OpticElement::copy(const std::shared_ptr<OpticElement> obj)
{
    type = obj->type;
    if (parameters.size() == 0)
        parameters = obj->parameters;
};

void OpticElement::InsertErrors(const std::vector<double>& errors)
{
    alignment_errors = errors;
}

void OpticElement::GenerateErrors(const std::map<std::string, std::array<double, 8>>& errors,
                                  std::default_random_engine&       generator,
                                  std::normal_distribution<double>& distribution)
{
    auto fillup = [&](const std::string& key) {
        alignment_errors.resize(8);
        for (size_t i = 0; i < 8; i++)
        {
            alignment_errors[i] =
                ((errors.find(key)->second[i] / (2.0 * 3.0)) * distribution(generator));
        }
    };

    if (type == "RBEND")
    {
        fillup("rbend");
    }
    else if (type == "QUADRUPOLE")
    {
        fillup("quadrupole");
    }
    else if (type == "SBEND")
    {
        fillup("sbend");
    }
}
void OpticElement::GenerateErrors(double RelError, std::default_random_engine& generator,
                                  std::normal_distribution<double>& distribution)
{
    if (type == "RBEND" || type == "QUADRUPOLE" || type == "SBEND")
    {

        alignment_errors.resize(8);
        for (auto& er : alignment_errors)
        {
            er = ((RelError / (2.0 * 3.0)) * distribution(generator));
        }

        /*for (auto it = parameters.begin(); it != parameters.end(); ++it)
        {
        if (0 == it->second)
        alignment_errors.push_back((RelError / (2.0 * 3.0)) * distribution(generator));
        else
        alignment_errors.push_back((it->second * RelError / (2.0 * 3.0)) * distribution(generator));
        }*/
    }
};

OpticElement* OpticElement::copy() const
{
    OpticElement* result = new OpticElement();
    *result              = *this;
    return result;
};

std::string OpticElementsSequence::GetType(int i) const
{
    return elements[i].GetType();
};
std::string OpticElementsSequence::GetLabel(int i) const
{
    return elements[i].GetLabel();
};
void OpticElementsSequence::SetLabel(const int i, const std::string& label)
{
	elements[i].SetLabel(label);
}


int OpticElementsSequence::insertelem(std::vector<std::shared_ptr<OpticElement>>& elementsBase,
                                      std::vector<std::shared_ptr<OpticElementsSequence>>& tmpSeq,
                                      const std::string&                                   name)
{
    int k = findOpticsElem(elementsBase, name);
    if (k != -1)
        elements.push_back(*(elementsBase[k]));
    else
    {
        if (name.size() && name != "\n")
        {
            int m = findline(tmpSeq, name);
            if (m != -1)
            {
                for (int i = 0; i < tmpSeq[m]->elements.size(); i++)
                    elements.push_back(tmpSeq[m]->elements[i]);
            }
            else
                return 0;
        }
    }
    return 1;
}
void OpticElementsSequence::insertMonitors(int nmonitors)
{
    elements.erase(std::remove_if(std::begin(elements), std::end(elements),
                                  [](auto el) { return el.GetType() == "MONITOR"; }),
                   elements.end());

    int dn = elements.size() / nmonitors;

    std::vector<int> dns(nmonitors);
    std::fill(std::begin(dns), std::end(dns), dn);

    for (int k = 0; k < elements.size() - dn * nmonitors; k++)
    {
        int ind = rand() % dns.size();
        dns[ind]++;
    };

    int insPos = 0;
    for (int n = 0; n < nmonitors; n++)
    {
        insPos = insPos + dns[n];
        if (insPos == elements.size())
            break;
        auto el = elements.insert(elements.begin() + insPos, OpticElement("MONITOR"));
        el->at  = elements[insPos - 1].at;
        el->parameters.insert(std::make_pair("L", 0.0));
        insPos++;
    }
}
void OpticElementsSequence::SetFlags()
{
    std::string tmp;
    elements.insert(elements.begin(), OpticElement({name + "$START", "MARKER", "L", "0"}, tmp));
    elements.push_back(OpticElement({name + "$END", "MARKER", "L", "0"}, tmp));
    setParameter(0, "at", 0);
    setParameter(elements.size() - 1, "at", L);
}
void OpticElementsSequence::SetL(double Lin)
{
    L = Lin;
}
double OpticElementsSequence::GetL() const
{
    return L;
}
int OpticElementsSequence::setParameter(int i, std::string marker, double value)
{
    if (marker == "at")
        elements[i].at = value;
    else
        elements[i].parameters[marker] = value;
    return 1;
};
std::vector<int> OpticElementsSequence::findType(const std::string& name)
{
    std::vector<int> result;
    for (int i = 0; i < elements.size(); i++)
    {
        if (elements[i].GetType() == name)
            result.push_back(i);
    };
    return result;
};

double OpticElementsSequence::getParameter(int i, std::string marker) const
{
    if (marker == "at")
        return elements[i].at;
    else
    {
        /*auto tmp =c;

        if (elements[i].GetErrors().size()!=0)
        {
                size_t k = 0;
                for (auto it = tmp.begin(); it != tmp.end(); ++it)
                {
                        it->second = it->second + elements[i].GetErrors()[k];
                        k++;
                }
        }*/
        return elements[i].parameters.find(marker)->second;
        // return elements[i].getParameter(marker);

        // return val;
    }
};
size_t OpticElementsSequence::length() const
{
    return elements.size();
};
OpticElementsSequence::OpticElementsSequence(const std::string& name)
{
    L          = 0;
    this->name = name;
};
OpticElementsSequence::OpticElementsSequence()
{
    distribution_loc = std::normal_distribution<double>(0, 1);

    for (size_t i = 0; i < 1e5; i++)
    {
        distribution_loc(generator_loc);
    }
};

void OpticElementsSequence::InsertErrors(
    const std::vector<std::pair<std::string, std::array<double, 8>>>& errors)
{
    volatile int ind = 0;
    for (int i = 0; i < length(); i++)
    {

        auto it = std::find_if(errors.begin(), errors.end(), [&](const auto& elem) {
            return elem.first == elements[i].GetLabel();
        });

        if (it != errors.end())
        {
            std::vector<double> cur_err(it->second.begin(), it->second.end());

            elements[i].InsertErrors(cur_err);
        }

        /*	size_t s_cur = elements[i].GetErrors().size();
                std::vector<double> cur_err(errors.begin()+ind, errors.begin() + ind + s_cur);
                elements[i].InsertErrors(cur_err);
                ind = ind + s_cur;*/
    }
}

void OpticElementsSequence::SaveMADObsCommands(const std::string& name)
{
    FILE* fid = fopen(name.c_str(), "w");
    for (int i = 0; i < length(); i++)
    {
        fprintf(fid, "ptc_observe, place = %s;\n", GetLabel(i).c_str());
    }
    fclose(fid);
}

void OpticElementsSequence::SaveMADAlignmentCommands(const std::string& name)
{
    FILE* fid = fopen(name.c_str(), "w");

    for (int i = 0; i < length(); i++)
    {
        auto err = elements[i].GetErrors();

        if (GetType(i) == "QUADRUPOLE") // TODO remove
        {
            err.resize(6);
            err[0] = 0.1;
            err[1] = 0.1;
            err[2] = 0.1;
            err[3] = 0.1;
            err[4] = 0.1;
            err[5] = 0.1;
        }

        if (err.size() != 6)
        {
            continue;
        }

        fprintf(fid, "SELECT, FLAG=ERROR, PATTERN = \"%s\";\n", GetLabel(i).c_str());
        // fprintf(fid, "SELECT, FLAG=ERROR, CLASS = quadrupole;\n");

        fprintf(fid,
                "EALIGN,  DS = %lf, DX = %lf, DY = %lf, DPSI = %lf, DTHETA = %lf, DPHI = %lf; \n",
                err[0], err[1], err[2], err[3], err[4], err[5]);
    }

    fclose(fid);
}

int findOpticsElem(std::vector<std::shared_ptr<OpticElement>>& optics, const std::string& elem)
{
    int i;
    for (i = 0; i < optics.size(); i++)
    {
        if (elem == optics[i]->GetLabel())
            return i;
    }
    i = -1;
    return i;
};
int findline(std::vector<std::shared_ptr<OpticElementsSequence>>& tmpSeq, const std::string& str)
{
    // ����� ������� �������� �� �����
    for (int i = 0; i < tmpSeq.size(); i++)
    {
        if (str == tmpSeq[i]->GetName())
            return i;
    }
    return -1;
}
arma::vec twiss(const arma::mat& M, const arma::mat& x)
{
    double m11 = M(0, 0);
    double m12 = M(0, 1);
    double m21 = M(1, 0);
    double m22 = M(1, 1);

    arma::mat Mt(3, 3);

    Mt(0, 0) = m11 * m11;
    Mt(0, 1) = -2 * m11 * m12;
    Mt(0, 2) = m12 * m12;
    Mt(1, 0) = -m11 * m21;
    Mt(1, 1) = m11 * m22 + m12 * m21;
    Mt(1, 2) = -m22 * m12;
    Mt(2, 0) = m21 * m21;
    Mt(2, 1) = -2 * m22 * m21;
    Mt(2, 2) = m22 * m22;

    return Mt * x;
};

void updateBeamPositions(const arma::mat& Mx, const arma::mat& My, LinacDynamicsAccel& dyn)
{
    double tmp;
    for (int i = 0; i < dyn.x.size(); i++)
    {
        tmp       = dyn.x[i];
        dyn.x[i]  = Mx(0, 0) * dyn.x[i] + Mx(0, 1) * dyn.dx[i];
        dyn.dx[i] = Mx(1, 0) * tmp + Mx(1, 1) * dyn.dx[i];

        tmp       = dyn.y[i];
        dyn.y[i]  = My(0, 0) * dyn.y[i] + My(0, 1) * dyn.dy[i];
        dyn.dy[i] = My(1, 0) * tmp + My(1, 1) * dyn.dy[i];
    };
};

void setOutputData(int i, int circle, const std::vector<int>& ind, const std::vector<int>& indP,
                   std::vector<void*>                                          dataIn,
                   std::vector<std::vector<std::vector<std::vector<double>>>>& data){

};
bool isEqual_el(const OpticElement* A, const OpticElement* B)
{
    if (A->type == B->label)
        return true;
    return false;
};

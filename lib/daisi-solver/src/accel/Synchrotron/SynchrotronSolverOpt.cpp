#define NOMINMAX
#include "SynchrotronSolver.h"
#include "../base/AccelFlow.h"
#include "Results.h"
#include "SynchrotronDevice.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
#include "SynchrotronToolsSim.h"
#include "Tools.h"

#include <notk/controller.hpp>
#include <common_tools/constants.h>

#include <fstream>
#include <memory>
#include <sstream>

#ifdef WIN32
#include <Windows.h>
#endif

#undef RELATIVE
#undef ABSOLUTE


std::vector<double> weights;

std::vector<double> extract(const std::vector<std::pair<std::string, std::array<double, 8>>>& data,
                            const std::pair<size_t, size_t>& inds)
{
    std::vector<double> result;

    for (size_t i = inds.first; i < inds.second; i++)
    {
        for (size_t ii = 0; ii < 8; ii++)
        {
            result.push_back(data[i].second[ii]);
        }
    }

    return result;
}

void insert(std::vector<std::pair<std::string, std::array<double, 8>>>& data,
            const std::pair<size_t, size_t>& inds, const std::vector<double>& errors)
{
    size_t ind = 0;
    for (size_t i = inds.first; i < inds.second; i++)
    {
        for (size_t ii = 0; ii < 8; ii++)
        {
            data[i].second[ii] = errors[ind];
            ind++;
        }
    }
}

double FitnessCenterMassRecostruction(const arma::vec&                        initialPos,
                                      const std::vector<std::vector<double>>& xcmS,
                                      std::shared_ptr<SynchrotronFlow>&       flow,
                                      std::pair<size_t, size_t> indexes,
                                      std::shared_ptr<OpticElementsSequence>& reconstructed_seq,
                                      const std::vector<int>&                 monitiors,
                                      const std::vector<double>&              errors)
{

    std::shared_ptr<OpticElementsSequence> new_seq =
        std::make_shared<OpticElementsSequence>(*reconstructed_seq);

    auto errors_current = new_seq->GetErrorsStruct();

    insert(errors_current, indexes, errors);

    new_seq->InsertErrors(errors_current);

    new_seq->mtrans();

    auto result = simulationTrace<double>(initialPos, new_seq, 1, monitiors, 1);

    double fit = 0;
    for (int i = 0; i < xcmS[0].size(); i++)
    {
        fit = fit + weights[i] * (result[1][i] - xcmS[1][i]) * (result[1][i] - xcmS[1][i]);
        fit = fit + weights[i] * (result[2][i] - xcmS[2][i]) * (result[2][i] - xcmS[2][i]);
    }
    return fit / double(xcmS[0].size());
}

double FitnessCenterMassRecostructionMatrix(
    const arma::vec& initialPos, const arma::mat& Rx, const arma::mat& Ry,
    std::shared_ptr<SynchrotronFlow>& flow, std::pair<size_t, size_t> indexes,
    std::shared_ptr<OpticElementsSequence>& reconstructed_seq, const std::vector<int>& monitiors,
    const std::vector<double>& errors)
{
    arma::mat RxS;
    arma::mat RyS;

    std::shared_ptr<OpticElementsSequence> new_seq =
        std::make_shared<OpticElementsSequence>(*reconstructed_seq);
    auto errors_current = new_seq->GetErrorsStruct();

    insert(errors_current, indexes, errors);

    new_seq->InsertErrors(errors_current);

    new_seq->mtrans();

    CorrMatrix(initialPos, new_seq, RxS, RyS, monitiors);

    double fit  = 0;
    double nind = 0;
    for (int i = 0; i < RxS.n_rows; i++)
    {
        for (int j = 0; j < RxS.n_cols; j++)
        {
            fit = fit + weights[i] * (Rx(i, j) - RxS(i, j)) * (Rx(i, j) - RxS(i, j));
            fit = fit + weights[i] * (Ry(i, j) - RyS(i, j)) * (Ry(i, j) - RyS(i, j));
            nind++;
        }
    }
    return fit / nind;
}

void SynchrotronSolver::SynchrotronOptimization(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
    std::vector<std::string>& status)
{
    /* status.clear();

     DynamicsSimulationCenterMass(device, progress, flagAbort, outputData, errorMsg);
     // Rosenbrock r(24);
     int nThreads = solverParameters[SynchrotronSolversOptimization]->find("Number of threads");
     int nIters =
             solverParameters[SynchrotronSolversOptimization]->find("Maximal number of iterations");
     double Range = solverParameters[SynchrotronSolversOptimization]->find(
             "Parameters search range (Relative)");
     int populationSize = solverParameters[SynchrotronSolversOptimization]->find("Population size");

     int optMethod = solverParametersFlags[SynchrotronSolversOptimization][0];
     int fitType   = solverParametersFlags[SynchrotronSolversOptimization][1];

     std::vector<double> result;
     double              start;
     double              end;

     if (nIters <= 0 || nThreads <= 0)
     {
             errorMsg = "Incorrect parameters";
             return;
     }

     std::unique_ptr<OptimizationBase> solver;
     switch (optMethod)
     {
     case 0:
             solver = std::unique_ptr<OptimizationBase>(new ParticlesSwarmOptimization());
             break;
     }

     solver->init(std::bind(FitnessBeam, device, fitType, 0, std::placeholders::_1),
                              std::bind(GeneratePopulation, device, Range, 0,
     std::placeholders::_1),
                              std::bind(CheckPopulation, device, Range, 0, std::placeholders::_1),
     nThreads);
     std::vector<float> bestValues;

     start  = omp_get_wtime();
     result = solver->SearchMin(populationSize, nIters, flagAbort, status, bestValues);
     end    = omp_get_wtime();

     status.push_back("Time elapsed " + std::to_string(end - start));

     device->GetOpticElementsSequence()->InsertErrors(0, result);
     device->GetOpticElementsSequence()->mtrans();
     DynamicsSimulationCenterMass(device, progress, flagAbort, outputData, errorMsg);*/
}

void SynchrotronSolver::SynchrotronSolversTolerancesReconstruction(
    std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
    std::vector<std::string>& status)
{
   /* auto pt_error = readJSONFile(inputFileNames[SynchrotronSolversTolerancesNameReconstruction][3]);

    double step_delta;
    int    n_steps;

    if (!pt_error)
    {
        errorMsg = "An error occured during read error data file";
        return;
    }

    auto read_error_main = [&](const boost::property_tree::ptree& pt) -> bool {

        n_steps    = pt.get<int>("n_steps");
        step_delta = pt.get<double>("step_delta");
        return true;
    };

    if (!read_property_tree(*pt_error, "error read error data file", read_error_main))
    {
        errorMsg = "An error occured during read initial data file";
        return;
    }

    OpticElementsSequence final_seq = *device->GetOpticElementsSequence();

    for (int n_err_step = 0; n_err_step < n_steps; n_err_step++)
    {
        outputData.push_back(std::shared_ptr<SimulationDataAccel>(new SimulationDataAccel(
            SynchrotronRecCMFlags, SynchrotronSolversTolerancesNameReconstruction)));

        auto res0_CM = simulationCenterMass<float>(device->GetLinacFlows()[0],
                                                   device->GetOpticElementsSequence(), 1, "ALL", 0);
        auto res0_TW = DynamicsSimulationTwissF<float>(device->GetLinacFlows()[0]->GetTwissVector(),
                                                       device->GetOpticElementsSequence(), 1);
        outputData.back()->YData.push_back(res0_CM);
        outputData.back()->YData.push_back(res0_TW);

        // double Range =
        //	solverParameters[SynchrotronSolversTolerancesNameReconstruction]->find("Relative
        // error");

        int fitType = solverParametersFlags[SynchrotronSolversTolerancesNameReconstruction][0];

        auto pt_initial =
            readJSONFile(inputFileNames[SynchrotronSolversTolerancesNameReconstruction][2]);

        if (!pt_initial)
        {
            errorMsg = "An error occured during read initial data file";
            return;
        }

        std::vector<arma::vec> initial_pos;

        auto read_initial = [&](const boost::property_tree::ptree& pt) -> bool {
            bool result = true;
            for (auto& item : pt.get_child("initial_data"))
            {
                arma::vec pos(5);
                pos(0) = item.second.get<double>("x");
                pos(1) = item.second.get<double>("px");
                pos(2) = item.second.get<double>("y");
                pos(3) = item.second.get<double>("py");
                pos(4) = item.second.get<double>("dp");

                initial_pos.push_back(pos);
            }
            return result;
        };

        if (!read_property_tree(*pt_initial, "error read initial data file", read_initial))
        {
            errorMsg = "An error occured during read initial data file";
            return;
        }

        std::map<std::string, std::array<double, 8>> errors;

        auto read_error = [&](const boost::property_tree::ptree& pt) -> bool {
            auto _read_loc = [&](const std::string name) {
                auto pt_child = pt.get_child(name);
                std::array<double, 8> errors_;
                errors_[0] = pt_child.get<double>("z") + step_delta * n_err_step;
                errors_[1] = pt_child.get<double>("x") + step_delta * n_err_step;
                errors_[2] = pt_child.get<double>("y") + step_delta * n_err_step;
                errors_[3] = pt_child.get<double>("xy") + step_delta * n_err_step;
                errors_[4] = pt_child.get<double>("xz") + step_delta * n_err_step;
                errors_[5] = pt_child.get<double>("yz") + step_delta * n_err_step;
				errors_[6] = pt_child.get<double>("sigma") + step_delta * n_err_step;
				errors_[7] = pt_child.get<double>("L") + step_delta * n_err_step;

                errors.emplace(name, errors_);
            };
            _read_loc("quadrupole");
            _read_loc("rbend");
            _read_loc("sbend");

            return true;
        };

        if (!read_property_tree(*pt_error, "error read error data file", read_error))
        {
            errorMsg = "An error occured during read errors data file";
            return;
        }

        std::ifstream fS;
        fS.open(inputFileNames[SynchrotronSolversTolerancesNameReconstruction][0],
                std::ifstream::in);

        std::vector<double> mon_flags;
        std::string         str;

        std::vector<int> mon_nums = device->GetOpticElementsSequence()->findType("MONITOR");

        mon_nums.insert(mon_nums.begin(), 0);

        weights.resize(mon_nums.size());
        std::fill(std::begin(weights), std::end(weights), 1);

        mon_flags.resize(mon_nums.size());
        mon_flags[0] = 1;
        weights[0]   = 1;

        int ii = 1;

        while (fS)
        {
            fS >> str;
            fS >> weights[ii];
            fS >> mon_flags[ii];
            ii++;
			if (ii == mon_flags.size())
			{
				break;
			}
        }

        std::shared_ptr<OpticElementsSequence> reconstructed_seq =
            std::make_shared<OpticElementsSequence>(*device->GetOpticElementsSequence());

        auto seq_err = *device->GetOpticElementsSequence();
        seq_err.GenerateErrors(errors);
        auto errors_expected = seq_err.GetErrorsStruct();
        auto errors_borders  = seq_err.GetErrorsBorders(errors);

        seq_err.mtrans();

        auto resEr_CM = simulationCenterMass<float>(
            device->GetLinacFlows()[0], std::make_shared<OpticElementsSequence>(seq_err), 1, "ALL",
            0);
        auto resEr_TW =
            DynamicsSimulationTwissF<float>(device->GetLinacFlows()[0]->GetTwissVector(),
                                            std::make_shared<OpticElementsSequence>(seq_err), 1);

        outputData.back()->YData.push_back(resEr_CM);
        outputData.back()->YData.push_back(resEr_TW);

        std::vector<float> bestValues;

        std::vector<std::vector<float>> arg_seqs(errors_expected.size() * 8);

        size_t pos_i = 0;

        auto initial_pos_curr = initial_pos;

        for (size_t mon_i = 0; mon_i < mon_flags.size(); mon_i++)
        {
			if (!flagAbort)
			{
				break;
			}
            if (mon_flags[mon_i] != 1)
            {
                continue;
            }
            for (size_t mon_ii = mon_i + 1; mon_ii < mon_flags.size(); mon_ii++)
            {
				if (!flagAbort)
				{
					break;
				}
                if (mon_flags[mon_ii] != 1)
                {
                    continue;
                }

                size_t ind_r = std::min(mon_ii + 1, mon_nums.size());

                std::vector<int> m_indexes(mon_nums.begin() + mon_i, mon_nums.begin() + ind_r);
                auto             err_ind = seq_err.GetErrIndexes(m_indexes);

                std::vector<double> x_left(errors_borders.first.begin() + err_ind.first * 8,
                                           errors_borders.first.begin() + err_ind.second * 8);

                std::vector<double> x_right(errors_borders.second.begin() + err_ind.first * 8,
                                            errors_borders.second.begin() + err_ind.second * 8);

                for (const auto initialPos : initial_pos_curr)
                {
					if (!flagAbort)
					{
						break;
					}

                    auto errors_current_ = reconstructed_seq->GetErrorsStruct();

                    std::vector<double> x_0 = extract(errors_current_, err_ind);

                    notk::NOTKController<double, double> NOTKController;

                    NOTKController.add_problem_config(
                        inputFileNames[SynchrotronSolversTolerancesNameReconstruction][1]);
                    NOTKController.set_borders(x_left, x_right);
                    NOTKController.set_x_0(x_0);

                    if (fitType == 0)
                    {
                        auto res = simulationTrace<double>(
                            initialPos, std::make_shared<OpticElementsSequence>(seq_err), 1,
                            m_indexes, 1);
                        NOTKController.set_fitness(
                            std::bind(FitnessCenterMassRecostruction, initialPos, res,
                                      device->GetLinacFlows()[0], err_ind, reconstructed_seq,
                                      m_indexes, std::placeholders::_1));
                    }
                    if (fitType == 1)
                    {
                        arma::mat Rx;
                        arma::mat Ry;
                        CorrMatrix(initialPos,
                                   std::make_shared<OpticElementsSequence>(*reconstructed_seq), Rx,
                                   Ry, m_indexes);
                        NOTKController.set_fitness(
                            std::bind(FitnessCenterMassRecostructionMatrix, initialPos, Rx, Ry,
                                      device->GetLinacFlows()[0], err_ind, reconstructed_seq,
                                      m_indexes, std::placeholders::_1));
                    }

                    if (!NOTKController.check_configuration())
                    {
                        errorMsg = "An error occured during notk library initialization. Please "
                                   "see log files for additional info.";
                        return;
                    }

                    status.push_back("Start automatical reconstruction between  " +
                                     std::to_string(mon_i) + " and " + std::to_string(mon_ii));

                    auto result = NOTKController.process(flagAbort, status);

					if (!result)
					{
						errorMsg = "An error occured during notk library calcualtion. Please "
							"see log files for additional info.";
						return;
					}
                    auto errors_current = reconstructed_seq->GetErrorsStruct();

                    insert(errors_current, err_ind, result->get_last_argument());


                    reconstructed_seq->InsertErrors(errors_current);

                    //	reconstructed_seq->InsertErrors(errors_current);

                    reconstructed_seq->mtrans();

                    auto fit_seq = result->get_fit_array<float>();

                    auto arg_seq = result->get_arg_arr();

                    for (size_t pos_i = 0; pos_i < arg_seq.size(); pos_i++)
                    {
                        for (size_t pos_ii = 0; pos_ii < arg_seq[pos_i].size(); pos_ii++)
                        {
                            arg_seqs[err_ind.first * 8 + pos_i].push_back(arg_seq[pos_i][pos_ii]);
                        }
                    }

                    std::stringstream str_r;
                    str_r << std::scientific << fit_seq.back();

                    status.push_back("done, fitness = " + str_r.str());

                    bestValues.insert(bestValues.end(), fit_seq.begin(), fit_seq.end());
                }
                size_t pos_i = 0;
                for (const auto initialPos : initial_pos_curr)
                {
                    auto res =
                        simulationTrace<double>(initialPos, reconstructed_seq, 1, m_indexes, 0);

                    initial_pos_curr[pos_i](0) = res[1].back();
                    initial_pos_curr[pos_i](1) = res[2].back();
                    initial_pos_curr[pos_i](2) = res[3].back();
                    initial_pos_curr[pos_i](3) = res[4].back();
                    pos_i++;
                }
                break;
            }
            pos_i++;
        }

        auto res_best_CM =
            simulationCenterMass<float>(device->GetLinacFlows()[0], reconstructed_seq, 1, "ALL", 0);
        auto res_best_TW = DynamicsSimulationTwissF<float>(
            device->GetLinacFlows()[0]->GetTwissVector(), reconstructed_seq, 1);

        outputData.back()->YData.push_back(res_best_CM);
        outputData.back()->YData.push_back(res_best_TW);

        outputData.back()->YData1.push_back(arg_seqs);

        outputData.back()->XData.push_back(bestValues);

        auto rec_err = reconstructed_seq->GetErrorsStruct();
        auto err     = seq_err.GetErrorsStruct();

        save_errors(err, outpuFilePaths[SynchrotronSolversTolerancesNameReconstruction][1]);
        save_errors(rec_err, outpuFilePaths[SynchrotronSolversTolerancesNameReconstruction][2]);

        std::vector<std::vector<float>> tt(1);

        final_seq = seq_err;

        std::ofstream fS_out(outpuFilePaths[SynchrotronSolversTolerancesNameReconstruction][0] +
                                 std::to_string(n_err_step),
                             std::ofstream::out);
        for (size_t i = 0; i < rec_err.size(); i++)
        {
            for (size_t j = 0; j < 8; j++)
            {
                tt[0].push_back(err[i].second[j]);
                fS_out << err[i].second[j] << "\t";
                fS_out << rec_err[i].second[j] << "\t";

                fS_out << "\n";
            }
        }
        outputData.back()->YData1.push_back(tt);
    }*/
}
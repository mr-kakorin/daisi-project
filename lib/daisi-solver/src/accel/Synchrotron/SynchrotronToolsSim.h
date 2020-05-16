#ifndef SynchrotronTools_H
#define SynchrotronTools_H
#include <armadillo>
#include <memory>
#include <vector>

class SynchrotronDevice;
class OpticElementsSequence;
class SimulationDataAccel;
class SynchrotronFlow;

template <class Tout>
std::vector<std::vector<Tout>>
simulationCenterMass(std::shared_ptr<SynchrotronFlow>&             flow,
                     const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle,
                     const std::string& regime, const int saveType);

template <class Tout>
std::vector<std::vector<Tout>>
simulationCenterMass(
	const std::vector<double>& initialPos,
	const std::shared_ptr<OpticElementsSequence>& sequence,
	const int ncircle,
	const std::string& regime,
	const int saveType);

template <class T>
std::vector<std::vector<T>>  DynamicsSimulationTwissF(const arma::vec& x0,
	const std::shared_ptr<OpticElementsSequence>& sequence,
	const int ncircle);

template <class Tout>
void simulationBeam(std::vector<std::vector<float>>&              Xdata,
                    std::vector<std::vector<std::vector<float>>>& Ydata,
                    std::shared_ptr<SynchrotronFlow>&             flow,
                    const std::shared_ptr<OpticElementsSequence>& sequence, const int ncircle,
                    double& progress, bool& flagAbort, const std::string& regime,
                    const int saveType, const int nParticles);

void InitOrbit(std::shared_ptr<SynchrotronDevice>& device, std::vector<arma::mat>& x,
               std::vector<arma::mat>& y);

void CorrMatrix(const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& seq,
	arma::mat& Rx, arma::mat& Ry, std::vector<int> nump);

std::vector<float> calculateTransmission(const std::vector<std::vector<std::vector<float>>>& Ydata,
                                         float rx, float ry, float ChannelRel = 1.0);

template <class Tout>
std::vector<std::vector<Tout>>
simulationTrace(const arma::vec& x0, const std::shared_ptr<OpticElementsSequence>& sequence,
	const int ncircle, std::vector<int> m_indexes, const int saveType);

std::vector<int> getIndexes(const std::shared_ptr<OpticElementsSequence>& sequence,
	const std::string&                            regime);

std::pair<std::pair<std::vector<double>, std::vector<double>>, std::string>
calc_closed_orbit(const std::shared_ptr<OpticElementsSequence>& sequence,
	const std::string& notk_config,	const double x_max, const double y_max);

#endif
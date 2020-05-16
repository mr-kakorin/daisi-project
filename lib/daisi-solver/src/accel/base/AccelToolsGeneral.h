#include <memory>
#include <vector>
class ExternalAccelGrid
{
  public:
    std::vector<double> Ez;
    std::vector<double> EzCol;
    std::vector<double> Ex;
    std::vector<double> Ey;
    std::vector<double> z;
    std::vector<double> rho;
    void CreateRFQGrid(const std::vector<double>& u1cells, const std::vector<double>& L,
                       const std::vector<double>& MinumalRadii, int nCells, double lambda, double dz, double voltage);
    void CreateDTLGrid();
};
class SimulationDataAccel;
class AccelFlow;
void SimulateFlowDynamicsEnvelopes(std::shared_ptr<ExternalAccelGrid> grid, int stepsPerCell, int startCell,
                                   int flagType, std::shared_ptr<AccelFlow>& flow, double frequency,
                                   const std::vector<double>& L, const std::vector<double>& MinumalRadii,
                                   std::shared_ptr<SimulationDataAccel> outputData, int nsaveTr, double startPhase,
                                   double& progress);
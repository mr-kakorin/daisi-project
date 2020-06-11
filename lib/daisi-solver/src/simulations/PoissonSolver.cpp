#include "PoissonSolver.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "GridData.h"
#include "MeshContainer2d.h"
#include "MeshContainer3d.h"
//#include "SOR.h"

//#include "mkl_vml.h"


template<class PointType>
template<class Archive>
void PoissonSolver<PointType>::save( Archive &ar, const unsigned int ) const {
	ar & param;
}

template<class PointType>
template<class Archive>
void PoissonSolver<PointType>::load( Archive &ar, const unsigned int ) {
	ar & param;
	solverFlags.resize( 4 );
	SetParameters( param );
}

template<class PointType>
PoissonSolver<PointType>::PoissonSolver() {
	param.resize( 12 );
	param[0] = 1.90;
	param[1] = 0.00011;
	param[2] = 1.70;
	param[3] = 0.01;
	param[4] = 0.00001;
	param[5] = 0.00001;
	param[6] = 5.0;
	param[7] = 1.0;
	param[8] = 0.0; // System solver
	param[9] = 1.0; // time depending flag
	param[10] = 1.0; // space-charge flag
	param[11] = 1.0; // fast initialization
	solverFlags.resize( 4 );
	SetParameters( param );
}

template<class PointType>
void PoissonSolver<PointType>::SetParameters( const std::vector<double> &paramIn ) {
	param = paramIn;
	w = param[0];
	eps_tolerance = param[1];
	w_charge = param[2];
	eps_tolerance_charge = param[3];

	eps_p = param[4];
	eps_h = param[5];

	ChargeSpace = param[6];
	RecalculateParameter = param[7];

	// for (int i = 0; i < solverFlags.size(); i++)
	//	solverFlags[i] = param[8 + i];
}

template<class PointType>
std::vector<double> PoissonSolver<PointType>::GetParameters() {
	/*std::vector<double> result(8);

    result[0] = w;
    result[1] = eps_tolerance;

    result[2] = w_charge;
    result[3] = eps_tolerance_charge;

    result[4] = eps_p;
    result[5] = eps_h;

    result[6] = ChargeSpace;
    result[7] = RecalculateParameter;

    for (int i = 0; i < solverFlags.size(); i++)
            result.push_back(int(solverFlags[i]));*/

	return param;
};

template<class PointType>
bool compNum( BoundaryPoint<PointType> a, BoundaryPoint<PointType> b ) {
	return (a.pnum < b.pnum);
}

int sign( double val ) {
	if (val == 0)
		return 0;
	if (val < 0)
		return -1;
	return 1;
}

// hx1 - left, hx2 - right, hy1 - down, hy2 - up;
double up2dpolar( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hy2 * (hy1 + hy2) * x * x);
}

double down2dpolar( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hy1 * (hy2 + hy1) * x * x);
}

double left2dpolar( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hx1 * (hx1 + hx2));
}

double right2dpolar( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 1 / (x * hx2) + 2 / (hx2 * (hx1 + hx2));
}

double middle2dpolar( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return -1 / (x * hx2) - 2 / (hx1 * hx2) - 2 / (hy1 * hy2 * x * x);
}

double up2d( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hy2 * (hy1 + hy2)); // poisson
}

double down2d( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hy1 * (hy2 + hy1)); // poisson
}

double left2d( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hx1 * (hx1 + hx2)); // poisson
}

double right2d( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hx2 * (hx1 + hx2)); // poisson
}

double middle2d( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return -2 / (hx1 * hx2) - 2 / (hy1 * hy2); // poisson
}

double up2daxs( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hy2 * (hy1 + hy2)); // poisson
}

double down2daxs( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hy1 * (hy1 + hy2)); // poisson
}

double left2daxs( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 2 / (hx1 * (hx1 + hx2)); // poisson
}

double right2daxs( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return 1 / (x * hx2) + 2 / (hx2 * (hx1 + hx2));
}

double middle2daxs( double hx1, double hx2, double hy1, double hy2, double x, double y ) {
	return -1 / (x * hx2) - 2 / (hx1 * hx2) - 2 / (hy1 * hy2);
}

double deep3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
               double y, double z ) {
	return 2 / (hz2 * (hz1 + hz2)); // poisson
}

double outward3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
                  double y, double z ) {
	return 2 / (hz1 * (hz2 + hz1)); // poisson
}

double up3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
             double y, double z ) {
	return 2 / (hy2 * (hy1 + hy2)); // poisson
}

double down3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
               double y, double z ) {
	return 2 / (hy1 * (hy2 + hy1)); // poisson
}

double left3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
               double y, double z ) {
	return 2 / (hx1 * (hx1 + hx2)); // poisson
}

double right3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
                double y, double z ) {
	return 2 / (hx2 * (hx1 + hx2)); // poisson
}

double middle3d( double hx1, double hx2, double hy1, double hy2, double hz1, double hz2, double x,
                 double y, double z ) {
	return -2 / (hx1 * hx2) - 2 / (hy1 * hy2) - 2 / (hz1 * hz2); // poisson
}

template<class PointType>
void PoissonSolver<PointType>::searchBoundaryPointsIntersection(
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {
	boundaryPointsIntersectionCondition.clear();
	boundaryPointsIntersectionDefault.clear();
	DGeo::Edge<PointType> tmpEdge;
	tmpEdge.point1.z = 0;
	tmpEdge.point2.z = 0;
	std::vector<int> boundariesList;
	std::vector<DGeo::Point<PointType>> intersectionPointsTmp;
	std::vector<DGeo::Edge<PointType>> intersectionEdgesTmp;

	BoundaryPointIntersection<PointType> pointTmp;
	boundCondNum = boundaryConditions->PropertyConditionListSize();
	int counter = 0;
	for (int i = 0; i < mesh->meshData.size(); i++) {
		progress = double( counter ) / double( mesh->meshData.size() +
		                                       (mesh->getnVertexX() - 1) * (mesh->getnVertexY() - 1));
		counter++;

		tmpEdge.point1.x = mesh->meshData[i][0].x;
		tmpEdge.point2.x = mesh->meshData[i].back().x;

		tmpEdge.point1.y = mesh->meshData[i][0].y;
		tmpEdge.point2.y = mesh->meshData[i].back().y;

		if (solverFlags[3] == 1) {
			if (mesh->meshData[i][0].isOut == false && mesh->meshData[i].back().isOut == false)
				continue;
		}

		for (int j = 0; j < boundaryConditions->PropertyConditionListSize(); j++) {
			boundariesList = boundaryConditions->GetPropertyConditionsBoundariesList( j );
			for (int k = 0; k < boundariesList.size(); k++) {
				boundaries[boundariesList[k]]->FindAllIntersections(
						eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );
				for (int s = 0; s < intersectionPointsTmp.size(); s++) {
					pointTmp.point = intersectionPointsTmp[s];
					pointTmp.edge = intersectionEdgesTmp[s];
					pointTmp.type = 'd';
					pointTmp.value = boundaryConditions->GetPotentialOffset( j );
					pointTmp.bcnum = j;
					boundaryPointsIntersectionCondition.push_back( pointTmp );
				}
			}
		}

		boundariesList = boundaryConditions->GetDefaultConditionsList();
		for (int k = 0; k < boundariesList.size(); k++) {
			boundaries[boundariesList[k]]->FindAllIntersections(
					eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );

			for (int s = 0; s < intersectionPointsTmp.size(); s++) {
				pointTmp.point = intersectionPointsTmp[s];
				pointTmp.edge = intersectionEdgesTmp[s];
				pointTmp.type = 'n';
				pointTmp.value = 0;
				pointTmp.bcnum = -1;
				boundaryPointsIntersectionDefault.push_back( pointTmp );
			}
		}
	}

	for (int i = 0; i < mesh->getnVertexX() - 1; i++) {
		for (int j = 0; j < mesh->getnVertexY() - 1; j++) {
			progress =
					double( counter ) / double( mesh->meshData.size() +
					                            (mesh->getnVertexX() - 1) * (mesh->getnVertexY() - 1));
			counter++;
			if (mesh->templNumb( i, j ) != -1) {
				if (mesh->templNumb( i, j + 1 ) != -1) {
					tmpEdge.point1.x = mesh->serialMeshData[mesh->templNumb( i, j )].x;
					tmpEdge.point2.x = mesh->serialMeshData[mesh->templNumb( i, j + 1 )].x;

					tmpEdge.point1.y = mesh->serialMeshData[mesh->templNumb( i, j )].y;
					tmpEdge.point2.y = mesh->serialMeshData[mesh->templNumb( i, j + 1 )].y;

					if (solverFlags[3] == 1) {
						if (mesh->serialMeshData[mesh->templNumb( i, j )].isOut == false &&
						    mesh->serialMeshData[mesh->templNumb( i, j + 1 )].isOut == false)
							continue;
					}

					for (int j = 0; j < boundaryConditions->PropertyConditionListSize(); j++) {
						boundariesList = boundaryConditions->GetPropertyConditionsBoundariesList( j );
						for (int k = 0; k < boundariesList.size(); k++) {
							boundaries[boundariesList[k]]->FindAllIntersections(
									eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );

							for (int s = 0; s < intersectionPointsTmp.size(); s++) {
								pointTmp.point = intersectionPointsTmp[s];
								pointTmp.edge = intersectionEdgesTmp[s];
								pointTmp.type = 'd';
								pointTmp.value = boundaryConditions->GetPotentialOffset( j );
								pointTmp.bcnum = j;
								boundaryPointsIntersectionCondition.push_back( pointTmp );
							}
						}
					}

					boundariesList = boundaryConditions->GetDefaultConditionsList();
					for (int k = 0; k < boundariesList.size(); k++) {
						boundaries[boundariesList[k]]->FindAllIntersections(
								eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );

						for (int s = 0; s < intersectionPointsTmp.size(); s++) {
							pointTmp.point = intersectionPointsTmp[s];
							pointTmp.edge = intersectionEdgesTmp[s];
							pointTmp.type = 'n';
							pointTmp.value = 0;
							pointTmp.bcnum = -1;
							boundaryPointsIntersectionDefault.push_back( pointTmp );
						}
					}
				}
			}
		}
	}
}

template<class PointType>
void PoissonSolver<PointType>::searchBoundaryPointsIntersectionPolar(
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {
	boundaryPointsIntersectionCondition.clear();
	boundaryPointsIntersectionDefault.clear();
	DGeo::Edge<PointType> tmpEdge;
	tmpEdge.point1.z = 0;
	tmpEdge.point2.z = 0;
	std::vector<int> boundariesList;
	std::vector<DGeo::Point<PointType>> intersectionPointsTmp;
	std::vector<DGeo::Edge<PointType>> intersectionEdgesTmp;

	BoundaryPointIntersection<PointType> pointTmp;
	boundCondNum = boundaryConditions->PropertyConditionListSize();

	for (int i = 0; i < mesh->meshData.size(); i++) {

		tmpEdge.point1.x = mesh->meshData[i][0].x;
		tmpEdge.point2.x = mesh->meshData[i].back().x;

		tmpEdge.point1.y = mesh->meshData[i][0].y;
		tmpEdge.point2.y = mesh->meshData[i].back().y;

		for (int j = 0; j < boundaryConditions->PropertyConditionListSize(); j++) {
			boundariesList = boundaryConditions->GetPropertyConditionsBoundariesList( j );
			for (int k = 0; k < boundariesList.size(); k++) {
				boundaries[boundariesList[k]]->FindAllIntersections(
						eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );
				for (int s = 0; s < intersectionPointsTmp.size(); s++) {
					pointTmp.point = intersectionPointsTmp[s];
					pointTmp.edge = intersectionEdgesTmp[s];
					pointTmp.type = 'd';
					pointTmp.value = boundaryConditions->GetPotentialOffset( j );
					pointTmp.bcnum = j;
					boundaryPointsIntersectionCondition.push_back( pointTmp );
				}
			}
		}

		boundariesList = boundaryConditions->GetDefaultConditionsList();
		for (int k = 0; k < boundariesList.size(); k++) {
			boundaries[boundariesList[k]]->FindAllIntersections(
					eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );

			for (int s = 0; s < intersectionPointsTmp.size(); s++) {
				pointTmp.point = intersectionPointsTmp[s];
				pointTmp.edge = intersectionEdgesTmp[s];
				pointTmp.type = 'n';
				pointTmp.value = 0;
				pointTmp.bcnum = -1;
				boundaryPointsIntersectionDefault.push_back( pointTmp );
			};
		}
	}

	for (int i = 0; i < mesh->getnVertexX() - 1; i++) {
		for (int j = 0; j < mesh->getnVertexY() - 1; j++) {
			if (mesh->templNumb( i, j ) != -1) {
				if (mesh->templNumb( i, j + 1 ) != -1) {
					tmpEdge.point1.x = mesh->serialMeshData[mesh->templNumb( i, j )].x;
					tmpEdge.point2.x = mesh->serialMeshData[mesh->templNumb( i, j + 1 )].x;

					tmpEdge.point1.y = mesh->serialMeshData[mesh->templNumb( i, j )].y;
					tmpEdge.point2.y = mesh->serialMeshData[mesh->templNumb( i, j + 1 )].y;

					for (int j = 0; j < boundaryConditions->PropertyConditionListSize(); j++) {
						boundariesList = boundaryConditions->GetPropertyConditionsBoundariesList( j );
						for (int k = 0; k < boundariesList.size(); k++) {
							boundaries[boundariesList[k]]->FindAllIntersections(
									eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );

							for (int s = 0; s < intersectionPointsTmp.size(); s++) {
								pointTmp.point = intersectionPointsTmp[s];
								pointTmp.edge = intersectionEdgesTmp[s];
								pointTmp.type = 'd';
								pointTmp.value = boundaryConditions->GetPotentialOffset( j );
								pointTmp.bcnum = j;
								boundaryPointsIntersectionCondition.push_back( pointTmp );
							}
						}
					}

					boundariesList = boundaryConditions->GetDefaultConditionsList();
					for (int k = 0; k < boundariesList.size(); k++) {
						boundaries[boundariesList[k]]->FindAllIntersections(
								eps_p, tmpEdge, intersectionPointsTmp, intersectionEdgesTmp );

						for (int s = 0; s < intersectionPointsTmp.size(); s++) {
							pointTmp.point = intersectionPointsTmp[s];
							pointTmp.edge = intersectionEdgesTmp[s];
							pointTmp.type = 'n';
							pointTmp.value = 0;
							pointTmp.bcnum = -1;
							boundaryPointsIntersectionDefault.push_back( pointTmp );
						}
					}
				}
			}
		}
	}
}

template<class PointType>
DGeo::Point<int>
PoissonSolver<PointType>::closestVertex( const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
                                         DGeo::Point<PointType> point ) {
	DGeo::Point<int> result;
	int i, j;
	i = 0; // x
	j = 0; // y
	int flag = 0;

	// y
	while ((point.y > mesh->meshData[j][i].y) && (j + 1 < mesh->meshData.size())) {
		if (j + 1 < mesh->meshData.size()) {
			if (i < mesh->meshData[j + 1].size()) {
				if (point.y > mesh->meshData[j + 1][i].y) {
					++j;
				} else
					break;
			}
		} else {
			int last = mesh->meshData[j + 1].size() - (mesh->meshData[j].size() - i);
			if (point.y >= mesh->meshData[j + 1][last].y) {
				++j;
				i = last - 1;
			}
		}
	}
	result.y = j;

	if (j + 1 < mesh->meshData.size()) {
		if (i < mesh->meshData[j + 1].size()) {
			if (std::abs( mesh->meshData[j + 1][i].y - point.y ) < std::abs( mesh->meshData[j][i].y - point.y )) {
				result.y = j + 1;
				j++;
			}
		}
	}

	// x
	while ((point.x > mesh->meshData[j][i].x) && (i + 1 < mesh->meshData[j].size())) {
		if (i + 1 < mesh->meshData[j].size()) {
			if (point.x > mesh->meshData[j][i + 1].x) {
				++i;
			} else
				break;
		} else {
			i = mesh->meshData[j].size() - 1;
			break;
		}
	}

	result.x = i;
	if (i + 1 < mesh->meshData[j].size()) {
		if (std::abs( mesh->meshData[j][i + 1].x - point.x ) < std::abs( mesh->meshData[j][i].x - point.x )) {
			result.x = i + 1;
		}
	}

	return result;
}

template<class PointType>
double PoissonSolver<PointType>::dist( DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2 ) {

	double dx, dy;
	dx = vert2.x - vert1.x;
	dy = vert2.y - vert1.y;

	double result = sqrt( dx * dx + dy * dy );
	return result;
}

template<class PointType>
double PoissonSolver<PointType>::distToEdgeX( DGeo::Edge<PointType> edge,
                                              DGeo::Point<PointType> point ) {

	if (polar == false) {
		double y1, y2, x1, x2;
		y1 = edge.point1.y;
		y2 = edge.point2.y;
		x1 = edge.point1.x;
		x2 = edge.point2.x;

		if (edge.IsOnLine( point, eps_p )) {
			double h1 = sqrt((x1 - point.x) * (x1 - point.x) + (y1 - point.y) * (y1 - point.y));
			double h2 = sqrt((x2 - point.x) * (x2 - point.x) + (y2 - point.y) * (y2 - point.y));

			return PointType( std::min( h1, h2 ));
		}
		if (std::abs( x1 - x2 ) < eps_p)
			return std::abs( x1 - point.x );
		else {
			double k = (y2 - y1) / (x2 - x1);
			double b = y1 - k * x1;
			return std::abs( point.x - (point.y - b) / k );
		}
		return 0;
	} else {
		PointType pi = 3.14159265358979;
		PointType k1, k2;
		PointType b;

		PointType y1, y2, x1, x2;
		y1 = edge.point1.y;
		y2 = edge.point2.y;
		x1 = edge.point1.x;
		x2 = edge.point2.x;

		PointType x, y, dx, dy;
		Dmath::Polar2Cartesian( point.x, point.y, x, y );

		if (std::abs( y1 - y2 ) < eps_p) {
			dx = 0;
			dy = y - y1;
		} else if (std::abs( x1 - x2 ) < eps_p) {
			dx = x - x1;
			dy = 0;
		} else {
			k1 = (y2 - y1) / (x2 - x1);
			b = y1 - k1 * x1;

			k2 = y / x;
			dx = x - b / (k2 - k1);
			dy = y - k2 * b / (k2 - k1);
		}
		return std::abs( sqrt( dx * dx + dy * dy ));
	}
	return 0;
}

template<class PointType>
double PoissonSolver<PointType>::distToEdgeY( DGeo::Edge<PointType> edge,
                                              DGeo::Point<PointType> point ) {

	if (polar == false) {
		double y1, y2, x1, x2;
		y1 = edge.point1.y;
		y2 = edge.point2.y;
		x1 = edge.point1.x;
		x2 = edge.point2.x;

		if (edge.IsOnLine( point, eps_p )) {
			double h1 = sqrt((x1 - point.x) * (x1 - point.x) + (y1 - point.y) * (y1 - point.y));
			double h2 = sqrt((x2 - point.x) * (x2 - point.x) + (y2 - point.y) * (y2 - point.y));

			return PointType( std::min( h1, h2 ));
		}

		if (std::abs( y1 - y2 ) < eps_p)
			return std::abs( y1 - point.y );
		else {
			double k = (y2 - y1) / (x2 - x1);
			double b = y1 - k * x1;
			return std::abs( point.y - k * point.x - b );
		}
	} else {
		PointType pi = 3.14159265358979;
		PointType k;
		PointType b;

		PointType y1, y2, x1, x2;
		y1 = edge.point1.y;
		y2 = edge.point2.y;
		x1 = edge.point1.x;
		x2 = edge.point2.x;

		PointType x, y, r;
		PointType x_1, x_2, y_1, y_2;
		double dphi_1, dphi_2;
		Dmath::Polar2Cartesian( point.x, point.y, x, y );
		r = sqrt( x * x + y * y );
		if (std::abs( y1 - y2 ) < eps_p) {
			if (r * r < y1 * y1) {
				return 0;
			}
			x_1 = sqrt( r * r - y1 * y1 );
			x_2 = -sqrt( r * r - y1 * y1 );
			y_1 = y1;
			y_2 = y1;
		} else if (std::abs( x1 - x2 ) < eps_p) {
			if (r * r < x1 * x1) {
				return 0;
			}
			x_1 = x1;
			x_2 = x1;
			y_1 = sqrt( r * r - x1 * x1 );
			y_2 = -sqrt( r * r - x1 * x1 );
		} else {
			k = (y2 - y1) / (x2 - x1);
			b = y1 - k * x1;

			x_1 = (-2 * k * b + sqrt( 4 * k * k * b * b - (4 * (b * b - r * r)) * (k + 1))) /
			      (2 * (k + 1));
			x_2 = (-2 * k * b - sqrt( 4 * k * k * b * b - (4 * (b * b - r * r)) * (k + 1))) /
			      (2 * (k + 1));

			y_1 = k * x_1 + b;
			y_2 = k * x_2 + b;
		}
		PointType r_1, r_2, phi_1, phi_2;

		Dmath::Cartesian2Polar( x_1, y_1, r_1, phi_1 );
		Dmath::Cartesian2Polar( x_2, y_2, r_2, phi_2 );

		dphi_1 = std::abs( phi_1 - point.y );
		dphi_2 = std::abs( phi_2 - point.y );

		if (dphi_1 < dphi_2) {
			return dphi_1;
		} else {
			return dphi_2;
		}
	}
	return 0;
}

template<class PointType>
double PoissonSolver<PointType>::distP2PX( DGeo::Point<PointType> vert1,
                                           DGeo::Point<PointType> vert2 ) {
	return std::abs( vert2.x - vert1.x );
}

template<class PointType>
double PoissonSolver<PointType>::distP2PY( DGeo::Point<PointType> vert1,
                                           DGeo::Point<PointType> vert2 ) {
	return std::abs( vert2.y - vert1.y );
}

template<class PointType>
double PoissonSolver<PointType>::distP2PZ( DGeo::Point<PointType> vert1,
                                           DGeo::Point<PointType> vert2 ) {
	return std::abs( vert2.z - vert1.z );
}

template<class PointType>
double PoissonSolver<PointType>::distP2P( DGeo::Point<PointType> vert1, DGeo::Point<PointType> vert2,
                                          char type ) {
	if (type == 'x') {
		return std::abs( vert2.x - vert1.x );
	} else if (type == 'y') {
		return std::abs( vert2.y - vert1.y );
	} else {
		return std::abs( vert2.z - vert1.z );
	}
}

template<class PointType>
DGeo::Point<PointType> normalVectToEdge( DGeo::Edge<PointType> edge ) {

	DGeo::Point<PointType> result;
	double x0, y0, x, y;

	x0 = edge.point2.x - edge.point1.x;
	y0 = edge.point2.y - edge.point2.y;

	x = 1;
	y = -x0 / y0;

	double normk = 1 / (sqrt( x * x + y * y ));
	result.x = x * normk;
	result.y = y * normk;
	return result;
}

template<class PointType>
BoundaryPoint<PointType> PoissonSolver<PointType>::findBoundaryPoint( int pnumber ) {
	for (int i = 0; i < boundaryPoints.size(); ++i) {
		if (boundaryPoints[i].pnum == pnumber)
			return boundaryPoints[i];
	}
	BoundaryPoint<PointType> bp;
	bp.type = 'f';
	bp.value = 0;
	bp.pnum = -1;
	bp.bcnum = -2;
	return bp;
}

template<class PointType>
void PoissonSolver<PointType>::addNearBoundaryPoint_n( int n ) {
	if (nearBoundaryPoints_n.size() == 0) {
		nearBoundaryPoints_n.push_back( n );
		return;
	}
	int i = 0;
	for (i = 0; i < nearBoundaryPoints_n.size(); ++i) {
		if (n >= nearBoundaryPoints_n[i]) {
			break;
		}
	}
	if (nearBoundaryPoints_n[i] != n) {
		nearBoundaryPoints_n.insert( nearBoundaryPoints_n.begin() + i, n );
	}
}

template<class PointType>
void PoissonSolver<PointType>::addNearBoundaryPoint_all( int n ) {
	if (nearBoundaryPoints_all.size() == 0) {
		nearBoundaryPoints_all.push_back( n );
		return;
	}
	int i = 0;
	for (i = 0; i < nearBoundaryPoints_all.size(); ++i) {
		if (n >= nearBoundaryPoints_all[i]) {
			break;
		}
	}
	if (nearBoundaryPoints_all[i] != n) {
		nearBoundaryPoints_all.insert( nearBoundaryPoints_all.begin() + i, n );
	}
}

template<class PointType>
std::vector<int> PoissonSolver<PointType>::getNearBoundarypoints_n() {
	return nearBoundaryPoints_n;
}

template<class PointType>
std::vector<int> PoissonSolver<PointType>::getNearBoundarypoints_all() {
	return nearBoundaryPoints_all;
}

template<class PointType>
void PoissonSolver<PointType>::addBoundaryPoint( BoundaryPoint<PointType> bpoint ) {
	IBPoint<PointType> intersec_p;
	intersec_p.bpnum = -1;
	intersec_p.ipoint = bpoint.edge.point1;

	if (boundaryPoints.size() == 0) {
		boundaryPoints.push_back( bpoint );
		ibPoints.push_back( intersec_p );
		return;
	}
	int i = 0;
	for (i = 0; i < boundaryPoints.size(); ++i) {
		if (bpoint.pnum >= boundaryPoints[i].pnum) {
			break;
		}
	}

	if (i == boundaryPoints.size()) {
		boundaryPoints.push_back( bpoint );
		ibPoints.push_back( intersec_p );
		return;
	}

	if (boundaryPoints[i].pnum != bpoint.pnum) {
		boundaryPoints.insert( boundaryPoints.begin() + i, bpoint );
		ibPoints.insert( ibPoints.begin() + i, intersec_p );
	}
}

template<class PointType>
void PoissonSolver<PointType>::addBoundaryPoint( BoundaryPoint<PointType> bpoint,
                                                 DGeo::Point<PointType> ipoint ) {

	IBPoint<PointType> intersec_p;
	intersec_p.bpnum = bpoint.pnum;
	intersec_p.ipoint = ipoint;

	if (boundaryPoints.size() == 0) {
		boundaryPoints.push_back( bpoint );
		ibPoints.push_back( intersec_p );
		return;
	}
	int i = 0;
	for (i = 0; i < boundaryPoints.size(); ++i) {
		if (bpoint.pnum >= boundaryPoints[i].pnum) {
			break;
		}
	}

	if (i == boundaryPoints.size()) {
		boundaryPoints.push_back( bpoint );
		ibPoints.push_back( intersec_p );
		return;
	}

	if (boundaryPoints[i].pnum != bpoint.pnum) {
		boundaryPoints.insert( boundaryPoints.begin() + i, bpoint );
		ibPoints.insert( ibPoints.begin() + i, intersec_p );
	}
	if (boundaryPoints[i].pnum == bpoint.pnum) {
		if (((std::abs( ipoint.x - ibPoints[i].ipoint.x ) > eps_p) ||
		     (std::abs( ipoint.y - ibPoints[i].ipoint.y ) > eps_p)) &&
		    (intersec_p.bpnum != -1)) {

			if (polar) {
				PointType x, y;
				Dmath::Polar2Cartesian( ipoint.x, ipoint.y, x, y );
				boundaryPoints[i].edge.point1.x = x;
				boundaryPoints[i].edge.point1.y = y;

				Dmath::Polar2Cartesian( ibPoints[i].ipoint.x, ibPoints[i].ipoint.y, x, y );
				boundaryPoints[i].edge.point2.x = x;
				boundaryPoints[i].edge.point2.y = y;
			} else {
				boundaryPoints[i].edge.point1 = ipoint;
				boundaryPoints[i].edge.point2 = ibPoints[i].ipoint;
			}
		}
	}
}

template<class PointType>
void PoissonSolver<PointType>::noteBoundaryOnMesh(
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {

	if (polar == false) {
		searchBoundaryPointsIntersection( mesh, boundaries, boundaryConditions, progress );
	}

	int stepy, stepx;
	BoundaryPoint<PointType> bpoint, bpoint2;
	DGeo::Point<int> mpoint1, mpoint2;
	DGeo::Point<PointType> point1, point2;
	// for z = z0
	for (int k = 0; k < boundaryPointsIntersectionCondition.size(); ++k) {
		point1 = boundaryPointsIntersectionCondition[k].point;
		mpoint1 = closestVertex( mesh, point1 );
		point2 = mesh->meshData[mpoint1.y][mpoint1.x];
		// steps
		stepx = 0;
		stepy = 0;
		if (std::abs( point2.x - point1.x ) > eps_p) {
			stepx = sign( point1.x - point2.x );
		}
		if (std::abs( point2.y - point1.y ) > eps_p) {
			stepy = sign( point1.y - point2.y );
		}
		if (stepx != 0) {
			mpoint2 = mesh->doStepX( mpoint1, stepx );
		}
		if (stepy != 0) {
			mpoint2 = mesh->doStepY( mpoint1, stepy );
		}

		if ((stepy != 0) && (stepx != 0)) {
			double flags = 0;
		}

		if ((stepy == 0) && (stepx == 0)) {
			mpoint2 = mpoint1;
		}
		int points = 0;
		if ((mesh->meshData[mpoint1.y][mpoint1.x].isOut == true) &&
		    (mesh->meshData[mpoint2.y][mpoint2.x].isOut == true)) {
			DGeo::Edge<PointType> step_edge;
			DGeo::Point<PointType> point;
			step_edge.point1 = mesh->meshData[mpoint1.y][mpoint1.x];
			step_edge.point2 = mesh->meshData[mpoint2.y][mpoint2.x];
			points = 0;
			for (int j = 0; j < boundaryPointsIntersectionDefault.size(); j++) {
				point = boundaryPointsIntersectionDefault[j].point;

				// is the point lay on the edge point1-point2
				if (step_edge.IsOnEdge( point, eps_p ) != 2) {
					points++;
				}
			}
		}
		if ((points % 2 == 0) || (points == 1)) {
			if (mesh->meshData[mpoint1.y][mpoint1.x].isOut == true) {
				bpoint.edge = boundaryPointsIntersectionCondition[k].edge;
				bpoint.type = boundaryPointsIntersectionCondition[k].type;
				bpoint.value = boundaryPointsIntersectionCondition[k].value;
				bpoint.bcnum = boundaryPointsIntersectionCondition[k].bcnum;
				bpoint.pnum = mesh->meshData[mpoint1.y][mpoint1.x].Number;
			} else if (mesh->meshData[mpoint2.y][mpoint2.x].isOut == true) {
				bpoint.edge = boundaryPointsIntersectionCondition[k].edge;
				bpoint.type = boundaryPointsIntersectionCondition[k].type;
				bpoint.value = boundaryPointsIntersectionCondition[k].value;
				bpoint.bcnum = boundaryPointsIntersectionCondition[k].bcnum;
				bpoint.pnum = mesh->meshData[mpoint2.y][mpoint2.x].Number;
			} else {
				bpoint.edge = boundaryPointsIntersectionCondition[k].edge;
				bpoint.type = boundaryPointsIntersectionCondition[k].type;
				bpoint.value = boundaryPointsIntersectionCondition[k].value;
				bpoint.pnum = mesh->meshData[mpoint1.y][mpoint1.x].Number;
				bpoint.bcnum = boundaryPointsIntersectionCondition[k].bcnum;
				mesh->meshData[mpoint1.y][mpoint1.x].isOut = true;
			}
			addBoundaryPoint( bpoint, boundaryPointsIntersectionCondition[k].point );
		}
	}

	for (int k = 0; k < boundaryPointsIntersectionDefault.size(); ++k) {
		point1 = boundaryPointsIntersectionDefault[k].point;
		mpoint1 = closestVertex( mesh, point1 );
		point2 = mesh->meshData[mpoint1.y][mpoint1.x];
		// steps
		stepx = 0;
		stepy = 0;
		if (std::abs( point2.x - point1.x ) > eps_p) {
			stepx = sign( point1.x - point2.x );
		}
		if (std::abs( point2.y - point1.y ) > eps_p) {
			stepy = sign( point1.y - point2.y );
		}
		if (stepx != 0) {
			mpoint2 = mesh->doStepX( mpoint1, stepx );
		}
		if (stepy != 0) {
			mpoint2 = mesh->doStepY( mpoint1, stepy );
		}

		if ((stepy != 0) && (stepx != 0)) {
			double flags = 0;
		}

		if ((stepy == 0) && (stepx == 0)) {
			mpoint2 = mpoint1;
		}
		int points = 0;
		if ((mesh->meshData[mpoint1.y][mpoint1.x].isOut == true) &&
		    (mesh->meshData[mpoint2.y][mpoint2.x].isOut == true)) {
			DGeo::Edge<PointType> step_edge;
			DGeo::Point<PointType> point;
			step_edge.point1 = mesh->meshData[mpoint1.y][mpoint1.x];
			step_edge.point2 = mesh->meshData[mpoint2.y][mpoint2.x];
			points = 0;
			for (int j = 0; j < boundaryPointsIntersectionDefault.size(); j++) {
				point = boundaryPointsIntersectionDefault[j].point;

				// is the point lay on the edge point1-point2
				if (step_edge.IsOnEdge( point, eps_p ) != 2) {
					points++;
				}
			}
		}
		if ((points % 2 == 0) || (points == 1)) {
			if (mesh->meshData[mpoint1.y][mpoint1.x].isOut == true) {
				bpoint.edge = boundaryPointsIntersectionDefault[k].edge;
				bpoint.type = boundaryPointsIntersectionDefault[k].type;
				bpoint.value = boundaryPointsIntersectionDefault[k].value;
				bpoint.bcnum = boundaryPointsIntersectionDefault[k].bcnum;
				bpoint.pnum = mesh->meshData[mpoint1.y][mpoint1.x].Number;
			} else if (mesh->meshData[mpoint2.y][mpoint2.x].isOut == true) {
				bpoint.edge = boundaryPointsIntersectionDefault[k].edge;
				bpoint.type = boundaryPointsIntersectionDefault[k].type;
				bpoint.value = boundaryPointsIntersectionDefault[k].value;
				bpoint.bcnum = boundaryPointsIntersectionDefault[k].bcnum;
				bpoint.pnum = mesh->meshData[mpoint2.y][mpoint2.x].Number;
			} else {
				bpoint.edge = boundaryPointsIntersectionDefault[k].edge;
				bpoint.type = boundaryPointsIntersectionDefault[k].type;
				bpoint.value = boundaryPointsIntersectionDefault[k].value;
				bpoint.bcnum = boundaryPointsIntersectionDefault[k].bcnum;
				bpoint.pnum = mesh->meshData[mpoint1.y][mpoint1.x].Number;
				mesh->meshData[mpoint1.y][mpoint1.x].isOut = true;
			}
			addBoundaryPoint( bpoint, boundaryPointsIntersectionDefault[k].point );
		}
	}

	for (int j = 0; j < mesh->meshData.size(); ++j) {
		for (int i = 0; i < mesh->meshData[j].size(); ++i) {
			if (mesh->meshData[j][i].isOut == true) {

				mpoint1.x = i;
				mpoint1.y = j;
				mpoint1.Number = mesh->meshData[j][i].Number;
				bpoint2 = findBoundaryPoint( mpoint1.Number );
				if (bpoint2.type != 'f') {
					bpoint.edge = bpoint2.edge;
					bpoint.type = bpoint2.type;
					bpoint.value = bpoint2.value;
					bpoint.bcnum = bpoint2.bcnum;
					bpoint.pnum = mesh->meshData[j][i].Number;
				} else {
					mpoint2 = mesh->doStepX( mpoint1, 1 );
					bpoint2 = findBoundaryPoint( mpoint2.Number );
					if (bpoint2.type != 'f') {
						bpoint.edge = bpoint2.edge;
						bpoint.type = bpoint2.type;
						bpoint.value = bpoint2.value;
						bpoint.bcnum = bpoint2.bcnum;
						bpoint.pnum = mesh->meshData[j][i].Number;
					} else {
						mpoint2 = mesh->doStepX( mpoint1, -1 );
						bpoint2 = findBoundaryPoint( mpoint2.Number );
						if (bpoint2.type != 'f') {
							bpoint.edge = bpoint2.edge;
							bpoint.type = bpoint2.type;
							bpoint.value = bpoint2.value;
							bpoint.bcnum = bpoint2.bcnum;
							bpoint.pnum = mesh->meshData[j][i].Number;
						} else {
							mpoint2 = mesh->doStepY( mpoint1, 1 );
							bpoint2 = findBoundaryPoint( mpoint2.Number );
							if (bpoint2.type != 'f') {
								bpoint.edge = bpoint2.edge;
								bpoint.type = bpoint2.type;
								bpoint.value = bpoint2.value;
								bpoint.bcnum = bpoint2.bcnum;
								bpoint.pnum = mesh->meshData[j][i].Number;
							} else {
								mpoint2 = mesh->doStepY( mpoint1, -1 );
								bpoint2 = findBoundaryPoint( mpoint2.Number );
								if (bpoint2.type != 'f') {
									bpoint.edge = bpoint2.edge;
									bpoint.type = bpoint2.type;
									bpoint.value = bpoint2.value;
									bpoint.bcnum = bpoint2.bcnum;
									bpoint.pnum = mesh->meshData[j][i].Number;
								} else {
									bpoint.type = 'n';
									bpoint.value = 0;
									bpoint.bcnum = -1;
									bpoint.pnum = mesh->meshData[j][i].Number;
								}
							}
						}
					}
				}
				addBoundaryPoint( bpoint );
			}
		}
	}

	// i need to sort boundary points
	std::sort( boundaryPoints.begin(), boundaryPoints.end(), compNum<PointType> );
}

template<class PointType>
void PoissonSolver<PointType>::createMatrix( const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
                                             coeff left, coeff right, coeff up, coeff down,
                                             coeff middle ) {

	int size = mesh->serialMeshData.size();
	int m = 0;
	systemSize = size;
	int step;

	int tempj;

	xCharge.resize( mesh->serialMeshData.size());
	linVect.resize( mesh->serialMeshData.size());
	x.resize( mesh->serialMeshData.size());

	for (int j = 0; j < size; ++j) {
		linVect[j] = 0; // right hand side
		x[j] = 1; // start solution for iteration method
		xCharge[j] = 0;
	}
	templ.resize( mesh->serialMeshData.size());
	Vtmp.resize( mesh->serialMeshData.size());

	for (int j = 0; j < boundaryPoints.size(); ++j) {
		x[boundaryPoints[j].pnum] = boundaryPoints[j].value;
	}

	// loop for all points, building linear system
	DGeo::Point<int> tmp1, tmp2, tmp3, tmp4;
	int row;
	// int step, stepy, stepx;
	double h, hx, hy, a, b;
	int flag = 0;

	double h_up, h_down, h_left, h_right;
	int flag_h = 0;
	double h_bp_value, bp_bcnum;

	hx = mesh->h1;
	hy = mesh->h2;

	nonZeroLinVect.clear();
	NonZeroLinVect<PointType> tmp_nonZero;

	double x_, y_, hx_, hy_;
	DGeo::Point<PointType> norm;
	BoundaryPoint<PointType> bp;
	int c_i = -1; // current index

	int kk = 0;

	for (int j = 0; j < mesh->meshData.size(); ++j) {
		tmp1.y = j;
		for (int i = 0; i < mesh->meshData[j].size(); ++i) {
			tmp1.x = i;
			if (mesh->meshData[j][i].isOut == false) {
				row = mesh->meshData[tmp1.y][tmp1.x].Number;
				middles.push_back( mesh->meshData[tmp1.y][tmp1.x].Number );
				n2k.push_back( kk );
				kk++;

				c_i++;

				flag = 0;
				flag_h = 0;
				// find distance to edges
				// left
				tmp2 = mesh->doStepX( tmp1, -1 );
				if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
					bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
					h_left = distToEdgeX( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
					if (h_left < eps_h) {
						flag_h = 1;
						h_bp_value = bp.value;
						bp_bcnum = bp.bcnum;
					}
				} else {
					h_left =
							distP2PX( mesh->meshData[tmp1.y][tmp1.x], mesh->meshData[tmp2.y][tmp2.x] );
				}
				// right
				tmp2 = mesh->doStepX( tmp1, 1 );
				if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
					bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
					h_right = distToEdgeX( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
					if (h_right < eps_h) {
						flag_h = 1;
						h_bp_value = bp.value;
						bp_bcnum = bp.bcnum;
					}
				} else {
					h_right =
							distP2PX( mesh->meshData[tmp1.y][tmp1.x], mesh->meshData[tmp2.y][tmp2.x] );
				}
				// up
				tmp2 = mesh->doStepY( tmp1, 1 );
				if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
					bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
					h_up = distToEdgeY( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
					if (h_up < eps_h) {
						flag_h = 1;
						h_bp_value = bp.value;
						bp_bcnum = bp.bcnum;
					}
				} else {
					h_up = distP2PY( mesh->meshData[tmp1.y][tmp1.x], mesh->meshData[tmp2.y][tmp2.x] );
				}
				// down
				tmp2 = mesh->doStepY( tmp1, -1 );
				if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
					bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
					h_down = distToEdgeY( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
					if (h_down < eps_h) {
						flag_h = 1;
						h_bp_value = bp.value;
						bp_bcnum = bp.bcnum;
					}
				} else {
					h_down =
							distP2PY( mesh->meshData[tmp1.y][tmp1.x], mesh->meshData[tmp2.y][tmp2.x] );
				}

				if (flag_h == 1) {
					// sys.setElement(row, row, 1);
					c_middles.push_back( 1 );
					linVect[row] = linVect[row] + h_bp_value;
					tmp_nonZero.bcnum = bp.bcnum;
					tmp_nonZero.coef = -1;
					tmp_nonZero.row = row;
					nonZeroLinVect.push_back( tmp_nonZero );

					addNearBoundaryPoint_all( mesh->meshData[j][i].Number );
					c_rights.push_back( 0 );
					c_lefts.push_back( 0 );
					c_ups.push_back( 0 );
					c_downs.push_back( 0 );

					rights.push_back( 0 );
					lefts.push_back( 0 );
					ups.push_back( 0 );
					downs.push_back( 0 );
				} else {
					x_ = mesh->meshData[j][i].x;
					y_ = mesh->meshData[j][i].y;
					// sys.setElement(row, row, middle(h_left, h_right, h_down, h_up, x_, y_));
					c_middles.push_back( middle( h_left, h_right, h_down, h_up, x_, y_ ));
					// right
					step = 1;
					h = hx;
					tmp2 = mesh->doStepX( tmp1, step );
					// col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
					rights.push_back( mesh->meshData[tmp2.y][tmp2.x].Number );
					if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
						flag = 1;
						bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
						if (bp.type == 'n') {
							/*norm = normalVectToEdge(bp.edge);
                            stepy = sign(norm.y);
                            tmp3 = mesh->doStepY(tmp2, stepy);
                            tmp4 = mesh->doStepY(tmp1, stepy);
                            a = norm.x;
                            b = norm.y;
                            hx1 = distToEdgeX(bp.edge,mesh->meshData[tmp1.y][tmp1.x]);
                            col2 = numbers[mesh->meshData[tmp4.y][tmp4.x].Number];
                            //col3 = numbers[mesh->meshData[tmp3.y][tmp3.x].Number];
                            */
							// 0 tmp1 row_row			+right(hx1,hy,x_,y_)
							// 1 tmp2 is out				linVect +bp.value()
							// 5 tmp3 is out/row_col3
							// 2 tmp4 row_col2
							// sys.setElement(row, row, sys.getElement(row, row) + right(h_left,
							// h_right, h_down, h_up, x_, y_));
							c_middles[c_i] =
									c_middles[c_i] + right( h_left, h_right, h_down, h_up, x_, y_ );
							c_rights.push_back( 0 );
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] = linVect[row] -
								               h * bp.value * step *
								               right( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef =
									h * step * right( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							addNearBoundaryPoint_n( mesh->meshData[tmp1.y][tmp1.x].Number );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
						if (bp.type == 'd') {
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] =
										linVect[row] -
										bp.value * right( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef = right( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							c_rights.push_back( 0 );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
					} else {
						// sys.setElement(row, col, right(h_left, h_right, h_down, h_up, x_, y_));
						c_rights.push_back( right( h_left, h_right, h_down, h_up, x_, y_ ));
					}
					// left
					step = -1;
					h = hx;
					tmp2 = mesh->doStepX( tmp1, step );
					// col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
					lefts.push_back( mesh->meshData[tmp2.y][tmp2.x].Number );
					if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
						flag = 1;
						bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
						if (bp.type == 'n') {
							// sys.setElement(row, row, sys.getElement(row, row) + left(h_left,
							// h_right, h_down, h_up, x_, y_));
							c_lefts.push_back( 0 );
							c_middles[c_i] =
									c_middles[c_i] + left( h_left, h_right, h_down, h_up, x_, y_ );
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] = linVect[row] -
								               h * bp.value * step *
								               left( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef =
									h * step * left( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							addNearBoundaryPoint_n( mesh->meshData[tmp1.y][tmp1.x].Number );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}

						if (bp.type == 'd') {
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] =
										linVect[row] -
										bp.value * left( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef = left( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							c_lefts.push_back( 0 );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
					} else {
						// sys.setElement(row, col, left(h_left, h_right, h_down, h_up, x_, y_));
						c_lefts.push_back( left( h_left, h_right, h_down, h_up, x_, y_ ));
					}
					// up
					h = hy;
					step = 1;
					tmp2 = mesh->doStepY( tmp1, step );
					// col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
					ups.push_back( mesh->meshData[tmp2.y][tmp2.x].Number );
					if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
						flag = 1;
						bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
						// c_ups.push_back(0);
						if (bp.type == 'n') {
							// sys.setElement(row, row, sys.getElement(row, row) + up(h_left,
							// h_right, h_down, h_up, x_, y_));
							c_middles[c_i] =
									c_middles[c_i] + up( h_left, h_right, h_down, h_up, x_, y_ );
							c_ups.push_back( 0 );
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] =
										linVect[row] -
										h * bp.value * step * up( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef = h * step * up( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							addNearBoundaryPoint_n( mesh->meshData[tmp1.y][tmp1.x].Number );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
						if (bp.type == 'd') {
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] = linVect[row] -
								               bp.value * up( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef = up( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							c_ups.push_back( 0 );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
					} else {
						// sys.setElement(row, col, up(h_left, h_right, h_down, h_up, x_, y_));
						c_ups.push_back( up( h_left, h_right, h_down, h_up, x_, y_ ));
					}
					// down
					h = hy;
					step = -1;
					tmp2 = mesh->doStepY( tmp1, step );
					// col = numbers[mesh->meshData[tmp2.y][tmp2.x].Number];
					downs.push_back( mesh->meshData[tmp2.y][tmp2.x].Number );
					if (mesh->meshData[tmp2.y][tmp2.x].isOut == true) {
						flag = 1;
						bp = findBoundaryPoint( mesh->meshData[tmp2.y][tmp2.x].Number );
						// c_downs.push_back(0);
						if (bp.type == 'n') {
							// sys.setElement(row, row, sys.getElement(row, row) + down(h_left,
							// h_right, h_down, h_up, x_, y_));
							c_middles[c_i] =
									c_middles[c_i] + down( h_left, h_right, h_down, h_up, x_, y_ );
							c_downs.push_back( 0 );
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] = linVect[row] -
								               h * bp.value * step *
								               down( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef =
									h * step * down( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							addNearBoundaryPoint_n( mesh->meshData[tmp1.y][tmp1.x].Number );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
						if (bp.type == 'd') {
							if (std::abs( bp.value ) > 0.00001) {
								linVect[row] =
										linVect[row] -
										bp.value * down( h_left, h_right, h_down, h_up, x_, y_ );
							}
							tmp_nonZero.bcnum = bp.bcnum;
							tmp_nonZero.coef = down( h_left, h_right, h_down, h_up, x_, y_ );
							tmp_nonZero.row = row;
							nonZeroLinVect.push_back( tmp_nonZero );

							c_downs.push_back( 0 );
							addNearBoundaryPoint_all( mesh->meshData[tmp1.y][tmp1.x].Number );
						}
					} else {
						// sys.setElement(row, col, down(h_left, h_right, h_down, h_up, x_, y_));
						c_downs.push_back( down( h_left, h_right, h_down, h_up, x_, y_ ));
					}
				}
			} else
				n2k.push_back( -1 );
		}
	}

	std::sort( nearBoundaryPoints_all.begin(), nearBoundaryPoints_all.end());
	std::sort( nearBoundaryPoints_n.begin(), nearBoundaryPoints_n.end());

	kZeros.resize( middles.size());
	for (int k = 0; k < middles.size(); ++k) {
		int flag = 0;

		if (c_middles[k] == 0)
			flag++;

		if (c_lefts[k] == 0)
			flag++;

		if (c_rights[k] == 0)
			flag++;

		if (c_ups[k] == 0)
			flag++;

		if (c_downs[k] == 0)
			flag++;

		//	if (flag>1)
		//		kZeros[k] = -1;
		//	else
		//		kZeros[k] = 0;
	}
}

template<class PointType>
void PoissonSolver<PointType>::createMatrix3d(
		const std::shared_ptr<MeshContainer3d<PointType>> &mesh, coeff3d left, coeff3d right,
		coeff3d up, coeff3d down, coeff3d outward, coeff3d deep, coeff3d middle ) {

}

template<class PointType>
void PoissonSolver<PointType>::InitSolver(
		const std::shared_ptr<GridData2d<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {

	boundaryPoints.clear();
	nearBoundaryPoints_all.clear();
	nearBoundaryPoints_n.clear();

	middles.clear();
	ups.clear();
	downs.clear();
	c_rights.clear();
	c_middles.clear();
	lefts.clear();
	rights.clear();
	c_downs.clear();
	c_lefts.clear();
	c_ups.clear();
	boundaryPointsIntersectionCondition.clear();
	boundaryPointsIntersectionDefault.clear();

	double hx = mesh->h1;
	double hy = mesh->h2;
	// eps_p = mesh->h1*0.001;
	// eps_h = mesh->h1*0.01;

	int k = 0;

	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {
			mesh->meshData[i][j].Number = k;
			++k;
		}
	}
	polar = false;
	noteBoundaryOnMesh( mesh, boundaries, boundaryConditions, progress );
	createMatrix( mesh, left2d, right2d, up2d, down2d, middle2d );
}

template<class PointType>
void PoissonSolver<PointType>::InitSolver(
		const std::shared_ptr<GridData2daxs<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {

	boundaryPoints.clear();
	ibPoints.clear();

	nearBoundaryPoints_all.clear();
	nearBoundaryPoints_n.clear();
	middles.clear();

	ups.clear();
	downs.clear();
	c_rights.clear();
	c_middles.clear();
	lefts.clear();
	rights.clear();
	c_downs.clear();
	c_lefts.clear();
	c_ups.clear();

	gridData->ZeroingFields();

	boundaryPointsIntersectionCondition.clear();
	boundaryPointsIntersectionDefault.clear();

	double hr = mesh->h1;
	double hz = mesh->h2;
	// eps_p = mesh->h1*0.001;
	// eps_h = mesh->h1*0.01;

	int k = 0;

	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {
			mesh->meshData[i][j].Number = k;
			++k;
		}
	}
	polar = false;
	noteBoundaryOnMesh( mesh, boundaries, boundaryConditions, progress );
	createMatrix( mesh, left2daxs, right2daxs, up2daxs, down2daxs, middle2daxs );
}

template<class PointType>
void PoissonSolver<PointType>::InitSolver(
		const std::shared_ptr<GridData2dpolar<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {

	boundaryPoints.clear();
	nearBoundaryPoints_all.clear();
	nearBoundaryPoints_n.clear();

	middles.clear();
	ups.clear();
	downs.clear();
	c_rights.clear();
	c_middles.clear();
	lefts.clear();
	rights.clear();
	c_downs.clear();
	c_lefts.clear();
	c_ups.clear();
	boundaryPointsIntersectionCondition.clear();
	boundaryPointsIntersectionDefault.clear();

	mesh_copy = std::shared_ptr<MeshContainer2d<PointType>>( new MeshContainer2d<PointType>());

	mesh_copy->meshData.clear();
	mesh_copy->serialMeshData.clear();
	mesh_copy->neighbourContainer.clear();

	int k = 0;
	PointType r, phi;
	// eps_p = mesh->h1*0.001;
	// eps_h = mesh->h1*0.01;

	std::vector<DGeo::Point<PointType>> tmp;

	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {
			tmp.push_back( mesh->meshData[i][j] );
		}
		mesh_copy->meshData.push_back( tmp );
		tmp.clear();
	}

	for (int i = 0; i < mesh->serialMeshData.size(); ++i) {
		mesh_copy->serialMeshData.push_back( mesh->serialMeshData[i] );
		mesh_copy->neighbourContainer.push_back( mesh->neighbourContainer[i] );
	}
	mesh_copy->h1 = mesh->h1;
	mesh_copy->h2 = mesh->h2;

	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {
			mesh_copy->meshData[i][j].Number = k;
			mesh->meshData[i][j].Number = k;
			Dmath::Cartesian2Polar( mesh->meshData[i][j].x, mesh->meshData[i][j].y, r, phi );
			mesh_copy->meshData[i][j].x = r;
			mesh_copy->meshData[i][j].y = phi;
			++k;
		}
	}

	DGeo::Point<PointType> point;
	searchBoundaryPointsIntersectionPolar( mesh, boundaries, boundaryConditions, progress );
	for (int i = 0; i < boundaryPointsIntersectionCondition.size(); ++i) {
		point = boundaryPointsIntersectionCondition[i].point;
		Dmath::Cartesian2Polar( point.x, point.y, r, phi );
		boundaryPointsIntersectionCondition[i].point.x = r;
		boundaryPointsIntersectionCondition[i].point.y = phi;
	}

	for (int i = 0; i < boundaryPointsIntersectionDefault.size(); ++i) {
		point = boundaryPointsIntersectionDefault[i].point;
		Dmath::Cartesian2Polar( point.x, point.y, r, phi );
		boundaryPointsIntersectionDefault[i].point.x = r;
		boundaryPointsIntersectionDefault[i].point.y = phi;
	}

	polar = true;
	noteBoundaryOnMesh( mesh_copy, boundaries, boundaryConditions, progress );
	createMatrix( mesh_copy, left2dpolar, right2dpolar, up2dpolar, down2dpolar, middle2dpolar );
}

template<class PointType>
void PoissonSolver<PointType>::InitSolver(
		const std::shared_ptr<GridData3d<PointType>> &gridData,
		const std::shared_ptr<MeshContainer3d<PointType>> &mesh,
		const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> &boundaries,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double &progress ) {

	boundaryPoints.clear();
	nearBoundaryPoints_all.clear();
	nearBoundaryPoints_n.clear();

	middles.clear();
	ups.clear();
	downs.clear();
	c_rights.clear();
	c_middles.clear();
	lefts.clear();
	rights.clear();
	c_downs.clear();
	c_lefts.clear();
	c_ups.clear();
	c_outwards.clear();
	c_deeps.clear();
	outwards.clear();
	deeps.clear();
	//	nSliceBPoints.clear();
	boundaryPointsIntersectionCondition.clear();
	boundaryPointsIntersectionDefault.clear();

	// eps_p = mesh->h1*0.0001;
	// eps_h = mesh->h1*0.01;
	nSlicePoints = mesh->serialMeshData.size() / mesh->mesh.size();

	int k = 0;

	for (int b = 0; b < mesh->mesh.size(); ++b) {
		for (int i = 0; i < mesh->mesh[b].meshData.size(); ++i) {
			for (int j = 0; j < mesh->mesh[b].meshData[i].size(); ++j) {
				mesh->mesh[b].meshData[i][j].Number = k;
				++k;
			}
		}
	}
	polar = false;

	BoundaryPoint<PointType> bp;
	DGeo::Point<int> tmp, tmp1;
	int count = 0;
	noteBoundaryOnMesh( std::shared_ptr<MeshContainer2d<PointType>>( &mesh->mesh[0] ), boundaries,
	                    boundaryConditions, progress );

	for (int b = 1; b < mesh->mesh.size(); ++b) {
		tmp.z = b;
		count = 0;
		for (int i = 0; i < mesh->mesh[b].meshData.size(); ++i) {
			tmp.y = i;
			for (int j = 0; j < mesh->mesh[b].meshData[i].size(); ++j) {
				tmp.x = j;
				tmp1 = mesh->doStepZ( tmp, -1 );
				if (mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x].isOut == true) {
					count++;
					mesh->mesh[tmp.z].meshData[tmp.y][tmp.x].isOut = true;
					bp = findBoundaryPoint( mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x].Number );
					bp.pnum = mesh->mesh[tmp.z].meshData[tmp.y][tmp.x].Number;
					addBoundaryPoint( bp );
				}
			}
		}
		std::sort( boundaryPoints.begin(), boundaryPoints.end(), compNum<PointType> );
		//	nSliceBPoints.push_back(count);
	}

	BoundaryPoint<PointType> tmp_bp;
	tmp_bp.type = 'n';
	tmp_bp.value = 0;
	tmp_bp.bcnum = -1;
	// note front and back sides of mesh like neumann boundary
	for (int i = 0; i < mesh->mesh[0].meshData.size(); ++i) {
		for (int j = 0; j < mesh->mesh[0].meshData[i].size(); ++j) {
			mesh->mesh[0].meshData[i][j].isOut = true;
			tmp_bp.pnum = mesh->mesh[0].meshData[i][j].Number;
			addBoundaryPoint( tmp_bp );
			mesh->mesh[mesh->mesh.size() - 1].meshData[i][j].isOut = true;
			tmp_bp.pnum = mesh->mesh[mesh->mesh.size() - 1].meshData[i][j].Number;
			addBoundaryPoint( tmp_bp );
		}
	}

	createMatrix3d( mesh, left3d, right3d, up3d, down3d, outward3d, deep3d, middle3d );
}

template<class PointType>
void PoissonSolver<PointType>::FieldSimulate(
		const std::shared_ptr<GridData2d<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t, double &progress ) {

	solve( gridData->Getrho(), gridData->GetVA(), mesh, boundaryConditions, t, progress );

	fieldCalculation_Simple_FDE( gridData->GetVA(), gridData->GetVA(), mesh, 'x',
	                             gridData->Get_ExA());
	fieldCalculation_Simple_FDE( gridData->GetVA(), gridData->GetVA(), mesh, 'y',
	                             gridData->Get_EyA());
}

template<class PointType>
void PoissonSolver<PointType>::FieldSimulate(
		const std::shared_ptr<GridData2daxs<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t, double &progress ) {

	solve( gridData->Getrho(), gridData->GetVA(), mesh, boundaryConditions, t, progress );

	fieldCalculation_Simple_FDE( gridData->GetVA(), gridData->GetVA(), mesh, 'x',
	                             gridData->Get_ErA());
	fieldCalculation_Simple_FDE( gridData->GetVA(), gridData->GetVA(), mesh, 'y',
	                             gridData->Get_EzA());
}

template<class PointType>
void PoissonSolver<PointType>::CreateSolutionSpace(
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh, std::vector<int> &nonZeros ) {
	nonZerosV = nonZeros;
	indexes.resize( nonZeros.size());

	for (int i = 0; i < nonZeros.size(); i++) {
		indexes[i] = mesh->linearToTemplNumb[nonZeros[i]];
		templ[nonZeros[i]] = 1;
	}

	int s = ChargeSpace;

	for (int i = 0; i < nonZeros.size(); i++) {
		for (int j = -s; j < 0; j++) {
			int newIndex = mesh->templNumb.getElem( indexes[i], j, 0 );
			if (newIndex < 0)
				break;
			if (!templ[newIndex]) {
				nonZerosV.push_back( newIndex );
				templ[newIndex] = 1;
				indexes.push_back( indexes[i] + j );
			} else
				break;
		}

		for (int j = 1; j <= s; j++) {
			int newIndex = mesh->templNumb.getElem( indexes[i], j, 0 );
			if (newIndex < 0)
				break;
			if (!templ[newIndex]) {
				nonZerosV.push_back( newIndex );
				templ[newIndex] = 1;
				indexes.push_back( indexes[i] + j );
			} else
				break;
		}
	}

	for (int i = 0; i < indexes.size(); i++) {
		for (int j = -s; j < 0; j++) {
			int newIndex = mesh->templNumb.getElem( indexes[i], 0, j );
			if (newIndex < 0)
				break;
			if (!templ[newIndex]) {
				nonZerosV.push_back( newIndex );
				templ[newIndex] = 1;
			} else
				break;
		}

		for (int j = 1; j <= s; j++) {
			int newIndex = mesh->templNumb.getElem( indexes[i], 0, j );
			if (newIndex < 0)
				break;
			if (!templ[newIndex]) {
				nonZerosV.push_back( newIndex );
				templ[newIndex] = 1;
			} else
				break;
		}
	}

	Vtmp = xCharge;
	memset( &xCharge[0], 0, xCharge.size() * sizeof( xCharge[0] ));
	for (int i = 0; i < nonZerosV.size(); i++)
		xCharge[nonZerosV[i]] = Vtmp[nonZerosV[i]];
}

template<class PointType>
void PoissonSolver<PointType>::FieldSimulateCharge(
		const std::shared_ptr<GridData2daxs<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t,
		std::vector<int> &nonZeros, int step ) {
	if (step % RecalculateParameter != 0 || !nonZeros.size())
		return;

	// if (solverFlags[2] == 0)
	//	return;

	CreateSolutionSpace( mesh, nonZeros );

	if (xCharge.size() == 0)
		xCharge.resize( mesh->serialMeshData.size());

	solveCharge( gridData->Getrho(), gridData->GetVCharge(), mesh, boundaryConditions, t, nonZerosV );

	fieldCalculation_Simple_FDE( gridData->Getrho(), gridData->GetVCharge(), mesh, 'x',
	                             gridData->Get_ErCol());
	fieldCalculation_Simple_FDE( gridData->Getrho(), gridData->GetVCharge(), mesh, 'y',
	                             gridData->Get_EzCol());

	/*	for (int i = 0; i < gridData->GetV().size(); i++)
            {
                    if (std::isnan(gridData->ErCol[i]) || std::isinf(gridData->ErCol[i]))
                            gridData->ErCol[i] = 0.0;
            }


            for (int i = 0; i < gridData->GetV().size(); i++)
            {
                    if (std::isnan(gridData->EzCol[i]) || std::isinf(gridData->EzCol[i]))
                            gridData->EzCol[i] = 0.0;
            }*/

	for (int i = 0; i < nonZerosV.size(); i++)
		templ[nonZerosV[i]] = 0;
}

template<class PointType>
void PoissonSolver<PointType>::FieldSimulateCharge(
		const std::shared_ptr<GridData2d<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t,
		std::vector<int> &nonZeros, int step ) {
	if (step % RecalculateParameter != 0 || !nonZeros.size())
		return;

	// if (solverFlags[2] == 0)
	//	return;

	CreateSolutionSpace( mesh, nonZeros );

	if (xCharge.size() == 0)
		xCharge.resize( mesh->serialMeshData.size());

	solveCharge( gridData->Getrho(), gridData->GetVCharge(), mesh, boundaryConditions, t, nonZerosV );

	fieldCalculation_Simple_FDE( gridData->Getrho(), gridData->GetVCharge(), mesh, 'x',
	                             gridData->Get_ExCol());
	fieldCalculation_Simple_FDE( gridData->Getrho(), gridData->GetVCharge(), mesh, 'y',
	                             gridData->Get_EyCol());

	/*	for (int i = 0; i < gridData->GetV().size(); i++)
    {
    if (std::isnan(gridData->ErCol[i]) || std::isinf(gridData->ErCol[i]))
    gridData->ErCol[i] = 0.0;
    }


    for (int i = 0; i < gridData->GetV().size(); i++)
    {
    if (std::isnan(gridData->EzCol[i]) || std::isinf(gridData->EzCol[i]))
    gridData->EzCol[i] = 0.0;
    }*/

	for (int i = 0; i < nonZerosV.size(); i++)
		templ[nonZerosV[i]] = 0;
}

template<class PointType>
void PoissonSolver<PointType>::FieldSimulate(
		const std::shared_ptr<GridData2dpolar<PointType>> &gridData,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t, double &progress ) {

	solve( gridData->Getrho(), gridData->GetVA(), mesh_copy, boundaryConditions, t, progress );

	fieldCalculation_Simple_FDE( gridData->GetVA(), gridData->GetV(), mesh_copy, 'x',
	                             gridData->Get_Er());
	fieldCalculation_Simple_FDE( gridData->GetVA(), gridData->GetV(), mesh_copy, 'y',
	                             gridData->Get_Ephi());
}

template<class PointType>
void PoissonSolver<PointType>::FieldSimulate(
		const std::shared_ptr<GridData3d<PointType>> &gridData,
		const std::shared_ptr<MeshContainer3d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t, double &progress ) {
	solve_3d( gridData->Getrho(), gridData->GetV(), mesh, boundaryConditions, t );
	fieldCalculation_Simple_FDE( gridData->GetV(), gridData->GetV(), mesh, 'x', gridData->Get_Ex());
	fieldCalculation_Simple_FDE( gridData->GetV(), gridData->GetV(), mesh, 'y', gridData->Get_Ey());
	fieldCalculation_Simple_FDE( gridData->GetV(), gridData->GetV(), mesh, 'z', gridData->Get_Ez());
}

template<class PointType>
void PoissonSolver<PointType>::fieldCalculation_Simple_FDE(
		std::vector<PointType> &rho, std::vector<PointType> &V,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh, char type, std::vector<PointType> &E ) {

	double phi0, phi1, phi2, phi3;
	DGeo::Point<int> tmp, tmp1, tmp2, tmp3;
	BoundaryPoint<PointType> bp;
	double h_, h1, h2, h3, delta, H;

	int flagE = 0;
	int m = -1;
	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {

			++m;

			if (!rho[m]) {
				E[m] = 0;
				continue;
			}

			tmp.x = j;
			tmp.y = i;

			flagE = 0;
			tmp1 = doStep( mesh, type, tmp, -1 );
			tmp2 = doStep( mesh, type, tmp, 1 );

			phi1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
			phi0 = V[mesh->meshData[tmp.y][tmp.x].Number];
			phi2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];

			h1 = distP2P( mesh->meshData[tmp1.y][tmp1.x], mesh->meshData[tmp.y][tmp.x], type );
			h2 = distP2P( mesh->meshData[tmp2.y][tmp2.x], mesh->meshData[tmp.y][tmp.x], type );

			if ((tmp.x == tmp1.x) && (tmp1.y == tmp.y) ||
			    (mesh->meshData[tmp.y][tmp.x].isOut == true) &&
			    (mesh->meshData[tmp1.y][tmp1.x].isOut == true) &&
			    (mesh->meshData[tmp2.y][tmp2.x].isOut == false)) {
				// 1
				// tmp,tmp1 - 1, tmp2 - 2
				bp = findBoundaryPoint( mesh->meshData[tmp.y][tmp.x].Number );
				tmp3 = doStep( mesh, type, tmp2, 1 );
				phi3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				h3 = distP2P( mesh->meshData[tmp2.y][tmp2.x], mesh->meshData[tmp3.y][tmp3.x], type );
				H = h2 + h3;
				delta = h3 / h2;
				E[m] = -(-(2 + delta) * phi0 + (1 + delta) * (1 + delta) * phi2 / delta -
				         phi3 / delta) /
				       H;
			} else if ((tmp.x == tmp2.x) && (tmp.y == tmp2.y) ||
			           (mesh->meshData[tmp.y][tmp.x].isOut == true) &&
			           (mesh->meshData[tmp2.y][tmp2.x].isOut == true) &&
			           (mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
				// 5
				// tmp,tmp2 - 5, tmp1 - 4
				tmp3 = doStep( mesh, type, tmp1, -1 );
				phi3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				h3 = distP2P( mesh->meshData[tmp1.y][tmp1.x], mesh->meshData[tmp3.y][tmp3.x], type );
				H = h1 + h3;
				delta = h1 / h3;
				E[m] = -(delta * phi3 - (1 + delta) * (1 + delta) * phi1 / delta +
				         (2 + delta) * phi0 / delta) /
				       H;
			} else {
				// 3
				// tmp1 - 2, tmp2 - 4
				H = h1 + h2;
				delta = h2 / h1;
				E[m] = -(-delta * phi1 + (delta * delta - 1) * phi0 / delta + phi2 / delta) / H;
			}
			if (std::isnan( E[m] ))
				E[m] = 0;
		}
	}
}

template<class PointType>
void PoissonSolver<PointType>::fieldCalculation_Simple_FDE(
		std::vector<PointType> &rho, std::vector<PointType> &V,
		const std::shared_ptr<MeshContainer3d<PointType>> &mesh, char type, std::vector<PointType> &E ) {

	double phi0, phi1, phi2, phi3;
	DGeo::Point<int> tmp, tmp1, tmp2, tmp3;
	BoundaryPoint<PointType> bp;
	double h_, h1, h2, h3, delta, H;

	int flagE = 0;
	int m = -1;
	for (int b = 0; b < mesh->mesh.size(); ++b) {
		for (int i = 0; i < mesh->mesh[b].meshData.size(); ++i) {
			for (int j = 0; j < mesh->mesh[b].meshData[i].size(); ++j) {

				++m;

				if (!rho[m]) {
					E[m] = 0;
					continue;
				}

				tmp.x = j;
				tmp.y = i;
				tmp.z = b;

				flagE = 0;
				tmp1 = doStep( mesh, type, tmp, -1 );
				tmp2 = doStep( mesh, type, tmp, 1 );

				phi1 = V[mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x].Number];
				phi0 = V[mesh->mesh[tmp.z].meshData[tmp.y][tmp.x].Number];
				phi2 = V[mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].Number];

				h1 = distP2P( mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x],
				              mesh->mesh[tmp.z].meshData[tmp.y][tmp.x], type );
				h2 = distP2P( mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x],
				              mesh->mesh[tmp.z].meshData[tmp.y][tmp.x], type );

				if ((tmp.x == tmp1.x) && (tmp1.y == tmp.y) && (tmp1.z == tmp.z) ||
				    (mesh->mesh[tmp.z].meshData[tmp.y][tmp.x].isOut == true) &&
				    (mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x].isOut == true) &&
				    (mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].isOut == false)) {
					// 1
					// tmp,tmp1 - 1, tmp2 - 2
					bp = findBoundaryPoint( mesh->mesh[tmp.z].meshData[tmp.y][tmp.x].Number );
					tmp3 = doStep( mesh, type, tmp2, 1 );
					phi3 = V[mesh->mesh[tmp3.z].meshData[tmp3.y][tmp3.x].Number];
					h3 = distP2P( mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x],
					              mesh->mesh[tmp3.z].meshData[tmp3.y][tmp3.x], type );
					H = h2 + h3;
					delta = h3 / h2;
					E[m] = -(-(2 + delta) * phi0 + (1 + delta) * (1 + delta) * phi2 / delta -
					         phi3 / delta) /
					       H;
				} else if ((tmp.x == tmp2.x) && (tmp.y == tmp2.y) && (tmp.z == tmp2.z) ||
				           (mesh->mesh[tmp.z].meshData[tmp.y][tmp.x].isOut == true) &&
				           (mesh->mesh[tmp2.z].meshData[tmp2.y][tmp2.x].isOut == true) &&
				           (mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x].isOut == false)) {
					// 5
					// tmp,tmp2 - 5, tmp1 - 4
					tmp3 = doStep( mesh, type, tmp1, -1 );
					phi3 = V[mesh->mesh[tmp3.z].meshData[tmp3.y][tmp3.x].Number];
					h3 = distP2P( mesh->mesh[tmp1.z].meshData[tmp1.y][tmp1.x],
					              mesh->mesh[tmp3.z].meshData[tmp3.y][tmp3.x], type );
					H = h1 + h3;
					delta = h1 / h3;
					E[m] = -(delta * phi3 - (1 + delta) * (1 + delta) * phi1 / delta +
					         (2 + delta) * phi0 / delta) /
					       H;
				} else {
					// 3
					// tmp1 - 2, tmp2 - 4
					H = h1 + h2;
					delta = h2 / h1;
					E[m] = -(-delta * phi1 + (delta * delta - 1) * phi0 / delta + phi2 / delta) / H;
				}
			}
		}
	}
}

template<class PointType>
double PoissonSolver<PointType>::distToEdge( DGeo::Edge<PointType> edge,
                                             DGeo::Point<PointType> point, char type ) {
	if (type == 'x') {
		return distToEdgeX( edge, point );
	} else {
		return distToEdgeY( edge, point );
	}
}

template<class PointType>
DGeo::Point<int>
PoissonSolver<PointType>::doStep( const std::shared_ptr<MeshContainer2d<PointType>> &mesh, char type,
                                  DGeo::Point<int> point, int step ) {
	if (type == 'x') {
		return mesh->doStepX( point, step );
	} else {
		return mesh->doStepY( point, step );
	}
}

template<class PointType>
DGeo::Point<int>
PoissonSolver<PointType>::doStep( const std::shared_ptr<MeshContainer3d<PointType>> &mesh, char type,
                                  DGeo::Point<int> point, int step ) {
	if (type == 'x') {
		return mesh->doStepX( point, step );
	} else if (type == 'y') {
		return mesh->doStepY( point, step );
	} else {
		return mesh->doStepZ( point, step );
	}
}

template<typename T>
void write_out( const std::vector<T> &to_out, const std::string &out_name ) {
	ofstream out_stream;
	out_stream.open( out_name );
	for (int i = 0; i < to_out.size(); ++i)
		out_stream << std::setprecision(16) << to_out[i] << std::endl;
	out_stream.close();
}

template<class PointType>
void PoissonSolver<PointType>::solve( std::vector<PointType> &rho, std::vector<PointType> &V,
                                      const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
                                      const std::shared_ptr<BoundaryConditions> &boundaryConditions,
                                      double t, double &progress ) {
	auto meshData = mesh->meshData;
	int k = 0;
	auto s = meshData.size();
	auto s2 = meshData[0].size();
	int max = -1;
	int idx = 0;
//	for(int i =0; i < meshData.size(); ++i)
//	{
//		int s1 = meshData[i].size();
//		if (s1 > max){
//			max = s1;
//			idx = i;
//		}
//	}
//	//std::vector<std::vector<double>> X;
//	//X.resize(s);
//	std::cout << std::flush;
//	for (int i = 0; i < s; ++i){
//		//X[i].resize(max);
//		int kk = 0;
//		auto s1 = meshData[i].size();
//		for(int j = 0; j < s1; ++j){
//			while ( meshData[idx][kk++].x < meshData[i][j].x) {
//				std::cout << (kk==1? 0 : kk == 2? 1e-08: 1.11111e-7) <<'\t';
//			}
//			std::cout<<meshData[i][j].x  << '\t';
//		}
//		std::cout<<std::endl;
//	}


	//std::cout<<s<<" "<<s2<<std::endl;
	int nx = mesh->getnVertexX();
	int ny = mesh->getnVertexY();
	double pi = 3.14159265358979323;
	double tw = std::cos( pi / nx ) / std::cos( pi / ny );
	int size = systemSize;
	std::vector<int> strElements;
	double tempValue;
	double eps0 = VACUUM_PERMITTIVITY();
	int flag = 0;
	double tmp_sum = 0.;
	volatile double diff_sum = 0.;
	vect.resize( size );
	std::vector<PointType> boundConds;
	boundConds.clear();
	for (int i = 0; i < boundCondNum; ++i) {
		boundConds.push_back( boundaryConditions->GetPotential( i, t ));
	}

	for (int i = 0; i < size; ++i) {
		linVect[i] = 0.;
	}

	// for (int i = 0; i < size; ++i) {
	//	rho[i] = -1e-6;
	//}
//    std::ofstream out;
//    out.open("/home/mrkakorin/diploma/solver.txt", std::ios_base::app);
//    out <<"size nonZeroLinVect"<<nonZeroLinVect.size()<<std::endl;
//    out <<"size linVect"<<linVect.size()<<std::endl;
//    out <<"size"<<size<<std::endl;

	for (int i = 0; i < nonZeroLinVect.size(); ++i) {
		if (nonZeroLinVect[i].bcnum != -1) {
//            out<<nonZeroLinVect[i].row<<std::endl;
//            out<<linVect[nonZeroLinVect[i].row]<<std::endl;
//            out<<nonZeroLinVect[i].bcnum <<std::endl;
//            out<<boundConds[nonZeroLinVect[i].bcnum]<<std::endl;
			linVect[nonZeroLinVect[i].row] =
					linVect[nonZeroLinVect[i].row] - nonZeroLinVect[i].coef * boundConds[nonZeroLinVect[i].bcnum];
		}
	}
	//out.close();
	for (int j = 0; j < boundaryPoints.size(); ++j) {
		if (boundaryPoints[j].bcnum != -1)
			x[boundaryPoints[j].pnum] = boundConds[boundaryPoints[j].bcnum];
	}

	for (int i = 0; i < rho.size(); ++i) {
		vect[i] = linVect[i] - rho[i] / eps0;
	}

	if (solverFlags[0] == 1) {
		std::vector<double> r, r_, v_, p, s, t_;
		for (int i = 0; i < size; ++i) {
			r.push_back( 0. );
			r_.push_back( 0. );
			v_.push_back( 0. );
			p.push_back( 0. );
			s.push_back( 0. );
			t_.push_back( 0. );
		}

		double ro, alpha, omega, beta, tmp_sum2;

		for (int i = 0; i < middles.size(); ++i) {
			if (std::isnan( c_rights[i] ) || std::isinf( c_rights[i] ) || std::isnan( c_lefts[i] ) ||
			    std::isinf( c_lefts[i] ) || std::isnan( c_ups[i] ) || std::isinf( c_ups[i] ) ||
			    std::isnan( c_downs[i] ) ||
			    std::isinf( c_downs[i] ) || std::isnan( c_middles[i] ) || std::isinf( c_middles[i] )) {
				int tt = 0;
				continue;
			}

			v_[middles[i]] = 0;
			p[middles[i]] = 0;
			r[middles[i]] =
					vect[middles[i]] -
					(x[middles[i]] * c_middles[i] + x[rights[i]] * c_rights[i] +
					 x[lefts[i]] * c_lefts[i] + x[ups[i]] * c_ups[i] + x[downs[i]] * c_downs[i]);
			r_[middles[i]] = r[middles[i]];
		}
		ro = 1;
		alpha = 1;
		omega = 1;
		while (flag == 0) {

			beta = alpha / (ro * omega);
			ro = 0;
			for (int i = 0; i < middles.size(); ++i) {
				ro = ro + r_[middles[i]] * r[middles[i]];
			}

			beta = beta * ro;
			tmp_sum = 0;
			for (int i = 0; i < middles.size(); ++i) {
				p[middles[i]] = r[middles[i]] + beta * (p[middles[i]] - omega * v_[middles[i]]);
			}
			for (int i = 0; i < middles.size(); ++i) {
				if (std::isnan( c_rights[i] ) || std::isinf( c_rights[i] ) || std::isnan( c_lefts[i] ) ||
				    std::isinf( c_lefts[i] ) || std::isnan( c_ups[i] ) || std::isinf( c_ups[i] ) ||
				    std::isnan( c_downs[i] ) ||
				    std::isinf( c_downs[i] ) || std::isnan( c_middles[i] ) || std::isinf( c_middles[i] )) {
					int tt = 0;
					continue;
				}

				v_[middles[i]] =
						(p[middles[i]] * c_middles[i] + p[rights[i]] * c_rights[i] +
						 p[lefts[i]] * c_lefts[i] + p[ups[i]] * c_ups[i] + p[downs[i]] * c_downs[i]);
				tmp_sum = tmp_sum + r_[middles[i]] * v_[middles[i]];
			}

			alpha = ro / tmp_sum;
			tmp_sum = 0;
			for (int i = 0; i < middles.size(); ++i) {
				s[middles[i]] = r[middles[i]] - alpha * v_[middles[i]];
				tmp_sum = tmp_sum + s[middles[i]];
			}

			if (std::abs( tmp_sum ) < eps_tolerance) {
				flag = 1;
			}

			for (int i = 0; i < middles.size(); ++i) {
				if (std::isnan( c_rights[i] ) || std::isinf( c_rights[i] ) || std::isnan( c_lefts[i] ) ||
				    std::isinf( c_lefts[i] ) || std::isnan( c_ups[i] ) || std::isinf( c_ups[i] ) ||
				    std::isnan( c_downs[i] ) ||
				    std::isinf( c_downs[i] ) || std::isnan( c_middles[i] ) || std::isinf( c_middles[i] )) {
					int tt = 0;
					continue;
				}

				t_[middles[i]] =
						(s[middles[i]] * c_middles[i] + s[rights[i]] * c_rights[i] +
						 s[lefts[i]] * c_lefts[i] + s[ups[i]] * c_ups[i] + s[downs[i]] * c_downs[i]);
			}

			tmp_sum = 0;
			tmp_sum2 = 0;
			for (int i = 0; i < middles.size(); ++i) {
				tmp_sum = tmp_sum + t_[middles[i]] * s[middles[i]];
				tmp_sum2 = tmp_sum2 + t_[middles[i]] * t_[middles[i]];
			}
			omega = tmp_sum / tmp_sum2;

			for (int i = 0; i < middles.size(); ++i) {
				r[middles[i]] = s[middles[i]] - omega * t_[middles[i]];
			}
			diff_sum = 0;
			for (int i = 0; i < middles.size(); ++i) {
				diff_sum = diff_sum + std::abs( omega * s[middles[i]] + alpha * p[middles[i]] );
				x[middles[i]] = x[middles[i]] + omega * s[middles[i]] + alpha * p[middles[i]];
			}
			if (diff_sum < eps_tolerance) {
				flag = 1;
			}
		}
	}

	if (solverFlags[0] == 0) {
		int n = 0;
		double sum;
		/*while (flag == 0){
                diff_sum = 0;
                tmp_sum = 0;
                for (int k = 0; k < middles.size(); ++k){
                        tmp_sum = 0;
                        n = middles[k];
                        tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] +
        x[ups[k]] * c_ups[k] + x[downs[k]] * c_downs[k])*w / c_middles[k];

                        if (std::isnan(c_rights[k]) || std::isinf(c_rights[k]) || std::isnan(c_lefts[k]) ||
        std::isinf(c_lefts[k]) || std::isnan(c_ups[k]) || std::isinf(c_ups[k]) || std::isnan(c_downs[k]) ||
        std::isinf(c_downs[k]) || std::isnan(c_middles[k]) || std::isinf(c_middles[k]))
                        {
                                int tt = 0;
                                continue;
                        }

                        tmp_sum = tmp_sum + (1 - w)*x[n] + w*vect[n] / c_middles[k];
                        if ((x[n]> eps_tolerance*0.1) || (x[n] < -eps_tolerance*0.1)){
                                if (std::abs((x[n] - tmp_sum) / x[n])>diff_sum)
                                        diff_sum = std::abs((x[n] - tmp_sum) / x[n]);
                        //	diff_sum = diff_sum + std::abs(x[n] - tmp_sum) / std::abs(x[n]);
                        }
                        else {
                                if (std::abs((x[n] - tmp_sum)) > diff_sum)
                                        diff_sum = std::abs((x[n] - tmp_sum));
                        }
                        x[n] = tmp_sum;
                }
                if (std::abs(diff_sum) < eps_tolerance) flag = 1;
        }*/


//		write_out<PointType>( x, "./x_begin.csv" );
//
//		write_out<int>( middles, "./middles.csv" );
//		write_out<double>( c_middles, "./c_middles.csv" );
//
//		write_out<int>( rights, "./rights.csv" );
//		write_out<double>( c_rights, "./c_rights.csv" );
//
//		write_out<int>( lefts, "./lefts.csv" );
//		write_out<double>( c_lefts, "./c_lefts.csv" );
//
//		write_out<int>( downs, "./downs.csv" );
//		write_out<double>( c_downs, "./c_downs.csv" );
//
//		write_out<int>( ups, "./ups.csv" );
//		write_out<double>( c_ups, "./c_ups.csv" );
//
//		write_out<double>( vect, "./vect.csv" );

		while (flag == 0) {
			diff_sum = 0;
			tmp_sum = 0;

			for (int k = 0; k < middles.size(); ++k) {
				n = middles[k];
				tmp_sum = -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] +
				            x[ups[k]] * c_ups[k] + x[downs[k]] * c_downs[k]) *
				          w / c_middles[k];
				tmp_sum = tmp_sum + (1 - w) * x[n] + w * vect[n] / c_middles[k];
				if (std::abs(x[n]) > eps_tolerance * 0.1) {
					diff_sum = diff_sum + std::abs( x[n] - tmp_sum ) / x[n];
				} else {
					diff_sum = diff_sum + std::abs( x[n] - tmp_sum );
				}
				x[n] = tmp_sum;
			}
			progress = std::abs( diff_sum ) / eps_tolerance;
			if (std::abs( diff_sum ) < eps_tolerance)
				flag = 1;
		}

		//write_out<PointType>( x, "x.csv" );

	}
	// write solution in gridData
	int m = 0;
	for (int i = 0; i < size; ++i) {
		V[i] = x[i];
	}

	// extrapolation potencial
	DGeo::Point<int> tmp, tmp1, tmp2;
	PointType tmp_V, v, h_, v1, v2, h1, h2;
	BoundaryPoint<PointType> bp;
	std::vector<DGeo::Point<int>> strange_points;
	double another_h;
	int nm;

	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {
			if (mesh->meshData[i][j].isOut == true) {
				bp = findBoundaryPoint( mesh->meshData[i][j].Number );
				nm = 0;
				tmp_V = 0;
				tmp.x = j;
				tmp.y = i;

				if (bp.type == 'n') {
					tmp1 = mesh->doStepX( tmp, 1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepX( tmp1, 1 );
						h1 = distP2PX( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) { //потому что нет резаных ячеек для неймана
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							h2 = distP2PX( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}
					tmp1 = mesh->doStepX( tmp, -1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepX( tmp1, -1 );
						h1 = distP2PX( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) {
							h2 = distP2PX( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}
					tmp1 = mesh->doStepY( tmp, 1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepY( tmp1, 1 );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						h1 = distP2PY( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) {
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							h2 = distP2PY( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}
					tmp1 = mesh->doStepY( tmp, -1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepY( tmp1, -1 );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						h1 = distP2PY( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) {
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							h2 = distP2PY( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}

					if (nm != 0) {
						V[mesh->meshData[tmp.y][tmp.x].Number] = tmp_V / (1.0 * nm);
					} else {
						strange_points.push_back( tmp );
					}
				} else {
					DGeo::Edge<PointType> points_edge;
					DGeo::Edge<PointType> edge;
					points_edge.point1 = mesh->meshData[tmp.y][tmp.x];
					tmp1 = mesh->doStepX( tmp, 1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeX( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PX( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepX( tmp1, 1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PX( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}
					tmp1 = mesh->doStepX( tmp, -1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeX( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PX( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepX( tmp1, -1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PX( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}
					// h = mesh->h2;
					tmp1 = mesh->doStepY( tmp, 1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeY( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PY( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepY( tmp1, 1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PY( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}

					tmp1 = mesh->doStepY( tmp, -1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeY( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PY( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepY( tmp1, -1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PY( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}

					if (nm != 0) {
						V[mesh->meshData[tmp.y][tmp.x].Number] = tmp_V / (1.0 * nm);
					} else {
						strange_points.push_back( tmp );
					}
				}
			}
		}
	}

	PointType v3;
	DGeo::Point<int> tmp3;
	BoundaryPoint<PointType> bp_1, bp_2;
	for (int g = 0; g < 10; ++g) {
		for (int i = 0; i < strange_points.size(); ++i) {
			tmp_V = 0;
			nm = 0;
			tmp = strange_points[i];
			// 1
			tmp1 = mesh->doStepX( tmp, 1 );
			tmp2 = mesh->doStepY( tmp1, 1 );
			tmp3 = mesh->doStepY( tmp, 1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}
			// 2
			tmp1 = mesh->doStepX( tmp, 1 );
			tmp2 = mesh->doStepY( tmp1, -1 );
			tmp3 = mesh->doStepY( tmp, -1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}
			// 3
			tmp1 = mesh->doStepY( tmp, -1 );
			tmp2 = mesh->doStepX( tmp1, -1 );
			tmp3 = mesh->doStepX( tmp, -1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}
			// 4
			tmp1 = mesh->doStepX( tmp, -1 );
			tmp2 = mesh->doStepY( tmp1, 1 );
			tmp3 = mesh->doStepY( tmp, 1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}

			V[mesh->meshData[tmp.y][tmp.x].Number] = tmp_V / (1.0 * nm);
		}
	}

	/*	for (int i = 0; i < mesh->meshData[1].size(); i++)
            {
                    for (int j = 0; j < mesh->meshData[0].size(); j++)
                    {
                            if (std::abs(mesh->meshData[0][j].x -mesh->meshData[1][i].x)<1e-7)
                            {
                                    V[mesh->meshData[0][j].Number] = V[mesh->meshData[1][i].Number];
                            };
                    };
            };*/

	strange_points.clear();
	progress = 1;
}

template<class PointType>
void PoissonSolver<PointType>::solveCharge(
		std::vector<PointType> &rho, std::vector<PointType> &V,
		const std::shared_ptr<MeshContainer2d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t,
		std::vector<int> &nonZeros ) {

	int nx = mesh->getnVertexX();
	int ny = mesh->getnVertexY();
	double pi = 3.14159265358979323;
	double tw = std::cos( pi / nx ) / std::cos( pi / ny );
	int size = systemSize;
	std::vector<int> strElements;
	double tempValue;
	double eps0 = VACUUM_PERMITTIVITY();
	int flag = 0;
	double tmp_sum = 0;
	volatile double diff_sum = 0;
	vect.resize( size );

	for (int i = 0; i < size; ++i) {
		linVect[i] = 0;
	}

	for (int j = 0; j < boundaryPoints.size(); ++j) {
		if (boundaryPoints[j].bcnum != -1)
			xCharge[boundaryPoints[j].pnum] = 0;
	}

	for (int i = 0; i < rho.size(); ++i) {
		vect[i] = linVect[i] - rho[i] / eps0;
	}

	int s1;
	for (s1 = 0; s1 < rho.size(); s1++) {
		if (rho[s1] != 0)
			break;
	}

	int s2;
	for (s2 = rho.size() - 1; s2 >= 0; s2--) {
		if (rho[s2] != 0)
			break;
	}

	if (s1 == rho.size() || s2 == -1) {
		memset( &xCharge[0], 0, xCharge.size() * sizeof( xCharge[0] ));
		return;
	}

	if (solverFlags[0] == 1) {
		std::vector<double> r, r_, v_, p, s, t_;
		for (int i = 0; i < size; ++i) {
			r.push_back( 0 );
			r_.push_back( 0 );
			v_.push_back( 0 );
			p.push_back( 0 );
			s.push_back( 0 );
			t_.push_back( 0 );
		}

		double ro, alpha, omega, beta, tmp_sum2;

		for (int i = 0; i < middles.size(); ++i) {
			if (std::isnan( c_rights[i] ) || std::isinf( c_rights[i] ) || std::isnan( c_lefts[i] ) ||
			    std::isinf( c_lefts[i] ) || std::isnan( c_ups[i] ) || std::isinf( c_ups[i] ) ||
			    std::isnan( c_downs[i] ) ||
			    std::isinf( c_downs[i] ) || std::isnan( c_middles[i] ) || std::isinf( c_middles[i] )) {
				int tt = 0;
				continue;
			}

			v_[middles[i]] = 0;
			p[middles[i]] = 0;
			r[middles[i]] = vect[middles[i]] -
			                (xCharge[middles[i]] * c_middles[i] + xCharge[rights[i]] * c_rights[i] +
			                 xCharge[lefts[i]] * c_lefts[i] + xCharge[ups[i]] * c_ups[i] +
			                 xCharge[downs[i]] * c_downs[i]);
			r_[middles[i]] = r[middles[i]];
		}
		ro = 1;
		alpha = 1;
		omega = 1;
		while (flag == 0) {

			beta = alpha / (ro * omega);
			ro = 0;
			for (int i = 0; i < middles.size(); ++i) {
				ro = ro + r_[middles[i]] * r[middles[i]];
			}

			beta = beta * ro;
			tmp_sum = 0;
			for (int i = 0; i < middles.size(); ++i) {
				p[middles[i]] = r[middles[i]] + beta * (p[middles[i]] - omega * v_[middles[i]]);
			}
			for (int i = 0; i < middles.size(); ++i) {
				if (std::isnan( c_rights[i] ) || std::isinf( c_rights[i] ) || std::isnan( c_lefts[i] ) ||
				    std::isinf( c_lefts[i] ) || std::isnan( c_ups[i] ) || std::isinf( c_ups[i] ) ||
				    std::isnan( c_downs[i] ) ||
				    std::isinf( c_downs[i] ) || std::isnan( c_middles[i] ) || std::isinf( c_middles[i] )) {
					int tt = 0;
					continue;
				}

				v_[middles[i]] =
						(p[middles[i]] * c_middles[i] + p[rights[i]] * c_rights[i] +
						 p[lefts[i]] * c_lefts[i] + p[ups[i]] * c_ups[i] + p[downs[i]] * c_downs[i]);
				tmp_sum = tmp_sum + r_[middles[i]] * v_[middles[i]];
			}

			alpha = ro / tmp_sum;
			tmp_sum = 0;
			for (int i = 0; i < middles.size(); ++i) {
				s[middles[i]] = r[middles[i]] - alpha * v_[middles[i]];
				tmp_sum = tmp_sum + s[middles[i]];
			}

			if (std::abs( tmp_sum ) < eps_tolerance_charge) {
				flag = 1;
			}

			for (int i = 0; i < middles.size(); ++i) {
				if (std::isnan( c_rights[i] ) || std::isinf( c_rights[i] ) || std::isnan( c_lefts[i] ) ||
				    std::isinf( c_lefts[i] ) || std::isnan( c_ups[i] ) || std::isinf( c_ups[i] ) ||
				    std::isnan( c_downs[i] ) ||
				    std::isinf( c_downs[i] ) || std::isnan( c_middles[i] ) || std::isinf( c_middles[i] )) {
					int tt = 0;
					continue;
				}

				t_[middles[i]] =
						(s[middles[i]] * c_middles[i] + s[rights[i]] * c_rights[i] +
						 s[lefts[i]] * c_lefts[i] + s[ups[i]] * c_ups[i] + s[downs[i]] * c_downs[i]);
			}

			tmp_sum = 0;
			tmp_sum2 = 0;
			for (int i = 0; i < middles.size(); ++i) {
				tmp_sum = tmp_sum + t_[middles[i]] * s[middles[i]];
				tmp_sum2 = tmp_sum2 + t_[middles[i]] * t_[middles[i]];
			}
			omega = tmp_sum / tmp_sum2;

			for (int i = 0; i < middles.size(); ++i) {
				r[middles[i]] = s[middles[i]] - omega * t_[middles[i]];
			}
			diff_sum = 0;
			for (int i = 0; i < middles.size(); ++i) {
				diff_sum = diff_sum + std::abs( omega * s[middles[i]] + alpha * p[middles[i]] );
				xCharge[middles[i]] =
						xCharge[middles[i]] + omega * s[middles[i]] + alpha * p[middles[i]];
			}
			if (diff_sum < eps_tolerance_charge) {
				flag = 1;
			}
		}
	}
	s2 = s2 + ChargeSpace;
	s1 = s1 - ChargeSpace;
	int iter = 0;
	int s = middles.size();

	int n;
	for (int k = 0; k < s; ++k) {
		n = middles[k];
		if (n > s2 || n < s1)
			xCharge[n] = 0;
	}

	if (solverFlags[0] == 0) {
		int n = 0;
		double sum;
		while (flag == 0) {
			diff_sum = 0;
			tmp_sum = 0;
			for (int iii = 0; iii < nonZeros.size(); ++iii) {

				// int k = nonZeros[kk];

				// n = middles[k];

				int n = nonZeros[iii];
				int k = n2k[n];

				if (k == -1)
					continue;
				// if (n>s2 || n<s1 )
				//	continue;

				if (kZeros[k] == -1) {
					xCharge[n] = 0;
					continue;
				}

				tmp_sum = -(xCharge[rights[k]] * c_rights[k] + xCharge[lefts[k]] * c_lefts[k] +
				            xCharge[ups[k]] * c_ups[k] + xCharge[downs[k]] * c_downs[k]) *
				          w_charge / c_middles[k];

				/*	if (std::isnan(c_rights[k]) || std::isinf(c_rights[k]) || std::isnan(c_lefts[k]) ||
                   std::isinf(c_lefts[k]) || std::isnan(c_ups[k]) || std::isinf(c_ups[k]) || std::isnan(c_downs[k]) ||
                   std::isinf(c_downs[k]) || std::isnan(c_middles[k])
                   || std::isinf(c_middles[k]))
                        {
                                int tt = 0;
                                continue;
                        }*/

				/*if (std::isnan(tmp_sum) || std::isinf(tmp_sum))
                {
                        __debugbreak();

                        int tt = 0;
                }*/

				tmp_sum = tmp_sum + (1 - w_charge) * xCharge[n] + w_charge * vect[n] / c_middles[k];

				if (std::abs( xCharge[n] ) > eps_tolerance_charge * 0.1 &&
				    std::abs((xCharge[n] - tmp_sum) / xCharge[n] ) > diff_sum)
					diff_sum = std::abs((xCharge[n] - tmp_sum) / xCharge[n] );

				xCharge[n] = tmp_sum;
			}
			iter++;
			if (std::abs( diff_sum ) < eps_tolerance_charge)
				flag = 1;
		}
	}
	// write solution in gridData
	int m = 0;
	for (int i = 0; i < size; ++i) {
		V[i] = xCharge[i];
	}

	DGeo::Point<int> tmp, tmp1, tmp2;
	PointType tmp_V, v, h_, v1, v2, h1, h2;
	BoundaryPoint<PointType> bp;
	std::vector<DGeo::Point<int>> strange_points;
	double another_h;
	int nm;

	for (int i = 0; i < mesh->meshData.size(); ++i) {
		for (int j = 0; j < mesh->meshData[i].size(); ++j) {
			if (mesh->meshData[i][j].isOut == true) {
				bp = findBoundaryPoint( mesh->meshData[i][j].Number );
				nm = 0;
				tmp_V = 0;
				tmp.x = j;
				tmp.y = i;

				if (bp.type == 'n') {
					tmp1 = mesh->doStepX( tmp, 1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepX( tmp1, 1 );
						h1 = distP2PX( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) { //потому что нет резаных ячеек для неймана
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							h2 = distP2PX( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}
					tmp1 = mesh->doStepX( tmp, -1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepX( tmp1, -1 );
						h1 = distP2PX( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) {
							h2 = distP2PX( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}
					tmp1 = mesh->doStepY( tmp, 1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepY( tmp1, 1 );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						h1 = distP2PY( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) {
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							h2 = distP2PY( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}
					tmp1 = mesh->doStepY( tmp, -1 );
					if (mesh->meshData[tmp1.y][tmp1.x].isOut == false) {
						tmp2 = mesh->doStepY( tmp1, -1 );
						v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
						h1 = distP2PY( mesh->meshData[tmp.y][tmp.x], mesh->meshData[tmp1.y][tmp1.x] );
						if ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y)) {
							v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
							h2 = distP2PY( mesh->meshData[tmp1.y][tmp1.x],
							               mesh->meshData[tmp2.y][tmp2.x] );
							tmp_V = tmp_V + h1 * (v1 - v2) / h2 + v1;
							nm++;
						}
					}

					if (nm != 0) {
						V[mesh->meshData[tmp.y][tmp.x].Number] = tmp_V / (1.0 * nm);
					} else {
						strange_points.push_back( tmp );
					}
				} else {
					DGeo::Edge<PointType> points_edge;
					DGeo::Edge<PointType> edge;
					points_edge.point1 = mesh->meshData[tmp.y][tmp.x];
					tmp1 = mesh->doStepX( tmp, 1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeX( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PX( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepX( tmp1, 1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PX( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}
					tmp1 = mesh->doStepX( tmp, -1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeX( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PX( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepX( tmp1, -1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PX( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}
					// h = mesh->h2;
					tmp1 = mesh->doStepY( tmp, 1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeY( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PY( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepY( tmp1, 1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PY( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}

					tmp1 = mesh->doStepY( tmp, -1 );
					points_edge.point2 = mesh->meshData[tmp1.y][tmp1.x];
					if ((mesh->meshData[tmp1.y][tmp1.x].isOut == false)) {
						if (isIntersection( bp.edge, points_edge, eps_p )) {
							h_ = distToEdgeY( bp.edge, mesh->meshData[tmp1.y][tmp1.x] );
							h1 = distP2PY( mesh->meshData[tmp.y][tmp.x],
							               mesh->meshData[tmp1.y][tmp1.x] );
							v1 = V[mesh->meshData[tmp.y][tmp.x].Number];
							v2 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
							v = h1 * (v1 - v2) / h_ + v2;
							nm++;
							if (h_ < eps_h) {
								tmp2 = mesh->doStepY( tmp1, -1 );
								v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
								h2 = distP2PY( mesh->meshData[tmp2.y][tmp2.x],
								               mesh->meshData[tmp1.y][tmp1.x] );
								v = (h1 - h_) * (v1 - v2) / h2 + v1;
							}
							tmp_V = tmp_V + v;
						}
					}

					if (nm != 0) {
						V[mesh->meshData[tmp.y][tmp.x].Number] = tmp_V / (1.0 * nm);
					} else {
						strange_points.push_back( tmp );
					}
				}
			}
		}
	}

	PointType v3;
	DGeo::Point<int> tmp3;
	BoundaryPoint<PointType> bp_1, bp_2;
	for (int g = 0; g < 10; ++g) {
		for (int i = 0; i < strange_points.size(); ++i) {
			tmp_V = 0;
			nm = 0;
			tmp = strange_points[i];
			// 1
			tmp1 = mesh->doStepX( tmp, 1 );
			tmp2 = mesh->doStepY( tmp1, 1 );
			tmp3 = mesh->doStepY( tmp, 1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}
			// 2
			tmp1 = mesh->doStepX( tmp, 1 );
			tmp2 = mesh->doStepY( tmp1, -1 );
			tmp3 = mesh->doStepY( tmp, -1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}
			// 3
			tmp1 = mesh->doStepY( tmp, -1 );
			tmp2 = mesh->doStepX( tmp1, -1 );
			tmp3 = mesh->doStepX( tmp, -1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}
			// 4
			tmp1 = mesh->doStepX( tmp, -1 );
			tmp2 = mesh->doStepY( tmp1, 1 );
			tmp3 = mesh->doStepY( tmp, 1 );
			if (((tmp1.x != tmp.x) || (tmp1.y != tmp.y)) &&
			    ((tmp3.x != tmp.x) || (tmp3.y != tmp.y)) &&
			    ((tmp2.x != tmp1.x) || (tmp2.y != tmp1.y))) {
				v1 = V[mesh->meshData[tmp1.y][tmp1.x].Number];
				v2 = V[mesh->meshData[tmp2.y][tmp2.x].Number];
				v3 = V[mesh->meshData[tmp3.y][tmp3.x].Number];
				tmp_V = tmp_V + v1 + v3 - v2;
				nm++;
			}

			V[mesh->meshData[tmp.y][tmp.x].Number] = tmp_V / (1.0 * nm);
		}
	}

	/*for (int i = 0; i < mesh->meshData.size(); i++)
    {
            V[mesh->meshData[i][0].Number] = V[mesh->meshData[i][1].Number];
    };

    for (int i = 0; i < mesh->meshData[0].size(); i++)
    {
            V[mesh->meshData[0][i].Number] = V[mesh->meshData[1][i].Number];
    };*/
}

template<class PointType>
void PoissonSolver<PointType>::solve_3d(
		std::vector<PointType> &rho, std::vector<PointType> &V,
		const std::shared_ptr<MeshContainer3d<PointType>> &mesh,
		const std::shared_ptr<BoundaryConditions> &boundaryConditions, double t ) {

}

template<class PointType>
bool PoissonSolver<PointType>::isIntersection( DGeo::Edge<PointType> boundary_edge,
                                               DGeo::Edge<PointType> points_edge, double eps ) {
	if (polar == false) {
		return boundary_edge.IsIntersectionEdge( points_edge, eps );
	} else {
		DGeo::Edge<PointType> edge;
		Dmath::Polar2Cartesian( points_edge.point1.x, points_edge.point1.y, edge.point1.x,
		                        edge.point1.y );
		Dmath::Polar2Cartesian( points_edge.point2.x, points_edge.point2.y, edge.point2.x,
		                        edge.point2.y );
		return boundary_edge.IsIntersectionEdge( edge, eps );
	}
}

template class PoissonSolver<float>;
template class PoissonSolver<double>;

template void
PoissonSolver<float>::load<boost::archive::binary_iarchive>( boost::archive::binary_iarchive &ar,
                                                             const unsigned int file_version );

template void
PoissonSolver<double>::save<boost::archive::binary_oarchive>( boost::archive::binary_oarchive &ar,
                                                              const unsigned int file_version ) const;

template void
PoissonSolver<double>::load<boost::archive::binary_iarchive>( boost::archive::binary_iarchive &ar,
                                                              const unsigned int file_version );

template void
PoissonSolver<float>::save<boost::archive::binary_oarchive>( boost::archive::binary_oarchive &ar,
                                                             const unsigned int file_version ) const;
template <class PointType>
void updatePositionsVector(PointType* r, PointType* z, PointType* phi, PointType* pr, PointType* pz, PointType* pphi,
                           int size, PointType timeStep, PointType* tmp1, PointType* tmp2);

template <class PointType>
void updateMomentumsVector(PointType* r, PointType* z, PointType* phi, PointType* pr, PointType* pz, PointType* pphi,
                           int size, PointType timeStep, PointType* tmp1, PointType* tmp2);
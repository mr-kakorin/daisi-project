template <class PointType>
void SOR(PointType* x, int* rights, PointType* c_rights, int* lefts, PointType* c_lefts, int* ups, PointType* c_ups,
         int* downs, PointType* c_downs, int* middles, PointType* c_middles, PointType w, int size, PointType* vect,
         double& diff_sum, PointType* tmp);
class ICell2d
{
    int          nPoints;
    int*         pointArray;
    virtual void InCell(double x1, double x2) = 0;
};
class CellCIC2d : public ICell2d
{
};
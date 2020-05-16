#include <project.h>

#define VAL(str) #str
#define TOSTRING(str) VAL(str)

int main()
{
    Dproject::project currentProject;

    std::string errorMsgOpen;

    std::string path = TOSTRING(DATA_PATH);

    currentProject.LoadProject(path + "/transport_2_as/transport_2.dproj", "1.0", errorMsgOpen);
}
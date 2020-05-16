#include "SynchrotronTools.h"
#include "Tools.h"

void savetrack(const std::string& fileName, std::shared_ptr<OpticElementsSequence> seq,
               const std::vector<int>&                                           ind,
               const std::vector<std::vector<std::vector<std::vector<double>>>>& madXData)
{
    FILE* fid = fopen(fileName.c_str(), "w");

    std::vector<int> indWrite = ind;

    for (int i      = 1; i < indWrite.size(); i++)
        indWrite[i] = indWrite[i] + 1;

    int         nCircle = madXData.size();
    std::string time;
    std::string date;

    GetTime(time, date);

    fprintf(fid, "@ NAME             %4s \"TRACKONE\"\n", "%08s");
    fprintf(fid, "@ TYPE             %4s \"TRACKONE\"\n", "%08s");

    fprintf(fid, "@ TITLE            %4s \"no-title\"\n", "%08s");
    fprintf(fid, "@ ORIGIN           %4s \n", "%24s");

    fprintf(fid, "@ DATE             %4s %8s\n", "%08s", date.c_str());
    fprintf(fid, "@ TIME             %4s %8s\n", "%08s", time.c_str());

    // ��������� �������
    fprintf(fid, "*     NUMBER       TURN                  X                 PX                  Y "
                 "                PY  "
                 "                T                 PT                  S                  E\n");
    fprintf(fid, "$         %2s         %2s                %3s                %3s                "
                 "%3s                %3s     "
                 "           %3s                %3s                %3s                %3s\n",
            "%d", "%d", "%le", "%le", "%le", "%le", "%le", "%le", "%le", "%le");

    int i_segment = 0;

    for (int i = 0; i < nCircle; i++)
    {
        for (int j = 0; j < madXData[i].size(); j++)
        {
            i_segment = i_segment + 1;
            fprintf(fid, "#segment %7s %7s %7s %7s %s \n", std::to_string(i_segment).c_str(),
                    std::to_string(i + 1).c_str(), std::to_string(madXData[i][j].size()).c_str(),
                    std::to_string(indWrite[j]).c_str(), seq->GetLabel(ind[j]).c_str());

            for (int k = 0; k < madXData[i][j].size(); k++)
            {
                fprintf(fid, "%11s %10s %18s %18s %18s %18s %18s %18s %18s %18s\n",
                        std::to_string(k).c_str(), std::to_string(i).c_str(),
                        std::to_string(madXData[i][j][k][2]).c_str(),
                        std::to_string(madXData[i][j][k][3]).c_str(),
                        std::to_string(madXData[i][j][k][4]).c_str(),
                        std::to_string(madXData[i][j][k][5]).c_str(), std::to_string(0).c_str(),
                        std::to_string(0).c_str(), std::to_string(seq->getParameter(ind[j], "at") +
                                                                  seq->getParameter(ind[j], "L"))
                                                       .c_str(),
                        std::to_string(0).c_str());
            }
        };
    }
    fclose(fid);
};

void savetwiss(const std::string& fileName, const std::vector<arma::mat>& xtwiss,
               const std::vector<arma::mat>& ytwiss, const std::string& particleType,
               double GeVmass, std::shared_ptr<OpticElementsSequence> seq,
               const std::vector<float>& mu_x, const std::vector<float>& mu_y)
{
    FILE* fp = fopen(fileName.c_str(), "w");

    fprintf(fp, "@ NAME            %5s \"TWISS\"\n", "%05s");
    fprintf(fp, "@ TYPE            %5s \"TWISS\"\n", "%05s");
    fprintf(fp, "@ SEQUENCE        %5s \"NUCLOTRON\"\n", "%09s");

    fprintf(fp, "@ PARTICLE         %%0%ds \"%s\"\n", particleType.size(), particleType.c_str());
    fprintf(fp, "@ MASS             %3s         %lf\n", "%1e", GeVmass);
    fprintf(fp, "@ CHARGE           %3s \n", "%1e");
    fprintf(fp, "@ ENERGY           %3s \n", "%1e");
    fprintf(fp, "@ PC               %3s \n", "%1e");
    fprintf(fp, "@ GAMMA            %3s \n", "%1e");
    fprintf(fp, "@ KBUNCH           %3s \n", "%1e");
    fprintf(fp, "@ BCURRENT         %3s \n", "%1e");
    fprintf(fp, "@ SIGE             %3s \n", "%1e");
    fprintf(fp, "@ SIGT             %3s \n", "%1e");
    fprintf(fp, "@ NPART            %3s \n", "%1e");
    fprintf(fp, "@ EX               %3s \n", "%1e");
    fprintf(fp, "@ EY               %3s \n", "%1e");
    fprintf(fp, "@ ET               %3s \n", "%1e");

    fprintf(fp, "@ LENGTH           %3s         %lf\n", "%1e", seq->GetL());
    fprintf(fp, "@ ALFA             %3s \n", "%1e");
    fprintf(fp, "@ ORBIT5           %3s \n", "%1e");
    fprintf(fp, "@ GAMMATR          %3s \n", "%1e");
    fprintf(fp, "@ Q1               %3s \n", "%1e");
    fprintf(fp, "@ Q2               %3s \n", "%1e");
    fprintf(fp, "@ DQ1              %3s \n", "%1e");
    fprintf(fp, "@ DQ2              %3s \n", "%1e");
    fprintf(fp, "@ DXMAX            %3s \n",
            "%1e"); // �������� ���������� �������������� ��������� Dx
    fprintf(fp, "@ DYMAX            %3s \n",
            "%1e"); // �������� ���������� �������������� ��������� Dy
    fprintf(fp, "@ XCOMAX           %3s \n", "%1e");
    fprintf(fp, "@ YCOMAX           %3s \n", "%1e");

    double BETXm = xtwiss[0](0, 0);
    double BETYm = ytwiss[0](0, 0);

    for (int i = 0; i < xtwiss.size(); i++)
    {
        if (xtwiss[i](0, 0) > BETXm)
            BETXm = xtwiss[i](0, 0);

        if (ytwiss[i](0, 0) > BETYm)
            BETYm = ytwiss[i](0, 0);
    }

    fprintf(fp, "@ BETXMAX          %3s         %lf\n", "%1e", BETXm);
    fprintf(fp, "@ BETYMAX          %3s         %lf\n", "%1e", BETYm);

    fprintf(fp, "@ XCORMS           %3s \n", "%1e");
    fprintf(fp, "@ YCORMS           %3s \n", "%1e");
    fprintf(fp, "@ DXRMS            %3s \n", "%1e");
    fprintf(fp, "@ DYRMS            %3s \n", "%1e");
    fprintf(fp, "@ DELTAP           %3s \n", "%1e");
    fprintf(fp, "@ SYNCH_1          %3s \n", "%1e");
    fprintf(fp, "@ SYNCH_2          %3s \n", "%1e");
    fprintf(fp, "@ SYNCH_3          %3s \n", "%1e");
    fprintf(fp, "@ SYNCH_4          %3s \n", "%1e");
    fprintf(fp, "@ SYNCH_5          %3s \n", "%1e");
    fprintf(fp, "@ TITLE            %4s \"no-title\"\n", "%08s");
    fprintf(fp, "@ ORIGIN           %4s \n", "%24s");

    std::string time;
    std::string date;

    GetTime(time, date);

    fprintf(fp, "@ DATE             %4s %8s\n", "%08s", date.c_str());
    fprintf(fp, "@ TIME             %4s %8s\n", "%08s", time.c_str());

    fprintf(fp, "* NAME               KEYWORD                             S               BETX     "
                "          ALFX      "
                "          MUX               BETY               ALFY                MUY \n");
    fprintf(fp, "$ %2s                 %2s                                %3s                %3s   "
                "             %3s        "
                "        %3s                %3s                %3s                %3s \n",
            "%s", "%s", "%1e", "%1e", "%1e", "%1e", "%1e", "%1e", "%1e");
    for (int i = 0; i < seq->length(); i++)
    {
        fprintf(fp, "\"%-18s %-21s %16s %18s %18s %18s %18s %18s %18s\n",
                (seq->GetLabel(i) + "\"").c_str(), ("\"" + seq->GetType(i) + "\"").c_str(),
                std::to_string(seq->getParameter(i, "at") + seq->getParameter(i, "L")).c_str(),
                std::to_string(xtwiss[i](0, 0)).c_str(), std::to_string(xtwiss[i](1, 0)).c_str(),
                //	std::to_string(xtwiss[i](2, 0)).c_str(),
                std::to_string(mu_x[i]).c_str(), std::to_string(ytwiss[i](0, 0)).c_str(),
                std::to_string(ytwiss[i](1, 0)).c_str(), std::to_string(mu_y[i]).c_str());
        //	std::to_string(ytwiss[i](2, 0)).c_str()),
    }

    fclose(fp);
};

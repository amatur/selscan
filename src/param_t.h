#ifndef __PARAM_T_H__
#define __PARAM_T_H__

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <cctype>

using namespace std;

const string ARG_HELP = "--help";

struct param_main{
    string hapFilename, hapFilename2 ;
    string thapFilename, thapFilename2 ;
    string mapFilename;
    string tpedFilename, tpedFilename2;
    string vcfFilename, vcfFilename2 ;
    string jsonFilename;

    string outFilename;
    string query;
    string query_ehh12;

    bool TPED, VCF, THAP;

    int SCALE_PARAMETER, MAX_GAP;
    double EHH_CUTOFF;
    bool ALT,WAGH,TRUNC;
    double MAF;
    int MAX_EXTEND;
    int MAX_EXTEND_NSL;
    bool SKIP, WRITE_DETAILED_IHS;
    bool CALC_XPNSL, CALC_XP;
    bool SINGLE_EHH;
    bool SINGLE_EHH12;

    bool CALC_NSL, CALC_IHS, CALC_SOFT;

    int numThreads;

    bool UNPHASED;
    bool USE_PMAP;

    int EHH1K;
    int QWIN;

    bool CALC_PI;
    int PI_WIN;

    // bool MULTI_CHR;
    // string CHR_LIST;

    string benchmark_flag = "NO_XOR";//XOR
    string benchmark_flag2 = ""; //"FLIP";

    string MISSING_MODE = "ONE_IMPUTE";
    //"NO_IMPUTE"; //ONE_IMPUTE, ZERO_IMPUTE, NO_IMPUTE 

    //string benchmark_flag = "BITSET";
    //string benchmark_flag = "FLIP_ONLY";
    //string benchmark_flag = "BASIC";

    string DEBUG_FLAG = "VCF";
    //string DEBUG_FLAG = "";

    bool MISSING_ALLOWED = false;
    bool MULTI_MAF = false;
    bool MULTI_PARAMS = false;

};


class param_t
{
public:

    bool addFlag(string flag, bool value, string label, string description);
    bool addFlag(string flag, double value, string label, string description);
    bool addFlag(string flag, int value, string label, string description);
    bool addFlag(string flag, char value, string label, string description);
    bool addFlag(string flag, string value, string label, string description);
    bool addFlag(string flag, const char value[], string label, string description);

    bool addListFlag(string flag, string value, string label, string description);
    bool addListFlag(string flag, const char value[], string label, string description);
    bool addListFlag(string flag, int value, string label, string description);


    void printHelp();

    bool parseCommandLine(int argc, char *argv[]);

    bool getBoolFlag(string flag);
    double getDoubleFlag(string flag);
    int getIntFlag(string flag);
    char getCharFlag(string flag);
    string getStringFlag(string flag);

    vector<string> getStringListFlag(string flag);
    vector<int> getIntListFlag(string flag);

    void setPreamble(string str);

    param_t();


private:

    map<string, bool> argb;
    map<string, double> argd;
    map<string, int> argi;
    map<string, char> argch;
    map<string, string> args;

    map<string, vector< string > > listargs;
    map<string, vector< int > > listargi;    

    map<string, string> help;
    map<string, bool> isSet;

    map<string, string> labels;

    bool goodDouble(string str);
    bool goodInt(string str);
    bool goodChar(string str);

    bool flagExists(string flag);

    string preamble;
};

#endif

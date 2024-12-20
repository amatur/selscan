#ifndef __HAPDATA_H__
#define __HAPDATA_H__

#include <vector>
#include <string>
#include "bitset.h"
#include <queue>
#include <iostream>
#include <random>

using namespace std;

struct HapEntry
{
    bool flipped = false;
    vector <int> positions; // 0-based indexing of derived alleles ("1") (if flipped false)
    vector <int> xors;
    
    // for unphased
    vector <int> g[3];
    vector <int> positions2; // 0-based indexing of allele "2"

    
   

    // for LOW_MEM
    MyBitset* hapbitset;
    MyBitset* xorbitset;

    //for missing
    MyBitset* missbitset = xorbitset;

    

};

class HapData
{
public:
    void xor_for_phased_and_unphased(); //experimental

    struct HapEntry* hapEntries = NULL; //vector of haplotype entries
    int nloci;
    int nhaps;

    string benchmark_flag = "NO_XOR";//XOR
    string benchmark_flag2 = ""; //"FLIP";

    string MISSING_MODE = "NO_IMPUTE"; //ONE_IMPUTE, ZERO_IMPUTE, NO_IMPUTE 

    //string benchmark_flag = "BITSET";
    //string benchmark_flag = "FLIP_ONLY";
    //string benchmark_flag = "BASIC";

    string DEBUG_FLAG = "VCF";
    //string DEBUG_FLAG = "";

    bool unphased;
    double MAF;
    bool SKIP;
    int num_threads;
    bool LOW_MEM = true;
    bool MISSING_ALLOWED = false;
    bool MULTI_CHR = false;

    queue<int> skipQueue;

    int missing_count = 0;
    
    //ofstream* flog;

    ~HapData();

    //allocates the 2-d array and populated it with -9,
    /** Sets up structure according to nhaps and nloci
     * 
    */
    void initHapData(int nhaps, int nloci);
    void releaseHapData();
    /**
     * reads in haplotype data and also does basic checks on integrity of format
     * returns a populated HaplotypeData structure if successful
     * throws an exception otherwise
     */ 
    void readHapData(string filename);
    void readHapDataTPED(string filename);
    void readHapDataVCF(string filename);
    void readHapDataVCFMissing(string filename);

    // void readHapDataVCF_bitset(string filename);
    // void readHapData_bitset(string filename);


    
    void initParams(bool UNPHASED, bool SKIP, double MAF, int num_threads, bool LOW_MEM, bool MISSING){
        this->unphased = UNPHASED;
        this->SKIP = SKIP;
        this->MAF = MAF;
        this->num_threads = num_threads;
        this->LOW_MEM = LOW_MEM;
        this->MISSING_ALLOWED = MISSING;
    }

    pair<int, int> countFieldsAndOnes(const string &str)
    {
        string::const_iterator it;
        int ones = 0;
        int result;
        int numFields = 0;
        int seenChar = 0;
        for (it = str.begin() ; it < str.end(); it++)
        {
            if(*it == '1'){
                ones++;
            }
            result = isspace(*it);
            if (result == 0 && seenChar == 0)
            {
                numFields++;
                seenChar = 1;
            }
            else if (result != 0)
            {
                seenChar = 0;
            }
        }
        return make_pair(numFields, ones);
    }

    int countFields(const string &str)
    {
        string::const_iterator it;
        int result;
        int numFields = 0;
        int seenChar = 0;
        for (it = str.begin() ; it < str.end(); it++)
        {
            result = isspace(*it);
            if (result == 0 && seenChar == 0)
            {
                numFields++;
                seenChar = 1;
            }
            else if (result != 0)
            {
                seenChar = 0;
            }
        }
        return numFields;
    }

    char randomZeroOrOne() {
        // Create a random device to seed the generator
        std::random_device rd; 
        
        // Use a Mersenne Twister engine
        std::mt19937 gen(rd()); 
        
        // Define a uniform distribution for integers between 0 and 1
        std::uniform_int_distribution<int> dist(0, 1);
        
        // Return the random number
        int randomNum  = dist(gen);

         // Convert integer 0 or 1 to character '0' or '1'
        char randomChar = static_cast<char>(randomNum + '0');
        return randomChar;
    }



    void print(){
        for(int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
            cout<<locus_after_filter<<"::: ";
            for(int i = 0; i < hapEntries[locus_after_filter].positions.size(); i++){
                cout<<hapEntries[locus_after_filter].positions[i]<<" ";
            }
            cout<<endl;
            
        }

        /*
        for (int i = 0; i < 3; i++)
        {
            cout << "Locus: " << i << endl;
            cout << "Flipped: " << hapEntries[i].flipped << endl;
            cout << "Xors: ";
            for (int j = 0; j < hapEntries[i].xors.size(); j++)
            {
                cout << hapEntries[i].xors[j] << " ";
            }
            cout << endl;
            cout << "Positions: ";
            for (int j = 0; j < hapEntries[i].positions.size(); j++)
            {
                cout << hapEntries[i].positions[j] << " ";
            }
            cout << endl;

            // if(hapEntries[i].positions2.size()>0){
            //     cout << "Positions2: ";
            //     for (int j = 0; j < hapEntries[i].positions2.size(); j++)
            //     {
            //         cout << hapEntries[i].positions2[j] << " ";
            //     }
            //     cout << endl;
            // }
            
        }
        */
    }

    

    /// @brief assumes no missing
    /// @param locus 
    /// @return 
    double calcFreq(int locus)
    {
        //assuming no flip
        double total = 0;
        double freq = 0;
        if(this->unphased){
            if(LOW_MEM){
                //freq = hapEntries[locus].hapbitset->num_1s + (hapEntries[locus].xorbitset->num_1s-hapEntries[locus].hapbitset->num_1s)*2;
                freq = hapEntries[locus].hapbitset->num_1s + (hapEntries[locus].xorbitset->num_1s)*2;
            }else{
                freq = hapEntries[locus].positions.size() + hapEntries[locus].positions2.size()*2;
            }
            total = this->nhaps*2;
        }else{
            freq = hapEntries[locus].positions.size();
            if(LOW_MEM){
                //freq = hapEntries[locus].hapbitset->count_1s();
                freq = hapEntries[locus].hapbitset->num_1s;
            }
            total = this->nhaps;
        }
        return (freq / total);
    }

    double calcMissingFreq(int locus){
        if(!MISSING_ALLOWED){
            throw "missing not allowed";
        }
         if(!LOW_MEM){
            throw "not implemented";
         }
        return hapEntries[locus].xorbitset->num_1s / (double) nhaps;
    }

    inline int get_n_c2(int locus)
    {
        return LOW_MEM? nhaps - hapEntries[locus].xorbitset->num_1s: nhaps - hapEntries[locus].positions2.size();
        // // if flip was enabled
        // int non_flipped_1s = 0;
        // non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // int result = hapEntries[locus].flipped ? non_flipped_1s : nhaps - non_flipped_1s;
        // return result;
    }


    inline int get_n_c0(int locus)
    {
        if(MISSING_ALLOWED){
            return nhaps - hapEntries[locus].hapbitset->num_1s - hapEntries[locus].xorbitset->num_1s; //hapEntries[locus].xorbitset->num_1s is missing count
            throw "not implemented";
        }else{
            return LOW_MEM? nhaps - hapEntries[locus].hapbitset->num_1s: nhaps - hapEntries[locus].positions.size();

        }
        // // if flip was enabled
        // int non_flipped_1s = 0;
        // non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // int result = hapEntries[locus].flipped ? non_flipped_1s : nhaps - non_flipped_1s;
        // return result;
    }

    inline int get_n_c1(int locus)
    {
        return LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // // if flip was enabled
        // int non_flipped_1s = 0;
        // non_flipped_1s = LOW_MEM? hapEntries[locus].hapbitset->num_1s: hapEntries[locus].positions.size();
        // int result = hapEntries[locus].flipped ? nhaps - non_flipped_1s : non_flipped_1s;
        // return result;
    }

    inline int get_n_c_missing(int locus)
    {
        if(!LOW_MEM || unphased){
            throw "not implemented";
        }
        return hapEntries[locus].xorbitset->num_1s;
    }


    template<typename T>
    void getThreeUnphasedGroups(const std::vector<T>& new1, const std::vector<T>& old1, const std::vector<T>& new2, const std::vector<T>& old2, std::vector<T> g[]) {
        std::vector<T> symdif1;
        std::vector<bool> symdif1_isOld;
        std::vector<T> symdif2;
        std::vector<bool> symdif2_isOld;

        auto it1 = new1.begin();
        auto it2 = old1.begin();
        while (it1 != new1.end() && it2 != old1.end()) {
            if (*it1 < *it2) {
                symdif1.push_back(*it1);
                ++it1;
                symdif1_isOld.push_back(false);
                //in1 not in2  (new1)  
            } else if (*it2 < *it1) {
                symdif1.push_back(*it2);
                ++it2;
                symdif1_isOld.push_back(true);
                //in2 not in1   (old1)
            } else {
                // If the elements are equal, skip them in both vectors
                ++it1;
                ++it2;
                //in both
            }
        }

        while (it1 != new1.end()) { // Copy any remaining elements from vec1
            symdif1.push_back(*it1);
            ++it1;
            symdif1_isOld.push_back(false);

            //in1 not in2
        }
        while (it2 != old1.end()) {  // Copy any remaining elements from vec2
            symdif1.push_back(*it2);
            ++it2;
            symdif1_isOld.push_back(true);

            //in2 not in1
        }


        it1 = new2.begin();
        it2 = old2.begin();
        while (it1 != new2.end() && it2 != old2.end()) {
            if (*it1 < *it2) {
                symdif2.push_back(*it1);
                ++it1;
                symdif2_isOld.push_back(false);
                //in1 not in2  (new1)  
            } else if (*it2 < *it1) {
                symdif2.push_back(*it2);
                ++it2;
                symdif2_isOld.push_back(true);
                //in2 not in1   (old1)
            } else {
                // If the elements are equal, skip them in both vectors
                ++it1;
                ++it2;
                //in both
            }
        }
        while (it1 != new2.end()) { // Copy any remaining elements from vec1
            symdif2.push_back(*it1);
            ++it1;
            symdif2_isOld.push_back(false);
            //in1 not in2
        }
        while (it2 != old2.end()) {  // Copy any remaining elements from vec2
            symdif2.push_back(*it2);
            ++it2;
            symdif2_isOld.push_back(true);
            //in2 not in1
        }


        it1 = symdif1.begin();
        it2 = symdif2.begin();

        auto it1_bool = symdif1_isOld.begin();
        auto it2_bool = symdif2_isOld.begin();

        //std::vector<T> g[3];

        int txor[3][3]; //old->new
        txor[0][0] = 0;
        txor[0][1] = 2;
        txor[0][2] = 1;
        txor[1][0] = 1;
        txor[1][1] = 0;
        txor[1][2] = 2;
        txor[2][0] = 2;
        txor[2][1] = 1;
        txor[2][2] = 0;

        while (it1 != symdif1.end() && it2 != symdif2.end()) {
            if (*it1 < *it2) {
                //not in 2
                if((*it1_bool)==true){ //old 1
                    //new 0
                    g[txor[1][0]].push_back(*it1);
                }else{ //new 1
                    //old 0
                    g[txor[0][1]].push_back(*it1);
                }
                ++it1;
                ++it1_bool;
                //1-2  (new1)  
            } else if (*it2 < *it1) {
                //not in 1
                if((*it2_bool)==true){ //old 2
                    //new 0
                    g[txor[2][0]].push_back(*it2);
                }else{ //new 2
                    //old 0
                    g[txor[0][2]].push_back(*it2);
                }
                ++it2;
                ++it2_bool;
                //2-1   (old1)
            } else {
                // If the elements are equal, skip them in both vectors
                if((*it1_bool)==true){ //old 1
                    //new 2
                    g[txor[1][2]].push_back(*it1);
                }else{ //new 1
                    //old 2
                    g[txor[2][1]].push_back(*it1);
                }
                ++it1;
                ++it2;
                ++it1_bool;
                ++it2_bool;
                //1 and 2
            }
        }
        while (it1 != symdif1.end()) { // Copy any remaining elements from vec1
            //not in 2
            if((*it1_bool)==true){ //old 1
                //new 0
                g[txor[1][0]].push_back(*it1);
            }else{ //new 1
                //old 0
                g[txor[0][1]].push_back(*it1);
            }
            ++it1;
            ++it1_bool;
            //1-2
        }
        while (it2 != symdif2.end()) {  // Copy any remaining elements from vec2
            //not in 1
            if((*it2_bool)==true){ //old 2
                //new 0
                g[txor[2][0]].push_back(*it2);
            }else{ //new 2
                //old 0
                g[txor[0][2]].push_back(*it2);
            }
            ++it2;
            ++it2_bool;
            //2-1
        }
    }

    private:
        bool INIT_SUCCESS = false;

};


#endif
/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#include "selscan-data.h"
#include "selscan-cli.h"
#include <sstream>
#include <queue> 
#include<algorithm>
using namespace std;

void HapData::readHapData_bitset(string filename)
{
    //PHASE 1: Read VCF File to get "nloci", "nhaps" and "skiplist"
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    //int fileStart = fin.tellg();
    string line;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    
    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception

    queue<int> skiplist;
    vector<int> num_1s_per_loci;
    int nloci_before_filter = 0;
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        
        nloci_before_filter++;
        pair<int, int> fo = countFieldsAndOnes(line);
        current_nhaps = fo.first;
        num_1s_per_loci.push_back(fo.second);
        if( SKIP && (fo.second*1.0/current_nhaps < MAF || 1-(fo.second*1.0/current_nhaps) < MAF ) ) {
            skiplist.push(nloci_before_filter-1);
        }

        //cout << "nloci: " << current_nloci << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci_before_filter << " of " << filename << " has " << current_nhaps
                 << ", but the previous line has " << previous_nhaps << ".\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
        
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();


    //PHASE 2: Open VCF File To Load into Data Structure
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps << " haplotypes and " << nloci_before_filter << " loci...\n";
    

    if (unphased){
        initHapData(current_nhaps/2, nloci_before_filter-skiplist.size());
    }
    else{
        initHapData(current_nhaps, nloci_before_filter-skiplist.size());
    }

    this->skipQueue = skiplist; // make a copy
    

    char allele1;
    //int locus_after_filter = 0;
    for (int locus = 0; locus < this->nloci; locus++)
    {
        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();
                getline(fin, line);
                continue;
            }
        }

        vector<bool> current_haps(current_nhaps, false);
        for (int hap = 0; hap < current_nhaps; hap++)
        {
            if(unphased){
                fin >> allele1;
                if (allele1 != '0' && allele1 != '1'){
                    cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                    cerr << allele1 << endl;
                    throw 0;
                }
                current_haps[hap] = (allele1 == '1');
            }
            else{
                fin >> allele1;
                if (allele1 != '0' && allele1 != '1')
                {
                    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                    throw 0;
                }
                if(allele1=='1'){
                    this->hapEntries[locus].hapbitset->num_1s = num_1s_per_loci[locus];
                    this->hapEntries[locus].hapbitset->set_bit(hap);
                }
            }
        }

        if (unphased){
            if (current_nhaps % 2 != 0)
            {
                cerr << "ERROR:  Number of haplotypes must be even for unphased.\n";
                throw 0;
            }

            for (int hap = 0; hap < current_nhaps/2; hap++){ 
                if (hap % 2 == 1){
                    if (current_haps[hap]  && current_haps[hap*2]){ 
                        //data->data[(hap-1)/2][locus] = '2';
                        this->hapEntries[locus].positions2.push_back(hap);
                        //this->hapEntries[locus].positions.push_back(2*hap);
                    }
                    else if ( (current_haps[hap] && !current_haps[hap*2]) || (!current_haps[hap] && current_haps[hap*2]) ){
                        //data->data[(hap-1)/2][locus] = '1';
                        this->hapEntries[locus].positions.push_back(hap);
                        //this->hapEntries[locus].positions.push_back(2*hap+1);

                    }
                }
                // else{
                //     data->data[hap/2][locus] = allele1;
                // }
            }
        }
    }
    fin.close();



    //PHASE 3: XOR
    for(int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
        if(locus_after_filter==0){
            MyBitset* b1 =(hapEntries[locus_after_filter].hapbitset);
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[locus_after_filter].xorbitset->bits[k] = b1->bits[k] ;
            }
            hapEntries[locus_after_filter].xorbitset->num_1s = b1->num_1s;
        }else{
            MyBitset* b1 =(hapEntries[locus_after_filter].hapbitset);
            MyBitset* b2 = (hapEntries[locus_after_filter-1].hapbitset);

            int sum = 0;
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[locus_after_filter].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                sum += __builtin_popcountll(hapEntries[locus_after_filter].xorbitset->bits[k]);
            }
            hapEntries[locus_after_filter].xorbitset->num_1s = sum;
            
        }
    }

    //PHASE 4: FLIP
    for (int locus_after_filter = 0; locus_after_filter < this->nloci; locus_after_filter++){
        if(hapEntries[locus_after_filter].hapbitset->num_1s > nhaps/2){
            hapEntries[locus_after_filter].flipped = true;
            MyBitset* b1;
            b1 = hapEntries[locus_after_filter].hapbitset;
            for(int k = 0; k<b1->nwords; k++){
                b1->bits[k] = ~(b1->bits[k]);   // negate all bits
            }
            for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
                b1->clear_bit(i);       // clear the trailing bits
            }
        }else{
            hapEntries[locus_after_filter].flipped = false;
        }
    }
}


void HapData::readHapDataVCF_bitset(string filename)
{
    igzstream fin;

    vector<int> number_of_1s_per_loci;
    vector<int> number_of_2s_per_loci;
    queue<int> skiplist;

    // Pass 1: Counting so that inititalization is smooth
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    string line;
    unsigned int nloci_before_filtering = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;

    int skipcount = 0;

    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    int num_meta_data_lines = 0;
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            num_meta_data_lines++;
            continue;
        }
        nloci_before_filtering++;
        current_nhaps = countFields(line) - numMapCols;

        /********/
        string junk;
        char allele1, allele2, separator;
        std::stringstream ss(line);
        int number_of_1s = 0;
        int number_of_2s = 0;

        for (int i = 0; i < numMapCols; i++) {
            ss >> junk;
        }
        for (int field = 0; field < current_nhaps; field++)
        {
            ss >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                cerr << allele1 << " " << allele2 << endl;
                throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}

            if(unphased){
                char allele = '0';
                if (allele1 == '1' && allele2 == '1'){
                    //allele = '2';
                    number_of_2s++;
                }
                else if (allele1 == '1' || allele2 == '1'){
                    number_of_1s++;

                }else{
                    // allele = '0' implied;
                }
                //data->data[field][locus] = allele;
            }
            else{
                if(allele1 == '1'){
                    number_of_1s++;
                }
                if(allele2 == '1'){
                    number_of_1s++;
                }
                // data->data[2 * field][locus] = allele1;
                // data->data[2 * field + 1][locus] = allele2;
            }
        }

        int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);

        if ( SKIP && (derived_allele_count*1.0/(current_nhaps*2) < MAF || 1-(derived_allele_count*1.0/(current_nhaps*2)) < MAF ) ) {
            skiplist.push(nloci_before_filtering-1);
            skipcount++;
        } else {
            number_of_1s_per_loci.push_back(number_of_1s);
            number_of_2s_per_loci.push_back(number_of_2s);
        }

        /*********/

        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci_before_filtering << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();

    //Pass 2: Load according to first pass information
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;

    cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering << " loci...\n";
    if(SKIP){
        cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
        (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    }

    string junk;
    char allele1, allele2, separator;
    bool skipLine = false; // to skip metadata lines

    initHapData(nhaps, nloci_before_filtering-skipcount);

    skipQueue = skiplist; 
    unsigned int nloci_after_filtering = 0;

    for (unsigned int locus = 0; locus < nloci_before_filtering; locus++)
    {
        for (int i = 0; i < numMapCols; i++) {
            fin >> junk;
            if (i == 0 && junk[0] == '#') { // to skip metadata lines
                skipLine = true;
                break;
            }
        }

        if (skipLine) { // to skip metadata lines
            getline(fin, junk);
            skipLine = false;
            locus--;
            continue;
        }

        if(!skiplist.empty()){
            if(skiplist.front() == locus){
                skiplist.pop();    
                getline(fin, junk);
                continue;
            }
        }
        
        if(unphased){
            //TODO
        }else{
            hapEntries[nloci_after_filtering].hapbitset->num_1s = number_of_1s_per_loci[nloci_after_filtering];
            if(number_of_1s_per_loci[nloci_after_filtering] > nhaps/2){
                hapEntries[nloci_after_filtering].flipped = true;
            }else{
                hapEntries[nloci_after_filtering].flipped = false;
            }
        }

        for (int field = 0; field <  current_nhaps ; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            separator = junk[1];
            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
            {
                cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                cerr << allele1 << " " << allele2 << endl;
                throw 0;
            }

            //if(separator != '|'){
            //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
            //    throw 0;
            //}
            if(unphased){
                //TODO
            }
            else{
                if(allele1 == '1'){
                    (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field);
                }
                if(allele2 == '1'){
                    (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field + 1);;
                }
            }
        }

        if(nloci_after_filtering==0){
            MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
            int sum = 0;
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ;
            }
            hapEntries[nloci_after_filtering].xorbitset->num_1s = b1->num_1s;
        }else{
            MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
            MyBitset* b2 = (hapEntries[nloci_after_filtering-1].hapbitset);

            int sum = 0;
            for (int k = 0; k < b1->nwords; k++) {
                hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
                sum += __builtin_popcountll(hapEntries[nloci_after_filtering].xorbitset->bits[k]);
            }
            hapEntries[nloci_after_filtering].xorbitset->num_1s = sum;
            
        }
        nloci_after_filtering++;
    }


    //handle fliiped
    for (int locus = 0; locus < nloci_after_filtering; locus++){
        if(hapEntries[locus].flipped){
            MyBitset* b1;
            b1 = hapEntries[locus].hapbitset;
            for(int k = 0; k<b1->nwords; k++){
                    b1->bits[k] = ~(b1->bits[k]);
            }
            for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
                b1->clear_bit(i);
            }
        }
    }
    
    
   //BEGIN TESTING CODES

   //END TESTING CODES

    if(SKIP){
        cerr << "Removed " << skipcount << " low frequency variants.\n";
        (*flog) << "Removed " << skipcount << " low frequency variants.\n";
    }

    fin.close();
}


// void HapData::readHapDataVCF_bitset(string filename)
// {
//     igzstream fin;

//     vector<int> number_of_1s_per_loci;
//     vector<int> number_of_2s_per_loci;
//     queue<int> skiplist;

//     // Pass 1: Counting so that inititalization is smooth
//     cerr << "Opening " << filename << "...\n";
//     fin.open(filename.c_str());

//     if (fin.fail())
//     {
//         cerr << "ERROR: Failed to open " << filename << " for reading.\n";
//         throw 0;
//     }

//     int numMapCols = 9;
//     string line;
//     unsigned int nloci_before_filtering = 0;
//     int previous_nhaps = -1;
//     int current_nhaps = 0;

//     int skipcount = 0;

//     //Counts number of haps (cols) and number of loci (rows)
//     //if any lines differ, send an error message and throw an exception
//     int num_meta_data_lines = 0;
//     while (getline(fin, line))
//     {
//         if (line[0] == '#') {
//             num_meta_data_lines++;
//             continue;
//         }
//         nloci_before_filtering++;
//         current_nhaps = countFields(line) - numMapCols;

//         /********/
//         string junk;
//         char allele1, allele2, separator;
//         std::stringstream ss(line);
//         int number_of_1s = 0;
//         int number_of_2s = 0;

//         for (int i = 0; i < numMapCols; i++) {
//             ss >> junk;
//         }
//         for (int field = 0; field < current_nhaps; field++)
//         {
//             ss >> junk;
//             allele1 = junk[0];
//             separator = junk[1];
//             allele2 = junk[2];
//             if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
//             {
//                 cerr << "ERROR: Alleles must be coded 0/1 only.\n";
//                 cerr << allele1 << " " << allele2 << endl;
//                 throw 0;
//             }

//             //if(separator != '|'){
//             //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
//             //    throw 0;
//             //}

//             if(unphased){
//                 char allele = '0';
//                 if (allele1 == '1' && allele2 == '1'){
//                     //allele = '2';
//                     number_of_2s++;
//                 }
//                 else if (allele1 == '1' || allele2 == '1'){
//                     number_of_1s++;

//                 }else{
//                     // allele = '0' implied;
//                 }
//                 //data->data[field][locus] = allele;
//             }
//             else{
//                 if(allele1 == '1'){
//                     number_of_1s++;
//                 }
//                 if(allele2 == '1'){
//                     number_of_1s++;
//                 }
//                 // data->data[2 * field][locus] = allele1;
//                 // data->data[2 * field + 1][locus] = allele2;
//             }
//         }

//         int derived_allele_count = (unphased? (number_of_1s + number_of_2s*2) : number_of_1s);

//         if ( SKIP && (derived_allele_count*1.0/(current_nhaps*2) < MAF || 1-(derived_allele_count*1.0/(current_nhaps*2)) < MAF ) ) {
//             skiplist.push(nloci_before_filtering-1);
//             skipcount++;
//         } else {
//             number_of_1s_per_loci.push_back(number_of_1s);
//             number_of_2s_per_loci.push_back(number_of_2s);
//         }

//         /*********/

//         if (previous_nhaps < 0)
//         {
//             previous_nhaps = current_nhaps;
//             continue;
//         }
//         else if (previous_nhaps != current_nhaps)
//         {
//             cerr << "ERROR: line " << nloci_before_filtering << " of " << filename << " has " << current_nhaps
//                  << " fields, but the previous line has " << previous_nhaps << " fields.\n";
//             throw 0;
//         }
//         previous_nhaps = current_nhaps;
//     }

//     fin.clear(); // clear error flags
//     //fin.seekg(fileStart);
//     fin.close();

//     //Pass 2: Load according to first pass information
//     fin.open(filename.c_str());

//     if (fin.fail())
//     {
//         cerr << "ERROR: Failed to open " << filename << " for reading.\n";
//         throw 0;
//     }

//     int nhaps = unphased ? (current_nhaps ) : (current_nhaps ) * 2;

//     cerr << "Loading " << nhaps << " haplotypes and " << nloci_before_filtering << " loci...\n";
//     if(SKIP){
//         cerr << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
//         (*flog)  << ARG_SKIP << " set. Removing all variants < " << MAF << ".\n";
    
//     }

//     string junk;
//     char allele1, allele2, separator;
//     bool skipLine = false; // to skip metadata lines

//     initHapData(nhaps, nloci_before_filtering-skipcount);

//     skipQueue = skiplist; 
//     unsigned int nloci_after_filtering = 0;

//     string prev_loc_str = "";
//     string curr_loc_str = "";
//     for (unsigned int locus = 0; locus < nloci_before_filtering; locus++)
//     {
//         curr_loc_str = "";
//         for (int i = 0; i < numMapCols; i++) {
//             fin >> junk;
//             if (i == 0 && junk[0] == '#') { // to skip metadata lines
//                 skipLine = true;
//                 break;
//             }
//         }

//         if (skipLine) { // to skip metadata lines
//             getline(fin, junk);
//             skipLine = false;
//             locus--;
//             continue;
//         }

//         if(!skiplist.empty()){
//             if(skiplist.front() == locus){
//                 skiplist.pop();    
//                 getline(fin, junk);
//                 continue;
//             }
//         }
        
//         if(unphased){
//             if(benchmark_flag == "XOR"){
//                 hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]+number_of_2s_per_loci[nloci_after_filtering]);
//             }
//             //hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]+number_of_2s_per_loci[nloci_after_filtering]);

//             if(benchmark_flag != "BITSET"){
//                 hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
//                 hapEntries[nloci_after_filtering].positions2.reserve(number_of_2s_per_loci[nloci_after_filtering]);
//             }
//         }else{
//             if(benchmark_flag == "XOR"){
//                 hapEntries[nloci_after_filtering].xors.reserve(number_of_1s_per_loci[nloci_after_filtering]);
//             }

//             if(benchmark_flag != "BITSET"){
//                 hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
//             }
//             if(benchmark_flag == "BITSET"){
//                 hapEntries[nloci_after_filtering].hapbitset->num_1s = number_of_1s_per_loci[nloci_after_filtering];
//                 if(number_of_1s_per_loci[nloci_after_filtering] > nhaps/2){
//                     hapEntries[nloci_after_filtering].flipped = true;
//                 }else{
//                     hapEntries[nloci_after_filtering].flipped = false;
//                 }
//             }
//             // if(number_of_1s_per_loci[nloci_after_filtering] > nhaps/2){
//             //     hapEntries[nloci_after_filtering].flipped = true;
//             //     hapEntries[nloci_after_filtering].positions.reserve(nhaps - number_of_1s_per_loci[nloci_after_filtering]);

//             // }else{
//             //     hapEntries[nloci_after_filtering].flipped = false;
//             //     hapEntries[nloci_after_filtering].positions.reserve(number_of_1s_per_loci[nloci_after_filtering]);
//             // }
            
            
//         }

//         for (int field = 0; field <  current_nhaps ; field++)
//         {
//             fin >> junk;
//             allele1 = junk[0];
//             separator = junk[1];
//             allele2 = junk[2];
//             if ( (allele1 != '0' && allele1 != '1') || (allele2 != '0' && allele2 != '1') )
//             {
//                 cerr << "ERROR: Alleles must be coded 0/1 only.\n";
//                 cerr << allele1 << " " << allele2 << endl;
//                 throw 0;
//             }

//             //if(separator != '|'){
//             //    cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
//             //    throw 0;
//             //}
//             if(unphased){
//                 char allele = '0';
//                 if (allele1 == '1' && allele2 == '1'){
//                     hapEntries[nloci_after_filtering].positions2.push_back(field);
//                     //hapEntries[nloci_after_filtering].positions.push_back(2*field); //10
//                     hapEntries[nloci_after_filtering].count2++;
//                     curr_loc_str += "2";


//                 }
//                 else if (allele1 == '1' || allele2 == '1'){
//                      hapEntries[nloci_after_filtering].positions.push_back(field);
//                      //hapEntries[nloci_after_filtering].positions.push_back(2*field);
//                     hapEntries[nloci_after_filtering].count1++;
//                     //hapEntries[nloci_after_filtering].positions.push_back(2*field+1); //01
//                     curr_loc_str += "1";

//                 }else{
//                     // allele = '0' implied;
//                     curr_loc_str += "0";
//                 }
//                 //data->data[field][locus] = allele;
//             }
//             else{
//                 if(benchmark_flag == "BITSET"){
//                     if(allele1 == '1'){
//                         (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field);
//                     }
//                     if(allele2 == '1'){
//                         (hapEntries[nloci_after_filtering].hapbitset)->set_bit(2 * field + 1);;
//                     }
//                 }

//                 if(benchmark_flag != "BITSET"){
//                     if (hapEntries[nloci_after_filtering].flipped){
//                         if(allele1 == '0'){
//                             hapEntries[nloci_after_filtering].positions.push_back(2 * field);
//                         }
//                         if(allele2 == '0'){
//                             hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
//                         }
//                     }else{
//                         if(allele1 == '1'){
//                         hapEntries[nloci_after_filtering].positions.push_back(2 * field);
//                         }
//                         if(allele2 == '1'){
//                             hapEntries[nloci_after_filtering].positions.push_back(2 * field + 1);
//                         }
//                     }
//                 }
                
                
//                 // data->data[2 * field][locus] = allele1;
//                 // data->data[2 * field + 1][locus] = allele2;
//             }
//         }


//         if(nloci_after_filtering==0){
//             if(benchmark_flag == "XOR"){
//                 vector<unsigned int>& source = hapEntries[nloci_after_filtering].positions;
//                 vector<unsigned int>& destination = hapEntries[nloci_after_filtering].xors;
//                 std::copy(source.begin(), source.end(), destination.begin());

//                 //std::copy(hapEntries[nloci_after_filtering].positions.begin(), hapEntries[nloci_after_filtering].positions.end(), hapEntries[nloci_after_filtering].xors1.begin());
//                 //std::copy(hapEntries[nloci_after_filtering].positions2.begin(), hapEntries[nloci_after_filtering].positions2.end(), hapEntries[nloci_after_filtering].xors2.begin());

//                 hapEntries[nloci_after_filtering].xors1 = hapEntries[nloci_after_filtering].positions;
//                 hapEntries[nloci_after_filtering].xors2 = hapEntries[nloci_after_filtering].positions2;


//             }

//             if(benchmark_flag == "BITSET"){
//                 MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
//                 int sum = 0;
//                 for (int k = 0; k < b1->nwords; k++) {
//                     hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ;
//                 }
//                 hapEntries[nloci_after_filtering].xorbitset->num_1s = b1->num_1s;
//             }
           
           
//         }else{
//             if(benchmark_flag == "BITSET"){
//                 MyBitset* b1 =(hapEntries[nloci_after_filtering].hapbitset);
//                 MyBitset* b2 = (hapEntries[nloci_after_filtering-1].hapbitset);

//                 int sum = 0;
//                 for (int k = 0; k < b1->nwords; k++) {
//                     hapEntries[nloci_after_filtering].xorbitset->bits[k] = b1->bits[k] ^ b2->bits[k];
//                     sum += __builtin_popcountll(hapEntries[nloci_after_filtering].xorbitset->bits[k]);
//                 }
//                 hapEntries[nloci_after_filtering].xorbitset->num_1s = sum;
//             }
//             if(benchmark_flag == "XOR"){
//                 vector<unsigned int>& curr_xor = hapEntries[nloci_after_filtering].xors;
//                 vector<unsigned int>& curr_positions = hapEntries[nloci_after_filtering].positions;
//                 vector<unsigned int>& prev_positions = hapEntries[nloci_after_filtering-1].positions;
//                 std::set_symmetric_difference(curr_positions.begin(), curr_positions.end(),prev_positions.begin(), prev_positions.end(),
//                                 std::back_inserter(curr_xor));
            

                            
//                 //unphased
//                 for(int i=0; i<curr_loc_str.length(); i++){
//                     int sub = curr_loc_str[i]-prev_loc_str[i];
//                     if(sub<0){
//                         sub = 3+sub;
//                     }
//                     if(sub==1){
//                         hapEntries[nloci_after_filtering].xors1.push_back(i);
//                     }else if(sub==2){
//                         hapEntries[nloci_after_filtering].xors2.push_back(i);
//                     }
//                 }

//             }
//         }
//         prev_loc_str = curr_loc_str;
//         nloci_after_filtering++;
//     }

//     if(benchmark_flag!="BITSET"){
//         for (int locus = 0; locus < nloci_after_filtering; locus++){
//             if(hapEntries[locus].flipped){
//                 vector<unsigned int> zero_positions(nhaps - hapEntries[locus].positions.size());

//                 int j = 0;
//                 unsigned int front_one = hapEntries[locus].positions[j++];
//                 for(int i=0; i<nhaps; i++){
//                     if(i==front_one){
//                         front_one = hapEntries[locus].positions[j++];
//                     }else{
//                         zero_positions.push_back(i);
//                     }   
//                 }
//                 hapEntries[locus].positions = zero_positions;
//             }
//         }
//     }

//     if(benchmark_flag=="BITSET"){
//         for (int locus = 0; locus < nloci_after_filtering; locus++){
//             if(hapEntries[locus].flipped){
//                 MyBitset* b1;
//                 b1 = hapEntries[locus].hapbitset;
//                 for(int k = 0; k<b1->nwords; k++){
//                      b1->bits[k] = ~(b1->bits[k]);
//                 }
//                 for(int i = b1->nbits; i<b1->nwords*b1->WORDSZ; i++){
//                     b1->clear_bit(i);
//                 }
//             }
//         }
//     }
    
    
//     //nloci = nloci_before_filtering - skipcount*int(SKIP);
//    // nloci = nloci_after_filtering;

//    /*
//      for (int locus = 0; locus < 1; locus++){
//         hapEntries[locus].hapbitset->print_pos();
//         // for (int j = 0; j < 10; j++){
//         //     cout<<hapEntries[locus].positions[j]<<" ";
//         // }
//         // cout<<endl;
//      }
//     exit(1);
//     */


//     /*
//     //iterate over all_positions vector
//     int newlocus = 0;
//     for (int locus = 0; locus < nloci_before_filtering; locus++)
//     {
//         hapEntries[locus].positions.reserve(all_positions[locus].size());
//         hapEntries[locus].xors.reserve(all_positions[locus].size());

//         double AF1 = all_positions[locus].size()*1.0/nhaps;
//         if(unphased){
//             AF1 = (all_positions[locus].size() + all_positions2[locus].size()*2.0)/nhaps;
//         }

//         //if(!SKIP || min(AF1, 1-AF1) >= MAF){
//             for (unsigned int pos: all_positions[locus])
//             {
//                 hapEntries[newlocus].positions.push_back(pos);
//                 hapEntries[newlocus].xors.push_back(pos);
//                 ++newlocus;
//             }
//             if(unphased){
//                 for (unsigned int pos: all_positions2[locus])
//                 {
//                     hapEntries[newlocus].positions2.push_back(pos);
//                 }
//             }
//         //}


//         if(nloci!=newlocus){
//             cerr << "DEBUG ERROR: nloci != newlocus\n"<< nloci<<" "<<newlocus<<endl;
//             throw 0;
//         }
        
//     }
//     */



//    //BEGIN TESTING CODES

//    //END TESTING CODES

//     if(SKIP){
//         cerr << "Removed " << skipcount << " low frequency variants.\n";
//         (*flog) << "Removed " << skipcount << " low frequency variants.\n";
//     }

//     fin.close();

// }

/** Sets up structure according to nhaps and nloci
 * 
*/
void HapData::initHapData_bitset(int nhaps, unsigned int nloci)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    this->hapEntries = new struct HapEntry[nloci];
    this->nhaps = nhaps;
    this->nloci = nloci;

    if(benchmark_flag == "BITSET"){
        for (unsigned int j = 0; j < nloci; j++){
            hapEntries[j].hapbitset = new MyBitset(nhaps);
            hapEntries[j].xorbitset = new MyBitset(nhaps);

        }
    }

    //this->data = new char *[nhaps];
    // for (unsigned int i = 0; i < nhaps; i++)
    // {
    //     data->data[i] = new char[nloci];
    //     for (unsigned int j = 0; j < nloci; j++)
    //     {
    //         data->data[i][j] = MISSING_CHAR;
    //     }
    // }

    
    // for (unsigned int j = 0; j < nloci; j++)
    // {
    //     //hapEntries[j].xors.reserve(nhaps);
    //     //hapEntries[j].positions.reserve(nhaps);
    // }
    
}

void HapData::releaseHapData_bitset()
{
    if (hapEntries == NULL) return;
    for (int i = 0; i < this->nhaps; i++)
    {
        //delete [] hapEntries[i];
        //delete hapEntries[i];
        hapEntries[i].xors.clear();
        hapEntries[i].positions.clear();    
    }

    delete [] hapEntries;

    hapEntries = NULL;
    this->nhaps = -9;
    this->nloci = -9;
    //data = NULL;
    return;
}


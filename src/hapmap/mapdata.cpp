#include "mapdata.h"
#include <iostream>
#include <string>
#include "../gzstream.h"
#include <queue>
#include <sstream>
#include <limits>

using namespace std;

/**
 * reads in map data and also does basic checks on integrity of format
 * @returns a populated MapData structure if successful
 * @throws an exception if not successful
 */
void MapData::readMapData(string filename, size_t expected_loci, bool USE_PMAP, queue<size_t>& skip_queue)
{
    igzstream fin;
    cerr << "Opening " << filename << " to read map data...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    size_t nloci_before_filter = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci_before_filter++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci_before_filter << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            exit(EXIT_FAILURE);
        }
    }

    if (nloci_before_filter - (skip_queue.size()) != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci_before_filter-skip_queue.size() << ".\n";
        exit(EXIT_FAILURE);
    }

    fin.clear(); // clear error flags
    fin.close();


    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loading map data for " << nloci_before_filter-skip_queue.size() << " loci\n";
    initMapData(nloci_before_filter-skip_queue.size());

    string chr;
    size_t locus_after_filter = 0;
    for (size_t locus_before_filter = 0; locus_before_filter < nloci_before_filter; locus_before_filter++)
    {
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus_before_filter){
                skip_queue.pop();
                string junk;
                fin >> junk;
                fin >> junk;
                fin >>  junk;
                fin >>  junk;
                getline(fin, line);
                continue;
            }
        }

        fin >> mapEntries[locus_after_filter].chr;
        fin >> mapEntries[locus_after_filter].locusName;
        fin >> mapEntries[locus_after_filter].geneticPos;
        fin >> mapEntries[locus_after_filter].physicalPos;
        

        //locus_query_map[mapEntries[locus_after_filter].locusName] = locus_after_filter;
        mapEntries[locus_after_filter].locId = locus_before_filter;

        double Mb = 1000000.0;

       // mapEntries[locus_after_filter].geneticPos = (long double)(mapEntries[locus_after_filter].physicalPos/Mb);
        if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos/Mb);

        locus_after_filter++;
        getline(fin, line);
    }

    fin.close();
}

void MapData::readMapDataTPED(string filename, size_t expected_loci, int expected_haps, bool USE_PMAP, queue<size_t>& skip_queue)
{   
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    size_t nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols + expected_haps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols + expected_haps << ".\n";
            exit(EXIT_FAILURE);
        }
    }

    if (nloci - (skip_queue.size()) != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci - skip_queue.size()  << ".\n";
        exit(EXIT_FAILURE);
    }

    fin.clear(); // clear error flags
    fin.close();



    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loading map data for " << nloci-skip_queue.size() << " loci\n";

    initMapData(nloci-skip_queue.size());
    
    double Mb = 1000000.0;
    
    string chr;
    size_t locus_after_filter = 0;
    for (size_t locus = 0; locus < nloci; locus++)  // locus = locus_before_filter
    {
        if(!skip_queue.empty()){
            if(skip_queue.front()==locus){
                skip_queue.pop();
                getline(fin, line);
                continue;
            }
        }
        fin >> mapEntries[locus_after_filter].chr;
        fin >> mapEntries[locus_after_filter].locusName;
        fin >> mapEntries[locus_after_filter].geneticPos;
        fin >> mapEntries[locus_after_filter].physicalPos;
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip rest of the line

        if (USE_PMAP) mapEntries[locus_after_filter].geneticPos = double(mapEntries[locus_after_filter].physicalPos)/Mb;
    
        locus_after_filter++;
    }

    fin.close();
}

// void MapData::readMapDataVCF(string filename, size_t expected_loci, queue <size_t>& skip_queue) {
//     igzstream fin;
//     cerr << "Opening " << filename << " to read map data...\n";
//     fin.open(filename.c_str());

//     if (fin.fail())
//     {
//         cerr << "ERROR: Failed to open " << filename << " for reading.\n";
//         exit(EXIT_FAILURE);
//     }

//     string line;
//     size_t nloci_before_filter = 0;
//     int numCommentedLines = 0;
//     while (getline(fin, line))
//     {
//         if (line[0] == '#') {
//             numCommentedLines++;
//         }
//         else {
//             nloci_before_filter++;
//         }
//     }

//     if (nloci_before_filter-skip_queue.size() != expected_loci)
//     {
//         cerr << "ERROR: Expected " << expected_loci << " loci in file but found " << nloci_before_filter-skip_queue.size() << ".\n";
//         exit(EXIT_FAILURE);
//     }

//     fin.clear(); // clear error flags
//     //fin.seekg(fileStart);
//     fin.close();
//     fin.open(filename.c_str());

//     if (fin.fail())
//     {
//         cerr << "ERROR: Failed to open " << filename << " for reading.\n";
//         exit(EXIT_FAILURE);
//     }

//     cerr << "Loading map data for " << nloci_before_filter-skip_queue.size() << " loci\n";

//     for (int i = 0; i < numCommentedLines; i++) {
//         getline(fin, line);
//     }

//     this->initMapData(nloci_before_filter-skip_queue.size()); 
    
//     double Mb = 1000000.0;
    
//     string chr;

//     //int n_chromosomes_included = 0;

//     size_t locus_after_filter = 0;
//     for (size_t locus_before_filter = 0; locus_before_filter < nloci_before_filter; locus_before_filter++)
//     {
//         if(!skip_queue.empty()){
//             if(skip_queue.front()==locus_before_filter){
//                 skip_queue.pop();
//                 string junk;
//                 fin >> junk;
//                 fin >>  junk;
//                 fin >>  junk;
//                 getline(fin, line);
//                 continue;
//             }
//         }
        
//         fin >> mapEntries[locus_after_filter].chr;
//         fin >> mapEntries[locus_after_filter].physicalPos;
//         fin >> mapEntries[locus_after_filter].locusName;
//         //locus_query_map[mapEntries[locus].locusName] = locus;
//         //locus_query_map[to_string(mapEntries[locus_after_filter].physicalPos)] = locus_after_filter;
//         mapEntries[locus_after_filter].locId = locus_before_filter;

//         //if exists in map do nothing else insert
//         // if(chr_list.find(mapEntries[locus_after_filter].chr) == chr_list.end()){
//         //     chr_list[mapEntries[locus_after_filter].chr] = n_chromosomes_included++;
//         // }

//         mapEntries[locus_after_filter].geneticPos = (long double)(mapEntries[locus_after_filter].physicalPos)/Mb;
//         //cout<<mapEntries[locus_after_filter].geneticPos<<" "<<mapEntries[locus_after_filter].physicalPos<<endl;
//         getline(fin, line);
//         locus_after_filter++;
//     }

//     // if(chr_list.size() == 0){
//     //     cerr<<"ERROR: No chromosomes found in map file.\n";
//     //     throw 0;
//     // }
    
//     fin.close();
// }

void MapData::readMapDataVCF(string filename, size_t expected_loci, queue<size_t>& skip_queue) {
    igzstream fin;
    cerr << "Opening " << filename << " to read map data...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    string line;
    size_t nloci_before_filter = 0;

    // ------------------------------------------------------------
    // PASS 1:
    // count non-header loci only
    // ------------------------------------------------------------
    while (getline(fin, line))
    {
        if (!line.empty() && line[0] != '#') {
            nloci_before_filter++;
        }
    }

    if (nloci_before_filter - skip_queue.size() != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci
             << " loci in file but found "
             << nloci_before_filter - skip_queue.size() << ".\n";
        exit(EXIT_FAILURE);
    }

    fin.clear();
    fin.close();

    fin.open(filename.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        exit(EXIT_FAILURE);
    }

    cerr << "Loading map data for " << nloci_before_filter - skip_queue.size() << " loci\n";

    this->initMapData(nloci_before_filter - skip_queue.size());

    const long double Mb = 1000000.0L;

    // Use a local copy since we consume the queue during loading
    queue<size_t> local_skip_queue = skip_queue;

    // ------------------------------------------------------------
    // Small helpers for fast manual parsing
    // ------------------------------------------------------------
    auto skip_spaces = [](const char* p, const char* end) -> const char* {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        return p;
    };

    auto find_field_end = [](const char* p, const char* end) -> const char* {
        while (p < end && *p != ' ' && *p != '\t') ++p;
        return p;
    };

    auto parse_int_fast = [](const char* s, const char* e) -> int {
        int x = 0;
        for (const char* p = s; p < e; ++p) {
            x = x * 10 + (*p - '0');
        }
        return x;
    };

    size_t locus_before_filter = 0;
    size_t locus_after_filter = 0;

    // ------------------------------------------------------------
    // PASS 2:
    // read kept loci only
    // ------------------------------------------------------------
    while (getline(fin, line))
    {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        if (!local_skip_queue.empty() && local_skip_queue.front() == locus_before_filter) {
            local_skip_queue.pop();
            locus_before_filter++;
            continue;
        }

        const char* p = line.c_str();
        const char* end = p + line.size();

        // Field 1: CHROM
        p = skip_spaces(p, end);
        const char* f1s = p;
        const char* f1e = find_field_end(p, end);
        mapEntries[locus_after_filter].chr.assign(f1s, f1e - f1s);

        // Field 2: POS
        p = skip_spaces(f1e, end);
        const char* f2s = p;
        const char* f2e = find_field_end(p, end);
        mapEntries[locus_after_filter].physicalPos = parse_int_fast(f2s, f2e);

        // Field 3: ID
        p = skip_spaces(f2e, end);
        const char* f3s = p;
        const char* f3e = find_field_end(p, end);
        mapEntries[locus_after_filter].locusName.assign(f3s, f3e - f3s);

        // original unfiltered locus index
        mapEntries[locus_after_filter].locId = locus_before_filter;

        // simple genetic position proxy in Mb
        mapEntries[locus_after_filter].geneticPos =
            (long double)(mapEntries[locus_after_filter].physicalPos) / Mb;

        locus_before_filter++;
        locus_after_filter++;
    }

    fin.close();
}

void MapData::initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        exit(EXIT_FAILURE);
    }

    mapEntries = new struct MapEntry[nloci];
    this->nloci = nloci;
}

void MapData::releaseMapData()
{
    if (mapEntries == NULL) return;
    this->nloci = -9;

    delete [] mapEntries;
    mapEntries = NULL;
    return;
}
/* selscan -- a program to calculate EHH-based scans for positive selection in genomes
   Copyright (C) 2014-2024  Zachary A Szpiech
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
#ifndef __SELSCAN_CLI_H__
#define __SELSCAN_CLI_H__


#include <thread>


#include "param_t.h"
#include <string>

using namespace std;







#define LOG(msg) do { \
        std::cerr << msg << std::endl; \
        if (flog) (*flog) << msg << std::endl; \
    } while (0)
    
    
    #ifdef DEBUG
        #define DBG(x) do { std::cerr << "[DEBUG] " << x << std::endl; } while (0)
    #else
        #define DBG(x) do {} while (0)
    #endif
    
    #include <stdexcept>
    
    #ifdef DEBUG
        #define HANDLE_ERROR(msg) do { \
            std::cerr << "[DEBUG ERROR] " << msg << std::endl; \
            *flog << "[ERROR] " << msg << std::endl; \
            throw 1; \
        } while(0)
        #define HANDLE_ERROR_NOLOG(msg) do { \
            std::cerr << "[ERROR] " << msg << std::endl; \
            throw 1; \
        } while(0)
    #else
        #define HANDLE_ERROR(msg) do { \
            std::cerr << "[ERROR] " << msg << std::endl; \
            *flog << "[ERROR] " << msg << std::endl; \
            exit(EXIT_FAILURE); \
        } while(0)
        #define HANDLE_ERROR_NOLOG(msg) do { \
                std::cerr << "[ERROR] " << msg << std::endl; \
                exit(EXIT_FAILURE); \
            } while(0)
    #endif





const string VERSION = "2.1.0";

const string PREAMBLE = "\nselscan v" + VERSION + " -- a program to calculate EHH-based scans for positive selection in genomes.\n\
Source code and binaries can be found at <https://www.github.com/szpiech/selscan>.\n\
\n\
selscan currently implements EHH, iHS, XP-EHH, and nSL.\n\
\n\
Citations:\n\
\n\
selscan: ZA Szpiech and RD Hernandez (2014) MBE 31: 2824-2827.\n\
         ZA Szpiech (2024) Bioinformatics 40: btae006.\n\
iHH12: R Torres et al. (2018) PLoS Genetics 15: e1007898.\n\
       N Garud et al. (2015) PLoS Genetics 11: 1–32.\n\
nSL: A Ferrer-Admetlla et al. (2014) MBE 31: 1275-1291.\n\
XP-nSL: Szpiech et al. (2021) Evol Lett 5: 408-421.\n\
XP-EHH: PC Sabeti et al. (2007) Nature 449: 913–918.\n\
        K Wagh et al. (2012) PloS ONE 7: e44751.\n\
iHS: BF Voight et al. (2006) PLoS Biology 4: e72.\n\
EHH: PC Sabeti et al. (2002) Nature 419: 832–837.\n\
\n\
To calculate EHH:\n\
\n\
./selscan --ehh <locusID> --vcf <vcf> --map <mapfile> --out <outfile>\n\
\n\
To calculate iHS:\n\
\n\
./selscan --ihs --vcf <vcf> --map <mapfile> --out <outfile>\n\
\n\
To calculate nSL:\n\
\n\
./selscan --nsl --vcf <vcf> --out <outfile>\n\
\n\
To calculate XP-nSL:\n\
\n\
./selscan --xpnsl --vcf <vcf> --vcf-ref <vcf> --out <outfile>\n\
\n\
To calculate iHH12:\n\
\n\
./selscan --ihh12 --vcf <vcf> --map <mapfile> --out <outfile>\n\
\n\
To calculate XP-EHH:\n\
\n\
./selscan --xpehh --vcf <vcf> --vcf-ref <vcf> --map <mapfile> --out <outfile>\n";

const string ARG_THREAD = "--threads";
const int DEFAULT_THREAD = std::thread::hardware_concurrency();
const string HELP_THREAD = "The number of threads to spawn during the calculation.\n\
\tPartitions loci across threads.";

const string ARG_FILENAME_POP1_TPED = "--tped";
const string DEFAULT_FILENAME_POP1_TPED = "__hapfile1";
const string HELP_FILENAME_POP1_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1";

const string ARG_FILENAME_POP2_TPED = "--tped-ref";
const string DEFAULT_FILENAME_POP2_TPED = "__hapfile2";
const string HELP_FILENAME_POP2_TPED = "A TPED file containing haplotype and map data.\n\
\tVariants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP calculations and should contain the same number\n\
\tof loci as the query population. Ignored otherwise.";

const string ARG_FILENAME_POP1_VCF = "--vcf";
const string DEFAULT_FILENAME_POP1_VCF = "__hapfile1";
const string HELP_FILENAME_POP1_VCF = "A VCF file containing haplotype data.\n\
\tA map file must be specified with --map.";

const string ARG_FILENAME_POP2_VCF = "--vcf-ref";
const string DEFAULT_FILENAME_POP2_VCF = "__hapfile2";
const string HELP_FILENAME_POP2_VCF = "A VCF file containing haplotype and map data.\n\
\tVariants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP calculations and should contain the same number\n\
\tof loci as the query population. Ignored otherwise.";

const string ARG_FILENAME_POP1 = "--hap";
const string DEFAULT_FILENAME_POP1 = "__hapfile1";
const string HELP_FILENAME_POP1 = "A hapfile with one row per haplotype, and one column per\n\
\tvariant. Variants should be coded 0/1";

const string ARG_FILENAME_POP2 = "--ref";
const string DEFAULT_FILENAME_POP2 = "__hapfile2";
const string HELP_FILENAME_POP2 = "A hapfile with one row per haplotype, and one column per\n\
\tvariant. Variants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP calculations.  Ignored otherwise.";


const string ARG_FILENAME_POP1_THAP = "--thap";
const string DEFAULT_FILENAME_POP1_THAP = "__thapfile1";
const string HELP_FILENAME_POP1_THAP = "A hapfile in IMPUTE hap format with one column per haplotype, and one row per\n\
\tvariant. Variants should be coded 0/1";

const string ARG_FILENAME_POP2_THAP = "--thap-ref";
const string DEFAULT_FILENAME_POP2_THAP = "__thapfile2";
const string HELP_FILENAME_POP2_THAP = "A hapfile in IMPUTE hap format with one column per haplotype, and row per\n\
\tvariant. Variants should be coded 0/1. This is the 'reference'\n\
\tpopulation for XP calculations.  Ignored otherwise.";


const string ARG_FILENAME_MAP = "--map";
const string DEFAULT_FILENAME_MAP = "__mapfile";
const string HELP_FILENAME_MAP = "A mapfile with one row per variant site.\n\
\tFormatted <chr#> <locusID> <genetic pos> <physical pos>.";

const string ARG_PMAP = "--pmap";
const bool DEFAULT_PMAP = false;
const string HELP_PMAP = "Use physical map instead of a genetic map.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The basename for all output files.";

const string ARG_CUTOFF = "--cutoff";
const double DEFAULT_CUTOFF = 0.05;
const string HELP_CUTOFF = "The EHH decay cutoff.";

const string ARG_MAX_GAP = "--max-gap";
const int DEFAULT_MAX_GAP = 200000;
const string HELP_MAX_GAP = "Maximum allowed gap in bp between two snps.";

const string ARG_GAP_SCALE = "--gap-scale";
const int DEFAULT_GAP_SCALE = 20000;
const string HELP_GAP_SCALE = "Gap scale parameter in bp. If a gap is encountered between\n\
\ttwo snps > GAP_SCALE and < MAX_GAP, then the genetic distance is\n\
\tscaled by GAP_SCALE/GAP.";

const string ARG_IHS = "--ihs";
const bool DEFAULT_IHS = false;
const string HELP_IHS = "Set this flag to calculate iHS.";

const string ARG_IHS_DETAILED = "--ihs-detail";
const bool DEFAULT_IHS_DETAILED = false;
const string HELP_IHS_DETAILED = "Set this flag to write out left and right iHH scores for '1' and '0' in addition to iHS.";

const string ARG_UNPHASED = "--unphased";
const bool DEFAULT_UNPHASED = false;
const string HELP_UNPHASED = "Set this flag to use multilocus genotypes.";

const string ARG_XPNSL = "--xpnsl";
const bool DEFAULT_XPNSL = false;
const string HELP_XPNSL = "Set this flag to calculate XP-nSL.";

const string ARG_NSL = "--nsl";
const bool DEFAULT_NSL = false;
const string HELP_NSL = "Set this flag to calculate nSL.";

const string ARG_SOFT = "--ihh12";
const bool DEFAULT_SOFT = false;
const string HELP_SOFT = "Set this flag to calculate iHH12.";

// const string ARG_SOFT_K = "--k";
// const int DEFAULT_SOFT_K = 2;
// const string HELP_SOFT_K = "Specify K to compute for iHH1K.";

const string ARG_XP = "--xpehh";
const bool DEFAULT_XP = false;
const string HELP_XP = "Set this flag to calculate XP-EHH.";

const string ARG_ALT = "--alt";
const bool DEFAULT_ALT = false;
const string HELP_ALT = "Set this flag to calculate homozygosity based on the sum of the\n\
\tsquared haplotype frequencies in the observed data instead of using\n\
\tbinomial coefficients.";

const string ARG_WAGH = "--wagh";
const bool DEFAULT_WAGH = false;
const string HELP_WAGH = "Set this flag to calculate XP-EHH using definition of EHH which\n\
\tseparates core SNP alleles in the denominator.";

const string ARG_MAF = "--maf";
const double DEFAULT_MAF = 0.05;
const string HELP_MAF = "If a site has a MAF below this value, the program will not use\n\
\tit as a core snp.";

const string ARG_EHH = "--ehh";
const string DEFAULT_EHH = "__NO_LOCUS__";
const string HELP_EHH = "Calculate EHH of the '1' and '0' haplotypes at the specified\n\
\tlocus. Output: <physical dist> <genetic dist> <'1' EHH> <'0' EHH>";



// const string ARG_EHHS = "--ehhs";
// const bool DEFAULT_EHHS = false;
// const string HELP_EHHS = "Calculate EHH of the '1' and '0' haplotypes at all loci. Outputs separate files for different loci. \n Output: <physical dist> <genetic dist> <'1' EHH> <'0' EHH>";



const string ARG_EHH12 = "--ehh12";
const string DEFAULT_EHH12 = "__NO_LOCUS__";
const string HELP_EHH12 = "Calculate EHH12 of the '1' and '0' haplotypes at the specified\n\
\tlocus. Output: <physical dist> <genetic dist> <'1' EHH> <'0' EHH>";

const string ARG_QWIN = "--ehh-win";
const int DEFAULT_QWIN = 100000;
const string HELP_QWIN = "When calculating EHH, this is the length of the window in bp\n\
\tin each direction from the query locus.";

const string ARG_MAX_EXTEND = "--max-extend";
const int DEFAULT_MAX_EXTEND = 1000000;         //1000000
const string HELP_MAX_EXTEND = "The maximum distance an EHH decay curve is allowed to extend from the core.\n\
\tSet <= 0 for no restriction.";

const string ARG_MAX_EXTEND_NSL = "--max-extend-nsl";
const int DEFAULT_MAX_EXTEND_NSL = 100;
const string HELP_MAX_EXTEND_NSL = "The maximum distance an nSL haplotype is allowed to extend from the core.\n\
\tSet <= 0 for no restriction.";

const string ARG_SKIP = "--skip-low-freq";
const bool DEFAULT_SKIP = false;
const string HELP_SKIP = "**This flag is now on by default. If you want to include low frequency variants\n\
in the construction of your haplotypes please use the --keep-low-freq flag.";

const string ARG_KEEP = "--keep-low-freq";
const bool DEFAULT_KEEP = false;
const string HELP_KEEP = "Include low frequency variants in the construction of your haplotypes.";

const string ARG_TRUNC = "--trunc-ok";
const bool DEFAULT_TRUNC = false;
const string HELP_TRUNC = "If an EHH decay reaches the end of a sequence before reaching the cutoff,\n\
\tintegrate the curve anyway.\n\
\tNormal function is to disregard the score for that core.";

const string ARG_PI = "--pi";
const bool DEFAULT_PI = false;
const string HELP_PI = "Set this flag to calculate mean pairwise sequence difference in a sliding window.";

const string ARG_PI_WIN = "--pi-win";
const int DEFAULT_PI_WIN = 100;
const string HELP_PI_WIN = "Sliding window size in bp for calculating pi.";

// @UNCOMMENT
// const string ARG_MISSING_FLAG = "--missing";
// const bool DEFAULT_MISSING_FLAG = false;
// const string HELP_MISSING_FLAG = "Set this flag to allow missing data.";

// const string ARG_MULTI_CHR = "--multi-chr";
// const string DEFAULT_MULTI_CHR = "__NO_CHR__";
// const string HELP_MULTI_CHR = "Comma-separated list of chromosomes to include in the analysis.";

// @UNCOMMENT
// const string ARG_IMPUTE_FLAG = "--impute";
// const bool DEFAULT_IMPUTE_FLAG = false;
// const string HELP_IMPUTE_FLAG = "Set this flag to allow imputing missing data using new algorithm.";

const string ARG_MULTI_PARAMS = "--multi-param";
const string DEFAULT_MULTI_PARAMS = "__jsonFile";
const string HELP_MULTI_PARAMS = "Specify a JSON file with multiple parameter sets.\n\
\tEach set should match the structure of command-line arguments.\n\
\tThe program will run the analysis for each set, generating separate outputs.\n\
\tUseful for batch processing and exploring different configurations.\n";

void initalizeParameters(param_t &params,int argc, char *argv[]);
void checkParameters(param_main &params);
void getBaseParamFromCmdLine(param_t& params, param_main &p);
void jsonParse(param_main &base_p, vector<param_main> &multi_params);
//bool checkParameters(param_t &params,int argc, char *argv[]);

#endif
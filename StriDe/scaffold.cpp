#include "SGUtil.h"
#include "scaffold.h"
#include "SGACommon.h"
#include "SequenceProcessFramework.h"
#include "ScaffoldProcess.h"
#include "SGWalk.h"
#define SUBPROGRAM "scaffold"

// Getopt
static const char *SCAFFOLD_USAGE_MESSAGE =
"Usage: StriDe/stride scaffold {option} MatePairRead.fa \n\n"
"      option:\n"
"      -g, --asqg=NAME                  ASQGFILE.(skip mapping) \n"
"      -i, --insertSize=N               The insert size of matepair. (default: 3000)\n"
// "      -k, --kmer-length=N              The length of the kmer to use. (default: 31)\n"
"      -n  --npair-threshold=N          The threshold of npair. (default: 1)\n"
"      -p, --prefix=NAME                prefix of FM-index of contigs. (bwt, rbwt, sai, rsai)\n"
// "      -q, --prefix2=NAME               prefix of FM-index of pair-end reads. (bwt, rbwt, sai, rsai)\n"
"      -t, --threads=NUM                use NUM threads. (default: 1)\n"
"      -o, --out-prefix=NAME            use NAME as the prefix of the output files. \n"
"      --help                           display this help and exit\n\n";

namespace opt
{
    static std::string asqgFile;
    static size_t insertSize = 3000;
	// static size_t kmerLength = 31;
	static std::string prefixContigs;
    static std::string prefixPair_end;
    static size_t npairThreshold = 1;
    static int threads = 1;
    static std::string outFile;    
    static std::string readsFile;
    
    static BWT* pContigsBWT =NULL;
    static BWT* pContigsRBWT =NULL;
    static SampledSuffixArray* pContigsSSA = NULL;
    
    // static BWT* pPair_endBWT =NULL;
    // static BWT* pPair_endRBWT =NULL;
    // static SampledSuffixArray* pPair_endSSA = NULL;
//---------------------------------------------
    static size_t algo = 1;
    static size_t repeatCutoff = 1;
    static bool skipMapping = false;
}

static const char* shortopts = "g:i:k:n:o:p:t:a:r:"; // q:
enum {OPT_HELP = 1};
static const struct option longopts[] = {    
    { "asqg",               required_argument, NULL, 'g' },
    { "insertSize",         required_argument, NULL, 'i' },    
	// { "kmer-length",        required_argument, NULL, 'k' },
    { "npairThreshold-threshold",    required_argument, NULL, 'n' },
    { "out-prefix",         required_argument, NULL, 'o' },
	{ "prefix1",            required_argument, NULL, 'p' },
    // { "prefix2",            required_argument, NULL, 'q' },
    { "threads",            required_argument, NULL, 't' }, 
    { "help",               no_argument,       NULL, OPT_HELP },
//-------------------------------------------------------------
    // { "algo",               required_argument, NULL, 'a' },
    // { "repeat-threshold",   required_argument, NULL, 'r' },  
    // { "skip-mapping",       no_argument,       NULL, 's' },
	{ NULL, 0, NULL, 0 }
};

// Main
int scaffoldMain(int argc, char** argv)
{
	Timer* pTimer = new Timer("StriDe scaffold");
	parseScaffoldOptions(argc, argv); 
    scaffold();
	delete pTimer;
	return 0;
}
int scaffold()
{
    std::cout << "Start to scaffold.\n";
    StringGraph* pGraph; 
    if(opt::skipMapping == false) 
    {
        std::cout << "[ Loading contigs to form a no edge graph ]\n";
        pGraph = SGUtil::loadFASTA(opt::prefixContigs + ".fa");
        #pragma omp parallel
        {        		        
            #pragma omp single nowait
            {
                std::cout << "[ Loading contigs BWT ]\n";
                opt::pContigsBWT = new BWT(opt::prefixContigs + BWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
            }
            #pragma omp single nowait
            {
                std::cout << "[ Loading contigs RBWT ]\n";
                opt::pContigsRBWT = new BWT(opt::prefixContigs + RBWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
            }
            #pragma omp single nowait
            {
                std::cout << "[ Loading contigs SAI ]\n";
                opt::pContigsSSA = new SampledSuffixArray(opt::prefixContigs + SAI_EXT, SSA_FT_SAI);
            }
            // #pragma omp single nowait
            // {
                // std::cout << "[ Loading pair-end BWT ]\n";
                // opt::pPair_endBWT = new BWT(opt::prefixPair_end + BWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
            // }
            // #pragma omp single nowait
            // {
                // std::cout << "[ Loading pair-end RBWT ]\n";
                // opt::pPair_endRBWT = new BWT(opt::prefixPair_end + RBWT_EXT, BWT::DEFAULT_SAMPLE_RATE_SMALL);
            // }
            // #pragma omp single nowait
            // {
                // std::cout << "[ Loading pair-end SAI ]\n";
                // opt::pPair_endSSA = new SampledSuffixArray(opt::prefixPair_end + SAI_EXT, SSA_FT_SAI);
            // }
        }           
        ScaffoldParameters params;    
        BWTIndexSet indexContigsSet;
        indexContigsSet.pBWT = opt::pContigsBWT;
        indexContigsSet.pRBWT = opt::pContigsRBWT;
        indexContigsSet.pSSA = opt::pContigsSSA;
        params.ContigsIndices = indexContigsSet;        
        // BWTIndexSet indexPair_endSet;
        // indexPair_endSet.pBWT = opt::pPair_endBWT;
        // indexPair_endSet.pRBWT = opt::pPair_endRBWT;
        // indexPair_endSet.pSSA = opt::pPair_endSSA;
        // params.Pair_endIndices = indexPair_endSet;    
        // params.kmerLength = opt::kmerLength;
        
        // params.kmerMatchMap = new KmerMatchMap;        
        params.kmerMatchMap2 = new KmerMatchMap2();
        
        // Load contig ReadTable
        ReadInfoTable* pContigRIT = new ReadInfoTable(opt::prefixContigs + ".fa");    
        params.pContigRIT = pContigRIT;
        //--------------------------------------    
        // params.algo = opt::algo;
        // params.repeatCutoff = opt::repeatCutoff;
    
        std::cout << "{My Parameter}" << std::endl;
        // std::cout << "  algo = " << params.algo << std::endl;
        std::cout << "  prefixContigs = " << opt::prefixContigs << std::endl; 
        // std::cout << "  prefixPair_end = " << opt::prefixPair_end << std::endl;       
        // std::cout << "  kmerLength = " << params.kmerLength << std::endl;       
        std::cout << "  threads = " << opt::threads << std::endl;       
        // std::cout << "  repeatCutoff = " << params.repeatCutoff << std::endl;            

        ScaffoldPostProcess* pPostProcessor = new ScaffoldPostProcess(params, pGraph);    
        if(opt::threads <= 1)
        {
            // Serial mode
            ScaffoldProcess processor(params);
            SequenceProcessFramework::processSequencesSerial<SequenceWorkItemPair,
                                                             ScaffoldResult,
                                                             ScaffoldProcess, 
                                                             ScaffoldPostProcess>(opt::readsFile, &processor, pPostProcessor);
        }
        else
        {
            // Parallel mode
            std::vector<ScaffoldProcess*> processorVector;
            for(int i = 0; i < opt::threads; ++i)
            {
                ScaffoldProcess* pProcessor = new ScaffoldProcess(params);
                processorVector.push_back(pProcessor);
            }
            SequenceProcessFramework::processSequencesParallel<SequenceWorkItemPair, 
                                                               ScaffoldResult, 
                                                               ScaffoldProcess, 
                                                               ScaffoldPostProcess>(opt::readsFile, processorVector, pPostProcessor);
            for(int i = 0; i < opt::threads; ++i)
                delete processorVector[i];
        }    
        std::cout << "Mapping Done." << std::endl;
        delete pPostProcessor;              
        //Make a ASQG of scaffold graph
        pGraph->writeScaffoldASQG(opt::prefixContigs + "-scaffold.asqg.gz");
    }
    else
    {
        std::cout << "Skip mapping stage.\n";        
        pGraph = SGUtil::loadScaffoldASQG(opt::asqgFile); 
        std::cout << "[Stats] Input scaffold graph:\n";  
        std::cout << "  asqgFile = " << opt::asqgFile << std::endl;    
    }        
    std::cout << "  insertSize = " << opt::insertSize << std::endl;       
    std::cout << "  npairThreshold = " << opt::npairThreshold << std::endl;        
    SGGraphStatsVisitor statsVisit;
    pGraph->visit(statsVisit);        
    pGraph->writeScaffoldDot("1_Raw_Graph.dot", 0); 
    /*---Calculate distance between each two contigs of Every Edge---*/ 
    // SGCalcInsertSizeVisitor calv1(opt::insertSize, "1_Raw_mappingInfo.txt", !opt::skipMapping);
    // pGraph->visit(calv1);
    /*---Remove Edges by poisition---*/ 
    std::cout << "/*---Remove edges by poisition. \n";     
    SGRemoveByPosVisitor dv(opt::insertSize);
    pGraph->visit(dv); 
    pGraph->visit(statsVisit);
    pGraph->writeScaffoldDot("2_RemoveByPos_Graph.dot", 0);
       
    /*---Remove Circle Edges---*/
    std::cout << "/*---Removing Circle edges. \n";
    SGRemoveCircleEdgeVisitor rrev;
    pGraph->visit(rrev);
    pGraph->visit(statsVisit);    
    pGraph->writeScaffoldDot("3_RemoveCircleEdge_Graph.dot", 0); 
    
    /*---Remove Transitive Edges---*/ 
    std::cout << "/*---Removing transitive edges. \n";
    SGScaTransitiveReductionVisitor tr;
    pGraph->visit(tr);      
    pGraph->visit(statsVisit);  
    pGraph->writeScaffoldDot("4_RemoveTransitive_Graph.dot", 0); 
        
    /*---Remove Split Edges---*/     
    std::cout << "/*---Removing split edges. \n";
    SGRemoveSplitEdgeVisitor scarm;
    pGraph->visit(scarm);
    pGraph->visit(statsVisit);
    pGraph->writeScaffoldDot("5_RemoveSplitEdge_Graph.dot", 0); 
    

    /*---Remove Edges by Npair---*/  
    // std::cout << "/*---Remove small than Npair Threshold Edges. \n";
    // SGRemoveByNpairVisitor rbn(opt::npairThreshold);
    // pGraph->visit(rbn);
    // pGraph->visit(statsVisit);
    // pGraph->writeScaffoldDot("RemoveByNpair_Graph.dot", 0); 

    
    pGraph->simplifyScaffold(opt::insertSize);
    pGraph->visit(statsVisit);
    
    std::cout << "\n<Printing the scaffold-contig file> : " << opt::outFile << " \n" << std::endl;
	SGFastaVisitor2 av(opt::outFile);
    pGraph->visit(av);
    
    delete opt::pContigsBWT;
    delete opt::pContigsRBWT;
    delete opt::pContigsSSA;
    // delete opt::pPair_endBWT;
    // delete opt::pPair_endRBWT;
    // delete opt::pPair_endSSA;    
    return 0;
}
// Handle command line arguments
void parseScaffoldOptions(int argc, char** argv)
{
	optind=1;	//reset getopt
	// Set defaults
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c)
		{
            case '?': die = true; break;            
            case 'g': arg >> opt::asqgFile; opt::skipMapping = true; break; 
            case 'i': arg >> opt::insertSize; break;  	
            // case 'k': arg >> opt::kmerLength; break;  
            case 'n': arg >> opt::npairThreshold; break;     
            case 'o': arg >> opt::outFile; break;             
            case 'p': arg >> opt::prefixContigs; break;	
            // case 'q': arg >> opt::prefixPair_end; break;        
            case 't': arg >> opt::threads; break;            
            //-----------------------------------------
            case 'a': arg >> opt::algo; break; 
            case 'r': arg >> opt::repeatCutoff; break;
            // case 's': opt::skipMapping = true; break;
            case OPT_HELP:
                std::cout << SCAFFOLD_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);	
		}
	}
    if (argc - optind < 1)
	{
		std::cerr << SUBPROGRAM ": missing arguments\n";
		die = true;
	}
	else if (argc - optind > 1)
	{
		std::cerr << SUBPROGRAM ": too many arguments\n";
		die = true;
	}    
    if (die)
	{
		std::cout << "\n" << SCAFFOLD_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}
    if(opt::outFile.empty())
    {
        std::ostringstream ss;
        ss << opt::npairThreshold;
        opt::outFile = opt::prefixContigs + "-scaffold_n" + ss.str() + ".fa";
    }    
    // Parse the input filenames
    opt::readsFile = argv[optind++];    
}

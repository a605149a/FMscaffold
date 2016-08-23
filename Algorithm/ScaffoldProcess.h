#include "Util.h"
#include "SGUtil.h"
#include "BWTIndexSet.h"
#include "SequenceWorkItem.h"
#include "SGVisitors.h"

#include "tbb/concurrent_hash_map.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

/// typedef std::tr1::hash< int64_t> int64Hasher;
// typedef SparseHashMap<int64_t, std::pair<int64_t, size_t>, int64Hasher > KmerMatchMap;

typedef tbb::concurrent_hash_map<int64_t, std::pair<int64_t, size_t> > KmerMatchMap2;

// Parameters
struct ScaffoldParameters
{
    BWTIndexSet ContigsIndices;
    // BWTIndexSet Pair_endIndices;
    // size_t kmerLength;
    ReadInfoTable* pContigRIT;    
    
    // KmerMatchMap* kmerMatchMap;     
    KmerMatchMap2* kmerMatchMap2;   
    
    size_t repeatCutoff;
    size_t algo;
};
enum KmerComp
{
        KC_FORWARD = 0,
        KC_REVERSE = 1
};
class MAPPINGTABLE
{
    public:
        int64_t c1_idx;
        int64_t c2_idx;
        KmerComp k1_comp;
        KmerComp k2_comp;
        size_t npair;
        std::vector<size_t> c1_pos;
        std::vector<size_t> c2_pos;
        std::vector<SeqItem> r1_reads;
        std::vector<SeqItem> r2_reads;
};
class KMERMAPINFO
{
    public:
        int64_t c_idx;
        KmerComp k_comp;
        size_t pos;
};
// Results object
class ScaffoldResult
{
    public:
        std::vector<KMERMAPINFO> c1_indices;
        std::vector<KMERMAPINFO> c2_indices;
        std::string r1_id;
        std::string r1_seq;
        std::string r2_id;
        std::string r2_seq;
};
class ScaffoldProcess
{
    public:
        ScaffoldProcess(ScaffoldParameters params); 
        ~ScaffoldProcess();        
        ScaffoldResult process(const SequenceWorkItemPair& itemPair);
        
        KmerMatchMap2::accessor immAccessor;
        
    private:  
    
        const ScaffoldParameters m_params;
        void findStartSymbolIndex(int64_t idx);
        void pushUniIndex(std::vector<KMERMAPINFO> *c_indices, KMERMAPINFO temp);  
};
bool compare(const MAPPINGTABLE & s1, const MAPPINGTABLE & s2);
size_t CalcMedian(std::vector<size_t> poss);
class ScaffoldPostProcess
{
    public:
        ScaffoldPostProcess(ScaffoldParameters params, StringGraph* pGraph);
        ~ScaffoldPostProcess();
        void process(const SequenceWorkItemPair& itemPair, const ScaffoldResult& idx);
        
    private:    
        const ScaffoldParameters m_params;
        StringGraph* m_pGraph;
        size_t m_total_matepair;
        size_t m_align;
        size_t m_align_same;
        size_t m_align_diff;
        size_t m_unalign;  
        
        size_t m_npair;        
        size_t m_npair_ff;
        size_t m_npair_rr;
        size_t m_npair_fr;
        size_t m_npair_rf;
        //case 1,2: + + , - -
        std::vector<MAPPINGTABLE> same_map_table;
        //case 3: + -
        std::vector<MAPPINGTABLE> reverse_map_table1;
        //case 4: - +
        std::vector<MAPPINGTABLE> reverse_map_table2;
        void insert_same_table(int64_t c1_idx, KmerComp k1_comp,
                          int64_t c2_idx, KmerComp k2_comp,
                          size_t c1_pos, size_t c2_pos,
                          std::string r1_id, std::string r1_seq,
                          std::string r2_id, std::string r2_seq);
        void insert_reverse_table1(int64_t c1_idx, KmerComp k1_comp,
                          int64_t c2_idx, KmerComp k2_comp,
                          size_t c1_pos, size_t c2_pos,
                          std::string r1_id, std::string r1_seq,
                          std::string r2_id, std::string r2_seq);       
        void insert_reverse_table2(int64_t c1_idx, KmerComp k1_comp,
                          int64_t c2_idx, KmerComp k2_comp,
                          size_t c1_pos, size_t c2_pos,
                          std::string r1_id, std::string r1_seq,
                          std::string r2_id, std::string r2_seq); 
        void check_result();
        void connectSameEdge(int64_t c1_idx, int64_t c2_idx, EdgeComp comp, size_t npair,
                         size_t c1_pos_median, size_t c2_pos_median,                         
                         std::vector<SeqItem> r1_reads, std::vector<SeqItem> r2_reads);
        void connectReverseEdge1(int64_t c1_idx, int64_t c2_idx, EdgeComp comp, size_t npair,
                         size_t c1_pos_median, size_t c2_pos_median,                         
                         std::vector<SeqItem> r1_reads, std::vector<SeqItem> r2_reads);
        void connectReverseEdge2(int64_t c1_idx, int64_t c2_idx, EdgeComp comp, size_t npair,
                         size_t c1_pos_median, size_t c2_pos_median,                         
                         std::vector<SeqItem> r1_reads, std::vector<SeqItem> r2_reads);
};

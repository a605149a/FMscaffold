#include "ScaffoldProcess.h"
#include "BWTAlgorithms.h"
#include "HashMap.h"

using namespace std;
ScaffoldProcess::ScaffoldProcess(ScaffoldParameters params) : m_params(params)
{
    // (*m_params.kmerMatchMap).set_deleted_key(-1);
    // (*m_params.kmerMatchMap2).rehash();
}
//
ScaffoldProcess::~ScaffoldProcess()
{   
    // cout << "Sparse Hash Map Size: " << (*m_params.kmerMatchMap).size() << endl;
    // cout << "Hash Map Size: " << (*m_params.kmerMatchMap2).size() << endl;
}
//
void ScaffoldProcess::pushUniIndex(std::vector<KMERMAPINFO> *c_indices, KMERMAPINFO temp)
{
    bool is_exist = false;
    // std::cout <<     temp.c_idx << ", " << temp.pos << std::endl;
    for(size_t i = 0; i < c_indices->size(); i++)
    {
        if((*c_indices)[i].c_idx == temp.c_idx && (*c_indices)[i].k_comp == temp.k_comp )//&& (*c_indices)[i].pos == temp.pos)
        {
            is_exist = true;
            break;
        }
    }
    if(is_exist == false)
    {
        c_indices->push_back(temp);
    }
}
//
ScaffoldResult ScaffoldProcess::process(const SequenceWorkItemPair& workItemPair)
{
    assert(m_params.ContigsIndices.pBWT != NULL);
	assert(m_params.ContigsIndices.pSSA != NULL);    
    
    ScaffoldResult result;    
    
    string read1 = workItemPair.first.read.seq.toString();
    size_t read1_len = read1.length();
    string read2_t = workItemPair.second.read.seq.toString();
    string read2 = reverseComplement(read2_t);
    size_t read2_len = read1.length();
        for(size_t i = 0; i < read1.length()/3; i++)
        {
            string r1_kmer = read1.substr(i, read1_len);        
            string r1_kmer_rev = reverseComplement(r1_kmer);
            BWTInterval intervalR1 = BWTAlgorithms::findInterval(m_params.ContigsIndices, r1_kmer);
            BWTInterval intervalR1_rev = BWTAlgorithms::findInterval(m_params.ContigsIndices, r1_kmer_rev);
            if(intervalR1.size() + intervalR1_rev.size() == 0)
            {
                    continue;
            }
            if(intervalR1.size() + intervalR1_rev.size() == 1)
            {
                for(int64_t j = intervalR1.lower; j <= intervalR1.upper ; ++j)
                {                 
                    KMERMAPINFO temp;
                    temp.k_comp = KC_FORWARD;
                
                
                    if((*m_params.kmerMatchMap2).find(immAccessor,j))	//already has
                    {
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    else	//first time
                    {              
                        immAccessor.release();
                        findStartSymbolIndex(j);
                        (*m_params.kmerMatchMap2).find(immAccessor,j);
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    pushUniIndex(&result.c1_indices, temp);
                }      
                for(int64_t j = intervalR1_rev.lower; j <= intervalR1_rev.upper ; ++j)
                {
                    KMERMAPINFO temp;
                    temp.k_comp = KC_REVERSE;
                    if((*m_params.kmerMatchMap2).find(immAccessor,j))	//already has
                    {
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    else	//first time
                    {                  
                        immAccessor.release();
                        findStartSymbolIndex(j);
                        (*m_params.kmerMatchMap2).find(immAccessor,j);
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    pushUniIndex(&result.c1_indices, temp);
                }
                break;
            }
            else
                break;
        }
        for(size_t i = 0; i < read2.length()/3; i++)
        {
            string r2_kmer = read2.substr(i, read2_len);        
            string r2_kmer_rev = reverseComplement(r2_kmer);
            BWTInterval intervalR2 = BWTAlgorithms::findInterval(m_params.ContigsIndices, r2_kmer);
            BWTInterval intervalR2_rev = BWTAlgorithms::findInterval(m_params.ContigsIndices, r2_kmer_rev);
            if(intervalR2.size() + intervalR2_rev.size() == 0)
            {
                continue;
            }
            if(intervalR2.size() + intervalR2_rev.size() == 1)
            {
                for(int64_t j = intervalR2.lower; j <= intervalR2.upper ; ++j)
                {                 
                    KMERMAPINFO temp;
                    temp.k_comp = KC_FORWARD;
                    if((*m_params.kmerMatchMap2).find(immAccessor,j))	//already has
                    {
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    else	//first time
                    {                  
                        immAccessor.release();
                        findStartSymbolIndex(j);
                        (*m_params.kmerMatchMap2).find(immAccessor,j);
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    pushUniIndex(&result.c2_indices, temp);
                }     
                for(int64_t j = intervalR2_rev.lower; j <= intervalR2_rev.upper ; ++j)
                {
                    KMERMAPINFO temp;
                    temp.k_comp = KC_REVERSE;
                    if((*m_params.kmerMatchMap2).find(immAccessor,j))	//already has
                    {
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    else	//first time
                    {                  
                        immAccessor.release();
                        findStartSymbolIndex(j);
                        (*m_params.kmerMatchMap2).find(immAccessor,j);
                        temp.c_idx = immAccessor->second.first;
                        temp.pos = immAccessor->second.second;
                        immAccessor.release();   
                    }
                    pushUniIndex(&result.c2_indices, temp);
                }
            }
            else
                break;
        }        
    result.r1_id = workItemPair.first.read.id;
    result.r1_seq = workItemPair.first.read.seq.toString();
    result.r2_id = workItemPair.second.read.id;
    result.r2_seq = workItemPair.second.read.seq.toString();
    return result;    
}
//loop ver, Backtrack the idx until we hit the starting symbol, building hashtable 
void ScaffoldProcess::findStartSymbolIndex(int64_t idx)
{   
    char b;
    int64_t new_idx;
    int64_t final_idx;
    std::vector<int64_t> idxTable;
    size_t distance = 0;
    while(1)
    {
        
        
        if((*m_params.kmerMatchMap2).find(immAccessor,idx))
        {            
            final_idx = immAccessor->second.first;
            distance += immAccessor->second.second;
            idxTable.push_back(idx);
            immAccessor.release();
            break;        
        }
        else
        {
            immAccessor.release();
            b = m_params.ContigsIndices.pBWT->getChar(idx);
            new_idx = m_params.ContigsIndices.pBWT->getPC(b) + m_params.ContigsIndices.pBWT->getOcc(b, idx - 1);        
            distance += 1;
            if(b == '$')
            {
                final_idx = m_params.ContigsIndices.pSSA->lookupLexoRank(new_idx);
                idxTable.push_back(idx);
                break;
            }
            else
            {
                idxTable.push_back(idx);
            }           
            idx = new_idx;
        }
    }
    // std::cout << "idxTable.size = " << idxTable.size() << std::endl;
    for(size_t i = 0; i < idxTable.size(); i++)
    {
        new_idx = idxTable[i];      
        (*m_params.kmerMatchMap2).insert(immAccessor,new_idx); 
        immAccessor->second.first = final_idx;
        immAccessor->second.second = distance - i;
        immAccessor.release(); 
    } 
    // if((*m_params.kmerMatchMap2).size() > 90000000)
    // (*m_params.kmerMatchMap2).rehash();
    
    idxTable.clear();        
}
ScaffoldPostProcess::ScaffoldPostProcess(ScaffoldParameters params, StringGraph* pGraph) :
    m_params(params), m_pGraph(pGraph), m_total_matepair(0), m_align_same(0), m_align_diff(0), m_unalign(0), m_npair(0)
{ 
    m_align = 0;
    m_npair_ff = 0;
    m_npair_rr = 0;
    m_npair_fr = 0;
    m_npair_rf = 0;
}
//
ScaffoldPostProcess::~ScaffoldPostProcess()
{    
    cout << "#Total_mate_pair: " << m_total_matepair << endl;
    cout << "#m_align: " << m_align;    
    cout << " Ratio: " << (double)m_align/m_total_matepair*100 << "%" << endl; 
    cout << "#m_align_same: " << m_align_same;    
    cout << " Ratio: " << (double)m_align_same/m_total_matepair*100 << "%" << endl; 
    cout << "#m_align_diff: " << m_align_diff;    
    cout << " Ratio: " << (double)m_align_diff/m_total_matepair*100 << "%" << endl; 
    cout << "#unalign: " << m_unalign;    
    cout << " Ratio: " << (double)m_unalign/m_total_matepair*100 << "%" << endl;
    cout << endl;
    cout << "#Total_npair: " << m_align_diff << endl;
    cout << "#Total_npair_ff: " << m_npair_ff;
    cout << " Ratio: " << (double)m_npair_ff/m_align_diff*100 << "%" << endl;
    cout << "#Total_npair_rr: " << m_npair_rr;
    cout << " Ratio: " << (double)m_npair_rr/m_align_diff*100 << "%" << endl;
    cout << "#Total_npair_fr: " << m_npair_fr;
    cout << " Ratio: " << (double)m_npair_fr/m_align_diff*100 << "%" << endl;
    cout << "#Total_npair_rf: " << m_npair_rf;
    cout << " Ratio: " << (double)m_npair_rf/m_align_diff*100 << "%" << endl;
    
    check_result();
}
//
void ScaffoldPostProcess::process(const SequenceWorkItemPair& /*itemPair*/, const ScaffoldResult& result)
{
    m_total_matepair += 1;
    if(result.c1_indices.size() > 0 && result.c2_indices.size() > 0)
        m_align += 1;
    else
        m_unalign += 1;
    
    // std::cout <<     result.c1_indices.size() << ", " << result.c2_indices.size() << std::endl;
    for(size_t i = 0; i < result.c1_indices.size(); i++)
    for(size_t j = 0; j < result.c2_indices.size(); j++)
    {
        
        if(result.c1_indices[i].c_idx != result.c2_indices[j].c_idx) 
        {
            
            m_align_diff += 1;                
            if(result.c1_indices[i].k_comp == KC_FORWARD 
            && result.c2_indices[j].k_comp == KC_FORWARD)
            {            
                insert_same_table(result.c1_indices[i].c_idx, result.c1_indices[i].k_comp,
                                result.c2_indices[j].c_idx, result.c2_indices[j].k_comp,
                                result.c1_indices[i].pos, result.c2_indices[j].pos,
                                result.r1_id, result.r1_seq,
                                result.r2_id, result.r2_seq);                    
                m_npair_ff += 1;
            }
            if(result.c1_indices[i].k_comp == KC_REVERSE 
            && result.c2_indices[j].k_comp == KC_REVERSE)
            {               
                insert_same_table(result.c2_indices[j].c_idx, result.c2_indices[j].k_comp,
                                    result.c1_indices[i].c_idx, result.c1_indices[i].k_comp,
                                    result.c2_indices[j].pos, result.c1_indices[i].pos,
                                    result.r2_id, result.r2_seq,
                                    result.r1_id, result.r1_seq);
                m_npair_rr += 1;
            }            
            if(result.c1_indices[i].k_comp == KC_FORWARD 
            && result.c2_indices[j].k_comp == KC_REVERSE)
            {
                insert_reverse_table1(result.c1_indices[i].c_idx, result.c1_indices[i].k_comp,
                                result.c2_indices[j].c_idx, result.c2_indices[j].k_comp,
                                result.c1_indices[i].pos, result.c2_indices[j].pos,
                                result.r1_id, result.r1_seq,
                                result.r2_id, result.r2_seq); 
                m_npair_fr += 1;            
            }
            if(result.c1_indices[i].k_comp == KC_REVERSE 
            && result.c2_indices[j].k_comp == KC_FORWARD)
            {
                insert_reverse_table2(result.c1_indices[i].c_idx, result.c1_indices[i].k_comp,
                                result.c2_indices[j].c_idx, result.c2_indices[j].k_comp,
                                result.c1_indices[i].pos, result.c2_indices[j].pos,
                                result.r1_id, result.r1_seq,
                                result.r2_id, result.r2_seq);
                m_npair_rf += 1;
            }                  
        }            
        else
            m_align_same += 1;
    
    }    
}
//
void ScaffoldPostProcess::insert_same_table(int64_t c1_idx, KmerComp k1_comp,
                                       int64_t c2_idx, KmerComp k2_comp,
                                       size_t c1_pos, size_t c2_pos,
                                       std::string r1_id, std::string r1_seq,
                                       std::string r2_id, std::string r2_seq)
{ 
    SeqItem r1_read;
    SeqItem r2_read;
    r1_read.id = r1_id;
    r1_read.seq = r1_seq;
    r2_read.id = r2_id;
    r2_read.seq = r2_seq;    
    bool is_exist = false;
    for(size_t i = 0; i < same_map_table.size(); i++)
    {
        if(same_map_table[i].c1_idx == c1_idx && same_map_table[i].c2_idx == c2_idx)
        // &&  same_map_table[i].k1_comp == k1_comp && same_map_table[i].k2_comp == k2_comp)
        {
            is_exist = true;
            same_map_table[i].npair += 1;
            same_map_table[i].c1_pos.push_back(c1_pos);
            same_map_table[i].c2_pos.push_back(c2_pos);
            same_map_table[i].r1_reads.push_back(r1_read);
            same_map_table[i].r2_reads.push_back(r2_read);            
            break;
        }
    }
    // for(size_t i = 0; i < same_map_table.size(); i++)
    // {
        // if(same_map_table[i].c1_idx == c2_idx && same_map_table[i].c2_idx == c1_idx)
        // {
            // is_exist = true;
            // same_map_table[i].npair += 1;
            // same_map_table[i].c1_pos.push_back(c1_pos);
            // same_map_table[i].c2_pos.push_back(c2_pos);
            // same_map_table[i].r1_reads.push_back(r1_read);
            // same_map_table[i].r2_reads.push_back(r2_read);            
            // break;
        // }
    // }
    if(is_exist == false)
    {        
        MAPPINGTABLE temp;
        temp.c1_idx = c1_idx;
        temp.c2_idx = c2_idx;
        temp.k1_comp = k1_comp;
        temp.k2_comp = k2_comp;
        temp.c1_pos.push_back(c1_pos);
        temp.c2_pos.push_back(c2_pos);
        temp.npair = 1;
        temp.r1_reads.push_back(r1_read);
        temp.r2_reads.push_back(r2_read);
        same_map_table.push_back(temp);
    }
}
void ScaffoldPostProcess::insert_reverse_table1(int64_t c1_idx, KmerComp k1_comp,
                                       int64_t c2_idx, KmerComp k2_comp,
                                       size_t c1_pos, size_t c2_pos,
                                       std::string r1_id, std::string r1_seq,
                                       std::string r2_id, std::string r2_seq)
{ 
    SeqItem r1_read;
    SeqItem r2_read;
    r1_read.id = r1_id;
    r1_read.seq = r1_seq;
    r2_read.id = r2_id;
    r2_read.seq = r2_seq;    
    bool is_exist = false;
    for(size_t i = 0; i < reverse_map_table1.size(); i++)
    {
        if(reverse_map_table1[i].c1_idx == c1_idx && reverse_map_table1[i].c2_idx == c2_idx)
        {
            is_exist = true;
            reverse_map_table1[i].npair += 1;
            reverse_map_table1[i].c1_pos.push_back(c1_pos);
            reverse_map_table1[i].c2_pos.push_back(c2_pos);
            reverse_map_table1[i].r1_reads.push_back(r1_read);
            reverse_map_table1[i].r2_reads.push_back(r2_read);            
            break;
        }
    }
    for(size_t i = 0; i < reverse_map_table1.size(); i++)
    {
        if(reverse_map_table1[i].c1_idx == c2_idx && reverse_map_table1[i].c2_idx == c1_idx)
        {
            is_exist = true;
            reverse_map_table1[i].npair += 1;
            reverse_map_table1[i].c1_pos.push_back(c2_pos);
            reverse_map_table1[i].c2_pos.push_back(c1_pos);
            reverse_map_table1[i].r1_reads.push_back(r2_read);
            reverse_map_table1[i].r2_reads.push_back(r1_read);            
            break;
        }
    }
    if(is_exist == false)
    {        
        MAPPINGTABLE temp;
        temp.c1_idx = c1_idx;
        temp.c2_idx = c2_idx;
        temp.k1_comp = k1_comp;
        temp.k2_comp = k2_comp;
        temp.c1_pos.push_back(c1_pos);
        temp.c2_pos.push_back(c2_pos);
        temp.npair = 1;
        temp.r1_reads.push_back(r1_read);
        temp.r2_reads.push_back(r2_read);
        reverse_map_table1.push_back(temp);
    }
}
void ScaffoldPostProcess::insert_reverse_table2(int64_t c1_idx, KmerComp k1_comp,
                                       int64_t c2_idx, KmerComp k2_comp,
                                       size_t c1_pos, size_t c2_pos,
                                       std::string r1_id, std::string r1_seq,
                                       std::string r2_id, std::string r2_seq)
{ 
    SeqItem r1_read;
    SeqItem r2_read;
    r1_read.id = r1_id;
    r1_read.seq = r1_seq;
    r2_read.id = r2_id;
    r2_read.seq = r2_seq;    
    bool is_exist = false;
    for(size_t i = 0; i < reverse_map_table2.size(); i++)
    {
        if(reverse_map_table2[i].c1_idx == c1_idx && reverse_map_table2[i].c2_idx == c2_idx)
        {
            is_exist = true;
            reverse_map_table2[i].npair += 1;
            reverse_map_table2[i].c1_pos.push_back(c1_pos);
            reverse_map_table2[i].c2_pos.push_back(c2_pos);
            reverse_map_table2[i].r1_reads.push_back(r1_read);
            reverse_map_table2[i].r2_reads.push_back(r2_read);            
            break;
        }
    }
    for(size_t i = 0; i < reverse_map_table2.size(); i++)
    {
        if(reverse_map_table2[i].c1_idx == c2_idx && reverse_map_table2[i].c2_idx == c1_idx)
        {
            is_exist = true;
            reverse_map_table2[i].npair += 1;
            reverse_map_table2[i].c1_pos.push_back(c2_pos);
            reverse_map_table2[i].c2_pos.push_back(c1_pos);
            reverse_map_table2[i].r1_reads.push_back(r2_read);
            reverse_map_table2[i].r2_reads.push_back(r1_read);            
            break;
        }
    }
    if(is_exist == false)
    {        
        MAPPINGTABLE temp;
        temp.c1_idx = c1_idx;
        temp.c2_idx = c2_idx;
        temp.k1_comp = k1_comp;
        temp.k2_comp = k2_comp;
        temp.c1_pos.push_back(c1_pos);
        temp.c2_pos.push_back(c2_pos);
        temp.npair = 1;
        temp.r1_reads.push_back(r1_read);
        temp.r2_reads.push_back(r2_read);
        reverse_map_table2.push_back(temp);
    }
}
size_t CalcMedian(vector<size_t> poss)
{
  size_t median;
  size_t size = poss.size();
  sort(poss.begin(), poss.end());
  if (size  % 2 == 0)
  {
      median = (poss[size / 2 - 1] + poss[size / 2]) / 2;
  }
  else 
  {
      median = poss[size / 2];
  }
  return median;
}
bool compare(const MAPPINGTABLE & s1, const MAPPINGTABLE & s2)
{
    if (s1.c1_idx != s2.c1_idx) 
        return s1.c1_idx < s2.c1_idx;
    return s1.npair > s2.npair;
}
void ScaffoldPostProcess::check_result()
{      
    for(size_t i = 0; i < same_map_table.size(); i++)
    {        
            size_t c1_pos_median = CalcMedian(same_map_table[i].c1_pos);
            size_t c2_pos_median = CalcMedian(same_map_table[i].c2_pos);
            connectSameEdge(same_map_table[i].c1_idx, same_map_table[i].c2_idx, EC_SAME, same_map_table[i].npair,
                        c1_pos_median, c2_pos_median,                        
                        same_map_table[i].r1_reads, same_map_table[i].r2_reads);       
    }    
    for(size_t i = 0; i < reverse_map_table1.size(); i++)
    {        
            size_t c1_pos_median = CalcMedian(reverse_map_table1[i].c1_pos);
            size_t c2_pos_median = CalcMedian(reverse_map_table1[i].c2_pos);
            connectReverseEdge1(reverse_map_table1[i].c1_idx, reverse_map_table1[i].c2_idx, EC_REVERSE, reverse_map_table1[i].npair,
                        c1_pos_median, c2_pos_median,                        
                        reverse_map_table1[i].r1_reads, reverse_map_table1[i].r2_reads);           
    }  
    for(size_t i = 0; i < reverse_map_table2.size(); i++)
    {        
            size_t c1_pos_median = CalcMedian(reverse_map_table2[i].c1_pos);
            size_t c2_pos_median = CalcMedian(reverse_map_table2[i].c2_pos);
            connectReverseEdge2(reverse_map_table2[i].c1_idx, reverse_map_table2[i].c2_idx, EC_REVERSE, reverse_map_table2[i].npair,
                        c1_pos_median, c2_pos_median,                        
                        reverse_map_table2[i].r1_reads, reverse_map_table2[i].r2_reads);           
    }  
}
//
void ScaffoldPostProcess::connectSameEdge(int64_t c1_idx, int64_t c2_idx, EdgeComp ec_comp, size_t npair,
                                      size_t c1_pos_median, size_t c2_pos_median,  
                                      std::vector<SeqItem> r1_reads, std::vector<SeqItem> r2_reads)
{    
    SeqItem scaffold_contig;
    string c1ContigID = m_params.pContigRIT->getReadInfo(c1_idx).id;    
    string c2ContigID = m_params.pContigRIT->getReadInfo(c2_idx).id;
    Vertex* pVL;
    pVL = m_pGraph->getVertex(c1ContigID);
    Vertex* pVR;
    pVR = m_pGraph->getVertex(c2ContigID);
    
	Edge* edgeLR=new Edge(pVR,ED_SENSE, ec_comp, npair, c1_pos_median, c2_pos_median, r1_reads, r2_reads);
	Edge* edgeRL=new Edge(pVL,ED_ANTISENSE, ec_comp, npair, c2_pos_median, c1_pos_median, r1_reads, r2_reads);
    edgeLR->setTwin(edgeRL);
	edgeRL->setTwin(edgeLR);
	m_pGraph->addEdge(pVL,edgeLR);
    m_pGraph->addEdge(pVR,edgeRL); 
}
void ScaffoldPostProcess::connectReverseEdge1(int64_t c1_idx, int64_t c2_idx, EdgeComp ec_comp, size_t npair,
                                      size_t c1_pos_median, size_t c2_pos_median,  
                                      std::vector<SeqItem> r1_reads, std::vector<SeqItem> r2_reads)
{    
    SeqItem scaffold_contig;
    string c1ContigID = m_params.pContigRIT->getReadInfo(c1_idx).id;    
    string c2ContigID = m_params.pContigRIT->getReadInfo(c2_idx).id;
    Vertex* pVL;
    pVL = m_pGraph->getVertex(c1ContigID);
    Vertex* pVR;
    pVR = m_pGraph->getVertex(c2ContigID);    
    
	Edge* edgeLR=new Edge(pVR,ED_SENSE, ec_comp, npair, c1_pos_median, c2_pos_median, r1_reads, r2_reads);
	Edge* edgeRL=new Edge(pVL,ED_SENSE, ec_comp, npair, c2_pos_median, c1_pos_median, r1_reads, r2_reads);
    edgeLR->setTwin(edgeRL);
	edgeRL->setTwin(edgeLR);
	m_pGraph->addEdge(pVL,edgeLR);
    m_pGraph->addEdge(pVR,edgeRL); 
}
void ScaffoldPostProcess::connectReverseEdge2(int64_t c1_idx, int64_t c2_idx, EdgeComp ec_comp, size_t npair,
                                      size_t c1_pos_median, size_t c2_pos_median,  
                                      std::vector<SeqItem> r1_reads, std::vector<SeqItem> r2_reads)
{    
    SeqItem scaffold_contig;
    string c1ContigID = m_params.pContigRIT->getReadInfo(c1_idx).id;    
    string c2ContigID = m_params.pContigRIT->getReadInfo(c2_idx).id;
    Vertex* pVL;
    pVL = m_pGraph->getVertex(c1ContigID);
    Vertex* pVR;
    pVR = m_pGraph->getVertex(c2ContigID);    
    
	Edge* edgeLR=new Edge(pVR,ED_ANTISENSE, ec_comp, npair, c1_pos_median, c2_pos_median, r1_reads, r2_reads);
	Edge* edgeRL=new Edge(pVL,ED_ANTISENSE, ec_comp, npair, c2_pos_median, c1_pos_median, r1_reads, r2_reads);
    edgeLR->setTwin(edgeRL);
	edgeRL->setTwin(edgeLR);
	m_pGraph->addEdge(pVL,edgeLR);
    m_pGraph->addEdge(pVR,edgeRL); 
}
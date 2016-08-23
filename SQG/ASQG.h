//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ASQG - Definitions and functions
// for handling ASGQ files
//
#ifndef ASQG_H
#define ASQG_H

#include "SQG.h"
#include "Match.h"

namespace ASQG
{
    enum RecordType
    {
        RT_HEADER = 0,
        RT_VERTEX,
        RT_EDGE
    };

    // A header record is just a tag:value pairs
    struct HeaderRecord
    {
        public:
            HeaderRecord();
            HeaderRecord(const std::string& recordLine);
            
            void setOverlapTag(int overlapLen);
            void setInputFileTag(const std::string& name);
            void setErrorRateTag(float errorRate);
            void setContainmentTag(int v);
            void setTransitiveTag(int v);

            const SQG::IntTag& getVersionTag() const { return m_versionTag; }
            const SQG::FloatTag& getErrorRateTag() const { return m_errorRateTag; }
            const SQG::StringTag& getInfileTag() const { return m_infileTag; }
            const SQG::IntTag& getOverlapTag() const { return m_overlapTag; }
            const SQG::IntTag& getContainmentTag() const { return m_containmentTag; };
            const SQG::IntTag& getTransitiveTag() const { return m_transitiveTag; };

            void write(std::ostream& out);
            void parse(const std::string& record);

        private:

            void setVersionTag(int version);
            
            SQG::IntTag m_versionTag;
            SQG::FloatTag m_errorRateTag;
            SQG::StringTag m_infileTag;
            SQG::IntTag m_overlapTag;
            SQG::IntTag m_containmentTag;
            SQG::IntTag m_transitiveTag;
    };

    // A vertex record is an id, sequence and an array of
    // tag:value 
    struct VertexRecord
    {
        public:
            VertexRecord() {}
            VertexRecord(const std::string& recordLine);
            VertexRecord(const std::string& i, const std::string& s) : m_id(i), m_seq(s) {}

            void setSubstringTag(bool b);
            
            const std::string& getID() const { return m_id; }
            const std::string& getSeq() const { return m_seq; }
            const SQG::IntTag& getSubstringTag() const { return m_substringTag; }

            void write(std::ostream& out);
            void parse(const std::string& record);

        private:

            std::string m_id;
            std::string m_seq;
            SQG::IntTag m_substringTag;
    };

    // An edge record is just an overlap object and tag:values
    struct EdgeRecord
    {
        public:
            EdgeRecord() {}
            EdgeRecord(const std::string& recordLine);
            EdgeRecord(const Overlap& o) : m_overlap(o) {}

            const Overlap& getOverlap() const { return m_overlap; }

            void write(std::ostream& out);
            void parse(const std::string& record);

        private:
        
            Overlap m_overlap;
    };

    struct ScaEdgeRecord
    {
        public:
            ScaEdgeRecord() {}
            ScaEdgeRecord(const std::string& recordLine);
            ScaEdgeRecord(const std::string StartID, const std::string EndID, const std::string EC_comp, const std::string ED_dir, const size_t npair, const size_t C1Pos, const size_t C2Pos) : 
                            m_StartID(StartID), m_EndID(EndID), m_EC_comp(EC_comp), m_ED_dir(ED_dir), m_npair(npair), m_C1Pos(C1Pos), m_C2Pos(C2Pos) {}

            const std::string& getStartID() const { return m_StartID; }
            const std::string& getEndID() const { return m_EndID; }
            const std::string& getECcomp() const { return m_EC_comp; }
            const std::string& getEDdir() const { return m_ED_dir; }
            const size_t& getNpair() const { return m_npair; }
            const size_t& getC1Pos() const { return m_C1Pos; }
            const size_t& getC2Pos() const { return m_C2Pos; }
            
            void write(std::ostream& out);
            void parse(const std::string& record);

        private:
            std::string m_StartID;
            std::string m_EndID;
            std::string m_EC_comp;
            std::string m_ED_dir;
            size_t m_npair;
            size_t m_C1Pos;
            size_t m_C2Pos;
    };
    
    // Parsing functions
    RecordType getRecordType(const std::string& record);
    HeaderRecord parseHeaderRecord(const std::string& record);
    VertexRecord parseVertexRecord(const std::string& record);
    EdgeRecord parseEdgeRecord(const std::string& record);
    ScaEdgeRecord parseScaEdgeRecord(const std::string& record);
    
    // Writing functions
    void writeFields(std::ostream& out, const StringVector& fields);
};

#endif

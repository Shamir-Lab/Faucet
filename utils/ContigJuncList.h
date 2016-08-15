#ifndef CONTIG_JUNC_LIST
#define CONTIG_JUNC_LIST

class JuncResult;

#include <deque>
#include "JuncPairs.h"

//Stores info about all interior junctions
//List of incremental distances and coverages
//Coverages include coverage at each end
class ContigJuncList{

public:

    typedef std::deque<unsigned char> junc_list;
    typedef junc_list::const_iterator const_iterator;


    ContigJuncList(std::string seq, junc_list dist, junc_list cov);
    ContigJuncList();

    const_iterator begin_distances() const{ return distances.begin();}
    const_iterator begin_coverages() const{ return coverages.begin();}
    const_iterator end_distances() const{ return distances.end();}
    const_iterator end_coverages() const{ return coverages.end();}
    void setSeq(std::string cont){seq = cont;}
    std::string getSeq(){ return seq;}
    bool isValidKmerPosition(int pos);//true iff getKmer makes sense on this
    kmer_type getKmer(int pos);//0 = first backward, 1 = first forward, 2 = second backward, etc.

    //Gets a list of JuncResults, specifying distance, coverage, and kmer for juncs, with reference to specified side 
    std::list<JuncResult> getJuncResults(bool startForward, int startDist, int maxDist);
    void printJuncResults(int side, int startDist, int maxDist);
    void printJuncResults(std::list<JuncResult> results);
    void printJuncValues();

    //Used for reversing a contig.  Simply reverses both lists
    void reverse();

    int length(); //returns length of sequence
    
    //Sums all distance values
    int getTotalDistance();

    //Concatenates this list of juncs with another 
    //Removes overlap of middle coverage and middle distance
    ContigJuncList concatenate(ContigJuncList otherList);

    //Averages all coverage values in list
    double getAvgCoverage();
    double getAvgCoverage(std::list<JuncResult> results);

    double getCoverageSampleVariance();
    double getCoverageSampleVariance(std::list<JuncResult> results);
    ContigJuncList getScaledContigJuncs(double scale_factor);


    //Prints distances then coverages to a string
    std::string getStringRep();

private: 
    junc_list distances;
    junc_list coverages;
    std::string seq; //string sequence, represented forward from side 1 to side 2
};

#endif

#include <sstream>
#include <deque>
#include <algorithm>   
#include <fstream>

#include "ContigJuncList.h"

using std::stringstream;

ContigJuncList::ContigJuncList(std::string sequence, junc_list dist, junc_list cov){
  seq = sequence;
  distances = dist;
  coverages = cov;   
} 

ContigJuncList::ContigJuncList(){
  distances.clear();
  coverages.clear();   
  seq = "";
}

int ContigJuncList::size(){
    return coverages.size();
}

void ContigJuncList::printJuncValues(){
    for(auto itC = coverages.begin(), itD = distances.begin(); itC != coverages.end(); itC++, itD++){
        std::cout << (unsigned int)*itD << ", " << (unsigned int)*itC <<" ; "; 
    }
    std::cout << "\n";
}

void ContigJuncList::printJuncResults(std::list<JuncResult> results){
    for(auto it = results.begin(); it != results.end(); ++it){
        std::cout << print_kmer(it->kmer) << " " << it->distance << " " 
            << it->coverage << " , ";
    }
    std::cout << "\n";
}

void ContigJuncList::printJuncResults(int side, int startDist, int maxDist){
    printJuncResults(getJuncResults(side, startDist, maxDist));
}

int ContigJuncList::length(){
    return seq.length();
}

bool ContigJuncList::isValidKmerPosition(int pos){
    return pos >= 0 && pos <= (getSeq().length() - sizeKmer + 1)*2;
}

//0 is the first backward kmer, 1 is the first forward kmer, 
//2 is the second backward kmer, etc.
kmer_type ContigJuncList::getKmer(int pos){
    if(pos % 2 == 1){
        return getKmerFromRead(getSeq(), pos/2);
    }
    else{
        return revcomp(getKmerFromRead(getSeq(), pos/2));
    }
}

//Gets a list of JuncResults, specifying distance, coverage, and kmer for juncs
//If startForward, assumes first junc faces forward (in towards contig)
//If !startForward, assumes first junc faces backward (away from contig)
//Always assumes no reverse needed- contig can reverse it if necessary before calling 
std::list<JuncResult> ContigJuncList::getJuncResults(bool startForward, int startDist, int maxDist){
    std::list<JuncResult> results = {};
    int startPos = 0;
    if(startForward){
        results.push_back(JuncResult(getKmer(3), 2+startDist,coverages.front()));
        startPos = 1;
    }
    int pos = startPos;
    for(auto itD = distances.begin(), itC = ++coverages.begin();
        itD != distances.end(); 
        itD++, itC++){
        pos += *itD;
        //real extension is 2 in front if the junc faces forward, or 2 behind if it faces backward
        int offset = 2;
        if(pos % 2 == 0){
            offset = -2;
        }
        if(pos + offset - startPos + startDist <= maxDist){
            if(isValidKmerPosition(pos+offset)){
                results.push_back(JuncResult(getKmer(pos + offset), pos + offset-startPos + startDist, *itC));
            }
        }
        else { break; }
        
    }
    return results;
}

//Used for reversing a contig.  Simply reverses both lists
void ContigJuncList::reverse(){
    std::reverse(coverages.begin(), coverages.end());
    std::reverse(distances.begin(), distances.end());
    seq = revcomp_string(seq);
}

//Concatenates this list of juncs with another 
//Removes overlap of middle coverage and middle distance
ContigJuncList ContigJuncList::concatenate(ContigJuncList otherList){

    junc_list newCov(coverages);
    newCov.pop_back();

    unsigned char lastCov = coverages.back();
    unsigned char firstCov = otherList.coverages.front();

    newCov.push_back((unsigned char) std::min((int)lastCov, (int)firstCov));

    newCov.insert(newCov.end(), ++otherList.coverages.begin(), otherList.coverages.end());

    junc_list newDist(distances);
    newDist.insert(newDist.end(), otherList.distances.begin(), otherList.distances.end());
    std::string newSeq  = getSeq().substr(0, getSeq().length()-sizeKmer) + otherList.getSeq();
    return ContigJuncList(newSeq, newDist, newCov);
}

ContigJuncList ContigJuncList::getShiftedCoverageContigJuncs(double shift){
    junc_list newCov(coverages);
    for(auto it = coverages.begin(); it != coverages.end(); ++it){
        double val = (double) *it+shift;   
        newCov.push_back((int) std::round((val > 255) ? 255: val));
        newCov.pop_front();
    }
    return ContigJuncList(seq,distances,newCov);
}


ContigJuncList ContigJuncList::getScaledContigJuncs(double scale_factor){
    junc_list newCov(coverages);
    for(auto it = coverages.begin(); it != coverages.end(); ++it){
        newCov.push_back((int) std::round(*it * scale_factor));
        newCov.pop_front();
    }
    return ContigJuncList(seq,distances,newCov);
}

//Averages all coverage values in list
double ContigJuncList::getAvgCoverage(){
    double covSum = 0;
    if(coverages.size()== 0){
        printf("ERROR: empty junctions list\n");
        return 0;
    }
    for(auto it = coverages.begin(); it != coverages.end(); ++it){
        covSum += (double) *it;
    }
    return covSum / coverages.size();
}

double ContigJuncList::getAvgCoverage(std::list<JuncResult> results){
    if (results.size()==0){return 0;}
    double covSum = 0;
    for(auto it = results.begin(); it != results.end(); ++it){
        covSum += (double) it->coverage; 
    }
    return covSum / results.size();
}

double ContigJuncList::getCoverageSampleVariance(){
    if(coverages.size() < 2){
        printf("ERROR: 1 or 0 values in junctions list\n");
        return 0;
    }
    double mean = getAvgCoverage();
    double sum_sqrs = 0;
    for(auto it = coverages.begin(); it != coverages.end(); ++it){
        sum_sqrs += pow(mean - (double) *it, 2); 
    }
    return pow(sum_sqrs / (coverages.size()-1), 0.5);
}

double ContigJuncList::getCoverageSampleVariance(std::list<JuncResult> results){
    if (results.size() < 2) {return 0;}
    double mean = getAvgCoverage(results);
    double sum_sqrs = 0;
    for(auto it = results.begin(); it != results.end(); ++it){
        sum_sqrs += pow(mean - (double) it->coverage, 2); 
    }
    return pow(sum_sqrs / (results.size()-1), 0.5);
}


//Sums all distance values
int ContigJuncList::getTotalDistance(){
    int totalDist = 0;
    for(auto it = distances.begin(); it != distances.end(); ++it){
        totalDist += (int) *it;
    }
    return totalDist;
}

//Prints distances then coverages to a string
std::string ContigJuncList::getStringRep(){
    stringstream stream;
    stream << getSeq() << "\n";
    stream << "Distances: ";
    for(auto it = distances.begin(); it != distances.end(); ++it){
        stream << (int)*it << " ";
    }
    stream << ". Coverages: ";
    for(auto it = coverages.begin(); it != coverages.end(); ++it){
        stream << (int)*it << " ";
    }
    return stream.str();
}

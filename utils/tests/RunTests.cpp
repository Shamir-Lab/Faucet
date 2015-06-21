#include "RunTests.h" 

//g++ ../Bloom.cpp ../Kmer.cpp ../Debloom.cpp JCheckTests.cpp KmerTests.cpp TestUtils.cpp RunTests.cpp -o RunTests

int main(int argc, char *argv[]){

    runKmerTests();
    runRollingHashTests();
    runBloomTests();
    runJCheckTests();
    //runFindNextJunctionTests();
    //runTraverseReadsTests();
    
    return 0;
}
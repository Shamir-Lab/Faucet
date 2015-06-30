#include "RunTests.h" 

//g++ ../Bloom.cpp ../Kmer.cpp ../Debloom.cpp JCheckTests.cpp KmerTests.cpp TestUtils.cpp RunTests.cpp -o RunTests

int main(int argc, char *argv[]){

    kmerTests::runKmerTests();
    rollingHashTests::runRollingHashTests();
    jCheckTests::runJCheckTests();
    runJunctionTests();
    runJunctionMapTests();
    bloomTests::runBloomTests();
    
    return 0;
}
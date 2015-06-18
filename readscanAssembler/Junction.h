#ifndef JUNCTION
#define JUNCTION

#include <iostream>
using std::ofstream;

class Junction{
public:
    unsigned char ext[4];
    unsigned char cov[4];
    unsigned char back;

    int numPathsOut();
    void writeToFile(ofstream* jFile);
    void update(int nucExt, int lengthFor, int lengthBack);
    Junction(int nucExt);
};

#endif
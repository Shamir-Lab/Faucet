#include "JuncPairs.h"

bool operator<(JuncResult a, JuncResult b){
    return a.distance < b.distance;
}

bool operator>(JuncResult a, JuncResult b){
    return a.distance > b.distance;
}
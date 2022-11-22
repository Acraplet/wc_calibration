#include "truth_alpha.hh"
#include <iostream>
#include <sstream>
#include <fstream>


int main(){
    //we need the true value of alpha
    std::cout << truth_alpha(wavelength=401.9, ABWFF=0.000486, RAYFF=10e10) std::endl;
    return 0;
}

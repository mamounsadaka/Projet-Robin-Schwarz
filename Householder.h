#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>



class Householder
{
private:

	std::vector<std::vector<double>> A_,Q_,R_;
	std::vector<double> b_;
    

public:

	Householder(std::vector<std::vector<double>> A ,std::vector<double> b) ;

public:

    void Decomposition_QR () ; // décomposition 

	std::vector<double> Solve(); // résolution  

	 // fonctions permettant de recupérer le Q et  le R  
	std::vector<std::vector<double>> get_Q() ;
	std::vector<std::vector<double>> get_R() ;
};
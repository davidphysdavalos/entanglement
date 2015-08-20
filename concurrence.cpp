//Plantilla lista pa chambear con itpp mas otras cosas
#include <iostream>
#include <cpp/dev_random.cpp>
#include <tclap/CmdLine.h>
#include <itpp/itbase.h>
#include <itpp/stat/histogram.h>
#include <cpp/itpp_ext_math.cpp>
#include <cpp/spinchain.cpp>
#include <itpp/stat/misc_stat.h>
#include <fstream>
#include <cpp/RMT.cpp>

using namespace std; 
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;
using namespace RMT;


TCLAP::CmdLine cmd("Programilla para calcular concurrencias", ' ', "0.1");
//TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"normalito", "string",cmd); // Para llamar strings
TCLAP::ValueArg<string> state("s","state", "Tipo de estado separable o bell" ,false,"bell", "string",cmd); // Para llamar strings
//TCLAP::ValueArg<unsigned int> seed("s","seed", "Random seed [0 for urandom]",false,
TCLAP::ValueArg<double> qubits("q","qubits", "qubits",false, 2,"int",cmd);
TCLAP::ValueArg<double> delta("d","delta", "delta del tiempo",false, 0.1,"double",cmd);
TCLAP::ValueArg<double> tiempo("t","tiempo", "tiempo",false, 2.0,"double",cmd);
//TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd); //llamando enteros


int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(12);

cvec init= zeros_c(4);
init(0)=1.0/sqrt(2);
cmat H = RandomGUE(4);
if(state.getValue()=="bell"){
init(3)=1.0/sqrt(2);
}
if(state.getValue()=="separable"){
init(0)=1.0;
}
int N = tiempo.getValue()/delta.getValue();

for(int i=0;i<N+1;i++){
	
cmat U = exponentiate_nonsym(-complex <double>(0,1)*(double)i*delta.getValue()*H);

init= U*init;

cout<<i*delta.getValue()<<" "<<Concurrence(init)<<endl;

}
cout<<exponentiate_nonsym(-complex <double>(0,1)*(double)0.0*delta.getValue()*H)<<endl;

}

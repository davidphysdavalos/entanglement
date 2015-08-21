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
TCLAP::ValueArg<unsigned int> seed("","seed", "Random seed [0 for urandom]",false, 0,"unsigned int",cmd);
TCLAP::ValueArg<int> members("m","members", "Ensemble members",false, 1,"int",cmd);
TCLAP::ValueArg<int> qubits("q","qubits", "qubits",false, 2,"int",cmd);
TCLAP::ValueArg<double> delta("d","delta", "delta del tiempo",false, 0.1,"double",cmd);
TCLAP::ValueArg<double> tiempo("t","tiempo", "tiempo",false, 2.0,"double",cmd);
//TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd); //llamando enteros


//double ConcurrenciaMia(cvec state){
	
	//conjugate(state)
	
	//return conc
//}

int main(int argc, char* argv[])
{

cmd.parse( argc, argv );
cout.precision(10);

// {{{ Set seed for random
unsigned int semilla=seed.getValue();
if (semilla == 0){
  Random semilla_uran; semilla=semilla_uran.strong();
} 
RNG_reset(semilla);
// }}}

cvec init= zeros_c(4);
init(0)=1.0/sqrt(2.0);

// state selection
if(state.getValue()=="bell"){
init(3)=1.0/sqrt(2.0);
}
if(state.getValue()=="separable"){
init(0)=1.0;
}

//cout<<ConcurrenceFromPure(init)<<endl;

//iterations for achieve time
int N = tiempo.getValue()/delta.getValue();


//ensemble loop
cvec psi;
//Concurrence list
vec list=zeros(N);
cmat H;
int i;
for(int s=1;s<members.getValue()+1;s++){
H = RandomGUE(4);
psi=init;
for(i=0;i<N+1;i++){

cmat U = exponentiate_nonsym(-complex <double>(0,1)*(double)i*delta.getValue()*H);

psi = U*psi;
//cout<<psi<<endl;
list(i)=list(i)*((double)s-1)/((double)s)+ConcurrenceFromPure(psi)/((double)s);

//cout<<i*delta.getValue()<<" "<<Concurrence(init)<<endl;

}


}
for(int i=0;i<N+1;i++){
	cout<<i*delta.getValue()<<" "<<list(i)<<endl;
}
//cmat H = RandomGUE(4);
//cout<<exponentiate_nonsym(-complex <double>(0,1)*(double)1.0*delta.getValue()*H)<<endl;
}

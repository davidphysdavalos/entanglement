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
//#include <itpp/base/vec.h>

using namespace std; 
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;
using namespace RMT;


TCLAP::CmdLine cmd("Programilla para calcular concurrencias", ' ', "0.9");
//TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"normalito", "string",cmd); // Para llamar strings
TCLAP::ValueArg<string> state("","state1", "Tipo de estado separable o bell" ,false,"separable","string",cmd); // Para llamar strings
TCLAP::ValueArg<string> state2("","state2", "Tipo de estado separable o bell" ,false,"separable","string",cmd); // Para llamar strings
TCLAP::ValueArg<string> ensemble("e","ensemble", "Ensamble a usar GOE o GUE" ,false,"GUE","string",cmd); // Para llamar strings
//TCLAP::ValueArg<string> interval("i","interval", "one instant or zero from the value given with -t" ,false,"int","string",cmd); // Para llamar strings
TCLAP::ValueArg<unsigned int> seed("","seed", "Random seed [0 for urandom]",false, 0,"unsigned int",cmd);
TCLAP::ValueArg<int> members("m","members", "Ensemble members",false, 1,"int",cmd);
TCLAP::ValueArg<int> qubits("q","qubits", "number of qubits",false, 2,"int",cmd);
TCLAP::ValueArg<double> delta("d","delta", "delta del tiempo",false, 0.1,"double",cmd);
TCLAP::ValueArg<double> theta("","theta", "theta de bell theta",false, 0.0,"double",cmd);
TCLAP::ValueArg<double> tiempo("t","tiempo", "tiempo",false, 2.0,"double",cmd);
//TCLAP::ValueArg<int> steps("","steps","steps",false, 100,"int",cmd); //llamando enteros


double ConcurrenceTri(cvec state){
	double M=length(state)-2.0;
	double sum=M;
	for(int i=1;i<M+1;i++){
		sum=sum-Purity(partial_trace_qubits(state,i));
	}
	return sqrt(sum)*(2/sqrt(M+2));
}

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


//{{{
cvec init, init2;
// state1 selection
if(qubits.getValue()==2){
init1=zeros_c(4);
init1(0)=1.0/sqrt(2.0);
if(state1.getValue()=="bell"){
init1(3)=1.0/sqrt(2.0);
}
if(state1.getValue()=="separable"){
init(0)=1.0;
}
if(state1.getValue()=="thetabell"){
	init(0)=cos(theta.getValue());
	init(3)=sin(theta.getValue());
}
}

if(qubits.getValue()==3){
init= zeros_c(8);
init(0)=1.0/sqrt(2.0);
if(state.getValue()=="ghz"){
init(7)=1.0/sqrt(2.0);
}
if(state1.getValue()=="separable"){
init(0)=1.0;
}
if(state1.getValue()=="bell"){
init(3)=1.0/sqrt(2.0);
}
if(state1.getValue()=="wstate"){
init(0)=0.0;
init(7)=0.0;
init(1)=1/sqrt(3.0);
init(2)=1/sqrt(3.0);
init(4)=1/sqrt(3.0);
}
if(state.getValue()=="thetaghz"){
	init(0)=cos(theta.getValue());
	init(7)=sin(theta.getValue());
}
if(state1.getValue()=="thetaw"){
init(0)=0.0;
init(7)=0.0;
init(1)=sin(theta.getValue())/sqrt(2.0);
init(2)=cos(theta.getValue());
init(4)=sin(theta.getValue())/sqrt(2.0);
}
if(state1.getValue()=="thetabell"){
	init(0)=cos(theta.getValue());
	init(3)=sin(theta.getValue());
}
}
int N = tiempo.getValue()/delta.getValue();


//ensemble loop
cvec psi;
cmat H;
if(interval.getValue()=="int"){
//Concurrence list
vec list=zeros(N);
int i;

if(qubits.getValue()==2){
for(int s=1;s<members.getValue()+1;s++){
if(ensemble.getValue()=="GUE"){
H = RandomGUE(4);
}
if(ensemble.getValue()=="GOE"){
H = to_cmat(RandomGOE(4));
}
for(i=0;i<N+1;i++){

cmat U = exponentiate_nonsym(-complex <double>(0,1)*(double)i*delta.getValue()*H);

psi = U*init;
//cout<<psi<<endl;
list(i)=list(i)*((double)s-1)/((double)s)+ConcurrenceFromPure(psi)/((double)s);

//cout<<i*delta.getValue()<<" "<<Concurrence(init)<<endl;

}


}
}

if(qubits.getValue()==3){
//for(int j=1;j<length(init)-2+1;j++){
//cout<<j<<"  "<<partial_trace_qubits(init,j)<<"    "<<Purity(partial_trace_qubits(init,j))<<endl;	
//}
//cout<<ConcurrenceTri(init)<<endl;
//cout<<(double)length(init)-2.0<<endl;
for(int s=1;s<members.getValue()+1;s++){
if(ensemble.getValue()=="GUE"){
H = RandomGUE(8);
}
if(ensemble.getValue()=="GOE"){
H = to_cmat(RandomGOE(8));
}
for(i=0;i<N+1;i++){

cmat U = exponentiate_nonsym(-complex <double>(0,1)*(double)i*delta.getValue()*H);

psi = U*init;
//cout<<psi<<endl;
list(i)=list(i)*((double)s-1)/((double)s)+ConcurrenceTri(psi)/((double)s);

//cout<<i*delta.getValue()<<" "<<Concurrence(init)<<endl;

}


}
}


for(int i=0;i<N+1;i++){
	cout<<i*delta.getValue()<<" "<<list(i)<<endl;
}

}

//cmat H = RandomGUE(4);
//cout<<exponentiate_nonsym(-complex <double>(0,1)*(double)1.0*delta.getValue()*H)<<endl;
}

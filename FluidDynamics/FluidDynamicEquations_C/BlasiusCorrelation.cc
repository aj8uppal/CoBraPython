#VERIFIED

#include "state.h"
#include <utility>
#include <cmath>
pair<double,State> BlasiusCorrelation(double Re) {
  
  if (Re < 2300){
    double Fd = 64/Re;
    State _State = State::Lam;
  } else if (Re < 20000) {
    Fd=0.3164/pow(Re,0.25);
    _State = State::Tur;
  } else {
    Fd=0.184/pow(Re,0.2);
    _State = State::Tur;
  }
  return std::make_pair(Fd, _State);
}

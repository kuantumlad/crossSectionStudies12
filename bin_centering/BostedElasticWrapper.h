#ifndef BOSTED_ELASTIC_WRAPPER_H
#define BOSTED_ELASTIC_WRAPPER_H


extern"C"{
  float elas_(float *beamEnergy, float *theta);
  float elasrad_(float *beamEnergy, float *theta, float *targetRadiationLengths, float *wCut);
}

float getRadiatedValue(float beamEnergy, float theta, float targetRadiationLengths, float wCut){
  return elasrad_(&beamEnergy, &theta, &targetRadiationLengths, &wCut);; 
}

float getValue(float beamEnergy, float theta){
  return elas_(&beamEnergy, &theta); 
}

#endif

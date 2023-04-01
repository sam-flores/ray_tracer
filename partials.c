


double uC_partialX(double xyz[3]){
  return 2*xyz[0];
}

double uC_partialY(double xyz[3]){
  return 2*xyz[1];
}

double uC_partialZ(double xyz[3]){
  return 2*xyz[2];
}

double hyper_partialX(double xyz[3]){
  return 2*xyz[0];
}

double hyper_partialY(double xyz[3]){
  return -2*xyz[1];
}

double hyper_partialZ(double xyz[3]){
  return 2*xyz[2];
}

double plane_partialX(double xyz[3]){
  return 0;
}

double plane_partialY(double xyz[3]){
  return 1;
}

double plane_partialZ(double xyz[3]){
  return 0;
}

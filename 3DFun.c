


//list of parametric functions that can be loaded into
//function variables, syntax:
/*


double (*f)(double t);
double (*g)(double t);

  f[objnum] = cos;
  g[objnum] = sin;

  then can call
  for(int t = tmin; t < tmax;t += deltaT){
  x = f(t);
  y = g(t);
}
*/

double cylX(double u, double v){
  return cos(u);
}
double cylY(double u, double v){
  return v;
}
double cylZ(double u, double v){
  return sin(u);
}

double sphX(double u, double v){
  return cos(v)*cos(u);
}
double sphY(double u, double v){
  return sin(v);
}
double sphZ(double u, double v){
  return cos(v)*sin(u);
}

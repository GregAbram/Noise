
void InitializeCatalyst(int numscripts, char *scripts);
void FinalizeCatalyst();
void CoProcessCatalyst(int extent[6], float spacing[3],
											 float *noise, float *gradient, 
											 double time, unsigned int timestep, bool lastTimeStep);

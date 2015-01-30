#include <mpi.h>

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <noise/noise.h>

#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"

using namespace noise;
using namespace std;








				vtkImageData *id = vtkImageData::New();
				id->Initialize();
				id->SetExtent(sz, ez, sy, ey, sx, ex);
				id->SetSpacing(d, d, d);
				id->SetOrigin(0, 0, 0);

				vtkFloatArray *scalars = vtkFloatArray::New();
				scalars->SetNumberOfComponents(1);
				scalars->SetArray(noise, np, 1);
				scalars->SetName("noise");
				id->GetPointData()->SetScalars(scalars);
				scalars->Delete();

				vtkFloatArray *vectors = vtkFloatArray::New();
				vectors->SetNumberOfComponents(3);
				vectors->SetArray(gradient, 3*np, 1);
				vectors->SetName("gradient");
				id->GetPointData()->SetVectors(vectors);
				vectors->Delete();

				id->Delete();

void
factor(int n, int *factors)
{
    int k = (int)(pow((double)n, 0.333) + 1);

    int i;
    for (i = k; i > 1; i--)
        if (n%i == 0) break;

    n /= i;
    k = (int)(pow((double)n, 0.5) + 1);

    int j;
    for (j = k; j > 1; j--)
        if (n%j == 0) break;

    factors[0] = i;
    factors[1] = j;
    factors[2] = n/j;
}

void
syntax(char *a)
{
    cerr << "syntax: " << a << " [options] (to stdout)\n";
    cerr << "options:\n";
    cerr << "  -r xres yres zres  overall grid resolution (256x256x256)\n";
    cerr << "  -O octave          noise octave (4)\n";
    cerr << "  -F frequency       noise frequency (8)\n";
    cerr << "  -P persistence     noise persistence (0.5)\n";
    cerr << "  -t dt nt           time series delta, number of timesteps (0, 1)\n";
    cerr << "  -m r s             set rank, size (for testing)\n";
    exit(1);
}

static module::Perlin myModule;

int main(int argc, char *argv[])
{
  int xsz = 256, ysz = 256, zsz = 256;
  float t = 3.1415926;
  int octave = 4;
  float freq = 8.0;
  float pers = 0.5;
  float delta_t = 0;
  int   nt = 1;
  int   psize = 100000000;
  int   ascii = 0;
  int		mpir, mpis;

	MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpis);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpir);

  for (int i = 1; i < argc; i++)
    if (argv[i][0] == '-') 
      switch(argv[i][1])
      {
				case 'a': ascii = 1; break;
				case 'p': psize = atoi(argv[++i]); break;
				case 'm': mpir = atoi(argv[++i]);
									mpis = atoi(argv[++i]); break;
				case 'r': xsz = atoi(argv[++i]);
									ysz = atoi(argv[++i]);
									zsz = atoi(argv[++i]); break;
				case 'P': pers = atof(argv[++i]); break;
				case 'F': freq = atof(argv[++i]); break;
				case 'O': octave = atoi(argv[++i]); break;
				case 't': delta_t = atof(argv[++i]); nt = atoi(argv[++i]); break;
				default: 
					syntax(argv[0]);
      }
    else
      syntax(argv[0]);

  myModule.SetOctaveCount(octave);
  myModule.SetFrequency(freq);
  myModule.SetPersistence(pers);

  int sz = (xsz > ysz) ? xsz : ysz;
  sz = (zsz > sz) ? zsz : sz;

  int factors[3], nparts;
  factor(mpis, factors);
  nparts = factors[0]*factors[1]*factors[2];

  float d = 1.0 / (sz-1);

  int dx = ((xsz + factors[0])-1) / factors[0];
  int dy = ((ysz + factors[1])-1) / factors[1];
  int dz = ((zsz + factors[2])-1) / factors[2];

  if (mpir < nparts)
	{
		for (int t = 0; t < nt; t++)
		{
			float T = t*delta_t;

			int ix = mpir / factors[1]*factors[2];
			int iy = (mpir % factors[1]*factors[2]) / factors[2];
			int iz = (mpi % factors[2]);

			int sx = ix * dx,
					sy = iy * dy,
					sz = iz * dz;

			int ex = (ix == (factors[0]-1)) ? xsz-1 : sx + dx,
					ey = (iy == (factors[1]-1)) ? ysz-1 : sy + dy,
					ez = (iz == (factors[2]-1)) ? zsz-1 : sz + dz;

			sx = (ix == 0) ? sx : sx-1;
			sy = (iy == 0) ? sy : sy-1;
			sz = (iz == 0) ? sz : sz-1;

			ex = (ix == (factors[0]-1)) ? ex : ex+1;
			ey = (iy == (factors[1]-1)) ? ey : ey+1;
			ez = (iz == (factors[2]-1)) ? ez : ez+1;

			int lxsz = (ex - sx) + 1,
					lysz = (ey - sy) + 1,
					lzsz = (ez - sz) + 1;

			int np = lxsz*lysz*lzsz;
			float *noise = new float[np];

				float *p = noise;
				for (int i = 0; i < lxsz; i++)
				{
					float X = (i+sx)*d;
					for (int j = 0; j < lysz; j++)
					{
						float Y = (j+sy)*d;
						for (int k = 0; k < lzsz; k++)
							*p++ = myModule.GetValue(X, Y, (k+sz)*d, T);
					}
				}

			int xstep = lzsz*lysz,
					ystep = lzsz,
					zstep = 1;

			float *gradient = new float[3*np];
			float *g = gradient;
			p = noise;
			for (int i = 0; i < lxsz; i++)
				for (int j = 0; j < lysz; j++)
					for (int k = 0; k < lzsz; k++, p++)
					{
						if (i == 0)
							*g++ = (p[xstep] - p[0]) / d;
						else if (i == (lxsz-1))
							*g++ = (p[0] - p[-xstep]) / d;
						else
							*g++ = (p[xstep] - p[-xstep]) / (2*d);

						if (j == 0)
							*g++ = (p[ystep] - p[0]) / d;
						else if (j == (lysz-1))
							*g++ = (p[0] - p[-ystep]) / d;
						else
							*g++ = (p[ystep] - p[-ystep]) / (2*d);

						if (k == 0)
							*g++ = (p[zstep] - p[0]) / d;
						else if (k == (lzsz-1))
							*g++ = (p[0] - p[-zstep]) / d;
						else
							*g++ = (p[zstep] - p[-zstep]) / (2*d);
					}
		}
  }

  MPI_Finalize();
}

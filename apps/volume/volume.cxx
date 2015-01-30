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
#include "vtkXMLImageDataWriter.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkXMLMultiBlockDataWriter.h"

using namespace noise;
using namespace std;

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
  int ascii = 0;

  for (int i = 1; i < argc; i++)
    if (argv[i][0] == '-') 
      switch(argv[i][1])
      {
				case 'a': ascii = 1; break;
				case 'p': psize = atoi(argv[++i]); break;
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
  float d = 1.0 / (sz-1);

  for (int t = 0; t < nt; t++)
  {
		float T = t*delta_t;

		int np = xsz*ysz*zsz;

		float *noise = new float[np];
		float *p = noise;
		for (int i = 0; i < xsz; i++)
		{
			float X = i*d;
			for (int j = 0; j < ysz; j++)
			{
				float Y = j*d;
				for (int k = 0; k < zsz; k++)
					*p++ = myModule.GetValue(X, Y, k*d, T);
			}
		}

    int xstep = zsz*ysz,
				ystep = zsz,
				zstep = 1;

		int kk = 0;
		float *gradient = new float[3*np];
		float *g = gradient;
		p = noise;
		for (int i = 0; i < xsz; i++)
		{
			float X = i*d;
			for (int j = 0; j < ysz; j++)
			{
				float Y = j*d;
				for (int k = 0; k < zsz; k++, p++)
				{
					if (i == 0)
						*g++ = (p[xstep] - p[0]) / d;
					else if (i == (xsz-1))
						*g++ = (p[0] - p[-xstep]) / d;
					else
						*g++ = (p[xstep] - p[-xstep]) / (2*d);
						
					if (j == 0)
						*g++ = (p[ystep] - p[0]) / d;
					else if (j == (ysz-1))
						*g++ = (p[0] - p[-ystep]) / d;
					else
						*g++ = (p[ystep] - p[-ystep]) / (2*d);
						
					if (k == 0)
						*g++ = (p[zstep] - p[0]) / d;
					else if (k == (zsz-1))
						*g++ = (p[0] - p[-zstep]) / d;
					else
						*g++ = (p[zstep] - p[-zstep]) / (2*d);
				}
			}
		}

	  vtkImageData *id = vtkImageData::New();
	  id->Initialize();
	  id->SetExtent(0, zsz-1, 0, ysz-1, 0, xsz-1);
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

	  vtkXMLImageDataWriter *wr = vtkXMLImageDataWriter::New();
	  wr->SetInputData(id);

		char fn[256];
		if (nt == 1)
			sprintf(fn, "noise.vti");
		else
			sprintf(fn, "noise-%05d.vti", t);
			
	  wr->SetFileName(fn);

	  if (ascii)
	    wr->SetDataModeToAscii();
	  else
	    wr->SetDataModeToBinary();

	  wr->Update();
	  wr->Write();
	  wr->Delete();
	  id->Delete();
	}
}

#include <mpi.h>
#include <unistd.h>

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

#include "noise.h"

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkIntArray.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDataSetWriter.h>

#include "MyServerSocket.h"

using namespace noise;
using namespace std;

int  mpir, mpis;

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
    if (mpir == 0)
    {
      cerr << "syntax: " << a << " -l layoutfile [options]\n";
      cerr << "options:\n";
      cerr << "  -l layoutfile      list of IPs or hostnames and ports of vis servers\n";
      cerr << "  -r xres yres zres  overall grid resolution (256x256x256)\n";
      cerr << "  -O octave          noise octave (4)\n";
      cerr << "  -F frequency       noise frequency (3)\n";
      cerr << "  -P persistence     noise persistence (0.5)\n";
      cerr << "  -t dt nt           time series delta, number of timesteps (0, 1)\n";
      cerr << "  -m r s             set rank, size (for testing)\n";
      cerr << "  -W                 wait for attachment\n";
    }
    MPI_Finalize();
    std::cerr << "EXIT\n";
    exit(1);
}

static module::Perlin myModule;

int main(int argc, char *argv[])
{
  int xsz = 256, ysz = 256, zsz = 256;
  float t = 3.1415926;
  int octave = 4;
  float freq = 3.0;
  float pers = 0.5;
  float delta_t = 0;
  int   nt = 1;
  int   psize = 100000000;
  char  *layoutfile = NULL;
  int   lerr, gerr;
  bool  wait_for_vis = false;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpis);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpir);

  for (int i = 1; i < argc; i++)
    if (argv[i][0] == '-') 
      switch(argv[i][1])
      {
        case 'W': wait_for_vis = true; break;
        case 'l': layoutfile = argv[++i]; break;
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

  if (! layoutfile)
    syntax(argv[0]);

  lerr = 0;
  string server;
  int port;
  ifstream ifs(layoutfile);
  for (int i = 0; !lerr && i <= mpir; i++)
  {
    ifs >> server >> port;
    if (ifs.eof())
      lerr = 1;
  }

  MPI_Allreduce(&lerr, &gerr, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (gerr)
  {
    if (mpir == 0)
      cerr << "not enough vis servers\n";
    MPI_Finalize();
    exit(1);
  }

  myModule.SetOctaveCount(octave);
  myModule.SetFrequency(freq);
  myModule.SetPersistence(pers);

  int sz = (xsz > ysz) ? xsz : ysz;
  sz = (zsz > sz) ? zsz : sz;

  int factors[3], nparts;
  factor(mpis, factors);
  nparts = factors[0]*factors[1]*factors[2];

  float d = 2.0 / (sz-1);

  int dx = ((xsz + factors[0])-1) / factors[0];
  int dy = ((ysz + factors[1])-1) / factors[1];
  int dz = ((zsz + factors[2])-1) / factors[2];


  MyServerSocket *serverSocket = MyServerSocket::New();
  serverSocket->CreateServer(port);

  if (mpir == 0)
    std::cerr << "nt: " << nt << "\n";

  for (int t = 0; t < nt; t++)
  {
    if (mpir == 0)
      std::cerr << t << "\n";

    float T = t*delta_t;

    int ix = mpir / (factors[1] * factors[2]);
    int iy = (mpir % (factors[1] * factors[2])) / factors[2];
    int iz = (mpir % (factors[1] * factors[2])) % factors[2];

    int sx = ix * dx,
        sy = iy * dy,
        sz = iz * dz;

    int ex = (ix == (factors[0]-1)) ? xsz-1 : sx + dx,
        ey = (iy == (factors[1]-1)) ? ysz-1 : sy + dy,
        ez = (iz == (factors[2]-1)) ? zsz-1 : sz + dz;

    int lxsz = (ex - sx) + 1,
        lysz = (ey - sy) + 1,
        lzsz = (ez - sz) + 1;

    int np = lxsz*lysz*lzsz;
    float *noise = new float[np];

    float *p = noise;
    for (int i = 0; i < lzsz; i++)
    {
      float Z = (i+sz)*d;
      for (int j = 0; j < lysz; j++)
      {
        float Y = (j+sy)*d;
        for (int k = 0; k < lxsz; k++)
				{
					float X = (k+sx)*d;
          *p++ = myModule.GetValue(X, Y, Z, T);
				}
      }
    }

    // Here's where we connect to the vis server and send over the data  In this
    // example we connect anew for each timestep, this isn't necessarily so - we 
    // could also connect and hold the connection until the sim ends or the server
    // goes away

    vtkSocket *skt = NULL;
    int vis_ready, ok;
    if (mpir == 0)
    {
      do
      {
	if (serverSocket->ConnectionWaiting())
	  skt = (vtkSocket *)serverSocket->WaitForConnection(1);

	if (wait_for_vis && skt == NULL) 
	{
	  std::cerr << ".";
	  sleep(1);
	}

      } while (wait_for_vis && skt == NULL);

      vis_ready = skt == NULL ? 0 : 1;
      ok = vis_ready;
    }

    MPI_Bcast(&vis_ready, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (vis_ready)
    {
      if (mpir > 0)
      {
	skt = (vtkSocket *)serverSocket->WaitForConnection(10000);
	ok = skt == NULL ? 0 : 1;
      }
    }

    int global_ok;
    MPI_Allreduce(&ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    if (vis_ready && !global_ok)
    {
      if (mpir == 0)
        std::cerr << "Root node saw vis server, but the others could not connect\n";
	
      if (skt)
      {
	int err = -1;
	skt->Send((const void *)&err, sizeof(sz));
	skt->CloseSocket();
	skt->Delete();
	skt = NULL;
      }
    }

    if (skt)
    {
      vtkImageData *id = vtkImageData::New();
      id->Initialize();
      id->SetDimensions(lxsz, lysz, lzsz);
      id->SetSpacing(d, d, d);
      id->SetOrigin(-1, -1, -1);
      
      vtkFloatArray *scalars = vtkFloatArray::New();
      scalars->SetNumberOfComponents(1);
      scalars->SetArray(noise, np, 1);
      scalars->SetName("noise");
      id->GetPointData()->SetScalars(scalars);
      scalars->Delete();

			vtkIntArray *offset = vtkIntArray::New();
			offset->SetNumberOfComponents(3);
			offset->SetNumberOfTuples(1);
			offset->SetName("offset");
			int *o = (int *)offset->GetVoidPointer(0);
			o[0] = sx;
			o[1] = sy;
			o[2] = sz;

			vtkFieldData *fd = vtkFieldData::New();
			fd->AddArray(offset);
			offset->Delete();
			
			id->SetFieldData(fd);
			fd->Delete();

#if 0
{
vtkXMLImageDataWriter *w = vtkXMLImageDataWriter::New();;
char buf[246];
sprintf(buf, "sim-%d.vti", mpir);
w->SetInputData(id);
w->SetFileName(buf);
w->Write();
w->Delete();
}
#endif

      vtkDataSetWriter *wr = vtkDataSetWriter::New();
      wr->WriteToOutputStringOn();
      wr->SetInputData(id);
      id->Delete();

      wr->Update();
    
      int sz = wr->GetOutputStringLength();
      void *ptr = wr->GetOutputString();

      skt->Send((const void *)&sz, sizeof(sz));
      skt->Send((const void *)ptr, sz);

      wr->Delete();

      skt->CloseSocket();
      skt->Delete();
    }

    if (mpir == 0)
      std::cerr << "finished timestep " << t << "\n";

		sleep(1);
  }

  if (serverSocket)
  {
    serverSocket->CloseSocket();
    serverSocket->Delete();
  }

  MPI_Finalize();
}

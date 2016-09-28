#include <mpi.h>

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <noise/noise.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDataSetWriter.h>
#include <vtkClientSocket.h>

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
      cerr << "  -F frequency       noise frequency (8)\n";
      cerr << "  -P persistence     noise persistence (0.5)\n";
      cerr << "  -t dt nt           time series delta, number of timesteps (0, 1)\n";
      cerr << "  -m r s             set rank, size (for testing)\n";
      cerr << "  -S                 open server socket and check for connection (default)\n";
      cerr << "  -C                 check to see if a vis server is waiting\n";
    }
    MPI_Finalize();
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
  int   np = -1;
  char  *layoutfile = NULL;
  int   lerr, gerr;
  bool  server_socket = true;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpis);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpir);

  for (int i = 1; i < argc; i++)
    if (argv[i][0] == '-') 
      switch(argv[i][1])
      {
        case 'S': server_socket = true; break;
        case 'C': server_socket = false; break;
        case 'l': layoutfile = argv[++i]; break;
        case 'n': np = atoi(argv[++i]); break;
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

  if (np == -1)
    np = mpis;

  myModule.SetOctaveCount(octave);
  myModule.SetFrequency(freq);
  myModule.SetPersistence(pers);

  int sz = (xsz > ysz) ? xsz : ysz;
  sz = (zsz > sz) ? zsz : sz;

  int factors[3], nparts;
  factor(np, factors);
  nparts = factors[0]*factors[1]*factors[2];

  float d = 2.0 / (sz-1);

  int dx = ((xsz + factors[0])-1) / factors[0];
  int dy = ((ysz + factors[1])-1) / factors[1];
  int dz = ((zsz + factors[2])-1) / factors[2];

  MyServerSocket *serverSocket = NULL;
  if (server_socket)
  {
    serverSocket = MyServerSocket::New();
    serverSocket->CreateServer(port);
  }

  for (int t = 0; t < nt; t++)
  {
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

    // Here's where we connect to the vis server and send over the data  In this
    // example we connect anew for each timestep, this isn't necessarily so - we 
    // could also connect and hold the connection until the sim ends or the server
    // goes away

    vtkSocket *skt = NULL;
    if (serverSocket && serverSocket->ConnectionWaiting())
    {
      skt = (vtkSocket *)serverSocket->WaitForConnection(1);
    }
    else if (! serverSocket)
    {
      vtkClientSocket *clientSocket = vtkClientSocket::New();
      if (clientSocket->ConnectToServer(server.c_str(), port))
        clientSocket->Delete();
      else
        skt = (vtkSocket *)clientSocket;
    }

    if (skt)
    {
      vtkImageData *id = vtkImageData::New();
      id->Initialize();
      id->SetExtent(sx, ex, sy, ey, sz, ez);
      id->SetSpacing(d, d, d);
      id->SetOrigin(-1, -1, -1);
      
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
  }

  if (serverSocket)
  {
    serverSocket->CloseSocket();
    serverSocket->Delete();
  }

  MPI_Finalize();
}

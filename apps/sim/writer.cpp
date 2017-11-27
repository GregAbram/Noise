#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <iostream>
#include <string>
#include <vtkSmartPointer.h>
#include <vtkCharArray.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkDataSetReader.h>
#include <vtkClientSocket.h>
#include <vtkServerSocket.h>
#include "vtkMPIController.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

using namespace vtk;
using namespace std;

int mpir, mpis;


void
syntax(char *a)
{
	if (mpir == 0)
	{
    cerr << "syntax: " << a << " -l layoutfile [options]\n";
    cerr << "options:\n";
    cerr << "  -l layoutfile      list of IPs or hostnames and ports of vis servers\n";
    cerr << "  -n knt             number of frames (1)\n";
    cerr << "  -S                 open server socket and check for connection\n";
    cerr << "  -C                 check to see if a vis server is waiting (default)\n";
	}
	MPI_Finalize();
	exit(1);
}

int
main(int argc, char *argv[])
{
	char *layoutfile = NULL;
	bool server_socket = false;
	int knt = 1;

	MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpis);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpir);

	for (int i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      switch(argv[i][1])
      {
        case 'S': server_socket = true; break;
        case 'C': server_socket = false; break;
        case 'n': knt = atoi(argv[++i]); break;
        case 'l': layoutfile = argv[++i]; break;
        default:
          syntax(argv[0]);
      }
    else
      syntax(argv[0]);

  if (! layoutfile)
    syntax(argv[0]);

	int gerr, lerr = 0;
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

  vtkServerSocket *serverSocket = NULL;
  if (server_socket)
  {
    serverSocket = vtkServerSocket::New();
    serverSocket->CreateServer(port);
  }

	int tstep = 0;
	for (int i = 0; i < knt; i++)
	{
    vtkSocket *skt = NULL;
    if (serverSocket)
    {
      skt = (vtkSocket *)serverSocket->WaitForConnection(4000000);
    }
    else
    {
      vtkClientSocket *clientSocket = vtkClientSocket::New();
      if (clientSocket->ConnectToServer(server.c_str(), port))
        clientSocket->Delete();
      else
        skt = (vtkSocket *)clientSocket;
    }

		if (! skt)
    {
      std::cerr << "socket error\n";
      break;
    }

		int sz;
		char *str;

		skt->Receive((void *)&sz, sizeof(sz));
		if (sz < 0) break;
	
		vtkSmartPointer<vtkCharArray> array = vtkSmartPointer<vtkCharArray>::New();
		array->Allocate(sz);
		array->SetNumberOfTuples(sz);
		void *ptr = array->GetVoidPointer(0);

		skt->Receive(ptr, sz);
		skt->CloseSocket();
		skt->Delete();
	
		vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
		reader->ReadFromInputStringOn();
		reader->SetInputArray(array.GetPointer());
		reader->Update();

		vtkSmartPointer<vtkXMLDataSetWriter> writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
		writer->SetInputConnection(reader->GetOutputPort());
		char filename[256];
		sprintf(filename, "noise_%04d_%04d.vti", tstep, mpir);
		writer->SetFileName(filename);
		writer->Update();

		tstep = tstep + 1;
		sleep(1);
	}

	serverSocket->Delete();

	MPI_Finalize();
	exit(0);
}

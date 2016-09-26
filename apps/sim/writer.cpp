#include <stdio.h>
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
#include "vtkCompositeRenderManager.h"
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
		cerr << "syntax: " << a << " -l layoutfile\n";
	MPI_Finalize();
	exit(1);
}

int
main(int argc, char *argv[])
{
	char *layoutfile = NULL;

	MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpis);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpir);

	for (int i = 1; i < argc; i++)
		if (! strcmp(argv[i], "-l")) { layoutfile = argv[++i]; break; }
		else syntax(argv[0]);

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

	vtkSmartPointer<vtkSocket> vtkSkt;

	vtkSmartPointer<vtkServerSocket> srvr = vtkSmartPointer<vtkServerSocket>::New();
	srvr->CreateServer(port);

	int tstep = 0;
	while (1 == 1)
	{
		vtkSmartPointer<vtkClientSocket> client = srvr->WaitForConnection(4000000);
		if (! client)
		{
			cerr << "timeout\n";
			exit(1);
		}

		vtkSkt = vtkSocket::SafeDownCast(client);

		int sz;
		char *str;

		vtkSkt->Receive((void *)&sz, sizeof(sz));
		if (sz < 0) break;
	
		vtkSmartPointer<vtkCharArray> array = vtkSmartPointer<vtkCharArray>::New();
		array->Allocate(sz);
		array->SetNumberOfTuples(sz);
		void *ptr = array->GetVoidPointer(0);

		vtkSkt->Receive(ptr, sz);
		vtkSkt->CloseSocket();
	
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
	}

	srvr->Delete();

	MPI_Finalize();
	exit(0);
}

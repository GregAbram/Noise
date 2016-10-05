#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vtkSmartPointer.h>
#include <vtkCharArray.h>
#include <vtkDataSetReader.h>
#include <vtkClientSocket.h>
#include <vtkServerSocket.h>
#include <vtkMPIController.h>
#include <vtkCompositeRenderManager.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPNGWriter.h>
#include <vtkProperty.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkFieldData.h>


using namespace vtk;
using namespace std;

int frame = 0;

vtkMPIController *controller = NULL;
vtkRenderer* renderer = NULL;
vtkRenderWindow *renderWindow = NULL;

int mpir, mpis;

void
process(vtkMultiProcessController *controller, void* arg)
{
  vtkCompositeRenderManager *renderManager = vtkCompositeRenderManager::New();
  renderManager->SetRenderWindow(renderWindow);

  renderWindow->Render();

  if (controller->GetLocalProcessId() == 0)
  {
    vtkWindowToImageFilter *window2image = vtkWindowToImageFilter::New();
    window2image->SetInput(renderWindow);

    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputConnection(window2image->GetOutputPort());
    window2image->Delete();

    char frameName[256];
    sprintf(frameName, "frame-%04d.png", frame);
    writer->SetFileName(frameName);

    writer->Write();
    writer->Delete();
  }

  renderManager->Delete();
}

void
syntax(char *a)
{
  if (mpir == 0)
  {
    cerr << "syntax: " << a << " -l layoutfile [options]\n";
    cerr << "options:\n";
		cerr << "  -n knt             number of frames to save (1)\n";
    cerr << "  -l layoutfile      list of IPs or hostnames and ports of vis servers\n";
    cerr << "  -S                 open server socket and check for connection\n";
    cerr << "  -C                 check to see if a vis server is waiting (default)\n";
  }

  controller->Finalize();
  exit(1);
}

int
main(int argc, char *argv[])
{
  char *layoutfile = NULL;
  bool  server_socket = false;
	int knt = 1;

  controller = vtkMPIController::New();
  controller->Initialize(&argc, &argv);

  mpis = controller->GetNumberOfProcesses();
  mpir = controller->GetLocalProcessId();

  for (int i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      switch(argv[i][1])
      {
        case 'n': knt = atoi(argv[++i]); break;
        case 'S': server_socket = true; break;
        case 'C': server_socket = false; break;
        case 'l': layoutfile = argv[++i]; break;
        default:
          syntax(argv[0]);
      }
    else
      syntax(argv[0]);

  if (! layoutfile)
    syntax(argv[0]);

  // Select the mpir'th line of the layout file.  That'll contain the 
  // IP address or hostname of the node this process is running on and 
  // the port it should listen on.

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

  controller->AllReduce(&lerr, &gerr, 1, vtkCommunicator::MAX_OP);
  if (gerr)
  {
    if (mpir == 0)
      cerr << "not enough vis servers\n";
    controller->Finalize();
    exit(1);
  }

  // Set up the parallel rendering pipeline

  controller->SetSingleMethod(process, (void *)NULL);
 
  renderer = vtkRenderer::New();

  vtkCamera *c = vtkCamera::New();
  c->SetPosition(3.0, 4.0, -5.0);
  c->SetFocalPoint(0.0, 0.0, 0.0);
  c->SetViewUp(0.0, 1.0, 0.0);
  c->SetClippingRange(2.0, 12.0);
  renderer->SetActiveCamera(c);
  c->Delete();

  renderWindow = vtkRenderWindow::New();
  renderWindow->SetSize(512, 512);
  renderWindow->SetOffScreenRendering(1);
  renderWindow->AddRenderer(renderer);

  // Only the root node writes

  // In this simple example, we expect there to be a 1-1 relationship between 
  // sim processes and vis processes. If told to, we create a single server socket and 
  // expect the sim to connect to it. Otherwise we connect to server on sim process

  vtkServerSocket *serverSocket = NULL;
  if (server_socket)
  {
    serverSocket = vtkServerSocket::New();
    serverSocket->CreateServer(port);
  }

  for (frame = 0; frame < knt; frame++)
  {
    renderer->RemoveAllViewProps();

    // Start with nothing to render.    Each input object will be "visualized" and one
    // or more actors will be added to the renderer for each.  In this example, a single
    // object is expected per time step

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

    // Got a connection.   Set up the reader

    int sz;
    char *str;
    skt->Receive((void *)&sz, sizeof(sz));
    if (sz < 0) break;

    std::cerr << mpir << " " << sz << " bytes\n";
  
    vtkCharArray *array = vtkCharArray::New();
    array->Allocate(sz);
    array->SetNumberOfTuples(sz);
    void *ptr = array->GetVoidPointer(0);

    skt->Receive(ptr, sz);
    skt->CloseSocket();
    skt->Delete();
  
    vtkDataSetReader *reader = vtkDataSetReader::New();
    reader->ReadFromInputStringOn();
    reader->SetInputArray(array);
    array->Delete();

		reader->Update();

		vtkImageData *id = vtkImageData::SafeDownCast(reader->GetOutput());
		if (! id)
		{
			std::cerr << "error -- received something other than an vtkImageData object\n";
			break;
		}

		vtkFieldData *fd = id->GetFieldData();
		if (fd)
		{
			int indx;
			vtkIntArray *oarray = vtkIntArray::SafeDownCast(fd->GetArray("offset", indx));
			if (oarray)
			{
				int *offsets = (int *)oarray->GetVoidPointer(0);
				double *spacing = id->GetSpacing();

				double origin[3];
				id->GetOrigin(origin);

				origin[0] += offsets[0] * spacing[0];
				origin[1] += offsets[1] * spacing[1];
				origin[2] += offsets[2] * spacing[2];

				id->SetOrigin(origin);
			}
		}

    // Here is where we set up the "visualization pipeline" - in this test,
    // we create a single (reader)->contour->mapper->actor and add it to
    // the renderer.  We could also create multiple pipelines starting 
    // with the same reader and resulting in more than one actor being     
    // added to the renderer

    vtkContourFilter *contourFilter = vtkContourFilter::New();
    contourFilter->SetValue(0, 0.00);
    contourFilter->SetInputConnection(reader->GetOutputPort());
    reader->Delete();

    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection(contourFilter->GetOutputPort());
    mapper->ScalarVisibilityOff();
    contourFilter->Delete();
 
    vtkActor *actor = vtkActor::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor((mpir & 0x1) ? 1.0 : 0.0, (mpir & 0x2) ? 1.0 : 0.5, (mpir & 0x4) ? 1.0 : 0.5);
    renderer->AddActor(actor);
    actor->Delete();

    // Render and write the result (in the process routine)

    controller->SingleMethodExecute();
  }

  if (serverSocket) serverSocket->Delete();
  renderer->Delete();
  renderWindow->Delete();

  controller->Finalize();
  exit(0);
}

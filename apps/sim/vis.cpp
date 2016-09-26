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

using namespace vtk;
using namespace std;

vtkMPIController *controller = NULL;
vtkRenderer* renderer = NULL;
vtkRenderWindow *renderWindow = NULL;
vtkPNGWriter *writer = NULL;

int mpir, mpis;
int tstep = 0;

void
process(vtkMultiProcessController *controller, void* arg)
{
	vtkCompositeRenderManager *renderManager = vtkCompositeRenderManager::New();
  renderManager->SetRenderWindow(renderWindow);

	if (controller->GetLocalProcessId() == 0)
    writer->Write();
  else
    renderWindow->Render();

	renderManager->Delete();
}

void
syntax(char *a)
{
	if (mpir == 0)
		cerr << "syntax: " << a << " -l layoutfile\n";
	controller->Finalize();
	exit(1);
}

int
main(int argc, char *argv[])
{
	char *layoutfile = NULL;

	controller = vtkMPIController::New();
  controller->Initialize(&argc, &argv);

	mpis = controller->GetNumberOfProcesses();
	mpir = controller->GetLocalProcessId();

	for (int i = 1; i < argc; i++)
		if (! strcmp(argv[i], "-l")) { layoutfile = argv[++i]; break; }
		else syntax(argv[0]);

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

	if (mpir == 0)
	{
		vtkWindowToImageFilter *window2image = vtkWindowToImageFilter::New();
    window2image->SetInput(renderWindow);

		writer = vtkPNGWriter::New();
    writer->SetInputConnection(window2image->GetOutputPort());
		window2image->Delete();
	}

	// In this simple example, we expect there to be a 1-1 relationship between 
	// sim processes and vis processes.   We create a single server socket and 
	// expect the sim to connect to it.

	vtkSmartPointer<vtkSocket> vtkSkt;
	vtkSmartPointer<vtkServerSocket> srvr = vtkSmartPointer<vtkServerSocket>::New();
	srvr->CreateServer(port);

	while (1 == 1)
	{
		// Start with nothing to render.    Each input object will be "visualized" and one
		// or more actors will be added to the renderer for each.  In this example, a single
		// object is expected per time step

		renderer->RemoveAllViewProps(); 

		vtkSmartPointer<vtkClientSocket> client = srvr->WaitForConnection(4000000);
		if (! client)
		{
			cerr << "timeout\n";
			exit(1);
		}

		vtkSkt = vtkSocket::SafeDownCast(client);

		// Got a connection.   Set up the reader

		int sz;
		char *str;
		vtkSkt->Receive((void *)&sz, sizeof(sz));
		if (sz < 0) break;

		std::cerr << mpir << " " << sz << " bytes\n";
	
		vtkCharArray *array = vtkCharArray::New();
		array->Allocate(sz);
		array->SetNumberOfTuples(sz);
		void *ptr = array->GetVoidPointer(0);

		vtkSkt->Receive(ptr, sz);
		vtkSkt->CloseSocket();
	
		vtkDataSetReader *reader = vtkDataSetReader::New();
		reader->ReadFromInputStringOn();
		reader->SetInputArray(array);
		array->Delete();

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
		contourFilter->Delete();
 
		vtkActor *actor = vtkActor::New();
		actor->SetMapper(mapper);
		renderer->AddActor(actor);
		actor->Delete();

		// This is the end of the setup of the pipeline.  Now we set up the 
		// writer with a filename and do the parallel rendering (which will 
		// execute the visualization pipeline(s)

		if (mpir == 0)
		{
			char frameName[256];
			sprintf(frameName, "frame-%04d.png", tstep);
			writer->SetFileName(frameName);
		}

		// Render and write the result (in the process routine)

		controller->SingleMethodExecute();
	}

	if (writer) writer->Delete();
	srvr->Delete();
	renderer->Delete();
	renderWindow->Delete();

	controller->Finalize();
	exit(0);
}

#include "catalyst.h"

#include "paraview-4.3/vtkCPProcessor.h"
#include "paraview-4.3/vtkCPDataDescription.h"
#include "paraview-4.3/vtkCPPythonScriptPython.h"

#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"

vtkCPProcessor *Processor = NULL;
vtkImageData *VTKGrid = NULL;

void InitializeCatalyst(int numscripts, char *scripts)
{
	if (Processor == NULL)
	{
		Processor = vtkCPProcessor::New();
		Processor->Initialize();
	}
	else
		Processor->RemoveAllPipelines();

	for (int i = 0; i < numscripts; i++)
	{
		vtkCPPythonScriptPipeline *pipeline = vtkCPPythonScriptPython::New();
		pipeline->Initialize(scripts[i]);
		Processor->AddPipeline(pipeline);
		pipeline->Delete();
	}
}

void FinalizeCatalyst()
{
	if (Processor)
	{
		Processor->Delete();
		Processor = NULL;
	}

  if (VTKGrid)
	{
		VTKGrid->Delete();
		VTKGrid = NULL;
	}
}

void CoProcessCatalyst(int extent[6], float spacing[3], 
											 float *noise, float *gradient, 
											 double time, unsigned int timestep, bool lastTimeStep)
{

	if (VTKGrid == NULL)
	{
		int np = ((extent[1] - extent[0]) + 1)
								* ((extent[3] - extent[2]) + 1)
								* ((extent[5] - extent[4]) + 1);

		VTKGrid = vtkImageData::New();
		id->Initialize();
		id->SetExtent(extent);
		id->SetSpacing(spacing);
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
	}

	vtkCPDataInputDescription *dataDesc = vtkCPDataInputDescription::New();
	dataDesc->AddInput("input");
	dataDesc->SetTimeData(time, timestep);
	if (lastTimeStep)
		dataDesc->ForceOutputOn();

  if (Processor->RequestDataInputDescription(dataDesc) != 0)
	{
		dataDesc->GetInputDescriptionByName("input")->SetGrid(VTKGrid);
		Processor->CoProcess(dataDesc);
	}

	dataDesc->Delete();
}

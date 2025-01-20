#include "vtkAutoInit.h" 
VTK_MODULE_INIT(vtkRenderingOpenGL2); // VTK was built with vtkRenderingOpenGL2
VTK_MODULE_INIT(vtkInteractionStyle);

#include "vtkTipsyReader.h"

#include "vtkActor.h"
#include "vtkGeometryFilter.h"
#include "vtkInformation.h"
#include "vtkLookupTable.h"
#include "vtkNew.h"
#include "vtkPartitionedDataSetCollection.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include <vtksys/CommandLineArguments.hxx>
#include <vtksys/SystemTools.hxx>
#include "vtkTestUtilities.h"
#include "vtkUnstructuredGrid.h"

#include <map>
#include <string>

using namespace std;

const vector<string> ptypes = {"gas", "dark", "star"};

int
TestSimpleTipsyReader(int argc, char* argv[])
{
  string filein;
  string varname, partname;
  bool vis = 0;

  double TimeStep = 0.0;
  int k, BlockIndex = 0;

  vtksys::CommandLineArguments args;
  args.Initialize(argc, argv);
  args.AddArgument(
    "-f", vtksys::CommandLineArguments::SPACE_ARGUMENT, &filein, "(the names of the Gadget (HDF5) files to read)");
  args.AddArgument(
    "-var", vtksys::CommandLineArguments::SPACE_ARGUMENT, &varname, "(the name of the SCALAR variable to display)");
  args.AddArgument(
    "-type", vtksys::CommandLineArguments::SPACE_ARGUMENT, &partname, "(the name of the SCALAR variable to display)");
  args.AddArgument(
    "-vis", vtksys::CommandLineArguments::NO_ARGUMENT, &vis, "(optional visualization) with display)");

  if ( !args.Parse() || argc == 1 || filein.empty())
    {
    cerr << "\nTestSimpleGadgetReader: Written by Jean M. Favre\n"
         << "options are:\n";
    cerr << args.GetHelp() << "\n";
    return EXIT_FAILURE;
    }

  if(!vtksys::SystemTools::FileExists(filein.c_str()))
    {
    cerr << "\nFile " << filein.c_str() << " does not exist\n\n";
    return EXIT_FAILURE;
    }

  vtkNew<vtkTipsyReader> reader;
  reader->DebugOff();
  reader->SetFileName(filein.c_str());
  reader->UpdateInformation();
  if(partname.empty())
    reader->SetParticleTypeToAll();
  else
    {
    auto it = find(ptypes.cbegin(), ptypes.cend(), partname);
    if (it != ptypes.cend())
      {
      BlockIndex = it - ptypes.cbegin();
      reader->SetParticleType(BlockIndex);
      cerr << "Partname GetParticleType" << reader->GetParticleType() << endl;
      }
    else
      {
      cerr << "Partname should be one of 'gas', 'dark', 'star'" << endl;
      return EXIT_FAILURE;
      }
    }
  for(auto i=0; i < reader->GetNumberOfPointArrays(); i++)
    cout << "found array (" << i << ") = " << reader->GetPointArrayName(i) << endl;
  
  if(varname.empty())
    reader->EnableAllPointArrays();
  else
    {
    reader->DisableAllPointArrays();
    reader->SetPointArrayStatus(varname.c_str(), 1);
    }
  if(vis)
    reader->GenerateVertexCellsOn();
  else
    reader->GenerateVertexCellsOn();
  reader->UpdateTimeStep(TimeStep); // time value
  reader->Update();

  double range[2];
      
  vtkUnstructuredGrid *FirstBlock = static_cast<vtkUnstructuredGrid *>(reader->GetOutput()->GetPartitionedDataSet(0)->GetPartition(BlockIndex));
  if(varname.size())
    {
    FirstBlock->GetPointData()->GetArray(0)->GetRange(range);
    cerr << varname.c_str() << ": scalar range = [" << range[0] << ", " << range[1] << "]\n";
    }
  
if(vis)
  {
  vtkNew<vtkLookupTable> lut;
  lut->SetHueRange(0.66,0.0);
  lut->SetNumberOfTableValues(256);
  lut->SetScaleToLog10();
  lut->Build();

  if(varname.size())
    {
    lut->SetTableRange(range[0], range[1]);
    lut->Build();
    }

  vtkNew<vtkGeometryFilter> geom1;
  geom1->SetInputData(FirstBlock);
  geom1->Update();
  
  vtkNew<vtkPolyDataMapper> mapper1;
  mapper1->SetInputConnection(geom1->GetOutputPort(0));
  mapper1->ScalarVisibilityOn();
  mapper1->SetScalarModeToUsePointFieldData();

  if(varname.size())
    {
    mapper1->SelectColorArray(varname.c_str());
    mapper1->SetLookupTable(lut);
    mapper1->UseLookupTableScalarRangeOn();
    }
  vtkNew<vtkActor> actor1;
  actor1->SetMapper(mapper1);

  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkRenderWindowInteractor> iren;

  iren->SetRenderWindow(renWin);
  renWin->AddRenderer(ren);
  ren->AddActor(actor1);

  renWin->SetSize(512, 512);
  renWin->Render();
  ren->ResetCamera();

  renWin->Render();

  iren->Start();
  }
  return EXIT_SUCCESS;
}

int
main(int argc, char* argv[])
{
  return TestSimpleTipsyReader(argc, argv);
}

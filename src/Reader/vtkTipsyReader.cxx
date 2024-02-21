/*=========================================================================

  Program:   ParaView
  Module:    vtkTipsyReader.cxx


=========================================================================*/

#include "vtkTipsyReader.h"

#include "vtkCellType.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vector>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include "vtkSmartPointer.h"

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkTipsyReader, Controller, vtkMultiProcessController);
#endif

#include <algorithm>
#include <functional>
#include <numeric>
#include "tipsy_file.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTipsyReader);
//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
  this->SetNumberOfInputPorts(0);
  this->Tipsyfile                = nullptr;
  this->SetParticleTypeToGas();
  this->TimeStep                 = 0;
  this->ActualTimeStep           = 0;
  this->GenerateVertexCells      = 1;
  this->FileName                 = nullptr;
  this->UpdatePiece              = 0;
  this->UpdateNumPieces          = 0;
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  //this->ParticleTypeSelection    = vtkDataArraySelection::New();
#ifdef PARAVIEW_USE_MPI
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}
//----------------------------------------------------------------------------
vtkTipsyReader::~vtkTipsyReader()
{
  this->CloseFile();
  delete [] this->FileName;
  this->FileName = nullptr;

  this->PointDataArraySelection->Delete();
  this->PointDataArraySelection = 0;
  
  //this->ParticleTypeSelection->Delete();
  //this->ParticleTypeSelection = 0;

#ifdef PARAVIEW_USE_MPI
  this->SetController(nullptr);
#endif
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void vtkTipsyReader::CloseFile()
{
  if (this->Tipsyfile != nullptr)
    {
    vtkDebugMacro(<< "vtkTipsyReader::CloseFile("<< this->FileName<<")\n");
    this->Tipsyfile->FileClose();
    delete this->Tipsyfile;
    this->Tipsyfile = nullptr;
    }
}
//----------------------------------------------------------------------------
int vtkTipsyReader::OpenFile()
{
  if (!this->FileName)
    {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
    }

  if (FileModifiedTime>FileOpenedTime)
    {
    this->CloseFile();
    }

  if (!this->Tipsyfile)
    {
    this->Tipsyfile = new TipsyFile(this->FileName);
    vtkDebugMacro(<< "vtkTipsyReader::OpenFile("<< this->FileName<<")\n");
    this->FileOpenedTime.Modified();
    }

  if (!this->Tipsyfile)
    {
    vtkErrorMacro(<< "Initialize: Could not open file " << this->FileName);
    return 0;
    }

  return 1;
}

//----------------------------------------------------------------------------
int vtkTipsyReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

#ifdef PARAVIEW_USE_MPI
  if (this->Controller)
    {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
    }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
#endif

  if (!this->OpenFile())
    {
    return 0;
    }

  this->Tipsyfile->read_header();
  bool CanReadFile = this->Tipsyfile->report_header();

  this->PointDataArraySelection->AddArray("mass");
  this->PointDataArraySelection->AddArray("vel");
  this->PointDataArraySelection->AddArray("rho");
  this->PointDataArraySelection->AddArray("temp");
  this->PointDataArraySelection->AddArray("hsmooth");
  this->PointDataArraySelection->AddArray("metals");
  this->PointDataArraySelection->AddArray("phi");
  
  double timeRange[2] = {this->Tipsyfile->h.time, this->Tipsyfile->h.time};

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &this->Tipsyfile->h.time, 1);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  this->CloseFile();
  return 1;
}

void  vtkTipsyReader::Read_Gas(vtkMultiBlockDataSet *mb, int N)
{
  int myType = TIPSY_TYPE_GAS;
  vtkPolyData *output = vtkPolyData::New();
  mb->SetBlock(myType, output);
  mb->GetMetaData(myType)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[myType]);
  output->Delete();
  if(N)
  {
  float *fp, *fp1;
  unsigned int fpoffset, phi_idx;

  fp = (float *)this->Tipsyfile->sph;
  fpoffset = sizeof(gas_particle)/sizeof(float);
  phi_idx = 11;
  
  vtkFloatArray *coords = vtkFloatArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(N);
  coords->SetName("coords");
  std::cerr << __LINE__ << ": Allocating coordinates array of size "<< N << "*" << sizeof(float) << " bytes = " << 3*N*sizeof(float) << " bytes\n";
  fp1 = fp;
  for(vtkIdType i=0; i < N; i++)
    {
    coords->SetTuple3(i, *(fp+1), *(fp+2), *(fp+3));
    fp += fpoffset;
    }

  vtkPointData *pd = output->GetPointData();
  if(this->GetPointArrayStatus("mass"))
    {
    fp = fp1;
    vtkFloatArray *mass = vtkFloatArray::New();
    mass->SetNumberOfComponents(1);
    mass->SetNumberOfTuples(N);
    mass->SetName("Mass");
    std::cerr << __LINE__ << ": Allocating Mass array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
    pd->AddArray(mass);
    mass->Delete();
  
    for(vtkIdType i=0; i < N; i++)
      {
      mass->SetTuple(i, fp);
      fp += fpoffset;
      }
    }
  if(this->GetPointArrayStatus("phi"))
    {
    fp = fp1;
    vtkFloatArray *phi = vtkFloatArray::New();
    phi->SetNumberOfComponents(1);
    phi->SetNumberOfTuples(N);
    phi->SetName("phi");
    std::cerr << __LINE__ << ": Allocating phi array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
    pd->AddArray(phi);
    phi->Delete();
  
    for(vtkIdType i=0; i < N; i++)
      {
      phi->SetTuple(i, fp+phi_idx);
      fp += fpoffset;
      }
    }
  if(this->GetPointArrayStatus("vel"))
    {
    fp = fp1;
    vtkFloatArray *vel = vtkFloatArray::New();
    vel->SetNumberOfComponents(3);
    vel->SetNumberOfTuples(N);
    vel->SetName("velocity");
    std::cerr << __LINE__ << ": Allocating velocity array of size "<< N << "*" << sizeof(float) << " bytes = " << 3*N*sizeof(float) << " bytes\n";
    for(vtkIdType i=0; i < N; i++)
      {
      vel->SetTuple3(i, *(fp+4), *(fp+5), *(fp+6));
      fp += fpoffset;
      }
    pd->AddArray(vel);
    vel->Delete();
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  output->SetPoints(points);
  coords->Delete();
  points->Delete();

  if(this->GetPointArrayStatus("rho"))
      {
      fp = fp1;
      vtkFloatArray *rho = vtkFloatArray::New();
      rho->SetNumberOfTuples(N);
      rho->SetName("rho");
      std::cerr << __LINE__ << ": Allocating rho array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
      fp = (float *)this->Tipsyfile->sph + 7; // the rho array;
      for(vtkIdType i=0; i < N; i++)
        {
        rho->SetTuple(i, fp);
        fp += fpoffset;
        }
      pd->AddArray(rho);
      rho->Delete();
      }

  if(this->GetPointArrayStatus("hsmooth"))
      {
      fp = fp1;
      vtkFloatArray *hsmooth = vtkFloatArray::New();
      hsmooth->SetNumberOfTuples(N);
      hsmooth->SetName("hsmooth");
      std::cerr << __LINE__ << ": Allocating hsmooth array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
      fp = (float *)this->Tipsyfile->sph + 9; // the hsmooth array;
      for(vtkIdType i=0; i < N; i++)
        {
        hsmooth->SetTuple(i, fp);
        fp += fpoffset;
        }
      pd->AddArray(hsmooth);
      hsmooth->Delete();
      }

  if (this->GenerateVertexCells)
    {
    vtkIdList *list = vtkIdList::New();
    std::cerr << __LINE__ << ": Allocating ID list of size "<< N << "*" << sizeof(vtkIdType) << " bytes = " << N*sizeof(vtkIdType) << " bytes\n";
    list->SetNumberOfIds(N);
    std::iota(list->GetPointer(0), list->GetPointer(N), 0);
    output->Allocate(1);
    output->InsertNextCell(VTK_POLY_VERTEX, list);
    list->Delete();
    }
  }
} // Read_Gas

void  vtkTipsyReader::Read_DarkMatter(vtkMultiBlockDataSet *mb, int N)
{
  int myType = TIPSY_TYPE_DARK;

  vtkPolyData *output = vtkPolyData::New();
  mb->SetBlock(myType, output);
  mb->GetMetaData(myType)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[myType]);
  output->Delete();
  if(N)
  {
  float *fp, *fp1;
  unsigned int fpoffset, phi_idx;

  fp = (float *)this->Tipsyfile->dark;
  fpoffset = sizeof(dark_particle)/sizeof(float);
  phi_idx = 8;

  vtkFloatArray *coords = vtkFloatArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(N);
  coords->SetName("coords");
  std::cerr << __LINE__ << ": Allocating coordinates array of size "<< N << "*" << sizeof(float) << " bytes = " << 3L*N*sizeof(float) << " bytes\n";
  fp1 = fp;
  for(vtkIdType i=0; i < N; i++)
    {
    coords->SetTuple3(i, *(fp+1), *(fp+2), *(fp+3));
    fp += fpoffset;
    }

  vtkPointData *pd = output->GetPointData();
  if(this->GetPointArrayStatus("mass"))
    {
    fp = fp1;
    vtkFloatArray *mass = vtkFloatArray::New();
    mass->SetNumberOfComponents(1);
    mass->SetNumberOfTuples(N);
    mass->SetName("Mass");
    std::cerr << __LINE__ << ": Allocating Mass array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
    pd->AddArray(mass);
    mass->Delete();
  
    for(vtkIdType i=0; i < N; i++)
      {
      mass->SetTuple(i, fp);
      fp += fpoffset;
      }
    }
  if(this->GetPointArrayStatus("phi"))
    {
    fp = fp1;
    vtkFloatArray *phi = vtkFloatArray::New();
    phi->SetNumberOfComponents(1);
    phi->SetNumberOfTuples(N);
    phi->SetName("phi");
    std::cerr << __LINE__ << ": Allocating phi array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
    pd->AddArray(phi);
    phi->Delete();
  
    for(vtkIdType i=0; i < N; i++)
      {
      phi->SetTuple(i, fp+phi_idx);
      fp += fpoffset;
      }
    }
  if(this->GetPointArrayStatus("vel"))
    {
    fp = fp1;
    vtkFloatArray *vel = vtkFloatArray::New();
    vel->SetNumberOfComponents(3);
    vel->SetNumberOfTuples(N);
    vel->SetName("velocity");
    std::cerr << __LINE__ << ": Allocating velocity array of size "<< N << "*" << sizeof(float) << " bytes = " << 3L*N*sizeof(float) << " bytes\n";
    for(vtkIdType i=0; i < N; i++)
      {
      vel->SetTuple3(i, *(fp+4), *(fp+5), *(fp+6));
      fp += fpoffset;
      }
    pd->AddArray(vel);
    vel->Delete();
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  output->SetPoints(points);
  coords->Delete();
  points->Delete();

  if (this->GenerateVertexCells)
    {
    vtkIdList *list = vtkIdList::New();
    std::cerr << __LINE__ << ": Allocating ID list of size "<< N << "*" << sizeof(vtkIdType) << " bytes = " << N*sizeof(vtkIdType) << " bytes\n";
    list->SetNumberOfIds(N);
    std::iota(list->GetPointer(0), list->GetPointer(N), 0);
    output->Allocate(1);
    output->InsertNextCell(VTK_POLY_VERTEX, list);
    list->Delete();
    }
  }
} // Read_DarkMatter

void  vtkTipsyReader::Read_Stars(vtkMultiBlockDataSet *mb, int N)
{
  int myType = TIPSY_TYPE_STAR;

  vtkPolyData *output = vtkPolyData::New();
  mb->SetBlock(myType, output);
  mb->GetMetaData(myType)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[myType]);
  output->Delete();
  if(N)
  {
  float *fp, *fp1;
  unsigned int fpoffset, phi_idx;

  fp = (float *)this->Tipsyfile->star;
  fpoffset = sizeof(star_particle)/sizeof(float);
  phi_idx = 10;
  
  vtkFloatArray *coords = vtkFloatArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(N);
  coords->SetName("coords");
  std::cerr << __LINE__ << ": Allocating coordinates array of size "<< N << "*" << sizeof(float) << " bytes = " << 3*N*sizeof(float) << " bytes\n";
  fp1 = fp;
  for(vtkIdType i=0; i < N; i++)
    {
    coords->SetTuple3(i, *(fp+1), *(fp+2), *(fp+3));
    fp += fpoffset;
    }

  vtkPointData *pd = output->GetPointData();
  if(this->GetPointArrayStatus("mass"))
    {
    fp = fp1;
    vtkFloatArray *mass = vtkFloatArray::New();
    mass->SetNumberOfComponents(1);
    mass->SetNumberOfTuples(N);
    mass->SetName("Mass");
    std::cerr << __LINE__ << ": Allocating Mass array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
    pd->AddArray(mass);
    mass->Delete();
  
    for(vtkIdType i=0; i < N; i++)
      {
      mass->SetTuple(i, fp);
      fp += fpoffset;
      }
    }
  if(this->GetPointArrayStatus("phi"))
    {
    fp = fp1;
    vtkFloatArray *phi = vtkFloatArray::New();
    phi->SetNumberOfComponents(1);
    phi->SetNumberOfTuples(N);
    phi->SetName("phi");
    std::cerr << __LINE__ << ": Allocating phi array of size "<< N << "*" << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";
    pd->AddArray(phi);
    phi->Delete();
  
    for(vtkIdType i=0; i < N; i++)
      {
      phi->SetTuple(i, fp+phi_idx);
      fp += fpoffset;
      }
    }
  if(this->GetPointArrayStatus("vel"))
    {
    fp = fp1;
    vtkFloatArray *vel = vtkFloatArray::New();
    vel->SetNumberOfComponents(3);
    vel->SetNumberOfTuples(N);
    vel->SetName("velocity");
    std::cerr << __LINE__ << ": Allocating velocity array of size "<< N << "*" << sizeof(float) << " bytes = " << 3*N*sizeof(float) << " bytes\n";
    for(vtkIdType i=0; i < N; i++)
      {
      vel->SetTuple3(i, *(fp+4), *(fp+5), *(fp+6));
      fp += fpoffset;
      }
    pd->AddArray(vel);
    vel->Delete();
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  output->SetPoints(points);
  coords->Delete();
  points->Delete();

  if (this->GenerateVertexCells)
    {
    vtkIdList *list = vtkIdList::New();
    std::cerr << __LINE__ << ": Allocating ID list of size "<< N << "*" << sizeof(vtkIdType) << " bytes = " << N*sizeof(vtkIdType) << " bytes\n";
    list->SetNumberOfIds(N);
    std::iota(list->GetPointer(0), list->GetPointer(N), 0);
    output->Allocate(1);
    output->InsertNextCell(VTK_POLY_VERTEX, list);
    list->Delete();
    }
  }
} // Read_DarkMatter

int vtkTipsyReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!mb)
    {
    return 0;
    }
#ifdef PARAVIEW_USE_MPI
  if (this->Controller &&
      (this->UpdatePiece != this->Controller->GetLocalProcessId() ||
       this->UpdateNumPieces != this->Controller->GetNumberOfProcesses()))
  {
    vtkDebugMacro(<< "Parallel failure, Id's not right (ignore)");
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
#endif

  if (!this->OpenFile())
    {
    return 0;
    }
  this->Tipsyfile->read_header();
  int n[3] = {0,0,0};

  if(this->ParticleType == TIPSY_TYPE_GAS)
    {
    n[0] = 1;
    this->Tipsyfile->read_gas_piece(this->UpdatePiece, this->UpdateNumPieces, n[0]);
    this->Read_Gas(mb, n[0]);
    this->Tipsyfile->Free_Sph_Buffer();
    }
 else if(this->ParticleType == TIPSY_TYPE_DARK)
    {
    n[1] = 1;
    this->Tipsyfile->read_dark_matter_piece(this->UpdatePiece, this->UpdateNumPieces, n[1]);
    this->Read_DarkMatter(mb, n[1]);
    this->Tipsyfile->Free_Dark_Buffer();
    }
  else if(this->ParticleType == TIPSY_TYPE_STAR)
    {
    n[2] = 1;
    this->Tipsyfile->read_star_piece(this->UpdatePiece, this->UpdateNumPieces, n[2]);
    this->Read_Stars(mb, n[2]);
    this->Tipsyfile->Free_Star_Buffer();
    }
  else if(this->ParticleType == TIPSY_TYPE_ALL)
    {
    n[0] = 1;
    this->Tipsyfile->read_gas_piece(this->UpdatePiece, this->UpdateNumPieces, n[0]);
    this->Read_Gas(mb, n[0]);
    this->Tipsyfile->Free_Sph_Buffer();
    n[1] = 1;
    this->Tipsyfile->read_dark_matter_piece(this->UpdatePiece, this->UpdateNumPieces, n[1]);
    this->Read_DarkMatter(mb, n[1]);
    this->Tipsyfile->Free_Dark_Buffer();
    n[2] = 1;
    this->Tipsyfile->read_star_piece(this->UpdatePiece, this->UpdateNumPieces, n[2]);
    this->Read_Stars(mb, n[2]);
    this->Tipsyfile->Free_Star_Buffer();
    }
  this->CloseFile();
  return 1;
}

//----------------------------------------------------------------------------


const char* vtkTipsyReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

int vtkTipsyReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

void vtkTipsyReader::SetPointArrayStatus(const char* name, int status)
{
  if (status != this->GetPointArrayStatus(name))
    {
    if (status)
      {
      this->PointDataArraySelection->EnableArray(name);
      }
    else
      {
      this->PointDataArraySelection->DisableArray(name);
      }
    this->Modified();
    }
}

void vtkTipsyReader::EnableAllPointArrays()
{
    this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkTipsyReader::DisableAllPointArrays()
{
    this->PointDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkTipsyReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
void vtkTipsyReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: " <<
    (this->FileName ? this->FileName : "(none)") << "\n";
}

int vtkTipsyReader::CanReadFile(const char* fname)
{
  FILE* fp;
  if ((fp = vtksys::SystemTools::Fopen(fname, "rb")) == nullptr)
  {
    return 0;
  }
  else
  {
    TipsyFile *file = new TipsyFile(fname);
    bool CanReadFile = file->read_header();
    //std::cerr << __LINE__ << ": CanReadFile  = " << CanReadFile << std::endl;
    file->FileClose();
    delete file;
    return CanReadFile;
  }
}



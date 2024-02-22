/*=========================================================================

  Program:   ParaView
  Module:    vtkTipsyReader.cxx


=========================================================================*/

#include "vtkTipsyReader.h"

#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPartitionedDataSetCollection.h"
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

//------------------------------------------------------------------------
int vtkTipsyReader::RequestDataObject(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkPartitionedDataSetCollection* output =
    vtkPartitionedDataSetCollection::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!output)
  {
    output = vtkPartitionedDataSetCollection::New();
    outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
    output->Delete();
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

vtkFloatArray *
vtkTipsyReader::GetVTKScalarArray(const char *name, unsigned int N, float *fp, unsigned int fpoffset)
{
  vtkFloatArray *vtkarray = vtkFloatArray::New();
  vtkarray->SetNumberOfComponents(1);
  vtkarray->SetNumberOfTuples(N);
  vtkarray->SetName(name);
  float *float_ptr = vtkarray->GetPointer(0);
  std::cout << __LINE__ << ": Allocating " << name << " array of size "<< N << "*"
            << sizeof(float) << " bytes = " << N*sizeof(float) << " bytes\n";

  for(vtkIdType i=0; i < N; i++)
    {
    //phi->SetTuple(i, fp+index);
    // avoid casting to double, and tuple intermediate copies
    *float_ptr++ = *fp;
    fp += fpoffset;
    }
  return vtkarray;
}

vtkFloatArray *
vtkTipsyReader::GetVTKVectorArray(const char *name, unsigned int N, float *fp, unsigned int fpoffset)
{
  vtkFloatArray *vtkarray = vtkFloatArray::New();
  vtkarray->SetNumberOfComponents(3);
  vtkarray->SetNumberOfTuples(N);
  vtkarray->SetName(name);
  float *float_ptr = vtkarray->GetPointer(0);
  std::cout << __LINE__ << ": Allocating " << name << " vector-array of size "<< N << "*"
            << sizeof(float) << " bytes = " << 3*N*sizeof(float) << " bytes\n";

  for(vtkIdType i=0; i < N; i++)
    {
    // vtkarray->SetTuple(i, fp+index);
    // avoid casting to double, and tuple intermediate copies
    std::copy(fp, fp+3, float_ptr); // source_first, source_last, destination
    fp += fpoffset; float_ptr +=3;
    }
  return vtkarray;
}

vtkPolyData *vtkTipsyReader::Read_Gas(int N)
{
  vtkPolyData *output = vtkPolyData::New();
  if(N)
  {
  unsigned int fpoffset = sizeof(gas_particle)/sizeof(float);
  float *fp = (float *)this->Tipsyfile->sph + 1;
  vtkFloatArray *coords = GetVTKVectorArray("coords", N, fp, fpoffset);

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  output->SetPoints(points);
  coords->Delete();
  points->Delete();

  vtkPointData *pd = output->GetPointData();
  if(this->GetPointArrayStatus("mass"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 0;
    vtkFloatArray *data = GetVTKScalarArray("mass", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if(this->GetPointArrayStatus("vel"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 4;
    vtkFloatArray *data = GetVTKVectorArray("vel", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if(this->GetPointArrayStatus("rho"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 7;
    vtkFloatArray *data = GetVTKScalarArray("rho", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if(this->GetPointArrayStatus("temp"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 8;
    vtkFloatArray *data = GetVTKScalarArray("temp", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if(this->GetPointArrayStatus("hsmooth"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 9;
    vtkFloatArray *data = GetVTKScalarArray("hsmooth", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }
      
  if(this->GetPointArrayStatus("metals"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 10;
    vtkFloatArray *data = GetVTKScalarArray("metals", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if(this->GetPointArrayStatus("phi"))
    {
    float *fp = (float *)this->Tipsyfile->sph + 11;
    vtkFloatArray *data = GetVTKScalarArray("metals", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if (this->GenerateVertexCells)
    {
    vtkIdList *list = vtkIdList::New();
    std::cout << __LINE__ << ": Allocating ID list of size "<< N << "*"
              << sizeof(vtkIdType) << " bytes = " << N*sizeof(vtkIdType) << " bytes\n";
    list->SetNumberOfIds(N);
    std::iota(list->GetPointer(0), list->GetPointer(N), 0);
    output->Allocate(1);
    output->InsertNextCell(VTK_POLY_VERTEX, list);
    list->Delete();
    }
  return output;
  }
  else
  {
    std::cerr << __LINE__ << "not carefully implemented. What to return?\n";
    return nullptr;
  }
} // Read_Gas

vtkPolyData *vtkTipsyReader::Read_DarkMatter(int N)
{
  vtkPolyData *output = vtkPolyData::New();
  if(N)
  {
  float *fp = (float *)this->Tipsyfile->dark + 1;
  unsigned int fpoffset = sizeof(dark_particle)/sizeof(float);

  vtkFloatArray *coords = GetVTKVectorArray("coords", N, fp, fpoffset);

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  output->SetPoints(points);
  coords->Delete();
  points->Delete();

  vtkPointData *pd = output->GetPointData();
  if(this->GetPointArrayStatus("mass"))
    {
    float *fp = (float *)this->Tipsyfile->dark + 0;
    vtkFloatArray *data = GetVTKScalarArray("mass", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }
  if(this->GetPointArrayStatus("phi"))
    {
    float *fp = (float *)this->Tipsyfile->dark + 8;
    vtkFloatArray *data = GetVTKScalarArray("phi", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }
  if(this->GetPointArrayStatus("vel"))
    {
    float *fp = (float *)this->Tipsyfile->dark + 4;
    vtkFloatArray *data = GetVTKVectorArray("vel", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if (this->GenerateVertexCells)
    {
    vtkIdList *list = vtkIdList::New();
    std::cout << __LINE__ << ": Allocating ID list of size "<< N << "*"
              << sizeof(vtkIdType) << " bytes = " << N*sizeof(vtkIdType) << " bytes\n";
    list->SetNumberOfIds(N);
    std::iota(list->GetPointer(0), list->GetPointer(N), 0);
    output->Allocate(1);
    output->InsertNextCell(VTK_POLY_VERTEX, list);
    list->Delete();
    }
  return output;
  }
  else
  {
    std::cerr << __LINE__ << "not carefully implemented. What to return?\n";
    return nullptr;
  }
} // Read_DarkMatter

vtkPolyData *vtkTipsyReader::Read_Stars(int N)
{
  vtkPolyData *output = vtkPolyData::New();
  if(N)
  {
  unsigned int fpoffset = sizeof(star_particle)/sizeof(float);
  float *fp = (float *)this->Tipsyfile->star + 1;
  vtkFloatArray *coords = GetVTKVectorArray("coords", N, fp, fpoffset);

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  output->SetPoints(points);
  coords->Delete();
  points->Delete();

  vtkPointData *pd = output->GetPointData();
  if(this->GetPointArrayStatus("mass"))
    {
    float *fp = (float *)this->Tipsyfile->star + 0;
    vtkFloatArray *data = GetVTKScalarArray("mass", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }
  if(this->GetPointArrayStatus("phi"))
    {
    float *fp = (float *)this->Tipsyfile->star + 10;
    vtkFloatArray *data = GetVTKScalarArray("phi", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }
  if(this->GetPointArrayStatus("vel"))
    {
    float *fp = (float *)this->Tipsyfile->star + 4;
    vtkFloatArray *data = GetVTKVectorArray("vel", N, fp, fpoffset);
    pd->AddArray(data); data->Delete();
    }

  if (this->GenerateVertexCells)
    {
    vtkIdList *list = vtkIdList::New();
    std::cerr << __LINE__ << ": Allocating ID list of size "<< N << "*"
              << sizeof(vtkIdType) << " bytes = " << N*sizeof(vtkIdType) << " bytes\n";
    list->SetNumberOfIds(N);
    std::iota(list->GetPointer(0), list->GetPointer(N), 0);
    output->Allocate(1);
    output->InsertNextCell(VTK_POLY_VERTEX, list);
    list->Delete();
    }
  return output;
  }
  else
  {
    std::cerr << __LINE__ << "not carefully implemented. What to return?\n";
    return nullptr;
  }
} // Read_Stars

int vtkTipsyReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  auto output =
    vtkPartitionedDataSetCollection::SafeDownCast(vtkDataObject::GetData(outputVector, 0));

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
  
  vtkPolyData *pd;
  if(this->ParticleType == TIPSY_TYPE_GAS)
    {
    output->SetNumberOfPartitionedDataSets(1);
    output->SetNumberOfPartitions(0, 1);
    n[0] = 1;
    this->Tipsyfile->read_gas_piece(this->UpdatePiece, this->UpdateNumPieces, n[0]);
    pd = this->Read_Gas(n[0]);
    output->SetPartition(0, 0, pd);
    pd->Delete();
    output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[TIPSY_TYPE_GAS]);
    this->Tipsyfile->Free_Sph_Buffer();
    }
 else if(this->ParticleType == TIPSY_TYPE_DARK)
    {
    output->SetNumberOfPartitionedDataSets(1);
    output->SetNumberOfPartitions(0, 1);
    n[1] = 1;
    this->Tipsyfile->read_dark_matter_piece(this->UpdatePiece, this->UpdateNumPieces, n[1]);
    pd = this->Read_DarkMatter(n[1]);
    output->SetPartition(0, 0, pd);
    pd->Delete();
    output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[TIPSY_TYPE_DARK]);
    this->Tipsyfile->Free_Dark_Buffer();
    }
  else if(this->ParticleType == TIPSY_TYPE_STAR)
    {
    output->SetNumberOfPartitionedDataSets(1);
    output->SetNumberOfPartitions(0, 1);
    n[2] = 1;
    this->Tipsyfile->read_star_piece(this->UpdatePiece, this->UpdateNumPieces, n[2]);
    pd = this->Read_Stars(n[2]);
    output->SetPartition(0, 0, pd);
    pd->Delete();
    output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[TIPSY_TYPE_STAR]);
    this->Tipsyfile->Free_Star_Buffer();
    }
  else if(this->ParticleType == TIPSY_TYPE_ALL)
    {
    output->SetNumberOfPartitionedDataSets(3);
    output->SetNumberOfPartitions(0, 1);
    output->SetNumberOfPartitions(1, 1);
    output->SetNumberOfPartitions(2, 1);
    n[0] = 1;
    this->Tipsyfile->read_gas_piece(this->UpdatePiece, this->UpdateNumPieces, n[0]);
    output->SetPartition(0, 0, this->Read_Gas(n[0]));
    output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[TIPSY_TYPE_GAS]);
    this->Tipsyfile->Free_Sph_Buffer();
    n[1] = 1;
    this->Tipsyfile->read_dark_matter_piece(this->UpdatePiece, this->UpdateNumPieces, n[1]);
    output->SetPartition(1, 0, this->Read_DarkMatter(n[1]));
    output->GetMetaData(1u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[TIPSY_TYPE_DARK]);
    this->Tipsyfile->Free_Dark_Buffer();
    n[2] = 1;
    this->Tipsyfile->read_star_piece(this->UpdatePiece, this->UpdateNumPieces, n[2]);
    output->SetPartition(2, 0, this->Read_Stars(n[2]));
    output->GetMetaData(2u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[TIPSY_TYPE_STAR]);
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


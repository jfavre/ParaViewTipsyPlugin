/*=========================================================================

  Program:   ParaView
  Module:    vtkTipsyReader.cxx

=========================================================================*/

#include "vtkTipsyReader.h"

#include "vtkCellArray.h"
#include "vtkCellType.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPartitionedDataSet.h"
#include "vtkPartitionedDataSetCollection.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkVersion.h"
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
#include <cstddef>
#include <cstring>
#include <functional>
#include <numeric>

#include "tipsy_file.h"

//----------------------------------------------------------------------------
// Field tables : one entry per readable field of each particle type. The
// offsets are computed with offsetof() so the GUI array list (RequestInformation)
// and the field extraction (BuildPolyData) can never disagree.
static_assert(offsetof(gas_particle, pos)  == sizeof(float),
              "gas pos must start at float 1 (coords = points + 1 offset)");
static_assert(offsetof(dark_particle, pos) == sizeof(float),
              "dark pos must start at float 1");
static_assert(offsetof(star_particle, pos) == sizeof(float),
              "star pos must start at float 1");

const std::vector<vtkTipsyReader::TipsyField>& vtkTipsyReader::GasParticleFields()
{
  static const std::vector<TipsyField> fields = {
    {"mass",    offsetof(gas_particle, mass)    / sizeof(float), 1},
    {"coords",  offsetof(gas_particle, pos)     / sizeof(float), 3},
    {"vel",     offsetof(gas_particle, vel)     / sizeof(float), 3},
    {"rho",     offsetof(gas_particle, rho)     / sizeof(float), 1},
    {"temp",    offsetof(gas_particle, temp)    / sizeof(float), 1},
    {"hsmooth", offsetof(gas_particle, hsmooth) / sizeof(float), 1},
    {"metals",  offsetof(gas_particle, metals)  / sizeof(float), 1},
    {"phi",     offsetof(gas_particle, phi)     / sizeof(float), 1},
  };
  return fields;
}

const std::vector<vtkTipsyReader::TipsyField>& vtkTipsyReader::DarkParticleFields()
{
  static const std::vector<TipsyField> fields = {
    {"mass",   offsetof(dark_particle, mass) / sizeof(float), 1},
    {"coords", offsetof(dark_particle, pos)  / sizeof(float), 3},
    {"vel",    offsetof(dark_particle, vel)  / sizeof(float), 3},
    {"eps",    offsetof(dark_particle, eps)  / sizeof(float), 1},
    {"phi",    offsetof(dark_particle, phi)  / sizeof(float), 1},
  };
  return fields;
}

const std::vector<vtkTipsyReader::TipsyField>& vtkTipsyReader::StarParticleFields()
{
  static const std::vector<TipsyField> fields = {
    {"mass",   offsetof(star_particle, mass)   / sizeof(float), 1},
    {"coords", offsetof(star_particle, pos)    / sizeof(float), 3},
    {"vel",    offsetof(star_particle, vel)    / sizeof(float), 3},
    {"metals", offsetof(star_particle, metals) / sizeof(float), 1},
    {"tform",  offsetof(star_particle, tform)  / sizeof(float), 1},
    {"eps",    offsetof(star_particle, eps)    / sizeof(float), 1},
    {"phi",    offsetof(star_particle, phi)    / sizeof(float), 1},
  };
  return fields;
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTipsyReader);
//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
  this->SetNumberOfInputPorts(0);
  this->Tipsyfile                = nullptr;
  this->SetParticleTypeToGas();
  this->FileOpenedMTime          = 0.0;
  this->VertexCells              = vertexCellType::PolyVertex;
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

  // reopen the file when its last-modification time changed on disk
  // since we opened it
  double mtime = 0.0;
  vtksys::SystemTools::Stat_t statBuf;
  if (vtksys::SystemTools::Stat(this->FileName, &statBuf) == 0)
    {
    mtime = static_cast<double>(statBuf.st_mtime);
    }
  if (this->Tipsyfile != nullptr && mtime != this->FileOpenedMTime)
    {
    vtkDebugMacro(<< "OpenFile: file " << this->FileName << " changed on disk, reopening");
    this->CloseFile();
    }

  if (!this->Tipsyfile)
    {
    this->Tipsyfile = new TipsyFile();
    }

  if (!this->Tipsyfile->FileOpen(this->FileName))
    {
    delete this->Tipsyfile;
    this->Tipsyfile = nullptr;
    vtkErrorMacro(<< "OpenFile: Could not open file " << this->FileName);
    return 0;
    }

  vtkDebugMacro(<< "vtkTipsyReader::OpenFile("<< this->FileName<<")\n");
  this->FileOpenedMTime = mtime;

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
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
#endif

  if (!this->OpenFile())
    {
    return 0;
    }

  if (!this->Tipsyfile->read_header() || !this->Tipsyfile->report_header())
    {
    vtkErrorMacro(<< "RequestInformation: invalid Tipsy header in " << this->FileName);
    this->CloseFile();
    return 0;
    }

  // the union of all readable fields over all particle types; 'coords' is
  // the points and is not selectable
  for (const auto* fields : {&GasParticleFields(), &DarkParticleFields(), &StarParticleFields()})
    {
    for (const auto& field : *fields)
      {
      if (strcmp(field.Name, "coords") != 0)
        {
        this->PointDataArraySelection->AddArray(field.Name);
        }
      }
    }
  
  double timeRange[2] = {this->Tipsyfile->h.time, this->Tipsyfile->h.time};

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &this->Tipsyfile->h.time, 1);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  this->CloseFile();
  return 1;
}

vtkSmartPointer<vtkFloatArray>
vtkTipsyReader::GetVTKScalarArray(const char *name, vtkIdType N, const float *fp, std::size_t fpoffset)
{
  vtkNew<vtkFloatArray> vtkarray;
  vtkarray->SetNumberOfComponents(1);
  vtkarray->SetNumberOfTuples(N);
  vtkarray->SetName(name);
  float *float_ptr = vtkarray->GetPointer(0);
  vtkDebugMacro(<< "Allocating " << name << " array of " << N
                << " floats = " << N*sizeof(float)/1024 << " Kbytes");

  for(vtkIdType i=0; i < N; i++)
    {
    // avoid casting to double, and tuple intermediate copies
    *float_ptr++ = *fp;
    fp += fpoffset;
    }
  return vtkSmartPointer<vtkFloatArray>(vtkarray.GetPointer());

}

vtkSmartPointer<vtkFloatArray>
vtkTipsyReader::GetVTKVectorArray(const char *name, vtkIdType N, const float *fp, std::size_t fpoffset)
{
  vtkNew<vtkFloatArray> vtkarray;
  vtkarray->SetNumberOfComponents(3);
  vtkarray->SetNumberOfTuples(N);
  vtkarray->SetName(name);
  float *float_ptr = vtkarray->GetPointer(0);
  vtkDebugMacro(<< "Allocating " << name << " 3-component array of " << 3*N
                << " floats = " << 3*N*sizeof(float)/1024 << " Kbytes");

  for(vtkIdType i=0; i < N; i++)
    {
    // vtkarray->SetTuple(i, fp+index);
    // avoid casting to double, and tuple intermediate copies
    std::copy(fp, fp+3, float_ptr); // source_first, source_last, destination
    fp += fpoffset; float_ptr +=3;
    }
  return vtkSmartPointer<vtkFloatArray>(vtkarray.GetPointer());
}

void vtkTipsyReader::GenerateCells(vtkPolyData *output, vtkIdType N)
{
  switch (this->VertexCells)
  {
    case vertexCellType::CellVertex:
    // one VTK_VERTEX cell per point, offsets generated internally as an
    // affine array by vtkCellArray::SetData(cellSize, connectivity)
    {
    vtkDebugMacro(<< "Generating " << N << " vertex cells");
    vtkNew<vtkIdTypeArray> polyVertex;
    polyVertex->SetNumberOfTuples(N);
    std::iota(polyVertex->GetPointer(0), polyVertex->GetPointer(0) + N, vtkIdType{0});
    vtkNew<vtkCellArray> verts;
    verts->SetData(1, polyVertex);
    output->SetVerts(verts);
    }
    break;
    case vertexCellType::PolyVertex:
    // identical to vtkConvertToPointCloud::POLYVERTEX_CELL
    {
    vtkDebugMacro(<< "Generating one poly-vertex cell over " << N << " points");
    std::vector<vtkIdType> polyVertex(N);
    std::iota(polyVertex.begin(), polyVertex.end(), 0);
    vtkNew<vtkCellArray> verts;
    verts->InsertNextCell(N, polyVertex.data());
    output->SetVerts(verts);
    }
    break;
    case vertexCellType::NoCells:
    default:
    break;
  }
}

vtkPolyData* vtkTipsyReader::BuildPolyData(const float *fp, std::size_t strideInFloats,
                                           const std::vector<TipsyField>& fields, vtkIdType N)
{
  if(N < 1)
    return nullptr;

  vtkPolyData *output = vtkPolyData::New();
  vtkPointData *pd = output->GetPointData();

  for (const auto& field : fields)
    {
    if (strcmp(field.Name, "coords") == 0)
      {
      // 'coords' becomes the vtkPoints of the dataset
      vtkSmartPointer<vtkFloatArray> coords =
        GetVTKVectorArray("coords", N, fp + field.OffsetInFloats, strideInFloats);
      vtkNew<vtkPoints> points;
      points->SetData(coords);
      output->SetPoints(points);
      }
    else if (this->GetPointArrayStatus(field.Name))
      {
      vtkSmartPointer<vtkFloatArray> data = field.Components == 3
        ? GetVTKVectorArray(field.Name, N, fp + field.OffsetInFloats, strideInFloats)
        : GetVTKScalarArray(field.Name, N, fp + field.OffsetInFloats, strideInFloats);
      pd->AddArray(data);
      }
    }

  if (this->VertexCells != vertexCellType::NoCells)
    this->GenerateCells(output, N);
  return output;
}

vtkPolyData* vtkTipsyReader::Read_Gas(vtkIdType N)
{
  return BuildPolyData(this->Tipsyfile->gas_ptr(),
                       TipsyFile::stride_of_gas_particle, GasParticleFields(), N);
} // Read_Gas

vtkPolyData *vtkTipsyReader::Read_DarkMatter(vtkIdType N)
{
  return BuildPolyData(this->Tipsyfile->dark_ptr(),
                       TipsyFile::stride_of_dark_particle, DarkParticleFields(), N);
} // Read_DarkMatter

vtkPolyData *vtkTipsyReader::Read_Stars(vtkIdType N)
{
  return BuildPolyData(this->Tipsyfile->star_ptr(),
                       TipsyFile::stride_of_star_particle, StarParticleFields(), N);
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
  if (!this->Tipsyfile->read_header())
    {
    vtkErrorMacro(<< "RequestData: invalid Tipsy header in " << this->FileName);
    this->CloseFile();
    return 0;
    }
  int n[3] = {0,0,0};

  if(this->ParticleType == particleType::Gas)
    {
    n[0] = 1;
    if (!this->Tipsyfile->read_gas_piece(this->UpdatePiece, this->UpdateNumPieces, n[0]))
      {
      vtkErrorMacro(<< "RequestData: error reading gas particles");
      }
    vtkPolyData *gas = this->Read_Gas(n[0]);
    if (gas)
      {
      output->SetPartition(0, 0, gas);
      output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[particleType::Gas]);
      gas->Delete();
      }
    this->Tipsyfile->Free_Sph_Buffer();
    }
 else if(this->ParticleType == particleType::Dark)
    {
    n[1] = 1;
    if (!this->Tipsyfile->read_dark_matter_piece(this->UpdatePiece, this->UpdateNumPieces, n[1]))
      {
      vtkErrorMacro(<< "RequestData: error reading dark matter particles");
      }
    vtkPolyData *darkmatter = this->Read_DarkMatter(n[1]);
    if (darkmatter)
      {
      output->SetPartition(0, 0, darkmatter);
      output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[particleType::Dark]);
      darkmatter->Delete();
      }
    this->Tipsyfile->Free_Dark_Buffer();
    }
  else if(this->ParticleType == particleType::Star)
    {
    n[2] = 1;
    if (!this->Tipsyfile->read_star_piece(this->UpdatePiece, this->UpdateNumPieces, n[2]))
      {
      vtkErrorMacro(<< "RequestData: error reading star particles");
      }
    vtkPolyData *stars = this->Read_Stars(n[2]);
    if (stars)
      {
      output->SetPartition(0, 0, stars);
      output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[particleType::Star]);
      stars->Delete();
      }
    this->Tipsyfile->Free_Star_Buffer();
    }
  else if(this->ParticleType == particleType::All)
    {
    output->SetNumberOfPartitionedDataSets(3);
    output->SetNumberOfPartitions(0, 1);
    output->SetNumberOfPartitions(1, 1);
    output->SetNumberOfPartitions(2, 1);
    n[0] = 1;
    this->Tipsyfile->read_gas_piece(this->UpdatePiece, this->UpdateNumPieces, n[0]);
    vtkPolyData *gas = this->Read_Gas(n[0]);
    if(gas)
      {
      output->SetPartition(0, 0, gas);
      output->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[particleType::Gas]);
      gas->Delete();
      }
    this->Tipsyfile->Free_Sph_Buffer();
    
    n[1] = 1;
    this->Tipsyfile->read_dark_matter_piece(this->UpdatePiece, this->UpdateNumPieces, n[1]);
    vtkPolyData *darkmatter = this->Read_DarkMatter(n[1]);
    if(darkmatter)
      {
      output->SetPartition(1, 0, darkmatter);
      output->GetMetaData(1u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[particleType::Dark]);
      darkmatter->Delete();
      }
    this->Tipsyfile->Free_Dark_Buffer();
    
    n[2] = 1;
    this->Tipsyfile->read_star_piece(this->UpdatePiece, this->UpdateNumPieces, n[2]);
    vtkPolyData *stars = this->Read_Stars(n[2]);
    if(stars)
      {
      output->SetPartition(2, 0, stars);
      output->GetMetaData(2u)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[particleType::Star]);
      stars->Delete();
      }
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
    this->Modified();
}

//----------------------------------------------------------------------------
void vtkTipsyReader::DisableAllPointArrays()
{
    this->PointDataArraySelection->DisableAllArrays();
    this->Modified();
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
  os << indent << "ParticleType: " << ParticleTypes[this->ParticleType] << "\n";
  os << indent << "VertexCells: " << VertexCellTypes[this->VertexCells] << "\n";
}

int vtkTipsyReader::CanReadFile(const char* fname)
{
  // a file is acceptable if it opens and its header is consistent
  // (nbodies == nsph + ndark + nstar, with or without byte swapping)
  TipsyFile file;
  if (!file.FileOpen(fname))
    {
    return 0;
    }
  return file.read_header() ? 1 : 0;
}


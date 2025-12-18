// .NAME vtkTipsyReader
// .SECTION Description
// vtkTipsyReader reads 
// 
// .SECTION Thanks
// Jean M. Favre
// CSCS - Swiss National Supercomputing Centre for creating and contributing
// this class.

#ifndef vtkTipsyReader_h
#define vtkTipsyReader_h

#include "vtkIOTipsyModule.h" // for export macro

#include "vtkPartitionedDataSetCollectionAlgorithm.h"
#include <string>
#include <vector>
#include <sstream>

class TipsyFile;
class vtkDataArraySelection;
class vtkStdString;
class vtkFloatArray;
class vtkPolyData;
class vtkMultiProcessController;

enum particleType :int {Gas=0, Dark=1, Star=2, All=3};
static std::vector<std::string> ParticleTypes = {"Gas", "Dark", "Star"};

class VTKIOTIPSY_EXPORT vtkTipsyReader : public vtkPartitionedDataSetCollectionAlgorithm
{
public:
  static vtkTipsyReader *New();
  vtkTypeMacro(vtkTipsyReader, vtkPartitionedDataSetCollectionAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Set/Get the timestep to be read
  vtkSetMacro(TimeStep,int);
  vtkGetMacro(TimeStep,int);

  vtkSetClampMacro(ParticleType, int, particleType::Gas, particleType::All);
  vtkGetMacro(ParticleType, int);
  void SetParticleTypeToGas() { this->SetParticleType(particleType::Gas); }
  void SetParticleTypeToDark() { this->SetParticleType(particleType::Dark); }
  void SetParticleTypeToStar() { this->SetParticleType(particleType::Star); }
  void SetParticleTypeToAll() { this->SetParticleType(particleType::All); }
  
  // Description:
  // When set (default no), the reader will generate a vertex cell
  // for each point/particle read. When using the points directly
  // this is unnecessary and time can be saved by omitting cell generation
  // vtkPointSpriteMapper does not require them.
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable
  vtkSetMacro(GenerateVertexCells, int);
  vtkGetMacro(GenerateVertexCells, int);
  vtkBooleanMacro(GenerateVertexCells, int);

  int         GetNumberOfParticleTypeArrays() { return 3; }

  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAllPointArrays();
  void        EnableAllPointArrays();
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

#ifdef PARAVIEW_USE_MPI

    // Description:
    // Set/Get the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);
#endif
    int CanReadFile(const char* fname);

protected:
   vtkTipsyReader();
  ~vtkTipsyReader() override;
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)  override;
  int   OpenFile();
  void  CloseFile();
  vtkPolyData* Read_Gas(int N);
  vtkPolyData* Read_DarkMatter(int N);
  vtkPolyData* Read_Stars(int N);
  void GenerateCells(vtkPolyData *output, int N);
  vtkFloatArray* GetVTKScalarArray(const char *name, unsigned int N, float *fp, unsigned int poffset);
  vtkFloatArray* GetVTKVectorArray(const char *name, unsigned int N, float *fp, unsigned int poffset);
  // Internal Variables
  //
  char         *FileName;
  TipsyFile    *Tipsyfile;
  int           TimeStep;
  int           ActualTimeStep;
  int           GenerateVertexCells;
  vtkTimeStamp  FileModifiedTime;
  vtkTimeStamp  FileOpenedTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  int           ParticleType;

  typedef std::vector<std::string>  stringlist;
  std::vector<stringlist>           FieldArrays;

  // To allow ParaView GUI to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkTipsyReader(const vtkTipsyReader&)  = delete;
  void operator=(const vtkTipsyReader&)  = delete;
};

#endif

// .NAME vtkTipsyReader
// .SECTION Description
// vtkTipsyReader reads 
// 
// .SECTION Thanks
// Tim Dykes
// Jean Favre
// CSCS - Swiss National Supercomputing Centre for creating and contributing
// this class.

#ifndef vtkTipsyReader_h
#define vtkTipsyReader_h

#include "vtkMultiBlockDataSetAlgorithm.h"
#include <string>
#include <vector>
#include <sstream>

class TipsyFile;
class vtkDataArraySelection;
class vtkStdString;
class vtkMultiProcessController;

//enum class particleType {Gas=0, Dark=1, Star=2, All=3};
#define TIPSY_TYPE_GAS 0
#define TIPSY_TYPE_DARK 1
#define TIPSY_TYPE_STAR 2
#define TIPSY_TYPE_ALL 3
static std::vector<std::string> ParticleTypes = {"Gas", "Dark", "Star"};

class vtkTipsyReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkTipsyReader *New();
  vtkTypeMacro(vtkTipsyReader,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);   

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Set/Get the timestep to be read
  vtkSetMacro(TimeStep,int);
  vtkGetMacro(TimeStep,int);

  vtkSetClampMacro(ParticleType, int, TIPSY_TYPE_GAS, TIPSY_TYPE_ALL);
  vtkGetMacro(ParticleType, int);
  void SetParticleTypeToGas() { this->SetParticleType(TIPSY_TYPE_GAS); }
  void SetParticleTypeToDark() { this->SetParticleType(TIPSY_TYPE_DARK); }
  void SetParticleTypeToStar() { this->SetParticleType(TIPSY_TYPE_STAR); }
  void SetParticleTypeToAll() { this->SetParticleType(TIPSY_TYPE_ALL); }
  
  //vtkSetEnumMacro(ParticleType, particleType);
  //vtkGetEnumMacro(ParticleType, particleType);
  
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
  //const char* GetParticleTypeArrayName(int index);
  //int         GetParticleTypeArrayStatus(const char* name);
  //void        SetParticleTypeArrayStatus(const char* name, int status);
  void        EnableAllParticleTypes();
  void        DisableAllParticleTypes();

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
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   OpenFile();
  void  CloseFile();
  void  Read_Gas(vtkMultiBlockDataSet *mb, int N);
  void  Read_DarkMatter(vtkMultiBlockDataSet *mb, int N);
  void  Read_Stars(vtkMultiBlockDataSet *mb, int N);
  //
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

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;
  // To allow paraview gui to enable/disable block (particle type) selective reading
  // vtkDataArraySelection* ParticleTypeSelection;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkTipsyReader(const vtkTipsyReader&)  = delete;
  void operator=(const vtkTipsyReader&)  = delete;
};

#endif

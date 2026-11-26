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
#include "vtkSmartPointer.h"
#include <cstddef>
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
static std::vector<std::string> ParticleTypes = {"Gas", "Dark", "Star", "All"};

enum vertexCellType :int {NoCells=0, CellVertex=1, PolyVertex=2};
static std::vector<std::string> VertexCellTypes = {"None", "Cell Vertex", "Poly Vertex"};

class VTKIOTIPSY_EXPORT vtkTipsyReader : public vtkPartitionedDataSetCollectionAlgorithm
{
public:
  static vtkTipsyReader *New();
  vtkTypeMacro(vtkTipsyReader, vtkPartitionedDataSetCollectionAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Set/Get the particle type to be read. The Tipsy format stores
  // gas (SPH), dark matter and star particles in three consecutive sections.
  vtkSetClampMacro(ParticleType, int, particleType::Gas, particleType::All);
  vtkGetMacro(ParticleType, int);
  void SetParticleTypeToGas() { this->SetParticleType(particleType::Gas); }
  void SetParticleTypeToDark() { this->SetParticleType(particleType::Dark); }
  void SetParticleTypeToStar() { this->SetParticleType(particleType::Star); }
  void SetParticleTypeToAll() { this->SetParticleType(particleType::All); }
  
  // Description:
  // Controls the cell generation of the particles/points read:
  // NoCells    -- no cells are generated;
  // CellVertex -- one VTK_VERTEX cell per particle/point;
  // PolyVertex -- one VTK_POLY_VERTEX cell containing all particles/points.
  // When using the points directly (vtkPointSpriteMapper) cell generation
  // is unnecessary and time can be saved by omitting it ('NoCells').
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable. 'CellVertex' is the most
  // common choice, 'PolyVertex' generates a single cell referencing
  // all points (identical to vtkConvertToPointCloud::POLYVERTEX_CELL).
  vtkSetClampMacro(VertexCells, int, vertexCellType::NoCells, vertexCellType::PolyVertex);
  vtkGetMacro(VertexCells, int);
  void SetVertexCellsToNone() { this->SetVertexCells(vertexCellType::NoCells); }
  void SetVertexCellsToCellVertex() { this->SetVertexCells(vertexCellType::CellVertex); }
  void SetVertexCellsToPolyVertex() { this->SetVertexCells(vertexCellType::PolyVertex); }

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

  // Description:
  // Description of one floating point field inside a particle struct:
  // its VTK array name, its offset inside the struct (in units of float)
  // and its number of components. The per-particle-type tables guarantee
  // that the GUI array list and the field extraction cannot disagree.
  struct TipsyField
  {
    const char*  Name;
    std::size_t  OffsetInFloats;
    int          Components;
  };

protected:
  vtkTipsyReader();
  ~vtkTipsyReader() override;
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *)  override;
  int   OpenFile();
  void  CloseFile();

  // field tables, one per particle type, offsets computed with offsetof()
  static const std::vector<TipsyField>& GasParticleFields();
  static const std::vector<TipsyField>& DarkParticleFields();
  static const std::vector<TipsyField>& StarParticleFields();

  // build a vtkPolyData (points + selected fields + vertex cells) from
  // N particles packed in a struct of size strideInFloats floats
  vtkPolyData* BuildPolyData(const float *fp, std::size_t strideInFloats,
                             const std::vector<TipsyField>& fields, vtkIdType N);
  vtkPolyData* Read_Gas(vtkIdType N);
  vtkPolyData* Read_DarkMatter(vtkIdType N);
  vtkPolyData* Read_Stars(vtkIdType N);
  void GenerateCells(vtkPolyData *output, vtkIdType N);
  vtkSmartPointer<vtkFloatArray> GetVTKScalarArray(const char *name, vtkIdType N,
      const float *fp, std::size_t strideInFloats);
  vtkSmartPointer<vtkFloatArray> GetVTKVectorArray(const char *name, vtkIdType N,
      const float *fp, std::size_t strideInFloats);
  // Internal Variables
  //
  char         *FileName;
  TipsyFile    *Tipsyfile;
  // last-modification time of the file when it was opened, used to detect
  // that the file changed on disk and must be reopened
  double        FileOpenedMTime;
  int           UpdatePiece;
  int           UpdateNumPieces;
  int           ParticleType;
  int           VertexCells;

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

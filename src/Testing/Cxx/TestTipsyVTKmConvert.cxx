#include <stddef.h>
#include <iomanip>
#include "tipsy_file.h"
#include <vtkm/cont/Timer.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/ArrayHandleExtractComponent.h>
#include <vtkm/cont/ArrayHandleCompositeVector.h>
#include <vtkm/cont/ArrayHandleStride.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/FieldRangeCompute.h>
#include <vtkm/cont/BoundsCompute.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/io/VTKDataSetWriter.h>
#include <vtkm/filter/resampling/HistSampling.h>
#include <vtkm/filter/density_estimate/ParticleDensityCloudInCell.h>
#include <vtkm/filter/geometry_refinement/VertexClustering.h>
#include <vtkm/filter/density_estimate/NDHistogram.h>
#include <vtkm/filter/entity_extraction/ExtractPoints.h>
#include <vtkm/filter/entity_extraction/ThresholdPoints.h>

int
TestTipsyVTKmConvert(int argc, char* argv[])
{
  vtkm::cont::DataSet dataSet;

  std::cout << "VTK-m::Initialize" << std::endl;

  vtkm::cont::Initialize(argc, argv);

  TipsyFile *filein = new TipsyFile(argv[1]);
  filein->read_all();

  vtkm::cont::DataSetBuilderExplicit dataSetBuilder;
  auto AOS = vtkm::cont::make_ArrayHandle<vtkm::Float32>(filein->gas_ptr(),
                                             filein->h.nsph * filein->stride_of_gas_particle,
                                             vtkm::CopyFlag::Off);
                                             
  vtkm::Id pos_offset = offsetof(struct gas_particle, pos)/sizeof(float);
  vtkm::cont::ArrayHandleStride<vtkm::Float32> pos_x (AOS, filein->h.nsph, filein->stride_of_gas_particle, pos_offset);
  vtkm::cont::ArrayHandleStride<vtkm::Float32> pos_y (AOS, filein->h.nsph, filein->stride_of_gas_particle, pos_offset+1);
  vtkm::cont::ArrayHandleStride<vtkm::Float32> pos_z (AOS, filein->h.nsph, filein->stride_of_gas_particle, pos_offset+2);

  auto coordsArray2 = vtkm::cont::make_ArrayHandleCompositeVector(pos_x, pos_y, pos_z);
  std::cout << "COORDSARRAY2 HANDLE" << std::endl;
  vtkm::cont::printSummary_ArrayHandle(coordsArray2, std::cout);
  std::cout << "--------------------------" << std::endl;
  vtkm::cont::ArrayHandle<vtkm::Vec3f>      positions2;
  vtkm::cont::ArrayCopy(coordsArray2, positions2);
  std::cout << "POSITIONS2 HANDLE" << std::endl;
  vtkm::cont::printSummary_ArrayHandle(positions2, std::cout);
  std::cout << "--------------------------" << std::endl;

  vtkm::cont::ArrayHandle<vtkm::Id> connectivity;
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandleIndex(static_cast<vtkm::Id>(filein->h.nsph)), connectivity);
  
  vtkm::IdComponent numberOfPointsPerCell = 1;
  dataSet = dataSetBuilder.Create(positions2,
                                  vtkm::CellShapeTagVertex(),
                                  numberOfPointsPerCell,
                                  connectivity, "coords");

  vtkm::Bounds bounds1 = vtkm::cont::BoundsCompute(dataSet);
  std::cout << bounds1 << std::endl;
  /*
  vtkm::cont::ArrayHandleStride<vtkm::Float32> aos0(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, mass)/sizeof(float));
  dataSet.AddPointField("mass", aos0);
  */
  vtkm::cont::ArrayHandleStride<vtkm::Float32> aos7(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, rho)/sizeof(float));
  dataSet.AddPointField("rho", aos7);
  /*
  vtkm::cont::ArrayHandleStride<vtkm::Float32> aos8(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, temp)/sizeof(float));
  dataSet.AddPointField("temp", aos8);
  
  vtkm::cont::ArrayHandleStride<vtkm::Float32> aos9(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, hsmooth)/sizeof(float));
  dataSet.AddPointField("hsmooth", aos9);
  */
  /*
  vtkm::Id vel_offset = offsetof(struct gas_particle, vel)/sizeof(float);
  vtkm::cont::ArrayHandleStride<vtkm::Float32> vx (AOS, filein->h.nsph, filein->stride_of_gas_particle, vel_offset);
  vtkm::cont::ArrayHandleStride<vtkm::Float32> vy (AOS, filein->h.nsph, filein->stride_of_gas_particle, vel_offset+1);
  vtkm::cont::ArrayHandleStride<vtkm::Float32> vz (AOS, filein->h.nsph, filein->stride_of_gas_particle, vel_offset+2);
  auto velocity = vtkm::cont::make_ArrayHandleCompositeVector(vx, vy, vz);
  dataSet.AddPointField("velocity", velocity);
  */
  std::cout << "TIPSY dataSet summary--------------------------" << std::endl;
  dataSet.PrintSummary(std::cout);
  std::cout << "--------------------------" << std::endl;

  // Get the overall min/max of a field named "rho"
  vtkm::cont::ArrayHandle<vtkm::Range> rho = vtkm::cont::FieldRangeCompute(dataSet, "rho");
  std::cout << "range(rho) = " << rho.ReadPortal().Get(0) << std::endl;
  //vtkm::cont::ArrayHandle<vtkm::Range> temp = vtkm::cont::FieldRangeCompute(dataSet, "temp");
  //std::cout << "range(temp) = " << temp.ReadPortal().Get(0) << std::endl;
  
  /* dataset is fully build. Can now save it to disk */
  //vtkm::io::VTKDataSetWriter writer(argv[2]);
  //writer.SetFileTypeToBinary();
  //writer.WriteDataSet(dataSet);

  vtkm::cont::Timer timer;
#define TestingThresholdPoints 1
#ifdef TestingThresholdPoints
  timer.Start();
  vtkm::filter::entity_extraction::ThresholdPoints thresholdPoints;
  thresholdPoints.SetThresholdBetween(7.0f, 1000.0f);
  thresholdPoints.SetActiveField("rho");
  thresholdPoints.SetFieldsToPass("rho");
  thresholdPoints.SetCompactPoints(true);
  vtkm::cont::DataSet outputData = thresholdPoints.Execute(dataSet);
  std::cout << "thresholdPoints.Execute : " << timer.GetElapsedTime() << " seconds"<< std::endl;
  //vtkm::io::VTKDataSetWriter cspWriter("/dev/shm/thresholdPoints.vtk");
  //cspWriter.SetFileTypeToBinary();
  //cspWriter.WriteDataSet(outputData);
#endif

//#define TestingExtractPoints 1
#ifdef TestingExtractPoints
  /********** ndHistFilter *****************/
  
     // Implicit function
  vtkm::Vec3f minPoint(-2.6f, -2.9f, -1.88f);
  vtkm::Vec3f maxPoint(4.09f, 1.95f, 2.08f);
  vtkm::Box box(minPoint, maxPoint);
    
  // Setup and run filter to extract by volume of interest
  vtkm::filter::entity_extraction::ExtractPoints extractPoints;
  extractPoints.SetImplicitFunction(box);
  extractPoints.SetExtractInside(true);
  extractPoints.SetCompactPoints(true);
  vtkm::cont::DataSet outputData = extractPoints.Execute(dataSet);
  
  vtkm::io::VTKDataSetWriter cspWriter("/dev/shm/ExtractPoints.vtk");
  cspWriter.SetFileTypeToBinary();
  cspWriter.WriteDataSet(outputData);
#endif

//#define NDHISTOGRAM 1
#ifdef NDHISTOGRAM
  /********** ndHistFilter *****************/
  vtkm::filter::density_estimate::NDHistogram ndHistFilter;
  int nx = 128;
  int ny = 128;
  ndHistFilter.AddFieldAndBin("rho", nx);
  ndHistFilter.AddFieldAndBin("temp", ny);
  vtkm::cont::DataSet outputData = ndHistFilter.Execute(dataSet);
  // before writing, we must add a coordsystem
  /*
  terminate called after throwing an instance of 'vtkm::cont::ErrorBadValue'
  what():  DataSet has no coordinate system, which is not supported by VTK file format.
  */
  vtkm::Id3 dimensions(nx, ny,1);
  vtkm::Vec3f origin(0., 0., 0.);
  float spacing = 1.0/(nx-1.0);
  vtkm::Vec3f Spacing(spacing, spacing, spacing);
  vtkm::cont::ArrayHandleUniformPointCoordinates coords(dimensions, origin, Spacing);
  vtkm::cont::CoordinateSystem cs("coords", coords);
  outputData.AddCoordinateSystem(cs);
  std::cout << "dataSet summary--------------------------" << std::endl;
  outputData.PrintSummary(std::cout);
  std::cout << "--------------------------" << std::endl;
  
  vtkm::io::VTKDataSetWriter cspWriter("/dev/shm/ndHistFilter.vtk");
  cspWriter.SetFileTypeToBinary();
  cspWriter.WriteDataSet(outputData);
#endif

//#define HISTSAMPLING 1
#ifdef HISTSAMPLING
  /********** HistSampling *****************/
  using AssocType = vtkm::cont::Field::Association;
  vtkm::filter::resampling::HistSampling histsample;
  histsample.SetNumberOfBins(128);
  histsample.SetSamplePercent(0.1);
  histsample.SetActiveField("rho", AssocType::Points);
  auto histsampleDataSet = histsample.Execute(dataSet);
  
  vtkm::io::VTKDataSetWriter histsampleWriter("/dev/shm/histsample.vtk");
  histsampleWriter.SetFileTypeToBinary();
  histsampleWriter.WriteDataSet(histsampleDataSet);
#endif

//#define PARTICLEDENSITY 1
#ifdef PARTICLEDENSITY
  /********** ParticleDensityCloudInCell *****************/
  vtkm::Id3 cellDims = { 512,512,512 };

  vtkm::filter::density_estimate::ParticleDensityCloudInCell cic;
  cic.SetDimension(cellDims);
  vtkm::Bounds bounds0(vtkm::Range(-200, 180), vtkm::Range(-266, 230), vtkm::Range(-150, 150));
  cic.SetBounds(bounds0);
  cic.SetActiveField("rho");
  auto density = cic.Execute(dataSet);
  
  vtkm::io::VTKDataSetWriter cicWriter("/dev/shm/ParticleDensityCloudInCell.vtk");
  cicWriter.SetFileTypeToBinary();
  cicWriter.WriteDataSet(density);
#endif

//#define VERTEXCLUSTERING 1
#ifdef VERTEXCLUSTERING
  vtkm::filter::geometry_refinement::VertexClustering vertexClustering;
  vertexClustering.SetNumberOfDivisions(vtkm::Id3(128, 128, 128));

  auto simplifiedCloud = vertexClustering.Execute(dataSet);
  vtkm::io::VTKDataSetWriter vcWriter("/dev/shm/vertexClustering.vtk");
  vcWriter.SetFileTypeToBinary();
  vcWriter.WriteDataSet(simplifiedCloud);
#endif

    delete filein;
  return 1;
}

int
main(int argc, char* argv[])
{
  return TestTipsyVTKmConvert(argc, argv);
}

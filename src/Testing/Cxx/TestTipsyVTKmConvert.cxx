#include <stddef.h>
#include <iomanip>
#include "tipsy_file.h"
#include <viskores/cont/Timer.h>

#include <viskores/cont/DataSet.h>
#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/DataSetBuilderExplicit.h>
#include <viskores/cont/ArrayHandleExtractComponent.h>
#include <viskores/cont/ArrayHandleCompositeVector.h>
#include <viskores/cont/ArrayHandleStride.h>
#include <viskores/cont/ArrayHandleIndex.h>
#include <viskores/cont/FieldRangeCompute.h>
#include <viskores/cont/BoundsCompute.h>
#include <viskores/cont/Initialize.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/filter/resampling/HistSampling.h>
#include <viskores/filter/density_estimate/ParticleDensityCloudInCell.h>
#include <viskores/filter/geometry_refinement/VertexClustering.h>
#include <viskores/filter/density_estimate/NDHistogram.h>
#include <viskores/filter/entity_extraction/ExtractPoints.h>
#include <viskores/filter/entity_extraction/ThresholdPoints.h>

//#include <viskores/rendering/Actor.h>
//#include <viskores/rendering/CanvasRayTracer.h>
//#include <viskores/rendering/MapperRayTracer.h>
//#include <viskores/rendering/Scene.h>
//#include <viskores/rendering/View3D.h>

viskores::cont::DataSet *
TipsyToviskoresDataSet(TipsyFile *filein, bool write=false)
{
  viskores::cont::DataSet *dataSet = new viskores::cont::DataSet();

  viskores::cont::DataSetBuilderExplicit dataSetBuilder;
  auto AOS = viskores::cont::make_ArrayHandle<viskores::Float32>(filein->gas_ptr(),
                                             filein->h.nsph * filein->stride_of_gas_particle,
                                             viskores::CopyFlag::Off);
                                             
  viskores::Id pos_offset = offsetof(struct gas_particle, pos)/sizeof(float);
  viskores::cont::ArrayHandleStride<viskores::Float32> pos_x (AOS, filein->h.nsph, filein->stride_of_gas_particle, pos_offset);
  viskores::cont::ArrayHandleStride<viskores::Float32> pos_y (AOS, filein->h.nsph, filein->stride_of_gas_particle, pos_offset+1);
  viskores::cont::ArrayHandleStride<viskores::Float32> pos_z (AOS, filein->h.nsph, filein->stride_of_gas_particle, pos_offset+2);

  auto coordsArray2 = viskores::cont::make_ArrayHandleCompositeVector(pos_x, pos_y, pos_z);
  std::cout << "COORDSARRAY2 HANDLE" << std::endl;
  viskores::cont::printSummary_ArrayHandle(coordsArray2, std::cout);
  std::cout << "--------------------------" << std::endl;
  viskores::cont::ArrayHandle<viskores::Vec3f>      positions2;
  viskores::cont::ArrayCopy(coordsArray2, positions2);
  std::cout << "POSITIONS2 HANDLE" << std::endl;
  viskores::cont::printSummary_ArrayHandle(positions2, std::cout);
  std::cout << "--------------------------" << std::endl;

  viskores::cont::ArrayHandle<viskores::Id> connectivity;
  viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleIndex(static_cast<viskores::Id>(filein->h.nsph)), connectivity);
  
  viskores::IdComponent numberOfPointsPerCell = 1;
  *dataSet = dataSetBuilder.Create(positions2,
                                  viskores::CellShapeTagVertex(),
                                  numberOfPointsPerCell,
                                  connectivity, "coords");

  viskores::Bounds bounds1 = viskores::cont::BoundsCompute(*dataSet);
  std::cout << bounds1 << std::endl;
  /*
  viskores::cont::ArrayHandleStride<viskores::Float32> aos0(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, mass)/sizeof(float));
  dataSet.AddPointField("mass", aos0);
  */
  viskores::cont::ArrayHandleStride<viskores::Float32> aos7(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, rho)/sizeof(float));
  dataSet->AddPointField("rho", aos7);

  viskores::cont::ArrayHandleStride<viskores::Float32> aos8(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, temp)/sizeof(float));
  dataSet->AddPointField("temp", aos8);
  /*
  viskores::cont::ArrayHandleStride<viskores::Float32> aos9(AOS, filein->h.nsph,
                                                         filein->stride_of_gas_particle,
                                                         offsetof(struct gas_particle, hsmooth)/sizeof(float));
  dataSet.AddPointField("hsmooth", aos9);
  */
  /*
  viskores::Id vel_offset = offsetof(struct gas_particle, vel)/sizeof(float);
  viskores::cont::ArrayHandleStride<viskores::Float32> vx (AOS, filein->h.nsph, filein->stride_of_gas_particle, vel_offset);
  viskores::cont::ArrayHandleStride<viskores::Float32> vy (AOS, filein->h.nsph, filein->stride_of_gas_particle, vel_offset+1);
  viskores::cont::ArrayHandleStride<viskores::Float32> vz (AOS, filein->h.nsph, filein->stride_of_gas_particle, vel_offset+2);
  auto velocity = viskores::cont::make_ArrayHandleCompositeVector(vx, vy, vz);
  dataSet.AddPointField("velocity", velocity);
  */
  std::cout << "TIPSY dataSet summary--------------------------" << std::endl;
  dataSet->PrintSummary(std::cout);
  std::cout << "--------------------------" << std::endl;

  // Get the overall min/max of a field named "rho"
  viskores::cont::ArrayHandle<viskores::Range> rho = viskores::cont::FieldRangeCompute(*dataSet, "rho");
  std::cout << "range(rho) = " << rho.ReadPortal().Get(0) << std::endl << std::endl;
  //viskores::cont::ArrayHandle<viskores::Range> temp = viskores::cont::FieldRangeCompute(dataSet, "temp");
  //std::cout << "range(temp) = " << temp.ReadPortal().Get(0) << std::endl;
  
  if(write){
    std::string fname = "/dev/shm/dataSet.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    writer.SetFileTypeToBinary();
    writer.WriteDataSet(*dataSet);
  }
  return dataSet;
}
  
void TestingThresholdPoints(const viskores::cont::DataSet &dataSet, bool write=false)
{
  viskores::cont::Timer timer;
  timer.Start();
  viskores::filter::entity_extraction::ThresholdPoints thresholdPoints;
  // this threshold range for /local/data/Tipsy/hr8799_bol_bd1.017300
  thresholdPoints.SetThresholdBetween(1.15573e-7, 0.517376);
  thresholdPoints.SetActiveField("rho");
  thresholdPoints.SetFieldsToPass("rho");
  thresholdPoints.SetCompactPoints(true);
  auto output = thresholdPoints.Execute(dataSet);
  std::cout << "thresholdPoints.Execute :            " << timer.GetElapsedTime() << " seconds"<< std::endl;
  
  if(write) {
    std::string fname = "/dev/shm/thresholdPoints.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    writer.SetFileTypeToBinary();
    writer.WriteDataSet(output);
  }
}

void TestingExtractPoints(const viskores::cont::DataSet dataSet, bool write=false)
{
  viskores::cont::Timer timer;
  timer.Start();
  viskores::Vec3f minPoint(-2.6f, -2.9f, -1.88f);
  viskores::Vec3f maxPoint(4.09f, 1.95f, 2.08f);
  viskores::Box box(minPoint, maxPoint);
    
  // Setup and run filter to extract by volume of interest
  viskores::filter::entity_extraction::ExtractPoints extractPoints;
  extractPoints.SetImplicitFunction(box);
  extractPoints.SetExtractInside(true);
  extractPoints.SetCompactPoints(true);
  auto output = extractPoints.Execute(dataSet);
  std::cout << "ExtractPoints.Execute :              " << timer.GetElapsedTime() << " seconds"<< std::endl;
  
  if(write) {
    std::string fname = "/dev/shm/ExtractPoints.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    writer.SetFileTypeToBinary();
    writer.WriteDataSet(output);
  }
}

void TestingNDHistogram(const viskores::cont::DataSet &dataSet, bool write=false)
{
  viskores::cont::Timer timer;
  timer.Start();
  viskores::filter::density_estimate::NDHistogram ndHistFilter;
  int nx = 128;
  int ny = 128;
  ndHistFilter.AddFieldAndBin("rho", nx);
  ndHistFilter.AddFieldAndBin("temp", ny);
  auto output = ndHistFilter.Execute(dataSet);
  std::cout << "NDHistogram.Execute :                " << timer.GetElapsedTime() << " seconds"<< std::endl;
  std::cout << "dataSet summary--------------------------" << std::endl;
  output.PrintSummary(std::cout);
  std::cout << "--------------------------" << std::endl<< std::endl;
  /* does not make sense right now
  if(write) {
  // before writing, we must add a coordsystem
  //terminate called after throwing an instance of 'viskores::cont::ErrorBadValue'
  //what():  DataSet has no coordinate system, which is not supported by VTK file format.

  viskores::Id3 dimensions(nx, ny,1);
  viskores::Vec3f origin(0., 0., 0.);
  float spacing = 1.0/(nx-1.0);
  viskores::Vec3f Spacing(spacing, spacing, spacing);
  viskores::cont::ArrayHandleUniformPointCoordinates coords(dimensions, origin, Spacing);
  viskores::cont::CoordinateSystem cs("coords", coords);
  output.AddCoordinateSystem(cs);
  
  viskores::cont::CellSetStructured<2> cellSet;
  cellSet.SetPointDimensions(viskores::Id2(nx, ny));
  output.SetCellSet(cellSet);
  
  std::cout << "dataSet summary--------------------------" << std::endl;
  output.PrintSummary(std::cout);
  std::cout << "--------------------------" << std::endl;

    char *fname = "/dev/shm/ndHistFilter.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    //writer.SetFileTypeToBinary();
    writer.WriteDataSet(output);
  }*/
}

void TestingHistSampling(const viskores::cont::DataSet &dataSet, bool write=false)
{
  viskores::cont::Timer timer;
  timer.Start();
  using AssocType = viskores::cont::Field::Association;
  viskores::filter::resampling::HistSampling histsample;
  histsample.SetNumberOfBins(128);
  histsample.SetSampleFraction(0.1);
  histsample.SetActiveField("rho", AssocType::Points);
  auto output = histsample.Execute(dataSet);
  std::cout << "HistSampling.Execute :               " << timer.GetElapsedTime() << " seconds"<< std::endl;
  
  if(write) {
    std::string fname = "/dev/shm/histsample.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    writer.SetFileTypeToBinary();
    writer.WriteDataSet(output);
  }
}

void TestingParticleDensityCloudInCell(const viskores::cont::DataSet &dataSet, bool write=false, bool render=false)
{
  viskores::cont::Timer timer;
  timer.Start();
  viskores::Id3 cellDims = { 511,511,511 };

  viskores::filter::density_estimate::ParticleDensityCloudInCell cic;
  cic.SetDimension(cellDims);
  viskores::Bounds bounds0(viskores::Range(-200, 180), viskores::Range(-266, 230), viskores::Range(-150, 150));
  cic.SetBounds(bounds0);
  cic.SetActiveField("rho");
  cic.SetDivideByVolume(true);
  cic.SetComputeNumberDensity(true);
  auto output = cic.Execute(dataSet);
  std::cout << "ParticleDensityCloudInCell.Execute : " << timer.GetElapsedTime() << " seconds"<< std::endl;

  if(write) {
    std::string fname = "/dev/shm/ParticleDensityCloudInCell.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    writer.SetFileTypeToBinary();
    writer.WriteDataSet(output);
  }
  /*
  if(render) {
    std::string fname = "/dev/shm/ParticleDensityCloudInCell.png";
    std::cout << "Rendering: " << fname << std::endl;
  //Creating Actor
  viskores::cont::ColorTable colorTable("viridis");
  viskores::rendering::Actor actor(output.GetCellSet(),
                               output.GetCoordinateSystem(),
                               output.GetField("density"),
                               colorTable);

  //Creating Scene and adding Actor
  viskores::rendering::Scene scene;
  scene.AddActor(std::move(actor));

  //Creating and initializing the View using the Canvas, Ray Tracer Mappers, and Scene
  viskores::rendering::MapperRayTracer mapper;
  viskores::rendering::CanvasRayTracer canvas(1080, 1080);
  viskores::rendering::View3D view(scene, mapper, canvas);

  //Setting the background and foreground colors; optional.
  view.SetBackgroundColor(viskores::rendering::Color(1.0f, 1.0f, 1.0f));
  view.SetForegroundColor(viskores::rendering::Color(0.0f, 0.0f, 0.0f));

  //Painting View
  view.Paint();

  //Saving View
  view.SaveAs(fname);
  }
  */
}

void TestingVertexClustering(const viskores::cont::DataSet &dataSet, bool write=false)
{
  viskores::cont::Timer timer;
  timer.Start();
  viskores::filter::geometry_refinement::VertexClustering vertexClustering;
  vertexClustering.SetNumberOfDivisions(viskores::Id3(128, 128, 128));

  auto output = vertexClustering.Execute(dataSet);
  std::cout << "VertexClustering.Execute :           " << timer.GetElapsedTime() << " seconds"<< std::endl;

  if(write) {
    std::string fname = "/dev/shm/vertexClustering.vtk";
    std::cout << "Writing: " << fname << std::endl;
    viskores::io::VTKDataSetWriter writer(fname);
    writer.SetFileTypeToBinary();
    writer.WriteDataSet(output);
  }
}

int
main(int argc, char* argv[])
{
  std::cout << "VTK-m::Initialize" << std::endl;
  viskores::cont::Initialize(argc, argv);

  TipsyFile *filein = new TipsyFile(argv[1]);
  filein->read_all();
  
  viskores::cont::DataSet *dataSet = TipsyToviskoresDataSet(filein, false);
  
  //TestingThresholdPoints(*dataSet, true);
  //TestingExtractPoints(*dataSet, false);
  //TestingNDHistogram(*dataSet, true);
  //TestingHistSampling(*dataSet, false);
  TestingParticleDensityCloudInCell(*dataSet, true, true);
  //TestingVertexClustering(*dataSet, false);
  std::cout << std::endl;

  delete dataSet;
  delete filein;

  return 1;
}

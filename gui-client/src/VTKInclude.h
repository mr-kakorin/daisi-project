#ifndef VTKINCLUDE_H
#define VTKINCLUDE_H

#ifdef WIN32
#define vtkRenderingCore_AUTOINIT                                                                                      \
    4(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingFreeTypeOpenGL, vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#endif

#include "vtkChartXYZ.h"
#include "vtkContext2D.h"
#include "vtkContextActor.h"
#include "vtkContextItem.h"
#include "vtkContextScene.h"
#include "vtkNew.h"
#include "vtkPicker.h"
#include "vtkPickingManager.h"
#include "vtkUnsignedCharArray.h"
#include "vtkVersion.h"
#include <QVTKWidget.h>
#include <vtkAbstractArray.h>
#include <vtkActor.h>
#include <vtkAnnotationLink.h>
#include <vtkAreaPicker.h>
#include <vtkAxis.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCellPicker.h>
#include <vtkCellTypes.h>
#include <vtkChartXY.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkCubeAxesActor.h>
#include <vtkCubeAxesActor2D.h>
#include <vtkDataArray.h>
#include <vtkDataObjectToTable.h>
#include <vtkDataRepresentation.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDoubleArray.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkExtractSelection.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkGL2PSExporter.h>
#include <vtkGraphLayoutView.h>
#include <vtkIOExportModule.h>
#include <vtkIdFilter.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkImageViewer.h>
#include <vtkInteractorStyleImage.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkInteractorStyleRubberBandZoom.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkJPEGReader.h>
#include <vtkLine.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkPen.h>
#include <vtkPlaneSource.h>
#include <vtkPlot.h>
#include <vtkPlotPoints.h>
#include <vtkPlotSurface.h>
#include <vtkPointData.h>
#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProp3DCollection.h>
#include <vtkPropPicker.h>
#include <vtkProperty.h>
#include <vtkQtTableView.h>
#include <vtkQuad.h>
#include <vtkRandomGraphSource.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSTLReader.h>
#include <vtkScalarBarActor.h>
#include <vtkSelectionNode.h> // for POINT and INDICES enum values
#include <vtkSelectionSource.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStringArray.h>
#include <vtkTable.h>
#include <vtkTextProperty.h>
#include <vtkTimerLog.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkViewUpdater.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkVersion.h>

#include "vtkPickingManager.h"
#include <vtkActor.h>
#include <vtkAreaPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkIdFilter.h>
#include <vtkIdTypeArray.h>
#include <vtkInteractorStyleRubberBandPick.h>
#include <vtkInteractorStyleTrackball.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkVertexGlyphFilter.h>
//#include <vtkTesting.h>
#include <vtkCommand.h>
#include <vtkLookupTable.h>
#include <vtkMutexLock.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#endif // VTKINCLUDE_H

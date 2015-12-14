#ifndef __EntryCluster_h
#define __EntryCluster_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared
#include "vtkIdList.h"
#include "vtkLookUpTable.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkExecutive.h"
#include "BaseCluster.h"
#include "HierarchicalCluster.h"
#include "KmeansCluster.h"
#include "SpectralCluster.h"
#include "vtkMath.h"
using namespace std;

class vtkPolyData;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;
class KmeansCluster;

class EntryCluster : public vtkPolyDataAlgorithm
{
 public:
	static EntryCluster *New();
	vtkTypeMacro(EntryCluster, vtkPolyDataAlgorithm);
	void PrintSelf(ostream &os, vtkIndent indent);
	//
	vtkSetMacro(NumOfLevel,int);
	vtkGetMacro(NumOfLevel,int);
	vtkSetMacro(BestClusterNum,int);
	vtkGetMacro(BestClusterNum,int);
	vtkGetMacro(ShowBestClusts,bool);
	vtkSetMacro(ShowBestClusts,bool);
	/**/
	vtkGetMacro(ExportClassifyInfo,bool);
	vtkSetMacro(ExportClassifyInfo,bool);
	/**/
	vtkSetVector3Macro(MeshSize,int);
	vtkGetVectorMacro(MeshSize,int,3);
	vtkSetMacro(CalType,int);
	vtkGetMacro(CalType,int);
	vtkSetMacro(ClusterType,int);
	vtkGetMacro(ClusterType,int);
	
	vtkSetMacro(ShowType,int);
	vtkGetMacro(ShowType,int);
	vtkSetMacro(MinNum,int);              
	vtkGetMacro(MinNum,int);
	vtkSetMacro(MaxNum,int);
	vtkGetMacro(MaxNum,int);

	//
protected:
	EntryCluster();
  	~EntryCluster();
  	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); 
	//
private:
	int MeshSize[3];
	int UsedMesh;
	int ClusterType;
	int CalType;
	int ShowType;
	int NumOfLevel;
	int BestClusterNum;
	int MaxNum;
	int MinNum;
	bool ShowBestClusts;
	bool ExportClassifyInfo;
	//int usedClusterType;
	unsigned int InitPoly;
	unsigned int InitType;
	BaseCluster *ClusterShow;
	SpectralCluster *SpecCluster;
	KmeansCluster *KMCluster;
	HierarchicalCluster *HierCluster;
	//
};
#endif


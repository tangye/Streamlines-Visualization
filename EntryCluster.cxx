#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPolyData.h"
#include "vtkIterativeClosestPointTransform.h"
#include "vtkLandmarkTransform.h"
#include "EntryCluster.h"
class KmeansCluster;
vtkStandardNewMacro(EntryCluster);

//-----------------------------------------------------------------------------0
EntryCluster::EntryCluster()
{
	ClusterShow=NULL;
	this->HierCluster=new HierarchicalCluster();;
	this->SpecCluster=new SpectralCluster();
	this->KMCluster=new KmeansCluster();
	//
	InitPoly=1;
	this->ShowBestClusts=false;
	//usedClusterType=-1;
}


//-----------------------------------------------------------------------------
EntryCluster::~EntryCluster()
{ 

}


//----------------------------------------------------------------------------
int EntryCluster::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector)
{
	//vtkOutputWindowDisplayWarningText("EntryCluster Begin-----------------------");
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);	
		// Check the size of the input.
		vtkIdType numPts = input->GetNumberOfPoints();
		if(numPts < 1)
		{
			vtkDebugMacro("No input!");
			return 1;
		}
		int numcell=input->GetNumberOfCells();
		bool tmp=false;
		switch(ClusterType)
		{
		case 1:
			tmp=true;
		case 2:
			if(tmp)
				this->HierCluster->SetHierType(1);
			else
				this->HierCluster->SetHierType(2);
			this->ClusterShow=this->HierCluster;
			break;
		case 3:
			this->ClusterShow=this->SpecCluster;
			break;
		case 4:
			this->ClusterShow=this->KMCluster;
			break;
		}
		//this->ClusterShow=this->KMCluster;

		if(InitPoly==1)
		{
			ClusterShow->InitPoly(input);
			InitPoly=100;
		}	
		ClusterShow->SetCurNum(NumOfLevel);
		//if((UsedMesh!=MeshSize[0]*MeshSize[1]*MeshSize[2])||(this->CalType!=(ClusterShow->GetCalType())))
		//{	
			//if(InitPoly>100)
			//	ClusterShow->Release();
		if((UsedMesh!=MeshSize[0]*MeshSize[1]*MeshSize[2])||(this->CalType!=(ClusterShow->GetCalType())))
		{
			ClusterShow->SetCalType(this->CalType);
			UsedMesh=MeshSize[0]*MeshSize[1]*MeshSize[2];
			ClusterShow->SetMeshSize(MeshSize);
				ClusterShow->CalDistInMesh();
			ClusterShow->CalAffiMatrix();
			//ClusterShow->CopyParaFromFile();
			ClusterShow->BestGroups();

			ClusterShow->AdjustCurClusts(2);

		}
			//this->usedClusterType=this->ClusterType;
		//}

		if(this->ShowBestClusts)
		{
			ClusterShow->ShowBestClusts();
			this->SetNumOfLevel(ClusterShow->GetCurNum());
		}
		else
		{
			ClusterShow->Main();
			ClusterShow->AdjustCurClusts(1);
		}
		if(MinNum>MaxNum)
		{
			int tmp=MinNum;
			MinNum=MaxNum;
			MaxNum=tmp;
		}
		if(MinNum>this->NumOfLevel)
		{
			MinNum=1;
			MaxNum=NumOfLevel;
		}
		if(MaxNum>this->NumOfLevel)
		{
				MaxNum=NumOfLevel;
		}
		switch(ShowType)
		{
			case 1:
				ClusterShow->GetCluster(output);
				
				//ClusterShow->GetSpectralLevel((unsigned int )NumOfLevel,output);
				//ClusterShow->GetKmeansLevel((unsigned int )NumOfLevel,output);
				break;
			case 2:
				ClusterShow->GetClusterOnly(MinNum,MaxNum,output);
				//ClusterShow->GetKmeansLevel((unsigned int )NumOfLevel,output);
				//ClusterShow->GetSpectralLevel((unsigned int )NumOfLevel,output);
				break;
			case 3:
				ClusterShow->GetCenterOnly(MinNum,MaxNum,output);
				break;
			case 4:
				ClusterShow->GetBundleOnly(MinNum,MaxNum,output);
				break;
			case 5:
				ClusterShow->GetStreamtapesOnly(MinNum,MaxNum,output);
				break;
			default:
				ClusterShow->GetCluster(output);
				break;
		}
		if(this->ExportClassifyInfo)
		{
			ClusterShow->ExportClassifyInfo(output);
		}
		//InitPoly++;
  return 1;
}

////////// External Operators /////////////

void EntryCluster::PrintSelf(ostream &os, vtkIndent indent)
{
}

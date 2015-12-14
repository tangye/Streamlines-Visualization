#include "SpectralCluster.h"
#include "KmeansCluster.h"
/**/
/**/
//template <class T> 
void MainKmeansCluster(const double * mat,int* cluster,int row,int k)
{
	int tmpi=0,tmpk=0,tmpj=0;
	int vari,varj,vark;
	const double* tmpMat;
	//const T* tmpMat;

	void* tmppointer;
	int* curCenter;
	int* curOldCenter;
	int* Center=new int[k];
	int* OldCenter=new int[k];
	int* LastPos=new int[k];
	int* ClusterIndex=new int[row];
	double MaxCondition;
	int time=10000;
	//the data distance to the cluster center
	//int * CenterDist=new int[row*k];
	//
	tmpMat=mat;
	tmppointer=NULL;
	//tmpvoid=NULL;
	curCenter=Center;
	curOldCenter=OldCenter;
	tmpj=0;
	for(tmpi=0,tmpk=row/(k);tmpi<k;tmpi++)
	{
		curCenter[tmpi]=tmpj;
		curOldCenter[tmpi]=tmpj;
		tmpj+=tmpk;
	}
	MaxCondition=100000;
	while(MaxCondition>_MaxDouble&&time>0)
	{
		/**/
		/**/
		for(tmpk=0;tmpk<k;tmpk++)
			LastPos[tmpk]=-1*(tmpk+1);
		/**/
		/**/
		for(tmpi=0;tmpi<row;tmpi++)
		{
			double minValue=MaxTValue;
			vari=tmpi*row;
			for(tmpk=0;tmpk<k;tmpk++)
				if(tmpMat[vari+curCenter[tmpk]]<minValue)
				{
					minValue=tmpMat[vari+curCenter[tmpk]];
					vark=tmpk;
				}
			ClusterIndex[tmpi]=LastPos[vark];
			LastPos[vark]=tmpi;
		}
		/**/
		/**/
		for(tmpk=0;tmpk<k;tmpk++)
		{
			vari=LastPos[tmpk];
			double sum=0;
			double min=MaxTValue;
			int idx=0;
			while(vari>=0)
			{
				sum=0;
				varj=LastPos[tmpk];
				while(varj>=0)
				{
					sum+=tmpMat[vari*row+varj];
					varj=ClusterIndex[varj];
				}
				if(sum<min)
				{
					min=sum;
					idx=vari;
				}
				vari=ClusterIndex[vari];
			}
			curCenter[tmpk]=idx;
		}
		/**/
		/**/
		MaxCondition=0.0;
		for(tmpk=0;tmpk<k;tmpk++)
		{
			//MaxCondition+=(curInterDist[tmpk]-curOldInterDist[tmpk])*(curInterDist[tmpk]-curOldInterDist[tmpk]);
			MaxCondition+=tmpMat[curCenter[tmpk]*row+curOldCenter[tmpk]];
		}
		/**/
		/**/
		tmppointer=curCenter;
		curCenter=curOldCenter;
		curOldCenter=(int*)tmppointer;
		/**/
		/**/
		time--;
	}
	/**/
	/**/
	for(tmpk=0;tmpk<k;tmpk++)
	{
			vari=LastPos[tmpk];
			while(vari>=0)
			{
				vark=ClusterIndex[vari];
				ClusterIndex[vari]=tmpk;
				vari=vark;
			}
	}
	for(tmpi=0;tmpi<row;tmpi++)
		cluster[tmpi]=ClusterIndex[tmpi];
	/**/
	/**/
	ofstream hFile("Kmeanscluster.index_3", ios::out, _SH_DENYRW); 
	for(tmpi=0;tmpi<row;tmpi++)
	{
		hFile<<ClusterIndex[tmpi]<<endl;
	}
	hFile.close();
	delete []Center;
	delete []OldCenter;
	delete []LastPos;
	delete []ClusterIndex;
	/**/
	/**/
};
/**/
/**/
KmeansCluster::KmeansCluster(void)
{
	this->caltype=-1;
};
KmeansCluster::~KmeansCluster(void)
{

};
/////**/
////void KmeansCluster::GetCluster(int levelnum,vtkPolyData *output)
////{
////	vtkPolyData *newclu=vtkPolyData::New();
////	vtkIntArray *newarray=vtkIntArray::New();
////	vtkCellArray *newlines=vtkCellArray::New();
////	vtkPoints *newpts=vtkPoints::New();
////
////	newpointnum=0;
////	newcellnum=0;
////	newclu->Allocate();
////	newclu->Initialize();
////	newarray->Initialize();
////	newarray->SetName("Cluster_kmeans");
////	int* idx=new int[this->cellnum];
////	/*for(int i=0;i<this->cellnum;i++)
////	{
////		for(int j=i+1;j<this->cellnum;j++)
////		{
////			tmpdou=this->alldist[i*this->cellnum+j];
////			if(tmpdou>threshold)
////			{
////					tmpDist[i*this->cellnum+j]=0.0;
////					tmpDist[j*this->cellnum+i]=0.0;
////			}
////			else
////			{
////				tmpDist[i*this->cellnum+j]=tmpdou;
////				tmpDist[j*this->cellnum+i]=tmpdou;
////			}
////		}
////	}
////	for(int i=0;i<this->cellnum;i++)
////		tmpDist[i*this->cellnum+i]=0.0;
////	/*	*/
////	//tmplist=level_entry->get(levelnum);
////	//clushow=list_entry(tmplist,cluster,embed);
////	//unsigned int leflag=level_entry->flag(levelnum);
////	//CopMatrixToFile(tmpDist,this->cellnum,this->cellnum,"test.my");
////	//NormalizedSpectralCluster(tmpDist,idx,this->cellnum,levelnum);
////
////	//ivar=1;
////	/*
////	while(level_entry->current<level_entry->size-1)
////	{
////		clushow=list_entry(tmplist,cluster,embed);
////		if(leflag>clushow->GetIndex())
////		{	
////				tmpidx=clushow->idxhead;
////				do
////				{
////				idxshow=list_entry((tmpidx),idxlist,embed);
////				CopyLine(idxshow->GetIndex(),newpts,newlines);
////				newarray->InsertNextTuple1(ivar);
////				tmpidx=tmpidx->right;
////				}while(tmpidx!=((clushow->idxtail)->right));
////			ivar++;
////		}
////		tmplist=*(level_entry->entry+level_entry->current);
////		level_entry->current++;
////	}
////	*/
////	for (int i=0;i<this->cellnum;i++)
////	{
////		CopyLine(i,newpts,newlines);
////		newarray->InsertNextTuple1(idx[i]);
////	}
////	newclu->SetPoints(newpts);
////	newclu->SetLines(newlines);
////	newclu->GetCellData()->SetScalars(newarray);
////	newclu->GetFieldData()->CopyAllOn();
////	output->DeepCopy(newclu);
////	newclu->Delete();
////	newarray->Delete();
////	newlines->Delete();
////	newpts->Delete();
////	delete []idx;
////	/*tmplist=level_entry->get(num);
////	clushow=list_entry(tmplist,cluster,embed);
////	unsigned int leflag=level_entry->flag(num);
////	ivar=0;
////	while(level_entry->current<level_entry->size-1)
////	{
////		clushow=list_entry(tmplist,cluster,embed);
////		tmpidx=clushow->idxhead;
////		if(leflag>clushow->GetIndex())
////		{
////			do
////			{
////			idxshow=list_entry((tmpidx),idxlist,embed);
////			justarray->SetTuple1(idxshow->GetIndex(),ivar);
////			tmpidx=tmpidx->right;
////			}while(tmpidx!=((clushow->idxtail)->right));	
////			ivar++;
////		}
////		tmplist=*(level_entry->entry+level_entry->current);
////		level_entry->current++;
////	}
////	lines->GetCellData()->SetScalars(justarray);
////	lines->GetFieldData()->CopyAllOn();
////	output->DeepCopy(lines);*/
////};
////void KmeansCluster::GetClusterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
////{
////
////};
////void KmeansCluster::GetCenterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
////{
////
////};
////void KmeansCluster::GetBundleOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
////{
////
////};
void KmeansCluster::Main(void)
{
		clock_t start,end; 
		double duration;

		start=clock();
		MainKmeansCluster((this->alldist),(this->cur_clusts),(int)this->cellnum,(this->cur_num));

		//
		end=clock();
		duration = (double)(end - start) / CLOCKS_PER_SEC; 
		writePara("Kmeans Cluster - all - num - time  ",(double)this->cellnum,(double)this->cur_num,duration,0.0);

};
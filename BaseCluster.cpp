
#include <vtkDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkStreamTracer.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkDataSetAttributes.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <Windows.h>
#include <math.h>
#include "HierarchicalCluster.h"
#include "SpectralCluster.h"
#include "KmeansCluster.h"
#include "BaseCluster.h"
#include "SelfTuningCluster.h"
#include "vtkFloatArray.h"
/**/
extern "C"
{
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"
int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *a,
	integer *lda, doublereal *wr, doublereal *wi, doublereal *vl,
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work,
	integer *lwork, integer *info);
};



/**/
using namespace std;
/**/
/**/
/**/
/**
points::points(int num)
{
	size=num;
	coor=new double[3*num];
};
/**/
/**
points::~points()
{
	delete [] coor;	
};
/**/
/**
void points::getpoint(int id,double *pt)
{
	*pt=*((this->coor)+3*id);
	*(pt+1)=*((this->coor)+3*id+1);
	*(pt+2)=*((this->coor)+3*id+2);
};
/**/
/**/
idxlist::idxlist(unsigned int num)
{
	this->size=num;
	this->idx=new double[num];
	this->length=VTK_DOUBLE_MAX;
};
/**/
/**/
idxlist::~idxlist()
{
	delete []idx;
};
/**/
/**/
double idxlist::GetLength(void)
{
	return this->length;
}
/**/
/**/
unsigned int idxlist::GetSize(void)
{
	return (this->size);
};
/**/
/**/
unsigned int idxlist::GetPtnum(void)
{
	return (this->size)/3;
};
/**/
/**/
void idxlist::SetSize(int num)
{
	this->size=num;
};
/**/
/**/
int idxlist::GetIndex()
{
	return (this->index);
};
/**/
/**/
void idxlist::SetIndex(int num)
{
	this->index=num;
};
/**/
/**/
double idxlist::GetValue(int valuenum)
{
	return *(this->idx+valuenum);
};
/**/
/**/
void idxlist::SetValue(int valuenum,double newvalue)
{
	*(this->idx+valuenum)=newvalue;
};
/**/	
/**/
void idxlist::GetPoint(unsigned int ptnum,double *ret)
{
	for(int i=0;i<3;i++)
		ret[i]=*(this->idx+ptnum*3+i);
		/**/ret[1]=*(this->idx+ptnum*3+1);
		/**/ret[2]=*(this->idx+ptnum*3+2);//?

};
/**/
/**/
void idxlist::SetPoint(unsigned int ptnum,double *newpt)
{
	for(int i=0;i<3;i++)
		*(this->idx+ptnum*3+i)=(newpt[i]);
};
/**/
/**/
void idxlist::CalLength(void)
{
	double fir[3];
	double sec[3];
	double *th;
	double *ot;
	double *var;
	double tmplength;
	double tmp;
	unsigned int i=0;
	tmp=0.0;
	this->GetPoint(0,fir);
	th=fir;
	ot=sec;
	for(i=3;i<(this->GetPtnum());i++)
	{
		this->GetPoint(i,ot);
		dist(th,ot,&tmplength);
		tmp+=tmplength;
		var=th;
		th=ot;
		ot=var;
	}
	this->length=tmp;
};
/**/
/**/
void idxlist::Release(void )
{
	delete []idx;
};
/**/
cluster::cluster(void)
{
	this->centerid=0;
	this->idxhead=NULL;
	this->idxtail=NULL;
	this->index=0;
}
/**/
cluster::~cluster(void)
{
}
/**/
/**/
int cluster::GetIndex(void)
{
	return this->index;
}
/**/
void cluster::SetIndex(int index)
{
	this->index=index;
};
/**/
/**/
unsigned int cluster::GetCenterId()
{
	return (this->centerid);
};
void cluster::SetCenterId(unsigned int it)
{
	(this->centerid)=it;
};
/**/
/**/
void list_add(list *add,list *sibling,list *next)
{
    next->sibling = add;
    add->left = next;
    add->sibling = sibling;
    sibling->right = add;
}
/**/
/**/
void list_del(list *sibling, list *next)
{
    next->sibling = sibling;
    sibling->left = next;
}
/**/
/**/

void list_init(list *head,list *tail)
{
	head->sibling=tail;
	tail->left=head;
	head->left=tail;
	tail->sibling=head;
};
/**/
/**/

level::level(unsigned int num)
{
	entry=new list*[2*num];
	for(unsigned int ivar=0;ivar<2*num;ivar++)
		*(entry+ivar)=NULL;
	size=2*num;
	current=0;
};
/**/
/**/
	/**/
	/**/
void linecontainer::push(list *add)
{
	*(entry+current)=add;//没有边界检查？？？
	current++;
};
/**/
/**/
/**/
/**/
void linecontainer::pop(list *set)
{
	if(current<this->size-1)
	{
	set=*(entry+current);
	cluster *tmp=list_entry(set,cluster,embed);
	current++;
	}
	else
	{
	current=VTK_INT_MAX;
	set=NULL;
	}
};
/**/
list* linecontainer::getvec(unsigned int vecnum)
{
	this->current=vecnum;
	return *(entry+current);
};
/**/
list* level::get(unsigned int level_num)
{
	this->current=size-2*level_num;
	return *(entry+current);
};
/**/
unsigned int  level::flag(unsigned int level_num)
{
	return ((size-level_num)*10-5);
};
/**/
polys::polys(unsigned int num)
{
	entry=new list*[num];
	for(unsigned int ivar=0;ivar<num;ivar++)
		*(entry+ivar)=NULL;
	size=num;
	current=0;
};
/**/
/**/
list* polys::getlist(unsigned int list_num)
{
	return *(entry+list_num);
};
/**
/**/
	vtkPolyData* BaseCluster::lines=NULL;
	//vtk线数据
	 double* BaseCluster::meshsize=NULL;
	//网格大小

	polys* BaseCluster::poly=NULL;
	//几何数据信息
	list* BaseCluster::idxentry=NULL;
	//原始线数据链表
	list* BaseCluster::lidxentry=NULL;
	//网格化线数据链表

	double* BaseCluster::mincoor=NULL;
	//三维最小坐标值
	double* BaseCluster::maxcoor=NULL;
	//三维最大坐标值
	double* BaseCluster::scale1=NULL;
	//坐标转化缩放值
	double* BaseCluster::alldist=NULL;
	double* BaseCluster::affinityegi=NULL;
	double* BaseCluster::affinitymat=NULL;

	//距离矩阵
	unsigned int BaseCluster::cellnum=0;
	//cell number
	unsigned int BaseCluster::newcellnum=0;
	unsigned int BaseCluster::newpointnum=0;
	unsigned int BaseCluster::caltype=1;
	bool BaseCluster::isinit=false;
	int* BaseCluster::best_clusts=NULL;
	unsigned int BaseCluster::best_num=0;
	int* BaseCluster::cur_clusts=NULL;
	int* BaseCluster::clusts_center=NULL;
	unsigned int BaseCluster::cur_num=0;

	//聚类类型

/**
/**/
BaseCluster::BaseCluster(void)
{
		/**
		meshsize=new double[3];
		level_entry=NULL;
		lines=vtkPolyData::New();
		cellnum=1000;
		poly=NULL;
		lentry=new list;
		idxentry=new list;
		lidxentry=new list;
		//idxlist *idxentry=new idxlist();
		alldist=NULL;
		mincoor=new double[3];
		maxcoor=new double[3];
		scale1=new double[3];
		//pts=NULL;
		meshsize[0]=40.0;
		meshsize[1]=40.0;
		meshsize[2]=40.0;
		/**/
}
/**
/**/
BaseCluster::~BaseCluster(void)
{
			

}
/**
/**/
void BaseCluster::SetCalType(unsigned int _type)
{
	this->caltype=_type;
};
//
unsigned int BaseCluster::GetCalType(void)
{
	return (this->caltype);
};
void BaseCluster::SetCurNum(unsigned int _type)
{
	this->cur_num=_type;
};
//
unsigned int BaseCluster::GetCurNum(void)
{
	return (this->cur_num);
};
//
//
void BaseCluster::SetMeshSize(int mesh[3])
{
	this->meshsize[0]=(double)mesh[0];
	this->meshsize[1]=(double)mesh[1];
	this->meshsize[2]=(double)mesh[2];
};
//
int BaseCluster::GetMeshSize(void)
{
	return (int)(this->meshsize[0]);
};
//
void BaseCluster::InitPoly(vtkPolyData *input)
{
	if(!this->isinit)
	{	/**/		
		double doupts2[1024];
		unsigned int itmp,ivar;
		list* tmplist;
		list* idxtmp;
		vtkCellData *vcd=vtkCellData::New();
		vtkCellArray* linecell=vtkCellArray::New();
		lines=vtkPolyData::New();
		linecell=vtkCellArray::New();//?
		mincoor=new double[3];
		maxcoor=new double[3];		
		this->idxentry=new list;

		this->lidxentry=new list;
		this->meshsize=new double[3];
		this->scale1=new double[3];
		meshsize[0]=40.0;
		meshsize[1]=40.0;
		meshsize[2]=40.0;
		lines->GetPointData()->PassData(input->GetPointData());
		lines->GetCellData()->PassData(input->GetCellData());
		lines->DeepCopy(input);
	
		linecell=lines->GetLines();
		linecell->InitTraversal();
		//pts=new points(lines->GetNumberOfPoints());
		this->poly=new polys(lines->GetNumberOfCells());
		double tmppts[3]={0.0,0.0,0.0};
		/**/
		for(int s=0;s<3;s++)
		{
			mincoor[s]=VTK_DOUBLE_MAX;
			maxcoor[s]=VTK_DOUBLE_MIN;
		}
		idxentry->left=NULL;
		idxentry->right=NULL;
		idxtmp=idxentry;
		cellnum=lines->GetNumberOfCells();
		for(int ii=0;ii<lines->GetNumberOfCells();ii++)
		{
			vtkIdList *ptidx=vtkIdList::New();
			linecell->GetNextCell(ptidx);
			idxlist *idxnew=new idxlist(ptidx->GetNumberOfIds()*3);
			for(int s=0;s<ptidx->GetNumberOfIds();s++)
			{
				lines->GetPoint(ptidx->GetId(s),tmppts);
				/*idxnew->SetValue(3*s,tmppts[0]);
				idxnew->SetValue(3*s+1,tmppts[1]);
				idxnew->SetValue(3*s+2,tmppts[2]);*/
				idxnew->SetPoint(s,tmppts);
				//idxnew->GetPoint(s,tmppts);
				for(int si=0;si<3;si++)
				{
					if(tmppts[si]>maxcoor[si])
						maxcoor[si]=tmppts[si];
					if(tmppts[si]<mincoor[si])
						mincoor[si]=tmppts[si];
				}
			}
			idxnew->embed.left=NULL;
			idxnew->embed.right=NULL;
			idxnew->SetIndex(ii);
			idxnew->CalLength();
			poly->push(&(idxnew->embed));
			idxtmp->sibling=&(idxnew->embed);
			idxtmp=idxtmp->sibling;
			ptidx->Delete();
		}
		

		idxtmp->sibling=NULL;
		this->alldist=new double[cellnum*cellnum];
		this->affinityegi=new double[cellnum*cellnum];
		this->affinitymat=new double[cellnum*cellnum];
		this->cur_clusts=new int[this->cellnum];
		this->best_clusts=new int[this->cellnum];
		this->clusts_center=new int[this->cellnum];

		for(itmp=0;itmp<cellnum;itmp++)
			for(ivar=0;ivar<cellnum;ivar++)
			{
				alldist[itmp*cellnum+ivar]=0.0;
			}
		lidxentry->left=NULL;
		lidxentry->right=NULL;
		for(itmp=0;itmp<1024;itmp++)
				doupts2[itmp]=VTK_DOUBLE_MAX;
		vcd->Delete();
		linecell->Delete();
		this->isinit=true;
		MatrixInitZeros(this->alldist,this->cellnum,this->cellnum);
		MatrixInitZeros(this->affinityegi,this->cellnum,this->cellnum);
		MatrixInitZeros(this->affinitymat,this->cellnum,this->cellnum);
		for(ivar=0;ivar<this->cellnum;ivar++)
		{
			this->cur_clusts[ivar]=0;
			this->best_clusts[ivar]=0;
			this->clusts_center[ivar]=-1;
		}
		/**/
	}
	else
	{
		return;
	}
};
void BaseCluster::CalDistInMesh( )
{
	//
	//this->Release();
	//
	if(this->isinit)
	{
		int ivar,itmp;
		idxlist* idxshow=NULL;
		list* lidxvar=NULL;

		double tmppts1[3];
		double doupts1[3];
		double doupts3[3];
		double doupts2[40960];
			//
		ivar=0;
		if(this->meshsize[0]*this->meshsize[1]*this->meshsize[2]>10)
		{
			for(int i=1;i<3;i++)
			{
				if(meshsize[ivar]>meshsize[i])
				{
					ivar=i;
				}
			}
		
			doupts3[ivar] = (maxcoor[ivar]-(mincoor[ivar]))/meshsize[ivar];
			tmppts1[0] = (maxcoor[ivar]-(mincoor[ivar]));
			for(int s=0;s<3;s++)
			{
			scale1[s]=((maxcoor[s]-(mincoor[s]))*doupts3[ivar])/(tmppts1[0]);
			}
		}
			//
			
		/**/
		list* idxtmp=idxentry->sibling;
		list* lidxtmp=lidxentry;
		idxlist *vartmp=NULL;
		while(idxtmp!=NULL)
		{
			ivar=0;
			idxshow=list_entry(idxtmp,idxlist,embed);
			if(this->meshsize[0]*this->meshsize[1]*this->meshsize[2]>10)
			{
				for(unsigned int s=0;s<idxshow->GetPtnum();s++)
				{
			
					idxshow->GetPoint(s,tmppts1);
					for(itmp=0;itmp<3;itmp++)
					{
						doupts1[itmp]=(int)((tmppts1[itmp]-mincoor[itmp])/scale1[itmp]);
						doupts3[itmp]=doupts1[itmp];
					}
					if((doupts2[ivar]!=doupts3[0])||(doupts2[ivar+1]!=doupts3[1])||(doupts2[ivar+2]!=doupts3[2]))
					{
						ivar+=3;
						doupts2[ivar]=doupts3[0]*scale1[0]+mincoor[0];
						doupts2[ivar+1]=doupts3[1]*scale1[1]+mincoor[1];
						doupts2[ivar+2]=doupts3[2]*scale1[2]+mincoor[2];
					}
				//s++;
				}
				idxlist *lidxnew=new idxlist(ivar);
				for(itmp=0;itmp<ivar;itmp++)
				{
					lidxnew->SetValue(itmp,doupts2[itmp+3]);
				}
				vartmp=lidxnew;
			}
			else
			{
				idxlist *llidxnew=new idxlist(idxshow->GetSize());
				for(unsigned int s=0;s<idxshow->GetSize();s++)
				{
					llidxnew->SetValue(s,idxshow->GetValue(s));
				}
				vartmp=llidxnew;
			}
			vartmp->embed.right=NULL;
			vartmp->embed.left=idxtmp;
			vartmp->SetIndex(idxshow->GetIndex());
			vartmp->CalLength();
			lidxtmp->sibling=&(vartmp->embed);
			lidxtmp=lidxtmp->sibling;
			idxtmp=idxtmp->sibling;
		}
		lidxtmp->sibling=NULL;
		lidxtmp=lidxentry->sibling;
		itmp=0;
		if(caltype==1)
		{
			while(lidxtmp!=NULL)
			{
				ivar=0;
				lidxvar=lidxentry->sibling;
				while(lidxvar!=NULL)
				{
					//lmindistance(lidxtmp,lidxvar,(alldist+itmp*cellnum+ivar));
					lchamferdistance(lidxtmp,lidxvar,(alldist+itmp*cellnum+ivar));
					lidxvar=lidxvar->sibling;
					ivar++;
				}
				lidxtmp=lidxtmp->sibling;
				itmp++;
			}
		}
		else if(caltype==2)
		{
			while(lidxtmp!=NULL)
			{
				ivar=0;
				lidxvar=lidxentry->sibling;
				while(lidxvar!=NULL)
				{
					lavgdistance(lidxtmp,lidxvar,(alldist+itmp*cellnum+ivar));
					lidxvar=lidxvar->sibling;
					ivar++;
				}
				lidxtmp=lidxtmp->sibling;
				itmp++;
			}	
		}
		//delete doupts2;
	}
 };
/**/
/**/
void BaseCluster::CalAffiMatrix(void)
{
	int coor[2];
	double* tmpDist=this->affinitymat;
	double threshold,tmpdou;	
	int nt=this->cellnum;
	int infot=0;
	int i,j,k;
	/*
	integer ldat =  nt;
    doublereal* wrt = (doublereal*)malloc( sizeof(doublereal) * nt) ;
    doublereal* wit = (doublereal*)malloc( sizeof(doublereal) * nt) ;   
    integer ldvrt = nt ;
    doublereal* vrt = (doublereal*)malloc( sizeof(doublereal) * nt * ldvrt) ;
    integer ldvlt = nt ;
    doublereal* vlt = (doublereal*)malloc( sizeof(doublereal) * nt * ldvlt) ;
    integer lworkt = nt * 4 ;
    doublereal *worktt = (doublereal*)malloc( sizeof(doublereal) * lworkt) ;
	dgeev_("V","V", (integer*)&nt,this->affinitymat, &ldat, wrt, wit, vlt , &ldvlt , vrt, &ldvrt, worktt, &lworkt, (integer*)&infot) ;
	*/
	int ldat =  nt;
    double* wrt = new double[nt] ;
    double* wit = new double[nt]  ;   
    int ldvrt = nt ;
    double* vrt = new double[nt * ldvrt] ;
    int ldvlt = nt ;
    double* vlt = new double[nt * ldvlt] ;
    int lworkt = nt * 4 ;
    double *worktt = new double[lworkt] ;
	MatrixInitZeros(wrt,nt,1);
	MatrixInitZeros(wit,nt,1);
	MatrixInitZeros(vrt,nt,nt);
	MatrixInitZeros(vlt,nt,nt);
	MatrixInitZeros(worktt,nt,4);
	MaxValue(this->alldist,nt,nt,coor);
	CopyMatrixToFile(this->alldist,nt,nt,"Dist_All");
	threshold=this->alldist[coor[0]*nt+coor[1]]*0.4;
	
	for(i=0;i<nt-1;i++)
	{
		for(j=i+1;j<nt;j++)
		{
			tmpdou=this->alldist[i*nt+j];
			if(tmpdou>threshold)
			{
				tmpDist[i*nt+j]=0.0;
				tmpDist[j*nt+i]=0.0;
			}
			else
			{
				tmpDist[i*nt+j]=0.0-tmpdou;
				tmpDist[j*nt+i]=0.0-tmpdou;
			}
		}
	}
	
	//CopyMatrixToFile(tmpDist,nt,nt,"tmpBefore");
	
	RowSum((tmpDist),wit,nt,nt);

	//CopyMatrixToFile(wit,1,nt,"wit");
	
	for(i=0;i<nt;i++)
	{
		wit[i]=0.0-wit[i];
		tmpDist[i*nt+i]+=wit[i];
	}
	
	for(i=0;i<nt;i++)
		wrt[i]=1.0/sqrt(wit[i]);
	//CopyMatrixToFile(wrt,1,nt,"wrt");
	
	for(i=0;i<nt;i++)
	{
		k=i*nt;
		for(j=0;j<nt;j++,k++)
			vrt[k]=tmpDist[k]*wrt[j];
	}
	//CopyMatrixToFile(vrt,nt,nt,"vrt");
	
	for(i=0;i<nt;i++)
	{
		k=i*nt;
		for(j=0;j<nt;j++,k++)
		{
			tmpDist[k]=vrt[k]*wrt[i];
		}
	}
	//CopyMatrixToFile(tmpDist,nt,nt,"Dist_L");
	
	dgeev_("V","V", (integer*)&nt,(doublereal*)(tmpDist),(integer*)&ldat,(doublereal*)wrt,(doublereal*)wit,(doublereal*)vlt,\
			(integer*)&ldvlt,(doublereal*)(this->affinityegi),(integer*)&ldvrt,(doublereal*)worktt,(integer*)&lworkt,(integer*)&infot) ;
	tmpDist=this->affinityegi;
	for(i=0;i<nt-1;i++)
	{
		for(j=i+1;j<nt;j++)
		{
			tmpdou=tmpDist[i*nt+j];
			tmpDist[i*nt+j]=tmpDist[j*nt+i];
			tmpDist[j*nt+i]=tmpdou;
		}
	}
	
	CopyMatrixToFile(this->affinityegi,nt,nt,"Dist_EgiV");
	//CopyMatrixToFile(wrt,1,nt,"Dist_EgiDR");
	//CopyMatrixToFile(wit,1,nt,"Dist_EgiDI");
	
		
		//
	
	delete []wrt;
	delete []wit;
	delete []vlt;
	delete []vrt;
	delete []worktt;

};
/**/
/**/
bool BaseCluster::BezerLine(unsigned int centerid,unsigned int ccid,vtkPoints *points,vtkCellArray *ln,vtkIntArray *tubearray,vtkIntArray *centerarray)
{
	idxlist *tmpcenter=NULL;
	idxlist *tmpcc=NULL;
	list *tmplist=NULL;
	list *centerlist=NULL;
	list *cclist=NULL;
	double centerhead[3]={0.0,0.0,0.0};
	double centertail[3]={0.0,0.0,0.0};
	double centercchead[3]={0.0,0.0,0.0};
	double centercctail[3]={0.0,0.0,0.0};
	double cchead[3]={0.0,0.0,0.0};
	double cctail[3]={0.0,0.0,0.0};	
	double pt[3]={0.0,0.0,0.0};
	double douhead=VTK_DOUBLE_MAX;
	double doutail=VTK_DOUBLE_MAX;
	double doutmp;
	int numpts;
	int numids;
	int numcellidx;
	int numhead;
	int numtail;
	list *tmpidx=NULL;
	if(centerid==ccid)
		return false;
	else
	{
		centerlist=poly->getlist(centerid);
		cclist=poly->getlist(ccid);
		tmpcenter=list_entry(centerlist,idxlist,embed);
		tmpcc=list_entry(cclist,idxlist,embed);
		tmpcc->GetPoint(0,cchead);
		tmpcc->GetPoint(tmpcc->GetPtnum()-1,cctail);
		for(unsigned int i=0;i<tmpcenter->GetPtnum();i++)
		{
			tmpcenter->GetPoint(i,pt);
			dist(pt,cchead,&doutmp);
			if(doutmp<douhead)
			{
				numpts=i;
				douhead=doutmp;
			}
			dist(pt,cctail,&doutmp);
			if(doutmp<doutail)
			{
			numids=i;
			doutail=doutmp;
			}
		}
		if(numpts>numids)
		{
			tmpcc->GetPoint(tmpcc->GetPtnum()-1,cchead);
			tmpcc->GetPoint(0,cctail);
			numcellidx=numpts;
			numpts=numids;
			numids=numcellidx;
		}
		tmpcenter->GetPoint(numpts,centerhead);
		tmpcenter->GetPoint(numids,centertail);
		//
		//
		dist(centerhead,cchead,&douhead);
		dist(centertail,cctail,&doutail);
		dist(cctail,cchead,&doutmp);
		if(doutmp>(douhead+doutail))
		{
			doutmp=0.0;
			numhead=3;
			numtail=3;
			{
				numcellidx=(numpts+numids)/2;
				while(numpts<(numcellidx)&&douhead>doutmp)
				{
				tmpcenter->GetPoint(++numpts,pt);
				dist(centerhead,pt,&doutmp);
				numhead++;
				}
				tmpcenter->GetPoint(numpts,centercchead);
				//
				numcellidx=centerarray->GetTuple1(numpts);
				numcellidx+=1;
				centerarray->InsertTuple1(numpts,numcellidx);
				//
				doutmp=0.0;
				while(numids>(numcellidx)&&doutail>doutmp)
				{
					tmpcenter->GetPoint(--numids,pt);
					dist(centertail,pt,&doutmp);
					numtail++;
				}
				tmpcenter->GetPoint(numids,centercctail);
				//
				numcellidx=centerarray->GetTuple1(numids);
				numcellidx--;
				centerarray->InsertTuple1(numids,numcellidx);
				//
			}
			double t=0.0;
			douhead=1.0/(double)(numhead);
			numcellidx=newpointnum;
			ln->InsertNextCell(numhead+1);
			points->InsertNextPoint(cchead);
			ln->InsertCellPoint(numcellidx++);
			tubearray->InsertNextTuple1(1);
			t=douhead;
			for(numpts=1;numpts<numhead;numpts++)
			{
				pt[0]=(1-t)*(1-t)*(cchead[0])+2*t*(1-t)*(centerhead[0])+t*t*(centercchead[0]);
				pt[1]=(1-t)*(1-t)*(cchead[1])+2*t*(1-t)*(centerhead[1])+t*t*(centercchead[1]);
				pt[2]=(1-t)*(1-t)*(cchead[2])+2*t*(1-t)*(centerhead[2])+t*t*(centercchead[2]);
				points->InsertNextPoint(pt);
				ln->InsertCellPoint(numcellidx++);
				tubearray->InsertNextTuple1(1);
				t+=douhead;
			}
			//
			points->InsertNextPoint(centercchead);
			ln->InsertCellPoint(numcellidx++);
			tubearray->InsertNextTuple1(1);
			//
			doutail=1.0/(double)(numtail);
			//
			ln->InsertNextCell(numtail+1);
			points->InsertNextPoint(cctail);
			ln->InsertCellPoint(numcellidx++);
			tubearray->InsertNextTuple1(1);
			//
			t=doutail;
			for(numpts=1;numpts<numtail;numpts++)
			{
				pt[0]=(1-t)*(1-t)*(cctail[0])+2*t*(1-t)*(centertail[0])+t*t*(centercctail[0]);
				pt[1]=(1-t)*(1-t)*(cctail[1])+2*t*(1-t)*(centertail[1])+t*t*(centercctail[1]);
				pt[2]=(1-t)*(1-t)*(cctail[2])+2*t*(1-t)*(centertail[2])+t*t*(centercctail[2]);
				//
				points->InsertNextPoint(pt);
				ln->InsertCellPoint(numcellidx++);
				tubearray->InsertNextTuple1(1);
				//
				t+=doutail;
			}
			//
			points->InsertNextPoint(centercctail);
			ln->InsertCellPoint(numcellidx++);
			tubearray->InsertNextTuple1(1);
			//
			newpointnum=numcellidx;
			return true;
		}
		else
		{
		return false;
		}
	}
};
//
/**/
/**/
/*
void BaseCluster::CalClusterCenter(cluster *clu)
{
	list* tmphead;
	list* tmplist;
	idxlist *tmpidx;
	idxlist *varidx;
	double tmpdistance=VTK_DOUBLE_MAX;
	unsigned int _tmpcenter;
	double tmpdouble=0.0;
	double clusize;
	tmplist=clu->idxhead;
	if(clu->idxhead==clu->idxtail)
	{
		_tmpcenter=clu->GetIndex();
	}
	else
	{
		do
		{
			tmphead=clu->idxhead;
			varidx=list_entry((tmplist),idxlist,embed);
			tmpdouble=0.0;
			clusize=0.0;
			do
			{
				tmpidx=list_entry((tmphead),idxlist,embed);
				tmpdouble+=*(this->alldist+(varidx->GetIndex())*cellnum+tmpidx->GetIndex());
				clusize++;
				tmphead=tmphead->right;
			}while(tmphead!=((clu->idxtail)->right));
			clusize--;
			if(clusize>0.0)
			{
				tmpdouble/=clusize;
			}
			if(tmpdouble<tmpdistance)
			{
				_tmpcenter=varidx->GetIndex();
				tmpdistance=tmpdouble;
			}
			tmplist=tmplist->right;
		}while(tmplist!=((clu->idxtail)->right));
	}
	clu->SetCenterId(_tmpcenter);
};
//
/**/
/**/
/*
void BaseCluster::CalClusterCentroid(cluster *clu)
{
	list* tmphead;
	list* tmplist;
	idxlist *tmpidx;
	idxlist *varidx;
	double tmpdistance=VTK_DOUBLE_MAX;
	unsigned int _tmpcenter;
	double tmpdouble=0.0;
	double avglength=0.0;
	tmplist=clu->idxhead;
	if((clu->idxhead==clu->idxtail)||clu->idxhead->right==clu->idxtail)
	{
		tmpidx=list_entry(clu->idxtail,idxlist,embed);
		varidx=list_entry(tmplist,idxlist,embed);
		_tmpcenter=varidx->GetIndex();
		if((tmpidx->GetLength())>(varidx->GetLength()))
			_tmpcenter=tmpidx->GetIndex();
	}
	else
	{
		_tmpcenter=0;
		do
		{
		varidx=list_entry((tmplist),idxlist,embed);
		avglength+=varidx->GetLength();
		tmplist=tmplist->right;
		_tmpcenter++;
		}while(tmplist!=((clu->idxtail)->right));
		avglength/=(double)_tmpcenter;
		tmplist=clu->idxhead;
		_tmpcenter=1;
		do
		{
			
			varidx=list_entry((tmplist),idxlist,embed);
			if(varidx->GetLength()>avglength)
			{
				tmphead=tmplist->right;
				tmpdouble=VTK_DOUBLE_MIN;
				do
				{
					tmpidx=list_entry((tmphead),idxlist,embed);
					if(tmpdouble<*(this->alldist+(varidx->GetIndex())*cellnum+tmpidx->GetIndex()))
						tmpdouble=*(this->alldist+(varidx->GetIndex())*cellnum+tmpidx->GetIndex());
					tmphead=tmphead->right;
				}while(tmphead!=((clu->idxtail)->right));
				if(varidx->GetLength()>tmpdouble)
				{
					tmpdouble/=(varidx->GetLength()*varidx->GetLength());
				}
				else
					tmpdouble=VTK_DOUBLE_MAX;
				if(tmpdouble<tmpdistance)
				{
					_tmpcenter=varidx->GetIndex();
					tmpdistance=tmpdouble;
				}
			}
			tmplist=tmplist->right;
		}while(tmplist!=((clu->idxtail)));
	}
	clu->SetCenterId(_tmpcenter);
};
//
/**/
void BaseCluster::CalClusterCentroid(int clust_num)
{
	list* tmphead;
	list* tmplist;
	idxlist *tmpidx;
	idxlist *varidx;
	double tmpdistance=VTK_DOUBLE_MAX;
	int _tmpcenter;
	double tmpdouble=0.0;
	double avglength=0.0;
	int vari,varj,vark;
	int i,j;
	int* int_tmp=new int[this->cellnum+1];
	for(i=0;i<this->cellnum+1;i++)
		int_tmp[i]=0;
	_tmpcenter=-1;
	for(i=0;i<this->cellnum;i++)
	{
		int_tmp[this->cur_clusts[i]]++;
		if(this->cur_clusts[i]==clust_num)
			this->clusts_center[i]=i;
	};
	if(int_tmp[clust_num]==1)
		return ;
	else
	{
		_tmpcenter=0;
		vari=0;
		for(i=0;i<this->cellnum;i++)
		{
			int_tmp[i]=-1;
			if(this->cur_clusts[i]==clust_num)
			{
				int_tmp[vari]=i;
				vari++;
			}
		}
		avglength=0.0;
		for(i=0;i<vari;i++)
		{
			tmplist=this->poly->getlist(int_tmp[i]);
			varidx=list_entry((tmplist),idxlist,embed);
			avglength+=varidx->GetLength();
		}
		avglength/=(double)(vari+1);
		_tmpcenter=0;
		for(i=0;i<vari;i++)
		{
			varj=int_tmp[i];
			tmplist=this->poly->getlist(varj);
			varidx=list_entry((tmplist),idxlist,embed);
			if(varidx->GetLength()>avglength)
			{
				tmpdouble=0.0;
				/**
				for(j=i+1;j<vari;j++)
				{
					vark=int_tmp[j];
					if(tmpdouble<this->alldist[varj*this->cellnum+vark])
						tmpdouble=this->alldist[varj*this->cellnum+vark];
				};
				*/
				for(j=0;j<vari;j++)
				{
					vark=int_tmp[j];
					tmpdouble+=this->alldist[varj*this->cellnum+vark];
				}
				tmpdouble/=(double)vari;
					/*
				if(varidx->GetLength()>tmpdouble*0.3)
					tmpdouble/=(varidx->GetLength());
				else
					tmpdouble=VTK_DOUBLE_MAX;
					/**/
				if(tmpdouble<tmpdistance)
				{
					_tmpcenter=varj;
					tmpdistance=tmpdouble;
				}
			}
		}
		for(i=0;i<vari;i++)
		{
			this->clusts_center[int_tmp[i]]=_tmpcenter;
		}
		
	}
};
/**/
//
unsigned int GetCenterid(cluster *clu)
{
	return clu->GetCenterId();
};
//
/**/
/**/
//
void BaseCluster::CopyLine(unsigned int linenum,vtkPoints *points,vtkCellArray *ln)
{
	idxlist *tmpidl=NULL;
	double pt[]={0.0,0.0,0.0};
	tmpidl=list_entry(this->poly->getlist(linenum),idxlist,embed);//?
	int numpts=this->newpointnum;
	int numids=(tmpidl->GetPtnum());
	int numcellidx=this->newpointnum;
	int numivar;
	ln->InsertNextCell(numids);
	for(numivar=0;numivar<numids;numivar++)
	{
		tmpidl->GetPoint(numivar,pt);
		points->InsertNextPoint(pt);
		ln->InsertCellPoint(numcellidx++);
	}
	this->newpointnum+=numids;
};
//
void BaseCluster::BestGroups()
{
	int* groups=new int[100];
	clock_t start,end; 
	double duration;
	for(int i=0;i<100;i++)
		groups[i]=i+2;
	double* rror=new double[this->cellnum*100];
	start=clock();
	double q=cluster_rotate((this->affinityegi),(this->cellnum),(this->cellnum),groups,6,(this->best_clusts),&(this->best_num),rror,1);
	end=clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC; 
	writePara("SelfTuneCluster - all - num - time  ",(double)this->cellnum,(double)this->cur_num,duration,0.0);
	delete []rror;
	delete []groups;
};
/**/
/**/
void BaseCluster::AdjustCurClusts(int type)
{
	int max,min;
	int i,j,le;
	int* clusts=NULL;	
	int* tmp=new int[this->cellnum];
	if (type==1)
		clusts=this->cur_clusts;
	else if(type==2)
		clusts=this->best_clusts;
	else
		return;
	max=clusts[0];
	min=clusts[0];
	for(i=1;i<this->cellnum;i++)
	{
		if(clusts[i]>max)
			max=clusts[i];
		else if(clusts[i]<min)
			min=clusts[i];
	}
	le=max-min+1;

	for(i=0;i<le;i++)
		tmp[i]=0;
	for(i=0;i<this->cellnum;i++)
		tmp[clusts[i]-min]++;	
	j=1;
	for(i=0;i<le;i++)
		if (tmp[i]>0)
		{
			tmp[i]=j;
			j++;
		}
	for(i=0;i<this->cellnum;i++)
		clusts[i]=tmp[clusts[i]-min];	
	if (type==1)
		this->cur_num=le;
	else if(type==2)
		this->best_num=le;
	delete []tmp;

};
//
inline void dist(double *fi,double *se,double *ret)
{
	//
	double dx=(*(fi)-*(se))*10;
	double dy=(*(fi+1)-*(se+1))*10;
	double dz=(*(fi+2)-*(se+2))*10;
	//
	*ret=dx*dx+dy*dy+dz*dz;
}
//
void lchamferdistance(list *thids,list *otids,double *ret)
{
	double thpt[3];
	double otpt[3];
	double var,douvar,doutmp;
	double tmp;	
	double gas[21]={0.0439,0.0796,0.1353,0.2163,0.3247,    
					0.4578,0.6065,0.7548,0.8825,0.9692,    
					1.0000,0.9692,0.8825,0.7548,0.6065,
					0.4578,0.3247,0.2163,0.1353,0.0796,    
					0.0439};
	int t,o;
	idxlist *thshow=NULL;
	idxlist *otshow=NULL;
	list *oth;
	list *the;
	list *listtmp;
	douvar=0.0;
	doutmp=VTK_DOUBLE_MIN;
	tmp=0.0;
	thshow=list_entry(thids,idxlist,embed);
	otshow=list_entry(otids,idxlist,embed);
	the=thids;
	oth=otids;
	int tsize=0;
	int osize=0;
	int gasindex=0;
	for(int i=0;i<2;i++)
	{
		douvar=0.0;
		thshow=list_entry(the,idxlist,embed);
		otshow=list_entry(oth,idxlist,embed);
		tsize=thshow->GetPtnum();
		osize=otshow->GetPtnum();
		for(t=0;t<(tsize);t++)
		{
			
			
			thshow->GetPoint(t,thpt);
			var=VTK_DOUBLE_MAX;
			{
				for(o=0;o<(osize);o++)
				{
				
					otshow->GetPoint(o,otpt);
					dist(thpt,otpt,&tmp);
					if(tmp<var)
							var=tmp;
				}
				tmp=(double)(t/tsize)*7;
				gasindex=tmp<21?(int)tmp:20;
				douvar+=(var*gas[gasindex]);
			}
		}
		douvar/=((double)t);
		if(douvar>doutmp)
			doutmp=douvar;
		listtmp=the;
		the=oth;
		oth=listtmp;
	}
	*ret=(doutmp);
};
/**/
//
void lmindistance(list *thids,list *otids,double *ret)
{
	double thpt[3];
	double otpt[3];
	double var,douvar;
	double tmp;	
	double gas[21]={0.0439,0.0796,0.1353,0.2163,0.3247,    
					0.4578,0.6065,0.7548,0.8825,0.9692,    
					1.0000,0.9692,0.8825,0.7548,0.6065,
					0.4578,0.3247,0.2163,0.1353,0.0796,    
					0.0439};
	unsigned int t,o;
	idxlist *thshow=NULL;
	idxlist *otshow=NULL;
	list *oth;
	list *the;
	douvar=0.0;
	tmp=0.0;
	thshow=list_entry(thids,idxlist,embed);
	otshow=list_entry(otids,idxlist,embed);
	if((thshow->GetSize())<(otshow->GetSize()))
	{
		the=thids;
		oth=otids;
	}
	else
	{
		the=otids;
		oth=thids;
	}
	{
		douvar=0.0;
		thshow=list_entry(the,idxlist,embed);
		otshow=list_entry(oth,idxlist,embed);
		for(t=0;t<(thshow->GetPtnum());t++)
		{	
		
			thshow->GetPoint(t,thpt);
			var=VTK_DOUBLE_MAX;
			{
				for(o=0;o<(otshow->GetPtnum());o++)
				{
			
					otshow->GetPoint(o,otpt);
					dist(thpt,otpt,&tmp);
					if(tmp<var)
						{
							var=tmp;
						}
				}
				douvar+=var;
			}
		}
	}
	*ret=(douvar/(double)t);
};
//
void lavgdistance(list *thids,list *otids,double *ret)
{
	double thpt[3];
	double otpt[3];
	double var=0.0;
	double tmp;	
	unsigned int t=0;
	unsigned int o=0;
	idxlist *thshow=NULL;
	idxlist *otshow=NULL;
	list *oth;
	list *the;
	var=0.0;
	tmp=0.0;
	oth=otids;
	the=thids;
	{
		thshow=list_entry(the,idxlist,embed);
		otshow=list_entry(oth,idxlist,embed);
		for(t=0;t<(thshow->GetPtnum());t++)
		{	

			thshow->GetPoint(t,thpt);
			{
				
				for(o=0;o<(otshow->GetPtnum());o++)
				{
			
					otshow->GetPoint(o,otpt);
					dist(thpt,otpt,&tmp);
					var+=(tmp);
				}

			}
		}
	}
	*ret=var/(double)(t*o);
};
//

//
void BaseCluster::ShowBestClusts()
{
	this->cur_num=this->best_num;
	for(int i=0;i<this->cellnum;i++)
		this->cur_clusts[i]=this->best_clusts[i];
}
/**/
/**/
void BaseCluster::GetCluster(vtkPolyData *output)
{
	vtkPolyData *newclu=vtkPolyData::New();
	vtkIntArray *newarray=vtkIntArray::New();
	vtkCellArray *newlines=vtkCellArray::New();
	vtkPoints *newpts=vtkPoints::New();
	this->newcellnum=0;
	this->newpointnum=0;
	newclu->Allocate();
	newclu->Initialize();
	newarray->Initialize();
	newarray->SetName("Cluster");
	int* idx=this->cur_clusts;
	//ivar=1;
	/*
	while(level_entry->current<level_entry->size-1)
	{
		clushow=list_entry(tmplist,cluster,embed);
		if(leflag>clushow->GetIndex())
		{	
				tmpidx=clushow->idxhead;
				do
				{
				idxshow=list_entry((tmpidx),idxlist,embed);
				CopyLine(idxshow->GetIndex(),newpts,newlines);
				newarray->InsertNextTuple1(ivar);
				tmpidx=tmpidx->right;
				}while(tmpidx!=((clushow->idxtail)->right));
			ivar++;
		}
		tmplist=*(level_entry->entry+level_entry->current);
		level_entry->current++;
	}
	*/
	for (unsigned int i=0;i<this->cellnum;i++)
	{
		CopyLine(i,newpts,newlines);
		newarray->InsertNextTuple1(idx[i]);
	}
	/**/
	newclu->SetPoints(newpts);
	newclu->SetLines(newlines);
	newclu->GetCellData()->SetScalars(newarray);
	newclu->GetFieldData()->CopyAllOn();
	output->DeepCopy(newclu);
	/**/
	newclu->Delete();
	newarray->Delete();
	newlines->Delete();
	newpts->Delete();

	/*tmplist=level_entry->get(num);
	clushow=list_entry(tmplist,cluster,embed);
	unsigned int leflag=level_entry->flag(num);
	ivar=0;
	while(level_entry->current<level_entry->size-1)
	{
		clushow=list_entry(tmplist,cluster,embed);
		tmpidx=clushow->idxhead;
		if(leflag>clushow->GetIndex())
		{
			do
			{
			idxshow=list_entry((tmpidx),idxlist,embed);
			justarray->SetTuple1(idxshow->GetIndex(),ivar);
			tmpidx=tmpidx->right;
			}while(tmpidx!=((clushow->idxtail)->right));	
			ivar++;
		}
		tmplist=*(level_entry->entry+level_entry->current);
		level_entry->current++;
	}
	lines->GetCellData()->SetScalars(justarray);
	lines->GetFieldData()->CopyAllOn();
	output->DeepCopy(lines);*/
};
/**
/**/
void BaseCluster::GetClusterOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
{
	vtkPolyData *newclu=vtkPolyData::New();
	vtkIntArray *newarray=vtkIntArray::New();
	vtkCellArray *newlines=vtkCellArray::New();
	vtkPoints *newpts=vtkPoints::New();
	this->newpointnum=0;
	this->newcellnum=0;
	newclu->Allocate();
	newclu->Initialize();
	newarray->Initialize();
	newarray->SetName("Cluster");
	int ivar=0;
	int* idx=this->cur_clusts;
	if(minnum>maxnum)
	{
		ivar=minnum;
		minnum=maxnum;
		maxnum=ivar;
	}
	if(minnum<1)
		minnum=1;
	if(maxnum>this->cur_num)
		maxnum=this->cur_num;
	for (unsigned int i=0;i<this->cellnum;i++)
	{	
		if(idx[i]>=(int)minnum&&idx[i]<=(int)maxnum)
		{
		CopyLine(i,newpts,newlines);
		newarray->InsertNextTuple1(idx[i]);
		}
	}
	/*ivar=1;
	while(level_entry->current<level_entry->size-1)
	{
		clushow=list_entry(tmplist,cluster,embed);
		if(leflag>clushow->GetIndex())
		{	
			if((ivar>=minnum)&&(ivar<=maxnum))
			{
				tmpidx=clushow->idxhead;
				do
				{
				idxshow=list_entry((tmpidx),idxlist,embed);
				CopyLine(idxshow->GetIndex(),newpts,newlines);
				newarray->InsertNextTuple1(ivar);
				tmpidx=tmpidx->right;
				}while(tmpidx!=((clushow->idxtail)->right));
			}
			if(ivar>maxnum)
			{
				leflag=VTK_INT_MAX;
				level_entry->current=level_entry->size;
			}
			ivar++;
		}
		tmplist=*(level_entry->entry+level_entry->current);
		level_entry->current++;
	}*/
	newclu->SetPoints(newpts);
	newclu->SetLines(newlines);
	newclu->GetCellData()->SetScalars(newarray);
	newclu->GetFieldData()->CopyAllOn();
	output->DeepCopy(newclu);
	newclu->Delete();
	newarray->Delete();
	newlines->Delete();
	newpts->Delete();
};

void BaseCluster::GetCenterOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
{
	vtkPolyData *newclu=vtkPolyData::New();
	vtkIntArray *newarray=vtkIntArray::New();
	vtkCellArray *newlines=vtkCellArray::New();
	vtkPoints *newpts=vtkPoints::New();
	this->newpointnum=0;
	this->newcellnum=0;
	newclu->Allocate();
	newclu->Initialize();
	newarray->Initialize();
	newarray->SetName("Cluster");
	int ivar;
	for(int i=0;i<this->cellnum;i++)
	{
		this->clusts_center[i]=-1;
	};
	if(minnum>maxnum)
	{
		ivar=minnum;
		minnum=maxnum;
		maxnum=ivar;
	}
	if(minnum<1)
		minnum=1;
	if(maxnum>this->cur_num)
		maxnum=this->cur_num;
	for(int i=0;i<this->cellnum;i++)
	{
		if(this->clusts_center[i]<0)
			CalClusterCentroid(this->cur_clusts[i]);
		if((i==this->clusts_center[i])&&(this->cur_clusts[i]>=minnum)&&(this->cur_clusts[i]<=maxnum))
		{
			CopyLine(i,newpts,newlines);
			newarray->InsertNextTuple1(this->cur_clusts[i]);
		}
	}
	newclu->SetPoints(newpts);
	newclu->SetLines(newlines);
	newclu->GetCellData()->SetScalars(newarray);
	newclu->GetFieldData()->CopyAllOn();
	output->DeepCopy(newclu);
	newclu->Delete();
	newarray->Delete();
	newlines->Delete();
	newpts->Delete();
};
//
void BaseCluster::GetBundleOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
{
	vtkPolyData *newclu=vtkPolyData::New();
	vtkIntArray *newarray=vtkIntArray::New();
	vtkIntArray *tubearray=vtkIntArray::New();
	vtkIntArray *centerarray=vtkIntArray::New();
	vtkCellArray *newlines=vtkCellArray::New();
	vtkPoints *newpts=vtkPoints::New();
	this->newpointnum=0;
	this->newcellnum=0;
	newclu->Allocate();
	newclu->Initialize();
	newarray->Initialize();
	newarray->SetName("Cluster");
	tubearray->Initialize();
	tubearray->SetName("Tube_");
	centerarray->Initialize();
	list* tmplist=NULL;
	idxlist *tmpcenter=NULL;
	int centerptnum;
	int tmps;
	int vars;
	int ivar;
	int tmpi,tmpj,tmpk;
	if(minnum>maxnum)
	{
		ivar=minnum;
		minnum=maxnum;
		maxnum=ivar;
	}
	if(minnum<1)
		minnum=1;
	if(maxnum>this->cur_num)
		maxnum=this->cur_num;
	for(int i=0;i<this->cellnum;i++)
	{
		this->clusts_center[i]=-1;
	};
	ivar=1;
	CopyMatrixToFile(this->cur_clusts,1,this->cellnum,"cluster_center_akk");
	for(tmpi=0;tmpi<this->cellnum;tmpi++)
	{
		if(this->cur_clusts[tmpi]>=minnum&&this->cur_clusts[tmpi]<=maxnum)
		{
			if(this->clusts_center[tmpi]<0)
				this->CalClusterCentroid(this->cur_clusts[tmpi]);
			CopyMatrixToFile(this->clusts_center,1,this->cellnum,"cluster_center");
			tmplist=this->poly->getlist((unsigned int)this->clusts_center[tmpi]);
			tmpcenter=list_entry(tmplist,idxlist,embed);
			//centerptnum=(list_entry(poly->getlist(this->clusts_center[tmpi]),idxlist,embed)->GetPtnum());
			centerptnum=(int)(tmpcenter->GetPtnum());
			for(tmps=0;tmps<centerptnum;tmps++)
				centerarray->InsertTuple1(tmps,0);
			centerarray->InsertTuple1(0,1);
				//
			if(BezerLine((this->clusts_center[tmpi]),tmpi,newpts,newlines,tubearray,centerarray))
			{
				newarray->InsertNextTuple1(this->cur_clusts[tmpi]);
				newarray->InsertNextTuple1(this->cur_clusts[tmpi]);
			}
			else
			{
						//CopyLine((intcenterid),newpts,newlines);
				CopyLine(this->clusts_center[tmpi],newpts,newlines);
				newarray->InsertNextTuple1(this->cur_clusts[tmpi]);
			}
			//
			vars=0;
			for(tmps=0;tmps<centerptnum;tmps++)
			{
				vars+=centerarray->GetTuple1(tmps);
				tubearray->InsertNextTuple1(vars);
				//centerarray->InsertTuple1(tmps,vars);
			}
				//
		}
	}

	newclu->SetPoints(newpts);
	newclu->SetLines(newlines);
	newclu->GetCellData()->AddArray(newarray);

	newclu->GetPointData()->AddArray(tubearray);
	//newclu->GetFieldData()->AddArray(newarray);
	newclu->GetFieldData()->CopyAllOn();
	output->DeepCopy(newclu);
	newclu->Delete();
	newarray->Delete();
	tubearray->Delete();
	centerarray->Delete();
	newlines->Delete();
	newpts->Delete();
};
void Minus(double *a,double *b,double *c){
	c[0]=a[0]-b[0];
	c[1]=a[1]-b[1];
	c[2]=a[2]-b[2];
}
void Multiple(double *a,double *b,double *c){
	//a × b= [a2b3-a3b2,a3b1-a1b3, a1b2-a2b1]
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
}
void add(double *a,double *v,double b,double *c){
	b/=10;
	c[0]=a[0]+b*v[0]*sqrt(v[0]*v[0]+v[1]*v[1])/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/sqrt(v[0]*v[0]+v[1]*v[1]);
	c[1]=a[1]+b*v[1]*sqrt(v[0]*v[0]+v[1]*v[1])/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/sqrt(v[0]*v[0]+v[1]*v[1]);
	c[2]=a[2]+b*v[2]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}
void BaseCluster::GetStreamtapesOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
{
	vtkPolyData *newclu=vtkPolyData::New();
	vtkIntArray *newarray=vtkIntArray::New();
	vtkIntArray *newptarray=vtkIntArray::New();
	vtkCellArray *newlines=vtkCellArray::New();
	vtkPoints *newpts=vtkPoints::New();
	this->newpointnum=0;
	this->newcellnum=0;
	newclu->Allocate();
	newclu->Initialize();
	newarray->Initialize();
	newarray->SetName("Cluster");
	newptarray->Initialize();
	newptarray->SetName("Point Cluster");
	list* tmplist=NULL;
	idxlist *tmpcenter=NULL;
	int tmps;
	int centerptnum;
	int ivar;
	for(int i=0;i<this->cellnum;i++)
	{
		this->clusts_center[i]=-1;
	};
	if(minnum>maxnum)
	{
		ivar=minnum;
		minnum=maxnum;
		maxnum=ivar;
	}
	if(minnum<1)
		minnum=1;
	if(maxnum>this->cur_num)
		maxnum=this->cur_num;
	double *T=new double[3];
	double *N=new double[3];
	double *B=new double[3];
	double *R=new double[3];
	double t=1.0;
	double *tmpt=new double[3];
	double *tmptl=new double[3];
	double *tmptll=new double[3];
	double *tmptlll=new double[3];
	double *tmptr=new double[3];
	double *tmptrr=new double[3];
	double *tmptrrr=new double[3];
	double *Tl=new double[3];
	double *Tr=new double[3];
	double *R1=new double[3];
	double *R2=new double[3];
	double *R3=new double[3];
	double *R4=new double[3];
	double v[3]={1.0,1.0,1.0};
	int pointid=0;
	int points=0;
	for(int i=0;i<this->cellnum;i++)
	{
		if(this->clusts_center[i]<0)
			CalClusterCentroid(this->cur_clusts[i]);
		if((i==this->clusts_center[i])&&(this->cur_clusts[i]>=minnum)&&(this->cur_clusts[i]<=maxnum))
		{  
			tmplist=this->poly->getlist(i);
			tmpcenter=list_entry(tmplist,idxlist,embed);
			centerptnum=(int)(tmpcenter->GetPtnum());
			if(centerptnum>7){
				points=pointid;
				tmpcenter->GetPoint(0,R1);
				add(R1,v,1.0,R2);
				newpts->InsertNextPoint(R1);
				newptarray->InsertNextTuple1(i);
				newpts->InsertNextPoint(R2);
				newptarray->InsertNextTuple1(this->cur_clusts[i]);
				pointid+=2;
				for(tmps=3;tmps<=centerptnum-4;tmps++){
					tmpcenter->GetPoint(tmps,tmpt);
					tmpcenter->GetPoint(tmps-1,tmptl);
					tmpcenter->GetPoint(tmps-2,tmptll);
					tmpcenter->GetPoint(tmps-3,tmptlll);
					tmpcenter->GetPoint(tmps+1,tmptr);
					tmpcenter->GetPoint(tmps+2,tmptrr);
					tmpcenter->GetPoint(tmps+3,tmptrrr);
					Minus(tmptr,tmptl,T);				
					Minus(tmpt,tmptll,Tl);
					Minus(tmptrr,tmpt,Tr);
					Minus(Tr,Tl,N);
					Minus(tmptl,tmptlll,Tl);
					Minus(tmptrrr,tmptr,Tr);
					Minus(Tr,Tl,R);
					Multiple(T,N,B);
					Multiple(T,R,Tl);
					t=T[0]*N[1]*R[2]+N[0]*R[1]*T[2]+R[0]*T[1]*N[2]-T[2]*N[1]*R[0]-N[0]*T[1]*R[2]-T[0]*N[2]*R[1];
					t/=sqrt(Tl[0]*Tl[0]+Tl[1]*Tl[1]+Tl[2]*Tl[2]);
					t=1-t;
					R3[0]=tmpt[0];
					R3[1]=tmpt[1];
					R3[2]=tmpt[2];
					add(R3,B,t,R4);
					newpts->InsertNextPoint(R3);
					newptarray->InsertNextTuple1(this->cur_clusts[i]);
					newpts->InsertNextPoint(R4);
					newptarray->InsertNextTuple1(this->cur_clusts[i]);
					pointid+=2;
					newlines->InsertNextCell(4);
					newlines->InsertCellPoint(pointid-4);
					newlines->InsertCellPoint(pointid-3);
					newlines->InsertCellPoint(pointid-1);
					newlines->InsertCellPoint(pointid-2);
					newarray->InsertNextTuple1(this->cur_clusts[i]);
				}
				R1[0]=(R3[0]+R4[0])/2;
				R1[1]=(R3[1]+R4[1])/2;
				R1[2]=(R3[2]+R4[2])/2;
				add(R1,B,t,R2);
				B[0]=0-B[0];
				B[1]=0-B[1];
				B[2]=0-B[2];
				add(R1,B,t,R3);
				add(R1,T,t,R4);
				newpts->InsertNextPoint(R2);
				newpts->InsertNextPoint(R3);
				newpts->InsertNextPoint(R4);
				newlines->InsertNextCell(3);
				newlines->InsertCellPoint(pointid++);
				newlines->InsertCellPoint(pointid++);
				newlines->InsertCellPoint(pointid++);
				newarray->InsertNextTuple1(this->cur_clusts[i]);
			}
		}
	}
	this->newpointnum+=pointid;
	//cout<<pointid<<endl;
	delete []T;
	delete []N;
	delete []B;
	delete []R;
	delete []tmpt;
	delete []tmptl;
	delete []tmptll;
	delete []tmptlll;
	delete []tmptr;
	delete []tmptrr;
	delete []tmptrrr;
	newclu->SetPoints(newpts);
	//newclu->SetLines(newlines);
	newclu->SetPolys(newlines);
	newclu->GetCellData()->SetScalars(newarray);
	newclu->GetPointData()->SetScalars(newptarray);
	newclu->GetFieldData()->CopyAllOn();
	output->DeepCopy(newclu);
	newclu->Delete();
	newarray->Delete();
	newlines->Delete();
	newpts->Delete();
};

void BaseCluster::CopyParaFromFile()
{
	CopyFileToMatrix(this->affinityegi,this->cellnum,this->cellnum,"Dist_EgiV");
	CopyFileToMatrix(this->alldist,this->cellnum,this->cellnum,"Dist_All");
	CopyFileToMatrix(this->affinitymat,this->cellnum,this->cellnum,"Dist_L");

};
/**/
void BaseCluster::ExportClassifyInfo(vtkPolyData* output)
{
		ofstream vtkfile;
		vtkfile.open ("streamLines.out", ios::out | ios::app | ios::binary);
		float tmppts_out[3]={0.0,0.0,0.0};
		double tmppts[3]={0.0,0.0,0.0};
		int stmp=0;
		vtkCellArray* linecell=vtkCellArray::New();
		stmp=output->GetNumberOfCells();
		vtkfile.write((char*)&stmp,sizeof(int));
		//vtkfile.write((char*)&space_out,1);

		if(vtkfile.is_open())
		{
			linecell=output->GetLines();
			linecell->InitTraversal();
			for(int ii=0;ii<output->GetNumberOfCells();ii++)
			{
				vtkIdList *ptidx=vtkIdList::New();
				linecell->GetNextCell(ptidx);
				stmp=ptidx->GetNumberOfIds();
				vtkfile.write((char*)&stmp,sizeof(int));
				//vtkfile.write((char*)&space_out,1);

				ptidx->Delete();
			}	
			linecell=output->GetLines();
			linecell->InitTraversal();
			for(int ii=0;ii<output->GetNumberOfCells();ii++)
			{
				vtkIdList *ptidx=vtkIdList::New();
				linecell->GetNextCell(ptidx);
				ptidx->GetNumberOfIds();
				for(int s=0;s<ptidx->GetNumberOfIds();s++)
				{
					output->GetPoint(ptidx->GetId(s),tmppts);
					tmppts_out[0]=tmppts[0];
					tmppts_out[1]=tmppts[1];
					tmppts_out[2]=tmppts[2];

					/*idxnew->SetValue(3*s,tmppts[0]);
					idxnew->SetValue(3*s+1,tmppts[1]);
					idxnew->SetValue(3*s+2,tmppts[2]);*/
					//idxnew->GetPoint(s,tmppts);
					vtkfile.write((char*)tmppts_out,sizeof(float)*3);
						//vtkfile.write((char*)&space_out,1);

						//<<tmppts[0]<<tmppts[1]<<tmppts[2];
				}
			ptidx->Delete();
			}
			ofstream vtkfile_2;
			vtkfile_2.open ("streamLines.CLASSIFY", ios::out | ios::app );
			vtkDataArray *classify=output->GetCellData()->GetScalars("Cluster");
			for(int ii=0;ii<classify->GetNumberOfTuples();ii++)
			{
				stmp=(int)classify->GetTuple1(ii);
				vtkfile.write((char*)&stmp,sizeof(int));
				vtkfile_2<<stmp<<endl;
			}
		}
		vtkfile.close();

};
/**
/**/
/**
/**/

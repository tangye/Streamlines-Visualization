
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
#include "HierarchicalCluster.h"
#include "SpectralCluster.h"
#include "KmeansCluster.h"
#include "BaseCluster.h"
//
using namespace std;
//
//

/**
points::points(int num)
{
	size=num;
	coor=new double[3*num];
};
//
points::~points()
{
	delete [] coor;	
};
//
void points::getpoint(int id,double *pt)
{
	*pt=*((this->coor)+3*id);
	*(pt+1)=*((this->coor)+3*id+1);
	*(pt+2)=*((this->coor)+3*id+2);
};
/**/
//
//-----------------------------------------------------------------------------
HierarchicalCluster::HierarchicalCluster(void)
{
		//meshsize=new double[3];
		//level_entry=NULL;
		tmplist=NULL;
		tmpidx=NULL;
		idxshow=NULL;
		//lines=vtkPolyData::New();
		ivar=0;
		itmp=0;
		dvar=0.0;
		dtmp=0.0;
		cellnum=1000;
		clushow=NULL;
		//poly=NULL;
		//linecell=vtkCellArray::New();
		clutmp=NULL;
		lstmp=NULL;
		//lentry=new list;
		up=NULL;
		//idxentry=new list;
		//lidxentry=new list;
		idxtmp=NULL;
		idxvar=NULL;
		lidxtmp=NULL;
		lidxvar=NULL;
		//idxlist *idxentry=new idxlist();
		s=0.0;
		//alldist=NULL;
		ps=new double[3];
	
		tmpdou=0.0;
		//mincoor=new double[3];
		//maxcoor=new double[3];
		//scale1=new double[3];
		tmppts1=new double[3];
		tmppts2=new double[3];
		////
		doupts1=new double[3];
		this->avg_lentry=new list;
		this->avg_level_entry=NULL;
		this->min_lentry=new list;
		this->min_level_entry=NULL;
		this->cur_level_entry=NULL;
		this->cur_lentry=NULL;
		doupts3=new double[3];
		this->hiertype=1;
		this->ismain=false;
		//pts=NULL;
		this->caltype=-1;

};
//----------------------------------------------------------------------------
HierarchicalCluster::~HierarchicalCluster()
{
		delete []ps;
		delete []mincoor;
		delete []maxcoor;
		delete []scale1;
		delete []tmppts1;
		delete []tmppts2;
		////
		delete []doupts1;
		delete []doupts2;
		delete []doupts3;
};
//
void HierarchicalCluster::Release(void)
{
	list *retmp=NULL;
	list *revar=NULL;
	idxlist *idxlistdel=NULL;
	cluster *clusterdel=NULL;
	retmp=lidxentry->sibling;
	/*while(retmp!=NULL)
	{
		revar=retmp->sibling;
		idxlistdel=list_entry(retmp,idxlist,embed);
		delete idxlistdel;
		retmp=revar;
	}*/
	for(int s=0;s<cur_level_entry->size-1;s++)
	{
		retmp=cur_level_entry->getvec(s);
		if(retmp!=NULL)
		{
		clusterdel=list_entry(retmp,cluster,embed);
		delete clusterdel;
		}
	}
}
/**
//
void HierarchicalCluster::InitPoly(vtkPolyData *input)
{
	vtkCellData *vcd=vtkCellData::New();
	
	lines->GetPointData()->PassData(input->GetPointData());
	lines->GetCellData()->PassData(input->GetCellData());
	lines->DeepCopy(input);
	
	linecell=lines->GetLines();
	linecell->InitTraversal();
	//pts=new points(lines->GetNumberOfPoints());
	poly=new polys(lines->GetNumberOfCells());
	double tmppts[3]={0.0,0.0,0.0};
	//
	
	for(int s=0;s<3;s++)
	{
		mincoor[s]=VTK_DOUBLE_MAX;
		maxcoor[s]=VTK_DOUBLE_MIN;
	}
	tmplist=lentry;
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
			//idxnew->SetValue(3*s,tmppts[0]);
			//idxnew->SetValue(3*s+1,tmppts[1]);
			//idxnew->SetValue(3*s+2,tmppts[2]);
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
		poly->push(&(idxnew->embed));
		idxtmp->sibling=&(idxnew->embed);
		idxtmp=idxtmp->sibling;
		ptidx->Delete();
	}
	idxtmp->sibling=NULL;
	alldist=new double[cellnum*cellnum];
	level_entry=new level(cellnum);
	for(itmp=0;itmp<cellnum;itmp++)
		for(ivar=0;ivar<cellnum;ivar++)
		{
			alldist[itmp*cellnum+ivar]=0.0;
		}
	lidxentry->left=NULL;
	lidxentry->right=NULL;
	for(itmp=0;itmp<1024;itmp++)
			doupts2[itmp]=VTK_DOUBLE_MAX;
};
/**/
/**
void HierarchicalCluster::CalDistInMesh(void)
{
	//
	//this->Release();
	//
	for(int s=0;s<3;s++)
	{
		scale1[s]=(maxcoor[s]-(mincoor[s]))/meshsize[s];
	}
	//
	idxtmp=idxentry->sibling;
	lidxtmp=lidxentry;
	idxlist *vartmp=NULL;
	while(idxtmp!=NULL)
	{
		ivar=0;
		idxshow=list_entry(idxtmp,idxlist,embed);
		if(this->meshsize!=0)
		{
			for(int s=0;s<idxshow->GetPtnum();s++)
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
					doupts2[ivar]=doupts3[0];
					doupts2[ivar+1]=doupts3[1];
					doupts2[ivar+2]=doupts3[2];
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
			for(int s=0;s<idxshow->GetSize();s++)
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
	if(type>0&&type<3)
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
	else if(type>2&&type<5)
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
};
/**/
/**/
void HierarchicalCluster::SetHierType(unsigned int type)
{
	this->hiertype=type;
};
unsigned int HierarchicalCluster::GetHierType(void)
{
	return this->hiertype;
};
void HierarchicalCluster::Main(void)
{	
	clock_t start,end; 
	double duration;
	start=clock();
	

	if(!this->ismain)
	{	
		avg_level_entry=new level(this->cellnum);
		min_level_entry=new level(this->cellnum);
		datainit(this->lidxentry,this->avg_lentry);
		datainit(this->lidxentry,this->min_lentry);
		minbottom2up(min_lentry,up,alldist,cellnum,min_level_entry);
		avgbottom2up(avg_lentry,up,alldist,cellnum,avg_level_entry);
		this->ismain=true;
	}
	if(this->hiertype==1)
	{
		this->cur_lentry=this->min_lentry;
		this->cur_level_entry=this->min_level_entry;
	}
	else if(this->hiertype==2)
	{
		this->cur_lentry=this->avg_lentry;
		this->cur_level_entry=this->avg_level_entry;	
	}
	tmplist=cur_level_entry->get(this->cur_num);
	clushow=list_entry(tmplist,cluster,embed);
	unsigned int leflag=cur_level_entry->flag(this->cur_num);
	ivar=0;
	while(cur_level_entry->current<cur_level_entry->size-1)
	{
		clushow=list_entry(tmplist,cluster,embed);
		if(leflag>clushow->GetIndex())
		{	
				tmpidx=clushow->idxhead;
				do
				{
					idxshow=list_entry((tmpidx->sibling),idxlist,embed);
					this->cur_clusts[idxshow->GetIndex()]=ivar;
					tmpidx=tmpidx->right;
				}while(tmpidx!=((clushow->idxtail)->right));
				ivar++;
		}
		tmplist=*(cur_level_entry->entry+cur_level_entry->current);
		cur_level_entry->current++;
	}
		//
	end=clock();
	duration = (double)(end - start) / CLOCKS_PER_SEC; 
	writePara("Hier Cluster - all - num - time  ",(double)this->cellnum,(double)this->cur_num,duration,0.0);

}
//
////void HierarchicalCluster::GetCluster(int levelnum,vtkPolyData *output)
////{
////	vtkPolyData *newclu=vtkPolyData::New();
////	vtkIntArray *newarray=vtkIntArray::New();
////	vtkCellArray *newlines=vtkCellArray::New();
////	vtkPoints *newpts=vtkPoints::New();
////	newpointnum=0;
////	newcellnum=0;
////	newclu->Allocate();
////	newclu->Initialize();
////	newarray->Initialize();
////	newarray->SetName("HierCluster");
////	tmplist=level_entry->get(levelnum);
////	clushow=list_entry(tmplist,cluster,embed);
////	unsigned int leflag=level_entry->flag(levelnum);
////	ivar=1;
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
////	newclu->SetPoints(newpts);
////	newclu->SetLines(newlines);
////	newclu->GetCellData()->SetScalars(newarray);
////	newclu->GetFieldData()->CopyAllOn();
////	output->DeepCopy(newclu);
////	newclu->Delete();
////	newarray->Delete();
////	newlines->Delete();
////	newpts->Delete();
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
/////**/
/////**/
/////**/
/////**/
////
////
////
////
////
/////**/
/////**/
//////
////void HierarchicalCluster::GetClusterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
////{
////	vtkPolyData *newclu=vtkPolyData::New();
////	vtkIntArray *newarray=vtkIntArray::New();
////	vtkCellArray *newlines=vtkCellArray::New();
////	vtkPoints *newpts=vtkPoints::New();
////	newpointnum=0;
////	newcellnum=0;
////	newclu->Allocate();
////	newclu->Initialize();
////	newarray->Initialize();
////	newarray->SetName("HierCluster");
////	tmplist=level_entry->get(levelnum);
////	clushow=list_entry(tmplist,cluster,embed);
////	unsigned int leflag=level_entry->flag(levelnum);
////	if(minnum>maxnum)
////	{
////		ivar=minnum;
////		minnum=maxnum;
////		maxnum=ivar;
////	}
////	if(minnum<1)
////		minnum=1;
////	if(maxnum>levelnum)
////		maxnum=levelnum;
////	ivar=1;
////	while(level_entry->current<level_entry->size-1)
////	{
////		clushow=list_entry(tmplist,cluster,embed);
////		if(leflag>clushow->GetIndex())
////		{	
////			if((ivar>=minnum)&&(ivar<=maxnum))
////			{
////				tmpidx=clushow->idxhead;
////				do
////				{
////				idxshow=list_entry((tmpidx),idxlist,embed);
////				CopyLine(idxshow->GetIndex(),newpts,newlines);
////				newarray->InsertNextTuple1(ivar);
////				tmpidx=tmpidx->right;
////				}while(tmpidx!=((clushow->idxtail)->right));
////			}
////			if(ivar>maxnum)
////			{
////				leflag=VTK_INT_MAX;
////				level_entry->current=level_entry->size;
////			}
////			ivar++;
////		}
////		tmplist=*(level_entry->entry+level_entry->current);
////		level_entry->current++;
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
////};
//////
////void HierarchicalCluster::GetCenterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
////{
////	vtkPolyData *newclu=vtkPolyData::New();
////	vtkIntArray *newarray=vtkIntArray::New();
////	vtkCellArray *newlines=vtkCellArray::New();
////	vtkPoints *newpts=vtkPoints::New();
////	this->newpointnum=0;
////	this->newcellnum=0;
////	newclu->Allocate();
////	newclu->Initialize();
////	newarray->Initialize();
////	newarray->SetName("HierCluster");
////	tmplist=level_entry->get(levelnum);
////	clushow=list_entry(tmplist,cluster,embed);
////	unsigned int leflag=level_entry->flag(levelnum);
////	if(minnum>maxnum)
////	{
////		ivar=minnum;
////		minnum=maxnum;
////		maxnum=ivar;
////	}
////	if(minnum<1)
////		minnum=1;
////	if(maxnum>levelnum)
////		maxnum=levelnum;
////	ivar=1;
////	while(level_entry->current<level_entry->size-1)
////	{
////		clushow=list_entry(tmplist,cluster,embed);
////		if(leflag>clushow->GetIndex())
////		{
////			if((ivar>=minnum)&&(ivar<=maxnum))
////			{
////				if(clushow->GetCenterId()<=0)
////					CalClusterCentroid(clushow);
////				CopyLine((clushow->GetCenterId()),newpts,newlines);
////				newarray->InsertNextTuple1(ivar);	
////			}
////			if(ivar>maxnum)
////			{
////				leflag=VTK_INT_MAX;
////				level_entry->current=level_entry->size;
////			}
////			ivar++;
////		}
////		tmplist=*(level_entry->entry+level_entry->current);
////		level_entry->current++;
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
////};
//////
////void HierarchicalCluster::GetBundleOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output)
////{
////	vtkPolyData *newclu=vtkPolyData::New();
////	vtkIntArray *newarray=vtkIntArray::New();
////	vtkIntArray *tubearray=vtkIntArray::New();
////	vtkIntArray *centerarray=vtkIntArray::New();
////	vtkCellArray *newlines=vtkCellArray::New();
////	vtkPoints *newpts=vtkPoints::New();
////	newpointnum=0;
////	newcellnum=0;
////	newclu->Allocate();
////	newclu->Initialize();
////	newarray->Initialize();
////	newarray->SetName("Cluster");
////	tubearray->Initialize();
////	tubearray->SetName("Tube_");
////	centerarray->Initialize();
////	tmplist=level_entry->get(levelnum);
////	clushow=list_entry(tmplist,cluster,embed);
////	unsigned int leflag=level_entry->flag(levelnum);
////	int intcenterid;
////	int centerptnum;
////	int tmps;
////	int vars;
////	int centervar;
////	if(minnum>maxnum)
////	{
////		ivar=minnum;
////		minnum=maxnum;
////		maxnum=ivar;
////	}
////	if(minnum<1)
////		minnum=1;
////	if(maxnum>levelnum)
////		maxnum=levelnum;
////	ivar=1;
////	while(level_entry->current<level_entry->size-1)
////	{
////		clushow=list_entry(tmplist,cluster,embed);
////		if(leflag>clushow->GetIndex())
////		{
////			if((ivar>=minnum)&&(ivar<=maxnum))
////			{
////				if(clushow->GetCenterId()<=0)
////					CalClusterCentroid(clushow);
////				intcenterid=clushow->GetCenterId();
////				//
////				centerptnum=(list_entry(poly->getlist(intcenterid),idxlist,embed)->GetPtnum());
////				for(tmps=0;tmps<centerptnum;tmps++)
////					centerarray->InsertTuple1(tmps,0);
////				centerarray->InsertTuple1(0,1);
////				//
////				tmpidx=clushow->idxhead;
////				do
////				{
////					idxshow=list_entry((tmpidx),idxlist,embed);
////					if(BezerLine(intcenterid,idxshow->GetIndex(),newpts,newlines,tubearray,centerarray))
////					{
////						newarray->InsertNextTuple1(ivar);
////						newarray->InsertNextTuple1(ivar);
////					}
////					else
////					{
////						//CopyLine((intcenterid),newpts,newlines);
////						centervar=ivar;
////					}
////				tmpidx=tmpidx->right;
////				}while(tmpidx!=((clushow->idxtail)->right));
////				CopyLine((intcenterid),newpts,newlines);
////				newarray->InsertNextTuple1(centervar);
////				//
////				vars=0;
////				for(tmps=0;tmps<centerptnum;tmps++)
////				{
////					vars+=centerarray->GetTuple1(tmps);
////					tubearray->InsertNextTuple1(vars);
////					//centerarray->InsertTuple1(tmps,vars);
////				}
////				//
////			}
////			if(ivar>maxnum)
////			{
////				leflag=VTK_INT_MAX;
////				level_entry->current=level_entry->size;
////			}
////			ivar++;
////		}
////		tmplist=*(level_entry->entry+level_entry->current);
////		level_entry->current++;
////	}
////	newclu->SetPoints(newpts);
////	newclu->SetLines(newlines);
////	newclu->GetCellData()->AddArray(newarray);
////
////	newclu->GetPointData()->AddArray(tubearray);
////	//newclu->GetFieldData()->AddArray(newarray);
////	newclu->GetFieldData()->CopyAllOn();
////	output->DeepCopy(newclu);
////	newclu->Delete();
////	newarray->Delete();
////	tubearray->Delete();
////	centerarray->Delete();
////	newlines->Delete();
////	newpts->Delete();
////};
//
//
//

//
/**
bool HierarchicalCluster::BezerLine(unsigned int centerid,unsigned int ccid,vtkPoints *points,vtkCellArray *ln)
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
		for(int i=0;i<tmpcenter->GetPtnum();i++)
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
				doutmp=0.0;
				while(numids>(numcellidx)&&doutail>doutmp)
				{
					tmpcenter->GetPoint(--numids,pt);
					dist(centertail,pt,&doutmp);
					numtail++;
				}
				tmpcenter->GetPoint(numids,centercctail);
			}
			double t=0.0;
			douhead=1.0/(double)(numhead);
			numcellidx=newpointnum;
			ln->InsertNextCell(numhead+1);
			points->InsertNextPoint(cchead);
			ln->InsertCellPoint(numcellidx++);
			t=douhead;
			for(numpts=1;numpts<numhead;numpts++)
			{
				pt[0]=(1-t)*(1-t)*(cchead[0])+2*t*(1-t)*(centerhead[0])+t*t*(centercchead[0]);
				pt[1]=(1-t)*(1-t)*(cchead[1])+2*t*(1-t)*(centerhead[1])+t*t*(centercchead[1]);
				pt[2]=(1-t)*(1-t)*(cchead[2])+2*t*(1-t)*(centerhead[2])+t*t*(centercchead[2]);
				points->InsertNextPoint(pt);
				ln->InsertCellPoint(numcellidx++);
				t+=douhead;
			}
			points->InsertNextPoint(centercchead);
			ln->InsertCellPoint(numcellidx++);
			doutail=1.0/(double)(numtail);
			ln->InsertNextCell(numtail+1);
			points->InsertNextPoint(cctail);
			ln->InsertCellPoint(numcellidx++);
			t=doutail;
			for(numpts=1;numpts<numtail;numpts++)
			{
				pt[0]=(1-t)*(1-t)*(cctail[0])+2*t*(1-t)*(centertail[0])+t*t*(centercctail[0]);
				pt[1]=(1-t)*(1-t)*(cctail[1])+2*t*(1-t)*(centertail[1])+t*t*(centercctail[1]);
				pt[2]=(1-t)*(1-t)*(cctail[2])+2*t*(1-t)*(centertail[2])+t*t*(centercctail[2]);
				points->InsertNextPoint(pt);
				ln->InsertCellPoint(numcellidx++);
				t+=doutail;
			}
			points->InsertNextPoint(centercctail);
			ln->InsertCellPoint(numcellidx++);
			newpointnum=numcellidx;
			return true;
		}
		else
		{
		return false;
		}
	}
};
/**/
/**/


void datainit(list *idxentry,list *lentry)
{
	int i=0,j=0;
	cluster *clu=new cluster();
	list *tmplist=NULL;
	list *idxtmp=NULL;
	tmplist=lentry;
	idxtmp=idxentry->sibling;
	tmplist->sibling=&(clu->embed);
	while(idxtmp!=NULL)
	{
		cluster *cluadd=new cluster();
		clu->SetIndex(i);
		tmplist=new list();
		tmplist->sibling=idxtmp;
		tmplist->left=NULL;
		tmplist->right=NULL;
		clu->idxhead=tmplist;
		clu->idxtail=tmplist;
		clu->embed.sibling=&(cluadd->embed);
		clu->embed.left=NULL;
		clu->embed.right=NULL;
		clu=cluadd;
		//list_entry(idxtmp,idxlist,embed)->Release();
		idxtmp=idxtmp->sibling;
		i++;
	}
	clu->embed.sibling=NULL;
	clu->idxhead=NULL;
};
list* mergerclu(list *left,list *right,unsigned int indextmp,level *level_entry)
{
	cluster *leftclu=NULL;
	cluster *rightclu=NULL;
	list *up=NULL;
	cluster *add =new cluster();
	leftclu=list_entry(left->sibling,cluster,embed);
	rightclu=list_entry(right->sibling,cluster,embed);
	(leftclu->idxtail)->right=rightclu->idxhead;
	add->idxhead=leftclu->idxhead;
	add->idxtail=rightclu->idxtail;
	add->SetIndex(indextmp);
	level_entry->push(left->sibling);
	level_entry->push(right->sibling);
	up=right->sibling->sibling;
	right->sibling=up;
	up=&(add->embed);
	up->left=left->sibling;
	up->right=right->sibling;
	up->sibling=left->sibling->sibling;
	left->sibling=up;
	return up;
}
//
void mergerclu_(list *left,list *right,unsigned int indextmp,level *level_entry)
{
	cluster *leftclu=NULL;
	cluster *rightclu=NULL;
	list *up=NULL;
	cluster *add =new cluster();
	leftclu=list_entry(left->sibling,cluster,embed);
	rightclu=list_entry(right->sibling,cluster,embed);
	(leftclu->idxtail)->right=rightclu->idxhead;
	add->idxhead=leftclu->idxhead;
	add->idxtail=rightclu->idxtail;
	add->SetIndex(indextmp);
	level_entry->push(left->sibling);
	level_entry->push(right->sibling);
	up=right->sibling->sibling;
	right->sibling=up;
	up=&(add->embed);
	up->left=left->sibling;
	up->right=right->sibling;
	up->sibling=left->sibling->sibling;
	left->sibling=up;

};
//
void minbottom2up(list *lentry,list *up,double *dist,unsigned cellnum,level *level_entry)
{
	cluster *cluth=NULL;
	cluster *cluot=NULL;
	cluster *clushow=NULL;
	cluster *clutmp=NULL;
	list *lsvar=NULL;
	list *lsth=NULL;
	list *lsot=NULL;
	list *lstmp=NULL;
	list *tmpclu=NULL;
	list *firclu=NULL;
	list *secclu=NULL;
	list *idlth=NULL;
	list *idlot=NULL;
	idxlist *idxth=NULL;
	idxlist *idxot=NULL;
	idxlist *idxshow=NULL;
	double mindist;
	double tmpdist;
	int ss=0;
	unsigned int intth;
	unsigned int intot;

	list *tmplist;
	//list *tmpidx;
	tmplist=lentry->sibling;
	//
	lsth=lentry;
	lstmp=lentry->sibling;
	unsigned int time=0;
	unsigned int indextmp=(cellnum*10);
	while(lstmp->sibling->sibling!=NULL)
	{
		lsth=lentry;
		mindist=VTK_DOUBLE_MAX;
		while(((lsth->sibling)->sibling)->sibling!=NULL)
		{
			lsot=lsth->sibling;
			cluth=list_entry(lsot,cluster,embed);
			while((lsot->sibling)->sibling!=NULL)
			{
					ss=0;
					tmpdist=VTK_DOUBLE_MAX;
					cluot=list_entry(lsot->sibling,cluster,embed);
					idlot=cluot->idxhead;
					idlth=cluth->idxhead;
					while(idlth!=NULL)
					{
						idxth=list_entry(idlth->sibling,idxlist,embed);
						intth=idxth->GetIndex();
						while(idlot!=NULL)
						{
							idxot=list_entry(idlot->sibling,idxlist,embed);
						intot=idxot->GetIndex();
						if(tmpdist>*(dist+intth*cellnum+intot))
							tmpdist=*(dist+intth*cellnum+intot);
						idlot=idlot->right;	
						}
						idlth=idlth->right;
					}
					if(mindist>tmpdist)
					{
					firclu=lsth;
					secclu=lsot;
					mindist=tmpdist;
					}
				lsot=lsot->sibling;
				}
			lsth=lsth->sibling;
			}
		tmpclu=mergerclu(firclu,secclu,indextmp,level_entry);	
		indextmp+=10;
		time++;
		lstmp=lentry->sibling;
	} 
	level_entry->push(lentry->sibling);
};
//
void avgbottom2up(list *lentry,list *up,double *dist,unsigned cellnum,level *level_entry)
{
	cluster *cluth=NULL;
	cluster *cluot=NULL;
	cluster *clushow=NULL;
	cluster *clutmp=NULL;
	list *lsvar=NULL;
	list *lsth=NULL;
	list *lsot=NULL;
	list *lstmp=NULL;
	list *tmpclu=NULL;
	list *firclu=NULL;
	list *secclu=NULL;
	list *idlth=NULL;
	list *idlot=NULL;
	idxlist *idxth=NULL;
	idxlist *idxot=NULL;
	idxlist *idxshow=NULL;
	double mindist;
	double tmpdist;
	unsigned int ss=0;
	unsigned int intth;
	unsigned int intot;

	list *tmplist;
	tmplist=lentry->sibling;
	
	//
	lsth=lentry;
	lstmp=lentry->sibling;
	unsigned int time=0;
	unsigned int indextmp=(cellnum*10);
	while(lstmp->sibling->sibling!=NULL)
	{
		lsth=lentry;
		mindist=VTK_DOUBLE_MAX;
		while(((lsth->sibling)->sibling)->sibling!=NULL)
		{
			lsot=lsth->sibling;
			cluth=list_entry(lsot,cluster,embed);
			while((lsot->sibling)->sibling!=NULL)
			{
					ss=0;
					tmpdist=0.0;
					cluot=list_entry(lsot->sibling,cluster,embed);
					idlot=cluot->idxhead;
					idlth=cluth->idxhead;
					while(idlth!=NULL)
					{
						idxth=list_entry(idlth->sibling,idxlist,embed);
						intth=idxth->GetIndex();
						while(idlot!=NULL)
						{
							idxot=list_entry(idlot->sibling,idxlist,embed);
						intot=idxot->GetIndex();
						tmpdist+=*(dist+intth*cellnum+intot);
						idlot=idlot->right;
						ss++;
						}
						idlth=idlth->right;
					} 
					tmpdist/=(double)ss;
					if(mindist>tmpdist)
					{
					firclu=lsth;
					secclu=lsot;
					mindist=tmpdist;
					}
				lsot=lsot->sibling;
				}
			lsth=lsth->sibling;
			}
		tmpclu=mergerclu(firclu,secclu,indextmp,level_entry);
		indextmp+=10;
		time++;
		lstmp=lentry->sibling;
	}
	level_entry->push(lentry->sibling);
};
//
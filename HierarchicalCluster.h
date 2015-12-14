
#include <vtkDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include "vtkIdList.h"
#include "vtkLookUpTable.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "BaseCluster.h"
#ifndef HierarchicalCluster_h
#define HierarchicalCluster_h
//
/*
	��������
	Min������̾���
	Avg����ƽ������
	��ϵ�һ������ȡ��������
	��ϵڶ����Ǿ�������
*/
//
/*class points
{
public:
	double *coor;
	int size;
	//points();
	points(int num);
	~points();
	void points::getpoint(int id,double *pt);
};*/
//
/*
	��Ƕ����listͨ���������ͣ�
	sibling ��ʾ�ֵܽڵ�
	left	��ʾ��ڵ�
	right	��ʾ�ҽڵ�

*/
/*
	�����װ���ͣ�
*/
class HierarchicalCluster : public BaseCluster
{
public:
	HierarchicalCluster(void);
  	~HierarchicalCluster(void);
	//void InitPoly(vtkPolyData *input);
	//void CalDistNoMesh(void);
	//void CalDistInMesh(void);
	//void Init(vtkPolyData *input);
	void Main(void);
	/*
	void GetCluster(int num,vtkPolyData *output);
	//void GetSpectralLevel(int num,vtkPolyData *output);
	//void GetKmeansLevel(int num,vtkPolyData *output);
	void GetClusterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetCenterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetBundleOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	//void HierarchicalCluster::CopyLine(unsigned int linenum,vtkPolyData *target);
	//bool BezerLine(unsigned int centerid,unsigned int cc,vtkPoints *points,vtkCellArray *ln);
	//bool BezerLine(unsigned int centerid,unsigned int ccid,vtkPoints *points,vtkCellArray *ln,vtkIntArray *tubearray,vtkIntArray *centerarray);
	//void CopyLine(unsigned int linenum,vtkPoints *points,vtkCellArray *line);
	//void CalClusterCenter(cluster *clu);
	//void CalClusterCentroid(cluster *clu);
	unsigned int GetCenterid(cluster *clu);*/
	void SetHierType(unsigned int type);
	unsigned int GetHierType(void);
	//void SetMeshSize(int mesh[3]);
	//int GetMeshSize(void);
	void Release(void);
private:
	//unsigned int centerid;
	//������ID
	//vtkPolyData* lines;
	//vtk������
	unsigned int hiertype;
	bool ismain;
	//��������
	//double *meshsize;
	//�����С
	//level* level_entry;
	//�����δ洢	
	//polys *poly;
	//����������Ϣ
	//list *idxentry;
	//ԭʼ����������
	//list *lidxentry;
	//��������������
	//list *lentry;
	//������������
	//double *mincoor;
	//��ά��С����ֵ
	//double *maxcoor;
	//��ά�������ֵ
	//double *scale1;
	//����ת������ֵ
	/***********-----------------*********/
		//���¾�Ϊ��ʱ���ݱ���
		//unsigned int cellnum;
		//unsigned int newcellnum;
		//unsigned int newpointnum;
		double dvar,dtmp;
		cluster *clushow;
		//vtkCellArray *linecell;
		cluster *clutmp;
		list *lstmp;
		list* tmplist;
		list* tmpidx;
		idxlist *idxshow;
		list *up;
		list *idxtmp;
		list *idxvar;
		list *lidxtmp;
		list *lidxvar;
		level* avg_level_entry;
		level* min_level_entry;
		level* cur_level_entry;
		//�����δ洢	
		list* avg_lentry;
		list* min_lentry;
		list* cur_lentry;
		//������������
		unsigned int ivar;
		unsigned int itmp;
		double s;
		double *ps;
		//double *alldist;

		double tmpdou;

		double *tmppts1;
		double *tmppts2;
		////
		
		double *doupts1;
		double *doupts2;
		double *doupts3;

		//points *pts;
	//
};

//
/*
	��������
*/
void list_add(list *add,list *prev,list *next);
//
/*
	����ɾ��
*/
void list_del(list *prev, list *next);
//
/**
	����ͷ��ʼ��
/**/
void list_init(list *head,list *tail);
//
/*void chamferdistance(list *thids,list *otids,points *pts,double *ret);
//
void mindistance(list *thids,list *otids,points *pts,double *ret);
//
void avgdistance(list *thids,list *otids,points *pts,double *ret);
*///

//
/*
	���������ֱ����
		��̡�ƽ���������
/**/
void minbottom2up(list *lentry,list *up,double *dist,unsigned int cellnum,level *level_entry);
/**
/**/
void avgbottom2up(list *lentry,list *up,double *dist,unsigned int cellnum,level *level_entry);
//
/**
	������ϲ�Ϊһ����
/**/
list* mergerclu(list *left,list *right,unsigned int indextmp,level *level_entry);
/**/
void mergerclu_(list *left,list *right,unsigned int indextmp,level *level_entry);
//
/**
	��������
/**/
inline void dist(double *fi,double *se,double *ret);
/**/
/**/
void datainit(list *idxentry,list *lentry);
/**/
#endif
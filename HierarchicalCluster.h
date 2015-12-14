
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
	聚类类型
	Min代表最短距离
	Avg代表平均距离
	组合第一个是求取距离类型
	组合第二个是聚类类型
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
	被嵌数据list通用链表类型，
	sibling 表示兄弟节点
	left	表示左节点
	right	表示右节点

*/
/*
	聚类封装类型，
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
	//中心线ID
	//vtkPolyData* lines;
	//vtk线数据
	unsigned int hiertype;
	bool ismain;
	//聚类类型
	//double *meshsize;
	//网格大小
	//level* level_entry;
	//聚类层次存储	
	//polys *poly;
	//几何数据信息
	//list *idxentry;
	//原始线数据链表
	//list *lidxentry;
	//网格化线数据链表
	//list *lentry;
	//聚类数据链表
	//double *mincoor;
	//三维最小坐标值
	//double *maxcoor;
	//三维最大坐标值
	//double *scale1;
	//坐标转化缩放值
	/***********-----------------*********/
		//以下均为临时数据变量
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
		//聚类层次存储	
		list* avg_lentry;
		list* min_lentry;
		list* cur_lentry;
		//聚类数据链表
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
	链表增加
*/
void list_add(list *add,list *prev,list *next);
//
/*
	链表删除
*/
void list_del(list *prev, list *next);
//
/**
	链表头初始化
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
	以下两个分别进行
		最短、平均距离聚类
/**/
void minbottom2up(list *lentry,list *up,double *dist,unsigned int cellnum,level *level_entry);
/**
/**/
void avgbottom2up(list *lentry,list *up,double *dist,unsigned int cellnum,level *level_entry);
//
/**
	两个类合并为一个类
/**/
list* mergerclu(list *left,list *right,unsigned int indextmp,level *level_entry);
/**/
void mergerclu_(list *left,list *right,unsigned int indextmp,level *level_entry);
//
/**
	两点间距离
/**/
inline void dist(double *fi,double *se,double *ret);
/**/
/**/
void datainit(list *idxentry,list *lentry);
/**/
#endif
#ifndef BaseCluster_h
#define BaseCluster_h

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
#include "SelfTuningCluster.h"
#include "myTools.h"
/*
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"
*/
using namespace std;

//
/*
	聚类类型
	Min代表最短距离
	Avg代表平均距离
	组合第一个是求取距离类型
	组合第二个是聚类类型
*/
/*
	定义嵌入指针提取宏
*/
#define list_entry(ptr, type, member) \
    ((type *)((char *)(ptr) - (unsigned long)(&((type *)0)->member)))
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
class list
{
public:
	list *sibling;
	list *left;
	list *right;
};
/*
	存储线数据的链表类型，
	embed	嵌入链表类型
	idx		线数据类型存储地址
	index	索引号
	length	线长
*/
class idxlist
{
public :
	list embed;
	//idxlist();
	idxlist(unsigned int num);
	~idxlist();
	unsigned int getidx(int id);
	void Release(void);
	void CalLength(void);
	double GetLength(void);
	unsigned int GetSize(void);
	unsigned int GetPtnum(void);
	void SetSize(int num);
	int GetIndex(void );
	void SetIndex(int num);
	double GetValue(int valuenum);
	void SetValue(int valuenum,double newvalue);
	void GetPoint(unsigned int ptnum,double *ret);
	void SetPoint(unsigned int idx,double *newpt);
private:
	unsigned int size;
	double *idx;
	unsigned int index;
	double length;
};
//
/*
	聚类用的链表类型，
	embed		嵌入链表类型
	idxhead		包含线数据链表首
	idxtail		包含线数据链表首
	index		索引编号
*/
class cluster
{
public:
		cluster(void);
		~cluster(void);
		list *idxhead;
		list *idxtail;
		list embed;
		int GetIndex(void );
		void SetIndex(int num);
		void SetCenterId(unsigned int value);
		unsigned int GetCenterId();
private:
	unsigned int centerid;
	unsigned int index;
};
//
/*
	线性指针存储类型，
	entry		链表指针表地址
	current		当前指针索引
	size		指针表大小
*/
class linecontainer
{
public:
	list **entry;
	unsigned int size;
	unsigned  int current;
	void push(list *add);
	void pop(list *del);
	list* getvec(unsigned int num);
};
//
/*
	聚类层次存储类型，继承与linecontainer，
	flag		用层次数生成的标识，用于从指针表剔除非此层的数据
*/
class level:public linecontainer
{
public:
	level();
	level(unsigned int num);
	~level();
	list *get(unsigned int level_num);
	unsigned int flag(unsigned int num);
};
//
/*
	集合信息存储类型，
	getlist	用于从指针表找到第n条线
*/
class polys:public linecontainer
{
public:
	polys(unsigned int num);
	~polys();
	list* getlist(unsigned int list_num);
private:
	//numofsize;
	//numofcell;
};
/**/
class BaseCluster
{
public:
	BaseCluster(void);
  	~BaseCluster(void);
	void InitPoly(vtkPolyData *input);
	void CalDistInMesh( );
	void CalAffiMatrix(void);
	void Init(vtkPolyData *input);
	virtual void Main(void)=0;
	void ShowBestClusts();
	void GetCluster(vtkPolyData *output);
	void AdjustCurClusts(int i);
	//void GetSpectralLevel(int num,vtkPolyData *output);
	//void GetKmeansLevel(int num,vtkPolyData *output);
	void GetClusterOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetCenterOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetBundleOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetStreamtapesOnly(unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	//void LevelOfCluster::CopyLine(unsigned int linenum,vtkPolyData *target);
	//bool BezerLine(unsigned int centerid,unsigned int cc,vtkPoints *points,vtkCellArray *ln);
	bool BezerLine(unsigned int centerid,unsigned int ccid,vtkPoints *points,vtkCellArray *ln,vtkIntArray *tubearray,vtkIntArray *centerarray);
	void CopyLine(unsigned int linenum,vtkPoints *points,vtkCellArray *line);
	//void CalClusterCenter(cluster *clu);
	//void CalClusterCentroid(cluster *clu);
	void CalClusterCentroid(int clust_num);
	//
	void ExportClassifyInfo(vtkPolyData* output);
	//
	void CopyParaFromFile();
	unsigned int GetCenterid(cluster *clu);
	void SetCalType(unsigned int type);
	unsigned int GetCalType(void);
	void SetCurNum(unsigned int type);
	unsigned int GetCurNum(void);
	void BestGroups();
	void SetMeshSize(int mesh[3]);
	int GetMeshSize(void);
	void Release(void);
protected:
	//unsigned int centerid;
	//中心线ID
	static vtkPolyData* lines;
	//vtk线数据
	static unsigned caltype;
	//

	//聚类类型
	static double *meshsize;
	//网格大小
	static polys *poly;
	//几何数据信息
	static list *idxentry;
	//原始线数据链表
	static list *lidxentry;
	//网格化线数据链表
	static double *mincoor;
	//三维最小坐标值
	static double *maxcoor;
	//三维最大坐标值
	static double *scale1;
	//坐标转化缩放值
	/***********-----------------*********/
	static unsigned int cellnum;
	static unsigned int newcellnum;
	static unsigned int newpointnum;
	static double* alldist;
	static double* affinitymat;
	static double* affinityegi;
	static unsigned int best_num;
	static int* best_clusts;
	static int* clusts_center;
	static unsigned int cur_num;
	static int* cur_clusts;
	static bool isinit;
		//points *pts;
	//
};
/**
	以下三个分别求取
		chamfer、最短、平均距离
**/
void lchamferdistance(list *thids,list *otids,double *ret);
/**/
void lmindistance(list *thids,list *otids,double *ret);
/**/
void lavgdistance(list *thids,list *otids,double *ret);
/**/
/**/
/**
	由线链表构造聚类链表

/**/
/**/
#endif
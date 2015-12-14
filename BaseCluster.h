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
	��������
	Min������̾���
	Avg����ƽ������
	��ϵ�һ������ȡ��������
	��ϵڶ����Ǿ�������
*/
/*
	����Ƕ��ָ����ȡ��
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
	��Ƕ����listͨ���������ͣ�
	sibling ��ʾ�ֵܽڵ�
	left	��ʾ��ڵ�
	right	��ʾ�ҽڵ�

*/
class list
{
public:
	list *sibling;
	list *left;
	list *right;
};
/*
	�洢�����ݵ��������ͣ�
	embed	Ƕ����������
	idx		���������ʹ洢��ַ
	index	������
	length	�߳�
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
	�����õ��������ͣ�
	embed		Ƕ����������
	idxhead		����������������
	idxtail		����������������
	index		�������
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
	����ָ��洢���ͣ�
	entry		����ָ����ַ
	current		��ǰָ������
	size		ָ����С
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
	�����δ洢���ͣ��̳���linecontainer��
	flag		�ò�������ɵı�ʶ�����ڴ�ָ����޳��Ǵ˲������
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
	������Ϣ�洢���ͣ�
	getlist	���ڴ�ָ����ҵ���n����
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
	//������ID
	static vtkPolyData* lines;
	//vtk������
	static unsigned caltype;
	//

	//��������
	static double *meshsize;
	//�����С
	static polys *poly;
	//����������Ϣ
	static list *idxentry;
	//ԭʼ����������
	static list *lidxentry;
	//��������������
	static double *mincoor;
	//��ά��С����ֵ
	static double *maxcoor;
	//��ά�������ֵ
	static double *scale1;
	//����ת������ֵ
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
	���������ֱ���ȡ
		chamfer����̡�ƽ������
**/
void lchamferdistance(list *thids,list *otids,double *ret);
/**/
void lmindistance(list *thids,list *otids,double *ret);
/**/
void lavgdistance(list *thids,list *otids,double *ret);
/**/
/**/
/**
	�����������������

/**/
/**/
#endif
#ifndef SpectralCluster_h
#define SpectralCluster_h    
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <Windows.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
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
using namespace std;

/**/
/**/
static double _MinDouble=-1.0E-3;
static double _MaxDouble=1.0E-3;
/**/
/**/
template <class T> void HeapSort(T* const num,int* const idx,int size);
/**/
/**/
template <class T> void BuildHeap(T* const num ,int* const idx,int size);  
/**/
/**/
template <class T> void PercolateDown(T* const num, int* const idx, int index,int size);
/**/
/**/
template <class T> void Swap(T* const num, int v, int u) ;
/**/
/**/
void Swap(int* const num, int v, int u);  
//**
/**/
template <class T> bool CopyMatrixToFile(const T* Mat,int row,int column,const char* const name);
//**
/**/
template <class T> bool CopyFileToMatrix(const T* Mat,int row,int column,const char* const name);
/**/
/**/
void RowMultiColumn(double* const Res,const double* const Fir,const double* const Sec,int row,int firnumofcolumn,int secnumofcolumn,int column);
/**/
/**/
void MatrixMultiMatrix(double* const Res,const double* const Fir,const double* const Sec,int firNumOfRow,int firNumOfColumn,int secNumOfRow,int secNumOfColumn);
/**/
/**/
bool MatrixIsSymmetry(const double* const Mat,int row,int column);
/**/
/**/
void MatrixSetOne(double* const Mat,int row,int column);
/**/
/**/
void OneMatrix(double* const Mat,int row,int column);
/**/
/**/
double MaxValue(const double* const Mat,int row,int column,int *res);
/**/
/**/
void MatrixInitZeros(double* const Mat,int row,int column);
/**/
/**/
void MatrixInitZeros(int* const Mat,int row,int column);
/**/
/**/
void MaxRowColumnValue(const double* const Mat,double* const Val,int* const Index,int row,int column,int type);
/**/
/**/
void MaxRowColumnAbsValue(const double* const Mat,double* const Val,int* const Index,int row,int column,int type);
/**/
/**/
void MinRowColumnValue(const double* const Mat,double* const Val,int* const Index,int row,int column,int type);
/**/
/**/
void Transpose(double* const Mat,double* const Tam,int row);
/**/
/**/
void DisplayMatrix(const double* const Mat,int row,int column);
/**/
/**/
void RowSum(const double* const Mat,double* const _Mat,int row,int column);
/**/
/**/
void MatrixSymmetryMulti(double* const Res,const double* Mat,const double alpha,const int row,const int column,int p,int q);
/**/
/**/
void MatrixRightRotate(double* const Res,const double* const Mat,const double alpha,const int row,const int column,const int p,const int q);
/**/
/**/
void MatrixLeftRotate(double* const Res,const double* const Mat,const double alpha,const int row,const int column,const int p,const int q);
/**/
/**/
void Jacobi(const double* const Mat,double* const matV,double* const matD,const int row);
/**/
/**/
/**/
/**/
double MatrixDistance(const double* const Fir,const double* const Sec,int FirRow,int SecRow,int Column);
/**/
/**/
/**/
/**/
void Kmeans(const double* const Mat,int* clusteridx,int row,int column,int k,char name[]);
/**/
/**/
template <class T> void AdjustVectorByIdx(T* matV,int* idx,int row,int column);
/**/
/**/
/**/
/**/
void MainSpectralCluster(double* Mat,int row,int k);

/**/
/**/
class SpectralCluster:public BaseCluster
{
public:
	SpectralCluster();
	//
	~SpectralCluster();
	//
	void Main(void);
	//
	void NormalizedSpectralCluster(int* const idx,int row,int k,int dk);
	//
	//
};
/**/
/**/
/**/
/**/
/**/




#endif
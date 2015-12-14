#ifndef KmeansCluster_h
#define KmeansCluster_h
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
using namespace std;
/**/
#define MaxTValue 100000;
//
//template <class T> 
void MainKmeansCluster(const double* mat,int* cluster,int row,int k);
/**/
/**/
class KmeansCluster:public BaseCluster
{
public:
	KmeansCluster(void);
	~KmeansCluster(void);
	void Main(void);
	/*
	void GetCluster(int num,vtkPolyData *output);
	void GetClusterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetCenterOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	void GetBundleOnly(unsigned int levelnum,unsigned int minnum,unsigned int maxnum,vtkPolyData *output);
	*/
};
/**/
/**/

#endif
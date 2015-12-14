#include"SpectralCluster.h"
template <class T> void HeapSort(T* const num,int* const idx,int size)
{
    int i;
    int iLength=size;
	for(i=0;i<size;i++)
		idx[i]=i;
    BuildHeap(num,idx,size);  
    
    for (i = iLength - 1; i >= 1; i--) 
	{   
        Swap(num, 0, i);
		Swap(idx, 0, i);
        size--;  
        PercolateDown(num,idx, 0,size);

    }
}

   
template <class T> void BuildHeap(T* const num ,int* const idx,int size) 
{ 
    int i; 
   
    for (i = size / 2 - 1; i >= 0; i--) 
	{   
        PercolateDown(num, idx,i,size);
   
    }   
}
    

template <class T> void PercolateDown(T* const num, int* const idx, int index,int size) 
{   
    int min;  
    while (index * 2 + 1<size) 
	{	   
        min = index * 2 + 1; 
        if (index * 2 + 2 < size) 
		{
            if (num[min] > num[index * 2 + 2]) 
			{
                min = index * 2 + 2;  
            }   
        }   
        if (num[index] < num[min]) 
		{   
            break; 
        } 
		else 
		{   
			Swap(idx, index, min);

			//printf("idx swap id is \t%d\t%d\t ,num is %d\t%d\n",index,min,idx[index],idx[min]);

			Swap(num, index, min);
			//printf("num swap id is \t%d\t%d\t ,num is %d\t%d\n",index,min,num[index],num[min]);

            index = min;
        }
    }  
}
    
template <class T> void Swap(T* const num, int v, int u) 
{  
    T temp = num[v];   
    num[v] = num[u];   
    num[u] = temp;   
}   
void Swap(int* const num, int v, int u) 
{  
    int temp = num[v];   
    num[v] = num[u];   
    num[u] = temp;   
}   
/**/
/**/
double _abs(double tmp)
{
	if(tmp>0.0)
		return tmp;
	else
		return (0.0-tmp);
}
//**
/**/
template <class T> bool CopyMatrixToFile(const T* Mat,int row,int column,const char* const name)
{
	ofstream hFile(name, ios::out, _SH_DENYNO);
	if(hFile)
	{
		for(int tmpi=0;tmpi<row;tmpi++)
		{
			for(int tmpj=0;tmpj<column;tmpj++)
			{
			hFile<<(double)Mat[tmpi*column+tmpj]<<'\t';
			}
			hFile<<endl;
		}
		hFile.close();
		return true;
	}
	else
	{
		hFile.close();	
		return false;
	}
};
/**/
/**/
template <class T> bool CopyFileToMatrix(const T* Mat,int row,int column,const char* const name)
{
	ifstream hFile(name, ios::in, _SH_DENYNO);
	if(hFile)
	{
		for(int tmpi=0;tmpi<row*column;tmpi++)
		{
			hFile>>(double)Mat[tmpi];
		}	
		hFile.close();
		return true;
	}
	else
	{
		hFile.close();
		return false;
	}
};
/**/
/**/
void RowMultiColumn(double* const Res,const double* const Fir,const double* const Sec,int row,int firnumofcolumn,int secnumofcolumn,int column)
{
	if(Fir==NULL||Sec==NULL||Res==NULL)
		return;
	else
	{
		*Res=0.0;
		for(int i=0;i<firnumofcolumn;i++)
		{
			*Res+=(*(Fir+row*firnumofcolumn+i))*(*(Sec+(i*secnumofcolumn)+column));
		}
	}
}
/**/
/**/
void MatrixMultiMatrix(double* const Res,const double* const Fir,const double* const Sec,int firNumOfRow,int firNumOfColumn,int secNumOfRow,int secNumOfColumn)
{
	if(secNumOfRow!=firNumOfColumn||Fir==NULL||Sec==NULL||Res==NULL)
		return;
	else
	{
		for(int i=0;i<firNumOfRow;i++)
			for(int j=0;j<secNumOfColumn;j++)
			{
				RowMultiColumn((Res+i*secNumOfColumn+j),Fir,Sec,i,firNumOfColumn,secNumOfColumn,j);
			}
	}
}
/**/
/**/
bool MatrixIsSymmetry(const double* const Mat,int row,int column)
{
	if(Mat==NULL||row!=column)
		return false;
	else
	{
		int i=0;
		int j=0;
		for(i=0;i<row-1;i++)
			for(j=i+1;j<row;j++)
			{
				if(_abs(Mat[i*column+j]-Mat[j*column+i])>_MaxDouble)
					return false;

			}
		return true;
	}
}
/**/
/**/
void MatrixSetOne(double* const Mat,int row,int column)
{
	if(Mat==NULL||row!=column)
		return ;
	else
	{
		int j=0;
		int i=0;
		for(i=0;i<row;i++)
			for(j=0;j<column;j++)
				Mat[i*column+j]=0.0f;
		for(i=0;i<row;i++)
			Mat[i*column+i]=1.0f;
	}
}
/**/
/**/
void OneMatrix(double* const Mat,int row,int column)
{
	if(Mat==NULL||row!=column)
		return ;
	else
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
			Mat[i*column+j]=1.0f;
	}
}
/**/
/**/
double MaxValue(const double* const Mat,int row,int column,int *res)
{
	if(Mat==NULL||res==NULL)
	{
		res[0]=-1;
		res[1]=-1;
		return 0.0;
	}
	else
	{
		double tmp=0.0f;
		res[0]=-1;
		res[1]=-1;
		int i=0;
		int j=0;
		for(i=0;i<row-1;i++)
			for(j=i+1;j<column;j++)
				if(_abs(Mat[i*column+j])>tmp)
				{
					tmp=_abs(Mat[i*column+j]);
					res[0]=i;
					res[1]=j;
				}
		return tmp;
	}
}
/**/
/**/
void MatrixInitZeros(double* const Mat,int row,int column)
{
	for(int i=0;i<row*column;i++)
		Mat[i]=0.0f;
}
/**/
/**/
void MatrixInitZeros(int* const Mat,int row,int column)
{
	for(int i=0;i<row*column;i++)
		Mat[i]=0;
}
/**/
/**/
void MaxRowColumnValue(const double* const Mat,double* const Val,int* const Index,int row,int column,int type)
{
	if(Mat==NULL||Val==NULL||Index==NULL)
		return ;
	else
	{
		double Max=-1000000.0f;
		int idx;
		int tmpi=0;
		double *doutmp;
		if(type==1)
		{
			for(int i=0;i<row;i++)
			{	
				doutmp=(double*)(Mat+i*column);
				Max=doutmp[0];
				idx=0;
				for(int j=1;j<column;j++)
					if(doutmp[j]>Max)
					{
					Max=doutmp[j];
					idx=j;
					}

				Val[i]=Max;
				Index[i]=idx;
			}
			return;
		}
		else if(type==2)
		{
			for(int i=0;i<column;i++)
			{
				Val[i]=Mat[i];
				Index[i]=0;
			};
			for(int i=1;i<row;i++)
			{	
				doutmp=(double*)(Mat+i*column);
				for(int j=0;j<column;j++)
					if(doutmp[j]>Val[j])
					{
						Val[j]=doutmp[j];
						Index[j]=i;
					}
			}
		}
		return;
	}
}
/**/
/**/
void MaxRowColumnAbsValue(const double* const Mat,double* const Val,int* const Index,int row,int column,int type)
{
	if(Mat==NULL||Val==NULL||Index==NULL)
		return ;
	else
	{
		double Max=-1000000.0f;
		int idx;
		int tmpi=0;
		double *doutmp;
		if(type==1)
		{
			for(int i=0;i<row;i++)
			{	
				doutmp=(double*)(Mat+i*column);
				Max=doutmp[0];
				idx=0;
				for(int j=1;j<column;j++)
					if((doutmp[j]*doutmp[j])>=(Max*Max))
					{
						Max=doutmp[j];
						idx=j;
					}
				Val[i]=Max;
				Index[i]=idx;
			}
			return;
		}
		else if(type==2)
		{
			for(int i=0;i<column;i++)
			{
				Val[i]=(Mat[i]);
			};
			for(int i=1;i<row;i++)
			{	
				doutmp=(double*)(Mat+i*column);
				for(int j=0;j<column;j++)
					if(doutmp[j]*doutmp[j]>=Val[j]*Val[j])
					{
						Val[j]=doutmp[j];
						Index[j]=i;
					}
			}
		}
		return;
	}
}
/**/
/**/
void MinRowColumnValue(const double* const Mat,double* const Val,int* const Index,int row,int column,int type)
{
	if(Mat==NULL||Val==NULL||Index==NULL)
		return ;
	else
	{
		double Max=-1000000.0f;
		int idx;
		int tmpi=0;
		double *doutmp;
		if(type==1)
		{
			for(int i=0;i<row;i++)
			{	
				doutmp=(double*)(Mat+i*column);
				Max=doutmp[0];
				idx=0;
				for(int j=1;j<column;j++)
					if(doutmp[j]<Max)
					{
					Max=doutmp[j];
					idx=j;
					}
				Val[i]=Max;
				Index[i]=idx;
			}
			return;
		}
		else if(type==2)
		{
			for(int i=0;i<column;i++)
			{
				Val[i]=Mat[i];
				Index[i]=0;
			};
			for(int i=1;i<row;i++)
			{	
				doutmp=(double*)(Mat+i*column);
				for(int j=0;j<column;j++)
					if(doutmp[j]<Val[j])
					{
						Val[j]=doutmp[j];
						Index[j]=i;
					}
			}
			
		}
		return;
	}
}
/**/
/**/
void Transpose(double* const Mat,double* const Tam,int row)
{
	for(int i=0;i<row;i++)
		for(int j=0;j<row;j++)
			Tam[j*row+i]=Mat[i*row+j];
}
/**/
/**/
void DisplayMatrix(const double* const Mat,int row,int column)
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
			cout<<Mat[i*column+j]<<"\t";
		cout<<endl;
	}
	cout<<"\t<-------------->\t"<<endl;
};
/**/
/**/
void RowSum(const double* const Mat,double* const _Mat,int row,int column)
{
	if(Mat==NULL||_Mat==NULL)
		return ;
	else
	{
		int i,j;
		for(i=0;i<row;i++)
		{		
			_Mat[i]=0.0f;	
			for(j=0;j<column;j++)
			{
			_Mat[i]+=Mat[i*column+j];
			}
		}
	}
}
void MatrixSymmetryMulti(double* const Res,const double* Mat,const double alpha,const int row,const int column,int p,int q)
{
	double* tmpmem=new double[2*row];
	double pp=cos(alpha);
	double qq=cos(alpha);
	double pq=0-sin(alpha);
	double qp=sin(alpha);
	int i=0;
	int j=0;
	for(i=0;i<row;i++)
		for(j=0;j<column;j++)
			Res[i*column+j]=Mat[i*column+j];
	for(i=0;i<row;i++)
	{
		Res[i*column+q]=Mat[i*column+p]*pq+Mat[i*column+q]*qq;
		Res[i*column+p]=Mat[i*column+p]*pp+Mat[i*column+q]*qp;
	}
	for(i=0,j=column;i<column;i++,j++)
	{
		tmpmem[i]=pp*Res[p*column+i]+qp*Res[q*column+i];
		tmpmem[j]=pq*Res[p*column+i]+qq*Res[q*column+i];
	}
	for(i=0,j=column;i<column;i++,j++)
	{
		Res[p*column+i]=tmpmem[i];
		Res[q*column+i]=tmpmem[j];
	}
}
/**/
void MatrixRightRotate(double* const Res,const double* const Mat,const double alpha,const int row,const int column,const int p,const int q)
{
	double pp=cos(alpha);
	double qq=cos(alpha);
	double pq=0-sin(alpha);
	double qp=sin(alpha);
	int i=0;
	int j=0;
	for(i=0;i<row;i++)
		for(j=0;j<column;j++)
			Res[i*column+j]=Mat[i*column+j];
	for(i=0;i<row;i++)
	{
		Res[i*column+q]=Mat[i*column+p]*pq+Mat[i*column+q]*qq;
		Res[i*column+p]=Mat[i*column+p]*pp+Mat[i*column+q]*qp;
	}
	/**
	for(i=0,j=column;i<column;i++,j++)
	{
		tmpmem[i]=pp*Res[p*column+i]+qp*Res[q*column+i];
		tmpmem[j]=pq*Res[p*column+i]+qq*Res[q*column+i];
	}
	for(i=0,j=column;i<column;i++,j++)
	{
		Res[p*column+i]=tmpmem[i];
		Res[q*column+i]=tmpmem[j];
	}
	/**/
}
/**/
/**/
void MatrixLeftRotate(double* const Res,const double* const Mat,const double alpha,const int row,const int column,const int p,const int q)
{
	double pp=cos(alpha);
	double qq=cos(alpha);
	double pq=0-sin(alpha);
	double qp=sin(alpha);
	int i=0;
	int j=0;
	for(i=0;i<row;i++)
		for(j=0;j<column;j++)
			Res[i*column+j]=Mat[i*column+j];
	for(i=0;i<column;i++)
	{
		Res[p*column+i]=pp*Mat[p*column+i]+qp*Mat[q*column+i];
		Res[q*column+i]=pq*Mat[p*column+i]+qq*Mat[q*column+i];
	}
}

/**/
/**/
/*
void DisplayMatrix(double* Mat,int row,int subrow)
{
	for(int i=0;i<row;i++)
	{
		int j=0;
		int k=row/subrow;
		int s=row-k*subrow;
		while(subrow>=-1)
		{
			subrow--;
		}
		cout<<endl;
	}
	cout<<"\t<-------------->\t"<<endl;
};*/
/**/
/**/
void Jacobi(const double* const Mat,double* const matV,double* const matD,const int row)
{
	if(Mat==NULL||matV==NULL||matD==NULL||!MatrixIsSymmetry(Mat,row,row))
		return ;
	else
	{
		SYSTEMTIME sys_time;
 
		//将变量值设置为本地时间
		GetLocalTime( &sys_time );
		cout<<"start time is\t"<<sys_time.wMinute<<":\t"<<sys_time.wSecond<<":\t"<<sys_time.wMilliseconds<<endl;
		double *tmpD=new double[row*row];
		double *tmpV=new double[row*row];
		double *tmpU=new double[row*row];
		
		MatrixSetOne(tmpD,row,row);
		MatrixSetOne(tmpV,row,row);
		MatrixSetOne(tmpU,row,row);
		MatrixSetOne(matD,row,row);
		MatrixSetOne(matV,row,row);

		int maxcoor[2]={0,0};

		double maxvalue=MaxValue(Mat,row,row,maxcoor);
		double phi=0.0f;
		double douvar;
		const double *tmpMat=Mat;
		double *tmpMatD=matD;
		double *tmpMatV=matV;
		double *tmpPtD=tmpD;
		double *tmpPtV=tmpV;

		double *tmp=NULL;
		int time=0;
		while(maxvalue>_MaxDouble)
		{
			douvar=tmpMat[maxcoor[0]*row+maxcoor[0]]-tmpMat[maxcoor[1]*row+maxcoor[1]];
			if(douvar>=0&&douvar<_MaxDouble)
				douvar=_MaxDouble;
			else if(douvar<0&&douvar>_MinDouble)
				douvar=_MinDouble;
			phi=atan((2*tmpMat[maxcoor[0]*row+maxcoor[1]])/(douvar))/2;
			/**/
			/**
			MatrixSetOne(tmpU,row,row);
			//DisplayMatrix(tmpU,row);
			tmpU[maxcoor[0]*row+maxcoor[0]]=cos(phi);
			tmpU[maxcoor[1]*row+maxcoor[1]]=cos(phi);
			tmpU[maxcoor[0]*row+maxcoor[1]]=0-sin(phi);
			tmpU[maxcoor[1]*row+maxcoor[0]]=sin(phi);
			
			//DisplayMatrix(tmpU,row);
			
			Transpose(tmpU,tmpV,row);

			//DisplayMatrix(tmpV,row);
			/**/
			/**
			MatrixMultiMatrix(tmpMatD,tmpMat,tmpU,row,row,row,row);
			MatrixMultiMatrix(tmpD,tmpV,tmpMatD,row,row,row,row);
			MatrixMultiMatrix(tmpV,tmpMatV,tmpU,row,row,row,row);
			/**/
			//MatrixSymmetryMulti(tmpD,tmpMat,phi,row,row,maxcoor[0],maxcoor[1]);
			MatrixRightRotate(tmpMatD,tmpMat,phi,row,row,maxcoor[0],maxcoor[1]);
			MatrixLeftRotate(tmpD,tmpMatD,phi,row,row,maxcoor[0],maxcoor[1]);
			MatrixRightRotate(tmpPtV,tmpMatV,phi,row,row,maxcoor[0],maxcoor[1]);			
			/**/
			tmp=tmpPtV;
			tmpPtV=tmpMatV;
			tmpMatV=tmp;
			/**/

			tmpMat=tmpD;
			maxvalue=MaxValue(tmpMatV,row,row,maxcoor);
			cout<<"Maxvalue is "<<maxvalue<<"\t Run Time is "<<time++<<endl;
			maxvalue=MaxValue(tmpMat,row,row,maxcoor);
			cout<<"Maxvalue is "<<maxvalue<<"\t Run Time is "<<time++<<endl;
			
			//getchar();
		}
		for(int i=0;i<row;i++)
		{			
			for(int j=0;j<row;j++)
			{
				if(tmpD[i*row+j]>_MinDouble&&tmpD[i*row+j]<_MaxDouble)
					matD[i*row+j]=0.0f;
				else
					matD[i*row+j]=tmpD[i*row+j];

				if(tmpMatV[i*row+j]>_MinDouble&&tmpMatV[i*row+j]<_MaxDouble)
					matV[i*row+j]=0.0f;
				else
					matV[i*row+j]=tmpMatV[i*row+j];
			}
		}
		for(int j=0;j<row;j++)
		{
			matD[j]=matD[j*row+j];
			matD[row+j]=matD[j];
			matD[(row*(row-1))+j]=matD[j];
		}
	
		GetLocalTime( &sys_time );
		cout<<"end time is\t"<<sys_time.wMinute<<":\t"<<sys_time.wSecond<<":\t"<<sys_time.wMilliseconds<<endl;
		CopyMatrixToFile(matV,row,row,"javobi_vector");
		CopyMatrixToFile(matD,row,row,"javobi_values");
			cout<<endl;
			//clear memory 
			delete []tmpD;
			delete []tmpV;
			delete []tmpU;

	}
}
/**/
/**/

/**/
/**/
double MatrixDistance(const double* const Fir,const double* const Sec,int FirRow,int SecRow,int Column)
{
	double tmp=0.0f;
	int tmpi=FirRow*Column;
	int tmpj=SecRow*Column;
	for(int i=0;i<Column;i++)
	{
		tmp+=((Fir[tmpi+i]-Sec[tmpj+i])*(Fir[tmpi+i]-Sec[tmpj+i]));
	}
	if(tmp<_MaxDouble)
		return 0.0f;
	else
		return sqrt(tmp);
}
/**/
/**/

/**/
/**/
void Kmeans(const double* const Mat,int* clusteridx,int row,int column,int k,char name[])
{
	int tmpi=0,tmpk=0,tmpj=0;
	double* tmpPointer;
	double* curCenter;
	double* curOldCenter;
	//
	tmpk=(row>column)?row:column;
	//the k Max value 
	double* _Max=new double[tmpk];
	//the k Max index
	int* _MaxIndex=new int[tmpk];
	//the sum weight
	double* MatSumWeight=new double[tmpk];
	//the k Min value
	double* _Min=new double[tmpk];
	//the k Min index
	int* _MinIndex=new int[tmpk];	
	//the k cluster center coordinate 
	double* Center=new double[k*column];
	//the old clust center coordinate
	double* OldCenter=new double[k*column];
	//the data distance to the cluster center
	double* CenterDist=new double[row*k];
	//the weight Matrix
	double* MatWeight=new double[row*k];
	//init all parameter to zero;
	MatrixInitZeros(_Max,tmpk,1);
	MatrixInitZeros(_MaxIndex,tmpk,1);
	MatrixInitZeros(_Min,tmpk,1);
	MatrixInitZeros(_MinIndex,tmpk,1);
	MatrixInitZeros(Center,k,column);
	MatrixInitZeros(OldCenter,k,column);
	MatrixInitZeros(CenterDist,row,k);
	MatrixInitZeros(MatWeight,row,k);

	//init k cluster center Max and Min value and index 
	/**/
	MaxRowColumnValue(Mat,_Max,_MaxIndex,row,column,2);
	MinRowColumnValue(Mat,_Min,_MinIndex,row,column,2);
	//init k cluster center 
	/**/
	for(tmpk=0;tmpk<k;tmpk++)
	{
		for(tmpj=0;tmpj<column;tmpj++)
		{
			OldCenter[tmpk*column+tmpj]=(Mat[tmpk*column+tmpj]);
			Center[tmpk*column+tmpj]=OldCenter[tmpk*column+tmpj];
			//cout<<Center[tmpk*column+tmpj];
		}
		cout<<endl;
	}
	/**/
	//calculate the distance to the center
	curCenter=Center;
	curOldCenter=OldCenter;
	double MaxCondition=100000.0f;
	bool over=true;
	cout<<"while "<<endl;
	int time=1000;
	while(MaxCondition>_MaxDouble&&time>0)
	{
		cout<<"go on 1"<<endl;
		for(tmpi=0;tmpi<row;tmpi++)
		{
			for(tmpk=0;tmpk<k;tmpk++)
			{
				CenterDist[tmpi*k+tmpk]=MatrixDistance(Mat,curOldCenter,tmpi,tmpk,column);
				//cout<<CenterDist[tmpi*k+tmpk]<<"\t";
			}
			//cout<<endl;
		}

		//select the min distance to the cluster center
		MinRowColumnValue(CenterDist,_Min,_MinIndex,row,k,1);
		//MatWeight update 
			/**
			for(tmpi=0;tmpi<row;tmpi++)
			{	
				MatWeight[tmpi*k+_MinIndex[tmpi]]=_Min[tmpi];	
			}
			/**/
			//calculate the k cluster center weight
			//ColumnSum(MatWeight,MatSumWeight,row,k,0);
			//update center position
			//zeros cluster center 
		MatrixInitZeros(curCenter,k,column);
		MatrixInitZeros(MatSumWeight,k,1);
		//sum all data to update cluster center 
		for(tmpi=0;tmpi<row;tmpi++)
		{
			for(tmpj=0;tmpj<column;tmpj++)
			{
				curCenter[_MinIndex[tmpi]*column+tmpj]+=(_Min[tmpi]*Mat[tmpi*column+tmpj]);
			}
			MatSumWeight[_MinIndex[tmpi]]+=_Min[tmpi];
		}
		cout<<"MatSumWeight"<<endl;
		for(tmpk=0;tmpk<k;tmpk++)
		{
			//cout<<MatSumWeight[tmpk]<<endl;
			if(MatSumWeight[tmpk]<_MaxDouble)
				MatSumWeight[tmpk]=1000000.0f;
		}
		cout<<"go on 2"<<endl;
		//norm the cluster center
		for(tmpk=0;tmpk<k;tmpk++)
			for(tmpj=0;tmpj<column;tmpj++)
			{
				curCenter[tmpk*column+tmpj]/=MatSumWeight[tmpk];
			}
			//
		//swap new and old center
		tmpPointer=curOldCenter;
		curOldCenter=curCenter;
		curCenter=tmpPointer;
		//condition calculate
		cout<<"go on 3"<<endl;
		MaxCondition=-1000000.0f;		
		for(tmpk=0;tmpk<k;tmpk++)
		{
			_Max[tmpk]=MatrixDistance(curOldCenter,curCenter,tmpk,tmpk,column);
			if(_Max[tmpk]>MaxCondition)
				MaxCondition=_Max[tmpk];
		}
		//cout<<"the cluster index id is "<<endl;
		//for(tmpi=0;tmpi<row;tmpi++)
		//	cout<<_MinIndex[tmpi]<<"\t";
		cout<<"the biggest values is \t"<<MaxCondition<<endl;
		//getchar();
		time--;
	}
	 ofstream hFile(name, ios::out, _SH_DENYRW); 
	for(tmpi=0;tmpi<row;tmpi++)
	{
		//cout<<_MinIndex[tmpi]<<endl;
		hFile<<_MinIndex[tmpi]<<endl;
		clusteridx[tmpi]=_MinIndex[tmpi];
	}
	hFile.close();
}
/**/
/**/
template <class T> void AdjustVectorByIdx(T* matV,int* idx,int row,int column)
{
	T tmp;
	int tmpi;
	int vari;
	int* tmpidx=new int[column];
	for(int i=0;i<column;i++)
	{
		tmpidx[idx[i]]=i;
		//cout<<idx[i]<<"\t";
	}
	//cout<<endl;
	/**
	cout<<"resize data :\t";
	for(int i=0;i<column;i++)
	{
		cout<<tmpidx[i]<<"\t";
	}
	cout<<endl;
	/**/
	/**/
	for(int tsi=column-1;tsi>=0;tsi--)
	{
		tmpi=idx[tsi];
		if(tmpi!=tsi)
		{
			idx[tmpidx[tsi]]=tmpi;
			tmpidx[tmpi]=tmpidx[tsi];
			vari=tsi;
			for(int tsj=0;tsj<row;tsj++)
			{
				tmp=matV[tmpi];
				matV[tmpi]=matV[vari];
				matV[vari]=tmp;
				tmpi+=column;
				vari+=column;
			}
			idx[tsi]=tsi;
		}
	}
	/**
	for(int tsi=column-1;tsi>=0;tsi--)
		for(int tsj=column-1;tsj>=0;tsj--)
		{
			
		}
	/**/
	delete []tmpidx;
}
/**/
/**/

/**/
/**/
void MainSpectralCluster(double* Mat,int row,int k)
{
	if(!MatrixIsSymmetry(Mat,row,row)||Mat==NULL||k>=row)
		return ;
	else
	{	
		double* matD=new double[row*row];
		double* matV=new double[row*row];
		double* matTmp=new double[row];
		int* sortidx=new int[row];
		int tmpi,tmpk;
		int vari,varj;
		int dk=30;
		double tmpDouble=0.0f;
		CopyMatrixToFile(Mat,row,row,"Mat.src");
		//
		cout<<"ColumnSum(Mat,matTmp,row,row,0)"<<endl;
		//
		RowSum(Mat,matTmp,row,row);
		CopyMatrixToFile(matTmp,row,1,"Mat.sum");
		//
		for(tmpi=0;tmpi<row*row;tmpi++)
		{
			Mat[tmpi]=0.0-Mat[tmpi];
		}
		//
		for(tmpi=0;tmpi<row;tmpi++)
		{
			Mat[tmpi*row+tmpi]+=matTmp[tmpi];
		}
		//
		CopyMatrixToFile(Mat,row,row,"Mat.L");
		//
		cout<<"Jacobi(Mat,matV,matD,row);"<<endl;
		if(CopyFileToMatrix(matV,row,row,"Mat.matv")&&CopyFileToMatrix(matD,row,row,"Mat.matdbefore"));
		else
		{
			Jacobi(Mat,matV,matD,row);
			CopyMatrixToFile(matV,row,row,"Mat.matv");
			CopyMatrixToFile(matD,row,row,"Mat.matdbefore");
		}
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<row;j++)
			cout<<matV[i*row+j]<<" ";
			cout<<endl;
		}
		//cout<<"The Matrix is :"<<endl;
		//DisplayMatrix(testdata,100);
		//cout<<"The Matrix V is :"<<endl;
		//DisplayMatrix(testV,100);
		//cout<<"The Matrix D is :"<<endl;
		//DisplayMatrix(testD,100);
		//cout<<"Random matrix "<<endl;
		/**/
		cout<<"Heap sort "<<endl;
		HeapSort(matD,sortidx,row);

		CopyMatrixToFile(sortidx,row,1,"Mat.matdsort");
		CopyMatrixToFile(matD,row,row,"Mat.matd");	
		CopyMatrixToFile(matV,row,row,"Mat.matV_");	
		//
		//
		cout<<"	AdjustVectorByIdx(matV,sortidx,row);"<<endl;
		//
		AdjustVectorByIdx(matV,sortidx,row,row);
		CopyMatrixToFile(matV,row,row,"MatV.matdafter");
		/**/
		//
		cout<<"	For select the min k ;"<<endl;
		vari=0;
		tmpDouble=matV[row-1];
		tmpDouble=(1.0/tmpDouble);
		for(tmpi=0;tmpi<row;tmpi++)
		{
			varj=tmpi*row+row-dk-1;
			for(tmpk=0;tmpk<dk;tmpk++)
			{
				matD[vari]=matV[varj]*tmpDouble;
				varj++;
				vari++;
			}
		}
		//
		CopyMatrixToFile(matD,row,dk,"Mat.matdafter");

		cout<<"	Kmeans(matV,row,k,k);;"<<endl;
		Kmeans(matD,sortidx,row,dk,k,"spectral_30_norm.cluster");
		//
		/**/
		//clear memory
		delete []matV;
		delete []matD;
		delete []sortidx;
		//
	}
}
/**/
/**/
void SpectralCluster::NormalizedSpectralCluster(int* const idx,int row,int cok,int dk)
{
	if(cok>=row)
		return ;
	else
	{	
		double* matTmp=new double[row*row];
		double* matSum=new double[row];
		int i,k,j;
		int vari,varj;
		int coor[2];
		double tmpdou=0.0f;
		double *Mat=this->alldist;
		double* matV=this->affinityegi;
		MatrixInitZeros(matSum,1,row);
		MatrixInitZeros(matTmp,row,row);
		MaxValue(Mat,row,row,coor);
		double threshold=Mat[coor[0]*row+coor[1]]*0.4;
		for(i=0;i<row;i++)
		{
			for(j=0;j<row;j++)
			{
			tmpdou=Mat[i*row+j];
			if(tmpdou<threshold)
				matSum[j]+=tmpdou;
			}
		}

		for(i=0;i<row;i++)
		{
			matSum[i]=sqrt(matSum[i]);
		}
		
		for(i=0;i<row;i++)
		{
			k=i*row;
			vari=dk*i;
			for(j=(row-dk-1);j<(row-1);j++,k++)
			{
				matTmp[vari]=matV[k]*matSum[i];
			}
		}
		/**/
		cout<<"	Kmeans(matV,row,k,k);;"<<endl;
		Kmeans(matTmp,idx,row,dk,cok,"Norm_spectral_30.cluster");
		//
		/**/
		//clear memory
		delete []matSum;
		delete []matTmp;

		//
	}
}
/**/
/**/
SpectralCluster::SpectralCluster()
{
		this->caltype=-1;
};
SpectralCluster::~SpectralCluster()
{

};


void SpectralCluster::Main()
{

		int* idx=this->cur_clusts;	
		clock_t start,end; 
		double duration;

		start=clock();
		NormalizedSpectralCluster(idx, (this->cellnum) , (this->cur_num) , 5);
		//
		end=clock();
		duration = (double)(end - start) / CLOCKS_PER_SEC; 
		writePara("Spectral Cluster - all - num - time  ",(double)this->cellnum,(double)this->cur_num,duration,0.0);


}; 
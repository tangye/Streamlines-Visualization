 /************
 *  mex interface to compute the gradient of the eigenvectors 
 *  alignment quality
 *
 *  To mexify:   
 *               mex  evrot.cpp;
 *   
 *  [clusters,Quality,Vrot] = evrot(V,method);
 *   
 *  Input:
 *    V = eigenvecors, each column is a vector
 *    method = 1   gradient descent
 *             2   approximate gradient descent
 *
 *  Output:
 *    clusts - Resulting cluster assignment
 *    Quality = The final quality
 *    Vr = The rotated eigenvectors
 *
 * 
 *  Lihi Zelnik (Caltech) March.2005
 * 
 * 
 ************/
#include <Windows.h>
#include <math.h>
/**
extern "C"
{
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"
};
	/**/
#include "SelfTuningCluster.h"
#include "SpectralCluster.h"
/**/
void buildA(double *X,double *U1,double *Vk,double *U2,double* ret,int row,int dim)
{
    //mxArray  *lhs[1], *rhs[2];
    double* lhs=new double[dim*dim];
    double* rhs=new double[dim*dim];
    MatrixInitZeros(lhs,dim,dim);
	MatrixInitZeros(rhs,dim,dim);
	MatrixInitZeros(ret,row,dim);                                                                                                         
    MatrixMultiMatrix(lhs,Vk,U2,dim,dim,dim,dim);
	MatrixMultiMatrix(rhs,U1,lhs,dim,dim,dim,dim);
    MatrixMultiMatrix(ret,X,rhs,row,dim,dim,dim);
	delete []lhs;
	delete []rhs;
};
/**/
/**/
void build_Uab(double *theta, int a, int b,const int* ik, const int* jk, int dim,double* ret)
{        
    double *p_Uab = ret;
    int ind, k,i,j,ind_ik,ind_jk;
	MatrixInitZeros(ret,dim,dim);
    /* set Uab to be an identity matrix */
    for( j=0; j<dim; j++ )
	{
        ind = dim*j + j;
        p_Uab[ind] = 1.0;
    }        
    
    if( b < a ) 
	{
        return ;
    }
    double *p_theta = theta;
    double tt,u_ik;
    for( k=a; k<=b; k++ )
	{
        tt = p_theta[k];
        for( i=0; i<dim; i++ )
		{
            ind_ik = ik[k]+i*dim;
            ind_jk = dim*i+jk[k];
            u_ik = p_Uab[ind_ik] * cos(tt) - p_Uab[ind_jk] * sin(tt);
            p_Uab[ind_jk] = p_Uab[ind_ik] * sin(tt) + p_Uab[ind_jk] * cos(tt);
            p_Uab[ind_ik] = u_ik;
        }                       
    }
};
/**/
/**/
void rotate_givens(double *X, double *theta, const int* ik, const int* jk, int angle_num,int row,int dim,double* ret)
{
	double* rot=new double[dim*dim];
	MatrixInitZeros(rot,dim,dim);
    build_Uab(theta, 0, angle_num-1,ik,jk,dim,rot);

    MatrixMultiMatrix(ret,X,rot,row,dim,dim,dim);

   	delete []rot;
};
/**/
/**/
void gradU(double *theta, int k,const int* ik,const int* jk, int dim,double* ret)
{
    
    double *p_theta =(theta);
   	double *p_V = ret; 
    p_V[ik[k]+dim*ik[k]] = -sin(p_theta[k]);
    p_V[jk[k]+dim*ik[k]] = cos(p_theta[k]);
    p_V[ik[k]+dim*jk[k]] = -cos(p_theta[k]);
    p_V[jk[k]+dim*jk[k]] = -sin(p_theta[k]);
};
/**/
/**/
double evqualitygrad(double *X, double* theta,const int *ik,const int *jk,int angle_num,int angle_index,int dim,int ndata)
{                   
    /* build V,U,A */
    double* matret=new double[dim*dim];
    double* U1=new double[dim*dim];
    double* U2=new double[dim*dim];
    double* A=new double[ndata*dim];
    double* matrot=new double[ndata*dim]; 
	double *max_values = new double[ndata];
    int *max_index = new int[ndata];
	double dJ=0, tmp1, tmp2;
	double* p_A;
	double *p_Y;
    /* find max of each row */
	MatrixInitZeros(matret,dim,dim);
	MatrixInitZeros(U1,dim,dim);
	MatrixInitZeros(U2,dim,dim);
	MatrixInitZeros(A,ndata,dim);
	MatrixInitZeros(matrot,ndata,dim);
	MatrixInitZeros(max_values,ndata,1);
	MatrixInitZeros(max_index,ndata,1);
    int i,j, ind = 0;
	
	//getchar();
    gradU(theta,angle_index,ik,jk,dim,matret);
    /**/

	//getchar();
    build_Uab(theta,0,angle_index-1,ik,jk,dim,U1);
	/**/

	/**/
    build_Uab(theta,angle_index+1,angle_num-1,ik,jk,dim,U2);


    buildA(X,U1,matret,U2,A,ndata,dim);
   
	
    /* rotate vecs according to current angles */   
    rotate_givens(X,theta,ik,jk,angle_num,ndata,dim,matrot);
	p_Y=matrot; 

	MaxRowColumnAbsValue(p_Y,max_values,max_index,ndata,dim,1);

    /* compute gradient */
    ind = 0;
	dJ=0;
	for( i=0; i<ndata; i++ )
	{ /* loop over all rows */
		p_A=(double*)(A+i*dim);
		p_Y=(double*)(matrot+i*dim);
		for( j=0; j<dim; j++ )
		{  /* loop over all columns */

            tmp1 = p_A[j] * p_Y[j] / (max_values[i]*max_values[i]);
            tmp2 = p_A[max_index[i]]*(p_Y[j]*p_Y[j])/(max_values[i]*max_values[i]*max_values[i]);
            dJ += tmp1-tmp2;
        }
    }
    dJ = 2*dJ/ndata/dim;
    delete []max_values;
    delete []max_index;
   	delete []matrot;
   	delete []matret;
   	delete []U2;
   	delete []U1;
    delete []A;
    return dJ;
};
/**/
/********** cluster assignments ************/  
/**/
void cluster_assign(double *X,int* clusts,int dim,int ndata)
{
    double *p_X = X;
    /* take the square of all entries and find max of each row */
    int i,j,index;
	double maxvalue,tmp;
	for(i=0;i<ndata;i++)
	{
		p_X=(double*)(X+i*dim);
		maxvalue=p_X[0]*p_X[0];
		index=0;
		for(j=1;j<dim;j++)
		{
			tmp=p_X[j]*p_X[j];
			if(tmp>=maxvalue)
			{
				index=j;
				maxvalue=tmp;
			}
		}
		clusts[i]=index;
	}
	
};

/**/
/**/
double evqual(double *X,int dim,int ndata)
{
    double *p_X = X;
    /* take the square of all entries and find max of each row */
    double *max_values = (double*)new double[ndata];
    int *max_index = (int*)new int[ndata];
    int i,j,ind = 0; 
	MatrixInitZeros(max_values,1,ndata);
	
 //   for( i=0; i<ndata; i++ )
	//{ /* loop over all rows */
	//	ind=i*dim;
	//	max_values[i]=p_X[ind]*p_X[ind];
	//	max_index[i]=0;
	//	for( j=1; j<dim; j++ )
	//	{  /* loop over all columns */
	//		ind++;
 //           if( max_values[i] <= p_X[ind]*p_X[ind]  )
	//		{                
 //               max_values[i] = p_X[ind]*p_X[ind];
 //               max_index[i] = j;
 //           }
 //       }
	//}
	//cout<<"evqual dim "<<dim<<"\t"<<ndata<<endl;
	/*
	for(int i=0;i<10;i++)
	{
			for(int j=0;j<dim;j++)
				//cout<<X[i*dim+j]<<"\t";
			//cout<<endl;
	}
	*/
	MaxRowColumnAbsValue(p_X,max_values,max_index,ndata,dim,1);
	for(i=0;i<ndata;i++)
		max_values[i]*=max_values[i];
	//{
	//	//cout<<"Max id \t"<<i<<"\t"<<max_values[i]<<"\t";
		//max_values[i]*=max_values[i];
		////cout<<max_values[i]<<"\t";
	//}
	/* compute cost *
	//cout<<"for max value grad\n"<<endl;
       for(int i=0;i<ndata;i++)
       {
           //cout<<max_values[i]<<"\t";  
		   if(i%5==4)       
                //cout<<endl;
        }
		/**/
    double J=0;
    for( i=0; i<ndata; i++ )
	{
		ind=i*dim;
		for( j=0; j<dim; j++ )
		{  /* loop over all columns */
			J += p_X[ind]*p_X[ind]/max_values[i];
			ind++;
		}
		
	}

    J = 1.0 - (J/ndata -1.0)/dim;

 	delete []max_values;
    delete []max_index;
    return J;
};
/**/
/**/
double evrot(double* mat,int dim,int ndata,int* out_clusts,double* out_xrot,int method)
{
		int tmpi=0;
		int tmpk=0;
		int tmpj=0;
		int vari=0;
		int varj=0;
		int dk=5;
		int max_iter = 200;    
		double dQ=0;
		double Q=0;
		double Q_new=0;
		double Q_old1=0;
		double Q_old2=0;
		double Q_up=0;
		double Q_down=0;
		double alpha=0;
		int iter=0;
		int d=0;
		double tmpDouble=0;
		int angle_num=0;
		angle_num=(int)(dim*(dim-1)/2);
		double* theta=new double[angle_num];
		double* theta_new=new double[angle_num];
		int* jk=new int[angle_num];
		int* ik=new int[angle_num];
		double* p_theta=theta;
		double* p_theta_new=theta_new;
		tmpi=0;
		for(vari=0;vari<dim-1;vari++)
			for(varj=vari+1;varj<=dim-1;varj++)
			{
				ik[tmpi]=vari;
				jk[tmpi]=varj;
				tmpi++;
			}

		Q = evqual(mat,dim,ndata);

				///getchar();
		Q_old1 = Q;
		Q_old2 = Q;
		iter = 0;

		MatrixInitZeros(out_xrot,dim,ndata);
		MatrixInitZeros(p_theta,1,angle_num);
		MatrixInitZeros(p_theta_new,1,angle_num);
		
		while( iter < max_iter )
		{/* iterate to refine quality */
			iter++;
			cout<<"iter "<<iter<<endl;
			for( d = 0; d < angle_num; d++ )
			{
				if( method == 2 )
            	{ /* descend through numerical drivative */
					alpha = 0.1;
               		/* move up */                
               		p_theta_new[d] = p_theta[d] + alpha;
               		rotate_givens(mat,theta_new,ik,jk,angle_num,ndata,dim,out_xrot);
					Q_up = evqual(out_xrot,dim,ndata);
				   	MatrixInitZeros(out_xrot,ndata,dim);
					 /* move down */
              		p_theta_new[d] = p_theta[d] - alpha;
              		rotate_givens(mat,theta_new,ik,jk,angle_num,ndata,dim,out_xrot);
              		Q_down = evqual(out_xrot,dim,ndata);
              		MatrixInitZeros(out_xrot,ndata,dim);

	                /* update only if at least one of them is better */
    	            if( Q_up > Q || Q_down > Q)
    	            {
     	               if( Q_up > Q_down )
    	                { 
    	                    p_theta[d] = p_theta[d] + alpha;
    	                    p_theta_new[d] = p_theta[d];
    	                    Q = Q_up;
    	                }
						else 
						{
							p_theta[d] = p_theta[d] - alpha;
							p_theta_new[d] = p_theta[d];
						    Q = Q_down;
						}
					}
				} 
				else 
				{ /* descend through true derivative */
					alpha = 1.0;

					dQ = evqualitygrad(mat,theta,ik,jk,angle_num,d,dim,ndata);
					//cout<<"the dQ is \t "<<dQ<<endl;
					p_theta_new[d] = p_theta[d] - alpha * dQ;
               		rotate_givens(mat,theta_new,ik,jk,angle_num,ndata,dim,out_xrot);
					Q_new = evqual(out_xrot,dim,ndata);                                    
				    if( Q_new > Q)
					{
						p_theta[d] = p_theta_new[d];
						Q = Q_new;
					 }
					else
					{
						p_theta_new[d] = p_theta[d];
					}
               		MatrixInitZeros(out_xrot,dim,ndata);
				}
			}        
			/* stopping criteria */
			if( iter > 2 )
			{
				if( Q - Q_old2 < 1e-3 )
				{
		
					break;
				}
			 }

			Q_old2 = Q_old1;
			Q_old1 = Q;
    }
    rotate_givens(mat,theta_new,ik,jk,angle_num,ndata,dim,out_xrot);
    cluster_assign(out_xrot,out_clusts,dim,ndata);

    /* free allocated memory */

    delete [] theta;
    delete [] theta_new;
	delete [] ik;
    delete [] jk;
	return Q;
};
/**/
/**/
/**/
/**/
double cluster_rotate(const double* mat,int row,int column,int* groups,int group_size,int* out_clusts,unsigned int* out_num,double* out_vr,int method)
{
	if(mat==NULL)
		return 0;
	else
	{
		double* cur_vr=new double[row*row];
		double* matV=new double[row*row];
		int* mycur_clusts=new int[group_size*row];
		double* out_Qs=new double[group_size];
		int tmpi=0;
		int tmpj=0;
		int tmpk=0;
		int vari=0;
		int vcur_column=0;
		double* doup_mat=NULL;
		double* doup_tmp=NULL;
		double* p_curvr=NULL;
		int* p_curclu=NULL;
		int* intp_tmp=NULL;
		double out_Q=0;
		MatrixInitZeros(cur_vr,row,row);
		MatrixInitZeros(matV,row,row);	
		vcur_column=groups[0];
		for(tmpj=0;tmpj<row;tmpj++)
		{
			mycur_clusts[tmpj]=0;
			doup_mat=(double*)(mat+tmpj*row);
			doup_tmp=(double*)(matV+tmpj*vcur_column);
			for(tmpi=0;tmpi<vcur_column;tmpi++)
			{
				doup_tmp[tmpi]=doup_mat[tmpi];
			}
		}

		p_curvr=out_vr;
		p_curclu=mycur_clusts;
		//cout<<"disp vcurr "<<endl;
		//DisplayMatrix(matV,10,vcur_column);
		out_Qs[0]=evrot(matV,vcur_column,row,p_curclu,p_curvr,method);
		out_Q=out_Qs[0];
		//cout<<"the "<<0<<"\tQuality is \t"<<out_Qs[tmpk]<<endl;
		for(tmpk=1;tmpk<group_size;tmpk++)
		{
				vari=groups[tmpk]-1;
				//cout<<"vari"<<vari<<endl;
				for(tmpj=0;tmpj<row;tmpj++)
				{
					doup_mat=(double*)(p_curvr+tmpj*vcur_column);
					doup_tmp=(double*)(matV+tmpj*(vcur_column+1));
					for(tmpi=0;tmpi<vcur_column;tmpi++)
						doup_tmp[tmpi]=doup_mat[tmpi];
					doup_tmp[tmpi]=mat[tmpj*row+vari];
				}
				vcur_column++;
				CopyFileToMatrix(doup_tmp,row,vcur_column,"vcur_vr");
				//cout<<"disp vcurr "<<endl;
				//DisplayMatrix(matV,10,vcur_column);
				p_curclu=(int*)(mycur_clusts+tmpk*row);
			 	out_Qs[tmpk]=evrot(matV,vcur_column,row,p_curclu,p_curvr,method);
				//cout<<"the "<<tmpk<<"\tQuality is \t"<<out_Qs[tmpk]<<endl;
		}
		p_curclu=mycur_clusts;
		CopyMatrixToFile(p_curclu,group_size,row,"all_out_clusts");
		tmpi=0;
		CopyMatrixToFile(out_Qs,group_size,1,"all_out_Quality");
		//cout<<"Out_Qs"<<endl;
		//cout<<out_Q<<"\t";
		for(tmpk=1;tmpk<group_size;tmpk++)
		{
			if(out_Q<out_Qs[tmpk])
			{
				out_Q=out_Qs[tmpk];
			}
			//cout<<out_Qs[tmpk]<<"\t";
		}
		//cout<<endl;
		for(tmpk=0;tmpk<group_size;tmpk++)
		{
			if((out_Q-out_Qs[tmpk])<=0.001)
			{
				tmpi=tmpk;
			}
		}
		/**/
		p_curclu=(int* )(mycur_clusts+tmpi*row);
		for(tmpk=0;tmpk<row;tmpk++)
			out_clusts[tmpk]=p_curclu[tmpk];
		*out_num=tmpi;
		/**/
		//clear memory
		delete []matV;
		delete []cur_vr;
		delete []mycur_clusts;
		delete []out_Qs;
		//
		return out_Q;
	}
};
/**/
/**/



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
#ifndef SELFTUNINGCLUSTER_h
#define SELFTUNINGCLUSTER_h
double cluster_rotate(const double* mat,int row,int column,int* groups,int group_size,int* out_clusts,unsigned int* out_num,double* out_vr,int method);
/**/
/**/
void buildA(double *X, double *U1, double *Vk, double *U2,double* ret,int row,int dim);
/**/
/**/
void build_Uab(double *theta, int a, int b,int* ik, int* jk, int dim,double* ret);
/**/
/**/
void gradU(double *theta, int k,int* ik, int* jk, int dim,double* ret);
/**/
/**/
double evqualitygrad(double *X, double* theta,const int *ik,const int *jk,int angle_num,int angle_index,int dim,int ndata);
/**/
/**/
void cluster_assign(double *X,int *ik, int *jk,int* clusts,int dim,int ndata);
/**/
/**/
void rotate_givens(double *X, double *theta, const int* ik, const int* jk, int angle_num,int row,int dim,double* ret);
/**/
/**/
double evqual(double *X,int *ik, int *jk,int dim,int ndata);
/**/
/**/
double evrot(double* mat,int row,int column,int* out_clusts,double* out_xrot,int method);
/**/
/**/
#endif
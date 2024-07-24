/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Decompositions

///////////////////////////

double cg_solve(double *mat, int length, double *mat2, double *mat3, int ncol, double tol)
{
//mat is matrix, size length x length
//mat2 stores the xks - at start contains initial guess, at end contains solution
//mat3 starts as rhs, will be used to store rs
//mat4 stores the pks
int i, j, count, count2;
double alpha, beta, value, diff, *mat4;
double *alphak, *betak, *rktrk, *rktrk2, *Apk;


if(length==0){return(0);}

//allocations
mat4=malloc(sizeof(double)*length*ncol);
alphak=malloc(sizeof(double)*ncol);
betak=malloc(sizeof(double)*ncol);
rktrk=malloc(sizeof(double)*ncol);
rktrk2=malloc(sizeof(double)*ncol);
Apk=malloc(sizeof(double)*length*ncol);

//initializations

//set rktrk2 to residual norms when guess is zero
for(j=0;j<ncol;j++)
{
rktrk2[j]=0;
for(i=0;i<length;i++){rktrk2[j]+=pow(mat3[i+j*length],2);}
}

//x0 already in mat2

//r0=b-Ax0 - saved in mat3 (which currently contains b)
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &length, &ncol, &length, &alpha, mat, &length, mat2, &length, &beta, mat3, &length);

//p0=r0 - saved in mat4
for(j=0;j<ncol;j++)
{
for(i=0;i<length;i++){mat4[i+j*length]=mat3[i+j*length];}
}

//get r0tr0s
for(j=0;j<ncol;j++)
{
rktrk[j]=0;
for(i=0;i<length;i++){rktrk[j]+=pow(mat3[i+j*length],2);}
}

//get number converged and difference
count2=0;for(j=0;j<ncol;j++){count2+=(fabs(rktrk[j]/rktrk2[j])<tol);}
for(j=0;j<ncol;j++)
{
if(j==0){diff=fabs(rktrk[j]/rktrk2[j]);}
if(fabs(rktrk[j]/rktrk2[j])>diff){diff=fabs(rktrk[j]/rktrk2[j]);}
}

for(count=0;count<length;count++)
{
//update Apk
alpha=1.0;beta=0.0;
dgemm_("N", "N", &length, &ncol, &length, &alpha, mat, &length, mat4, &length, &beta, Apk, &length);

//alphaks = rkTrk / pkT Apk
for(j=0;j<ncol;j++)
{
value=0;
for(i=0;i<length;i++){value+=mat4[i+j*length]*Apk[i+j*length];}
alphak[j]=rktrk[j]/value;
}

//xk+1 = xk + alphak pk
for(j=0;j<ncol;j++)
{
for(i=0;i<length;i++){mat2[i+j*length]+=alphak[j]*mat4[i+j*length];}
}

//rk+1 = rk - alphak Apk
for(j=0;j<ncol;j++)
{
for(i=0;i<length;i++){mat3[i+j*length]-=alphak[j]*Apk[i+j*length];}
}

//update rkTrk and get betas = rk+1Tr+1k / rkTrk
for(j=0;j<ncol;j++)
{
value=0;
for(i=0;i<length;i++){value+=pow(mat3[i+j*length],2);}
betak[j]=value/rktrk[j];
rktrk[j]=value;
}

//get number converged and difference
count2=0;for(j=0;j<ncol;j++){count2+=(fabs(rktrk[j]/rktrk2[j])<tol);}
for(j=0;j<ncol;j++)
{
if(j==0){diff=fabs(rktrk[j]/rktrk2[j]);}
if(fabs(rktrk[j]/rktrk2[j])>diff){diff=fabs(rktrk[j]/rktrk2[j]);}
}

if(count2==ncol){break;}

//pk+1 = rk + betak pk
for(j=0;j<ncol;j++)
{
for(i=0;i<length;i++){mat4[i+j*length]=mat3[i+j*length]+betak[j]*mat4[i+j*length];}
}
}	//end of while loop

printf("Took %d iterations, %d out of %d converged max %f\n", count+1, count2, ncol, diff);

return(0);
}

double eigen_invert(double *mat, int length, double *mat2, int ncol, double *mat3, int type)
{
//ncol=0 - get "cholesky", ncol=-1 - get inverse, ncol>0 - solve mat X = mat3
//ncol=0 - mat3 NULL, ncol=-1 - mat3 workspace size(mat), ncol>0 - mat3 = RHS
//type=0 - quiet, type=1 - complain
int i, j, count, count2, lwork, info;
double det, value, alpha, beta, wkopt, *work, *mat4;


if(length==0){return(0);}

lwork=-1;
dsyev_("V", "U", &length, mat, &length, mat2, &wkopt, &lwork, &info);
if(info!=0)
{printf("Error, eigen priming failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);

dsyev_("V", "U", &length, mat, &length, mat2, work, &lwork, &info);
if(info!=0)
{printf("Error, eigen decomp failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
free(work);

//get log determinant
det=0;count=0;
for(i=0;i<length;i++)
{
if(fabs(mat2[i])>=0.000001){det+=log(fabs(mat2[i]));}
else{det+=log(0.000001);}
if(mat2[i]<=-0.000001){count++;}
}
if(count>1&&type==1){printf("Warning, %d eigenvalues are negative\n\n", count);}

////////

if(ncol==0)	//solve XXT = mat, return in mat - X = UE^.5
{
for(j=0;j<length;j++)
{
if(mat2[j]>0.000001)
{
value=pow(mat2[j],.5);
for(i=0;i<length;i++){mat[(size_t)j*length+i]=mat[(size_t)j*length+i]*value;}
}
else
{
for(i=0;i<length;i++){mat[(size_t)j*length+i]=0;}
}
}
}

if(ncol==-1)	//solve X mat = I - therefore mat becomes UE^-1UT
{
//load U|E|^⁻.5 into mat3
for(j=0;j<length;j++)
{
if(fabs(mat2[j])>0.000001)
{
value=pow(fabs(mat2[j]),.5);
for(i=0;i<length;i++){mat3[(size_t)j*length+i]=mat[(size_t)j*length+i]/value;}
}
else
{
for(i=0;i<length;i++){mat3[(size_t)j*length+i]=0;}
}
}

//get number of positive eigenvalues (eigenvalues ranked in ascending order)
count2=length-count;

if(count2>0)	//deal with last count2 vectors (positive eigenvalues)
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &length, &length, &count2, &alpha, mat3+count*length, &length, mat3+count*length, &length, &beta, mat, &length);
}

if(count>0)	//deal with first count vectors (negative eigenvalues)
{
if(count2>0){alpha=-1.0;beta=1.0;}
else{alpha=-1.0;beta=0.0;}
dgemm_("N", "T", &length, &length, &count, &alpha, mat3, &length, mat3, &length, &beta, mat, &length);
}
}	//end of ncol=-1

if(ncol>0)	//solve mat X = mat3 - X = UE^-1UT mat3
{
mat4=malloc(sizeof(double)*length*ncol);
alpha=1.0;beta=0.0;
dgemm_("T", "N", &length, &ncol, &length, &alpha, mat, &length, mat3, &length, &beta, mat4, &length); 
for(i=0;i<length;i++)
{
for(j=0;j<ncol;j++)
{
if(fabs(mat2[i])>=0.000001){mat4[(size_t)j*length+i]=mat4[(size_t)j*length+i]/mat2[i];}
else{mat4[(size_t)j*length+i]=0;}
}
}
dgemm_("N", "N", &length, &ncol, &length, &alpha, mat, &length, mat4, &length, &beta, mat3, &length);
free(mat4);
}

return(det);
}	//end of eigen_invert

////////

/*
//the single version of the eigen decomposition seems unstable

float eigen_invert_single(float *mat, int length, float *mat2, int ncol, float *mat3, int type)
{
//ncol=0 - get "cholesky", ncol=-1 - get inverse, ncol>0 - solve mat X = mat3
//ncol=0 - mat3 NULL, ncol=-1 - mat3 workspace size(mat), ncol>0 - mat3 = RHS
//type=0 - quiet, type=1 - complain
int i, j, count, count2, lwork, info;
float det, value, alpha_single, beta_single, wkopt, *work, *mat4;


if(length==0){return(0);}

lwork=-1;
ssyev_("V", "U", &length, mat, &length, mat2, &wkopt, &lwork, &info);
if(info!=0)
{printf("Error, eigen priming failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(float)*lwork);

ssyev_("V", "U", &length, mat, &length, mat2, work, &lwork, &info);
if(info!=0)
{printf("Error, eigen decomp failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
free(work);

for(i=0;i<length;i++){printf("eigen %d is %f\n", i+1, mat2[i]);}

//get log determinant
det=0;count=0;
for(i=0;i<length;i++)
{
if(fabs(mat2[i])>=0.000001){det+=log(fabs(mat2[i]));}
else{det+=log(0.000001);}
if(mat2[i]<=-0.000001){count++;}
}
if(count>1&&type==1){printf("Warning, %d eigenvalues are negative\n\n", count);}

////////

if(ncol==0)	//solve XXT = mat, return in mat - X = UE^.5
{
for(j=0;j<length;j++)
{
if(mat2[j]>0.000001)
{
value=pow(mat2[j],.5);
for(i=0;i<length;i++){mat[(size_t)j*length+i]=mat[(size_t)j*length+i]*value;}
}
else
{
for(i=0;i<length;i++){mat[(size_t)j*length+i]=0;}
}
}
}

if(ncol==-1)	//solve X mat = I - therefore mat becomes UE^-1UT
{
//load U|E|^⁻.5 into mat3
for(j=0;j<length;j++)
{
if(fabs(mat2[j])>0.000001)
{
value=pow(fabs(mat2[j]),.5);
for(i=0;i<length;i++){mat3[(size_t)j*length+i]=mat[(size_t)j*length+i]/value;}
}
else
{
for(i=0;i<length;i++){mat3[(size_t)j*length+i]=0;}
}
}

//get number of positive eigenvalues (eigenvalues ranked in ascending order)
count2=length-count;

if(count2>0)	//deal with last count2 vectors (positive eigenvalues)
{
alpha_single=1.0;beta_single=0.0;
sgemm_("N", "T", &length, &length, &count2, &alpha_single, mat3+count*length, &length, mat3+count*length, &length, &beta_single, mat, &length);
}

if(count>0)	//deal with first count vectors (negative eigenvalues)
{
if(count2>0){alpha_single=-1.0;beta_single=1.0;}
else{alpha_single=-1.0;beta_single=0.0;}
sgemm_("N", "T", &length, &length, &count, &alpha_single, mat3, &length, mat3, &length, &beta_single, mat, &length);
}
}	//end of ncol=-1

if(ncol>0)	//solve mat X = mat3 - X = UE^-1UT mat3
{
mat4=malloc(sizeof(double)*length*ncol);
alpha_single=1.0;beta_single=0.0;
sgemm_("T", "N", &length, &ncol, &length, &alpha_single, mat, &length, mat3, &length, &beta_single, mat4, &length); 
for(i=0;i<length;i++)
{
for(j=0;j<ncol;j++)
{
if(fabs(mat2[i])>=0.000001){mat4[(size_t)j*length+i]=mat4[(size_t)j*length+i]/mat2[i];}
else{mat4[(size_t)j*length+i]=0;}
}
}
sgemm_("N", "N", &length, &ncol, &length, &alpha_single, mat, &length, mat4, &length, &beta_single, mat3, &length);
free(mat4);
}

return(det);
}	//end of eigen_invert_single
*/

////////

void eigen_strip(double *mat, int start, int end, int total, double strip)
{
int i, j, length, lwork, info;
double *U, *E, sum, sum2, value, alpha, beta, wkopt, *work;


length=end-start;
U=malloc(sizeof(double)*length*length);
E=malloc(sizeof(double)*length);
for(i=0;i<length;i++)
{
for(j=0;j<length;j++){U[(size_t)j*length+i]=mat[(size_t)(start+j)*total+(start+i)];}
}

lwork=-1;
dsyev_("V", "U", &length, U, &length, E, &wkopt, &lwork, &info);
if(info!=0)
{printf("Error, eigen priming failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &length, U, &length, E, work, &lwork, &info);
if(info!=0)
{printf("Error, eigen decomp failed; please tell Doug (info %d, length %d)\n\n", info, length);exit(1);}
free(work);

sum=0;for(i=0;i<length;i++){sum+=E[i];}
sum2=0;
for(i=length-1;i>=0;i--)
{
sum2+=E[i];
if(sum2/sum>1-strip){E[i]=0;}
}

for(j=0;j<length;j++)
{
value=pow(E[j],.5);
for(i=0;i<length;i++){U[(size_t)j*length+i]*=value;}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &length, &length, &length, &alpha, U, &length, U, &length, &beta, mat+start+start*total, &total);

free(U);free(E);
}

///////////////////////////

double ldlt_invert(double *mat, int length, int ncol, double *mat2, int *info, int type)
{
//ncol=-1 - get inverse, ncol>0 - solve mat X = mat2
//type=0 - quiet, type=1 - complain (and warn about eigendecomp)
int i, j, *ipiv, lwork;
double wkopt, *work, det;


if(length==0){return(0);}

//start with decomposition
ipiv=malloc(sizeof(int)*length);
lwork=-1;
dsytrf_("U", &length, mat, &length, ipiv, &wkopt, &lwork, info);
if(*info!=0)
{printf("Error, LDLT priming failed; please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsytrf_("U", &length, mat, &length, ipiv, work, &lwork, info);
free(work);

if(*info==0)	//decomposition worked
{
//get log determinant
det=0;
for(i=0;i<length;i++)
{
if(fabs(mat[(size_t)i*length+i])>=0.000001){det+=log(fabs(mat[(size_t)i*length+i]));}
else{det+=log(0.000001);}
}

if(ncol==-1)	//put inverse into mat
{
work=malloc(sizeof(double)*length);
dsytri_("U", &length, mat, &length, ipiv, work, info);
if(*info!=0)
{printf("Error LDLT invert failed, please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
for(i=0;i<length;i++)
{
for(j=0;j<i;j++){mat[(size_t)j*length+i]=mat[(size_t)i*length+j];}
}
free(work);
}

if(ncol>0)	//solve
{
dsytrs_("U", &length, &ncol, mat, &length, ipiv, mat2, &length, info);
if(*info!=0)
{printf("Error LDLT solve failed, please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
}
}
else	//decomposition failed
{
if(type==1)
{printf("Warning, LDLT failed (info %d, length %d); please tell Doug\n", *info, length);}
det=0;
}

free(ipiv);

return(det);
}

///////////////////////////

double cholesky_invert(double *mat, int length, int ncol, double *mat2, int *info, int type)
{
//ncol=-1 - get inverse, ncol>0 - solve mat X = mat2
//type=0 - quiet, type=1 - complain
int i, j;
double det;


if(length==0){return(0);}

//start with decomposition
dpotrf_("U", &length, mat, &length, info);

if(*info==0)	//decomposition worked
{
//get log determinant
det=0;
for(i=0;i<length;i++)
{
if(fabs(mat[(size_t)i*length+i])>=0.000001){det+=2*log(fabs(mat[(size_t)i*length+i]));}
else{det+=2*log(0.000001);}
}

if(ncol==-1)	//put inverse into mat
{
dpotri_("U", &length, mat, &length, info);
if(*info!=0)
{printf("Error Cholesky invert failed, please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
for(i=0;i<length;i++)
{
for(j=0;j<i;j++){mat[(size_t)j*length+i]=mat[(size_t)i*length+j];}
}
}

if(ncol>0)	//solve
{
dpotrs_("U", &length, &ncol, mat, &length, mat2, &length, info);
if(*info!=0)
{printf("Error Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", *info, length);exit(1);}
}
}
else	//decomposition failed
{
if(type==1){printf("Warning, Cholesky failed (info %d, length %d); please tell Doug\n", *info, length);}
det=0;
}

return(det);
}

///////////////////////////


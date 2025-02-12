/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Regression functions

///////////////////////////

double get_factor(double *Y, int ns, double prev, double ascer, char *outfile)
{
int i;
double sum, mean, value, factor;

char filename[500];
FILE *output;


if(ascer==-9999)	//get mean from phenotypes
{
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
mean=sum/ns;
}
else{mean=ascer;}

value=pow(2*M_PI,-.5)*exp(-.5*pow(normal_inv(1-prev),2));
factor=pow(prev*(1-prev),2)/mean/(1-mean)*pow(value,-2);
printf("Prevalence: %.6f, Ascertainment: %.2f, Scaling factor (obs->liab): %.2f\n\n", prev, mean, factor);

if(outfile!=NULL)
{
sprintf(filename,"%s.factor", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Prevalence %.6f\nAscertainment %.6f\nScaling_Factor_Observed->Liability %.6f\nScaling_Factor_Liability->Observed %.6f\n", prev, mean, factor, 1.0/factor);
fclose(output);
}

return(factor);
}	//end of get_factor

////////

double club_sandwich_one(double XTX, double *data, double *data2, int num_fams, int *famindex, int ns, int nf)
{
int i;
double sum, var;

double *mat;


mat=malloc(sizeof(double)*num_fams);

//elements of mat should contain t(es) Xs for each family
for(i=0;i<num_fams;i++){mat[i]=0;}
for(i=0;i<ns;i++){mat[famindex[i]]+=data[i]*data2[i];}

//can now get var
sum=0;for(i=0;i<num_fams;i++){sum+=pow(mat[i],2);}
var=sum*pow(XTX,-2)*(ns-1)/(ns-nf-1)*num_fams/(num_fams-1);

free(mat);

return(var);
}

////////

void club_sandwich_two(double *invXTX, double *data, double *data2, double *var, double *var2, int num_fams, int *famindex, int ns, int nf, int col2skip)
{
int i, j;
double sum, sum2;

double *mat, *mat2;


mat=malloc(sizeof(double)*num_fams*2);
mat2=malloc(sizeof(double)*num_fams*2);

//rows of mat should contain t(es) Xs for each family (for col 2 must include col2skip)
j=0;
for(i=0;i<num_fams;i++){mat[i+j*num_fams]=0;}
for(i=0;i<ns;i++){mat[famindex[i]+j*num_fams]+=data[i+j*ns]*data2[i];}
j=1;
for(i=0;i<num_fams;i++){mat[i+j*num_fams]=0;}
for(i=0;i<ns;i++){mat[famindex[i]+j*num_fams]+=data[i+(j+col2skip)*ns]*data2[i];}

//post-multiply mat by invXTX (do manually to avoid multi-threading)
for(j=0;j<2;j++)
{
for(i=0;i<num_fams;i++){mat2[i+j*num_fams]=0;}

for(i=0;i<num_fams;i++){mat2[i+j*num_fams]+=mat[i]*invXTX[0+2*j];}
for(i=0;i<num_fams;i++){mat2[i+j*num_fams]+=mat[i+num_fams]*invXTX[1+2*j];}
}

//get vars from the diagonal
sum=0;for(i=0;i<num_fams;i++){sum+=pow(mat2[i],2);}
*var=sum*(ns-1)/(ns-nf-2)*num_fams/(num_fams-1);

sum2=0;for(i=0;i<num_fams;i++){sum2+=pow(mat2[i+num_fams],2);}
*var2=sum2*(ns-1)/(ns-nf-2)*num_fams/(num_fams-1);


free(mat);free(mat2);
}

////////

void club_sandwich_three(double *invXTX, double *data, double *data2, double *var, double *var2, double *var3, int num_fams, int *famindex, int ns, int nf)
{
int i, j;
double sum, sum2, sum3;

double *mat, *mat2;


mat=malloc(sizeof(double)*num_fams*3);
mat2=malloc(sizeof(double)*num_fams*3);

//rows of mat should contain t(es) Xs for each family
for(j=0;j<3;j++)
{
for(i=0;i<num_fams;i++){mat[i+j*num_fams]=0;}

for(i=0;i<ns;i++){mat[famindex[i]+j*num_fams]+=data[i+j*ns]*data2[i];}
}

//post-multiply mat by invXTX (do manually to avoid multi-threading)
for(j=0;j<3;j++)
{
for(i=0;i<num_fams;i++){mat2[i+j*num_fams]=0;}

for(i=0;i<num_fams;i++){mat2[i+j*num_fams]+=mat[i]*invXTX[0+3*j];}
for(i=0;i<num_fams;i++){mat2[i+j*num_fams]+=mat[i+num_fams]*invXTX[1+3*j];}
for(i=0;i<num_fams;i++){mat2[i+j*num_fams]+=mat[i+2*num_fams]*invXTX[2+3*j];}
}

//get vars from the diagonal
sum=0;for(i=0;i<num_fams;i++){sum+=pow(mat2[i],2);}
*var=sum*(ns-1)/(ns-nf-3)*num_fams/(num_fams-1);

sum2=0;for(i=0;i<num_fams;i++){sum2+=pow(mat2[i+num_fams],2);}
*var2=sum2*(ns-1)/(ns-nf-3)*num_fams/(num_fams-1);

sum3=0;for(i=0;i<num_fams;i++){sum3+=pow(mat2[i+2*num_fams],2);}
*var3=sum3*(ns-1)/(ns-nf-3)*num_fams/(num_fams-1);

free(mat);free(mat2);
}

///////////////////////////

void copy_matrix(int ns, int length, double *X, double *X2, int type, double *weights)
//type=0, just copy; type=1, scale by weights; type=2, permute
//ok if X=X2, unless permuting
{
int i, j;
int *order;


if(X==X2&&type==2){printf("Error 300OP, please tell Doug\n\n");exit(1);}

if(type==0)	//just copy X into X2
{
#pragma omp parallel for private(j,i) schedule(static)
for(j=0;j<length;j++)
{
for(i=0;i<ns;i++){X2[(size_t)j*ns+i]=X[(size_t)j*ns+i];}
}
}
if(type==1)	//copy and pre-multiply by weights
{
#pragma omp parallel for private(j,i) schedule(static)
for(j=0;j<length;j++)
{
for(i=0;i<ns;i++){X2[(size_t)j*ns+i]=X[(size_t)j*ns+i]*weights[i];}
}
}
if(type==2)	//copy permuted version of X into X2
{
order=malloc(sizeof(int)*ns);
for(i=0;i<ns;i++){order[i]=i;}

for(j=0;j<length;j++)
{
permute_int(order, ns);
for(i=0;i<ns;i++){X2[(size_t)j*ns+i]=X[(size_t)j*ns+order[i]];}
}

free(order);
}
}

////////

void centre_matrix(double *X, int ns, int ns2, int length)
{
int i, j;
double sum, mean;


#pragma omp parallel for private(j, i, sum, mean) schedule(static)
for(j=0;j<length;j++)
{
sum=0;
for(i=0;i<ns;i++){sum+=X[(size_t)j*ns2+i];}
mean=sum/ns;
for(i=0;i<ns;i++){X[(size_t)j*ns2+i]-=mean;}
}
}

////////

void stand_matrix_nomiss(double *X, int ns, int ns2, int length)
{
int i, j;
double sum, sumsq, mean, var, value;


#pragma omp parallel for private(j, i, sum, sumsq, mean, var, value) schedule(static)
for(j=0;j<length;j++)
{
sum=0;sumsq=0;
for(i=0;i<ns;i++){sum+=X[(size_t)j*ns2+i];sumsq+=pow(X[(size_t)j*ns2+i],2);}
mean=sum/ns;
var=sumsq/ns-pow(mean,2);
if(var>0)
{
value=pow(var,-.5);
for(i=0;i<ns;i++){X[(size_t)j*ns2+i]=(X[(size_t)j*ns2+i]-mean)*value;}
}
else
{
for(i=0;i<ns;i++){X[(size_t)j*ns2+i]=0.0;}
}
}
}

////////

void stand_matrix_missing(double *X, int ns, int ns2, int length, double missingvalue)
{
int i, j, indcount;
double sum, sumsq, mean, var, value;


#pragma omp parallel for private(j, i, sum, sumsq, indcount, mean, var, value) schedule(static)
for(j=0;j<length;j++)
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<ns;i++)
{
if(X[(size_t)j*ns2+i]!=missingvalue){sum+=X[(size_t)j*ns2+i];sumsq+=pow(X[(size_t)j*ns2+i],2);indcount++;}
}
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);
if(var>0)
{
value=pow(var*indcount/ns,-.5);
for(i=0;i<ns;i++)
{
if(X[(size_t)j*ns2+i]!=missingvalue){X[(size_t)j*ns2+i]=(X[(size_t)j*ns2+i]-mean)*value;}
else{X[(size_t)j*ns2+i]=0;}
}
}
else
{
for(i=0;i<ns;i++){X[(size_t)j*ns2+i]=0.0;}
}
}
}

////////

void impute_matrix_missing(double *X, int ns, int ns2, int length, double missingvalue)
{
int i, j, indcount;
double sum, mean;


#pragma omp parallel for private(j, i, sum, indcount, mean) schedule(static)
for(j=0;j<length;j++)
{
sum=0;indcount=0;
for(i=0;i<ns;i++)
{
if(X[(size_t)j*ns2+i]!=missingvalue){sum+=X[(size_t)j*ns2+i];indcount++;}
}

if(indcount>0){mean=sum/indcount;}
else{mean=0;}

if(indcount<ns)	//replace missing values
{
for(i=0;i<ns;i++)
{
if(X[(size_t)j*ns2+i]==missingvalue){X[(size_t)j*ns2+i]=mean;}
}
}
}
}

///////////////////////////

void reg_covar_lin(double *Y, double *Z, int ns, int num_covars, int num_tops, double *thetas, double *thetasds, double *thetapvas, double *Yadj, int type, double *covher, double *topher)
//regress Y on Z linearly, then maybe fill Yadj with residuals and get fixed effect heritabilitiess
//type=0 - fill Yadj with regular residuals, type=1 - fill Yadj with standardized residuals
{
int i, j, num_fixed, one=1;
double sum, sumsq, mean, var, value, origss, fixedss, covarss, alpha, beta;
double *ZTZ, *ZTZ2, *ZTZ3, *ZTY, *thetas2;


num_fixed=num_covars+num_tops;

//allocate variables

ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTY=malloc(sizeof(double)*num_fixed);

//get original ss
sum=0;sumsq=0;for(i=0;i<ns;i++){sum+=Y[i];sumsq+=pow(Y[i],2);}
origss=sumsq-sum/ns*sum;

//solve for all fixed, and get fixed ss
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, thetas, &one);

sumsq=0;for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
fixedss=sumsq;for(j=0;j<num_fixed;j++){fixedss-=ZTY[j]*thetas[j];}

if(thetasds!=NULL)	//get sds and pvalues
{
for(j=0;j<num_fixed;j++)
{
thetasds[j]=pow(ZTZ[j+j*num_fixed]*fixedss/(ns-num_fixed),.5);
thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);
}
}

if(type==0||type==1)	//put residuals into Yadj
{
for(i=0;i<ns;i++){Yadj[i]=Y[i];}
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Yadj, &one);

if(type==1)	//standardize Yadj
{
sum=0;sumsq=0;for(i=0;i<ns;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/ns;var=sumsq/ns-pow(mean,2);
value=pow(var,-.5);
for(i=0;i<ns;i++){Yadj[i]=Yadj[i]*value;}
}
}

if(covher!=NULL)	//solve for just covariates, and get covher and topher
{
if(num_tops==0)	//already done
{covarss=fixedss;}
else
{
thetas2=malloc(sizeof(double)*num_covars);
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_covars);
(void)eigen_invert(ZTZ, num_covars, ZTZ2, -1, ZTZ3, 1);

alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_covars, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemv_("N", &num_covars, &num_covars, &alpha, ZTZ, &num_covars, ZTY, &one, &beta, thetas2, &one);

sumsq=0;for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
covarss=sumsq;for(j=0;j<num_covars;j++){fixedss-=ZTY[j]*thetas2[j];}
free(thetas2);
}

*covher=(origss-covarss)/origss;
*topher=(covarss-fixedss)/covarss;
}

free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTY);
}	//end of reg_covar_lin	

////////

void reg_covar_lin_missing(double *Y, double *Z, int ns, int num_fixed, double *thetas, double *thetasds, double *thetapvas, double *Yadj, double missingvalue)
//regress Y on Z linearly and fill Yadj with padded residuals
{
int i, j, indcount, one=1;
double sum, sumsq, mean, alpha, beta;

int *indexer;
double *Y2, *Z2, *ZTZ, *ZTZ2, *ZTZ3, *ZTY;


indexer=malloc(sizeof(int)*ns);
Y2=malloc(sizeof(double)*ns);
Z2=malloc(sizeof(double)*ns*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTY=malloc(sizeof(double)*num_fixed);

//index non missing
indcount=0;
for(i=0;i<ns;i++)
{
if(Y[i]!=missingvalue){indexer[indcount]=i;indcount++;}
}

//load values into Y2 and Z2
for(i=0;i<indcount;i++)
{
Y2[i]=Y[indexer[i]];
for(j=0;j<num_fixed;j++){Z2[i+j*indcount]=Z[indexer[i]+j*ns];}
}

//estimate fixed effects
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &indcount, &alpha, Z2, &indcount, Z2, &indcount, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

alpha=1.0;beta=0.0;
dgemv_("T", &indcount, &num_fixed, &alpha, Z2, &indcount, Y2, &one, &beta, ZTY, &one);
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, thetas, &one);

//get residual sumsq
sumsq=0;for(i=0;i<indcount;i++){sumsq+=pow(Y2[i],2);}
for(j=0;j<num_fixed;j++){sumsq-=ZTY[j]*thetas[j];}

//get sds and pvas
for(j=0;j<num_fixed;j++)
{
thetasds[j]=pow(ZTZ[j+j*num_fixed]*sumsq/(indcount-num_fixed),.5);
thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);
}

//get corrected phenotypes
alpha=-1.0;beta=1.0;
dgemv_("N", &indcount, &num_fixed, &alpha, Z2, &indcount, thetas, &one, &beta, Y2, &one);

//load up Yadj, setting missing values to the mean
sum=0;for(i=0;i<indcount;i++){sum+=Y2[i];}
mean=sum/indcount;
for(i=0;i<ns;i++){Yadj[i]=mean;}
for(i=0;i<indcount;i++){Yadj[indexer[i]]=Y2[i];}

free(indexer);free(Y2);free(Z2);free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTY);
}	//end of reg_covar_lin_missing

///////////////////////////

void get_log_probs(double *probs, double *Zthetas, double *Z, double *thetas, int ns, int num_fixed, double *offsets)
{
int i, one=1;
double alpha, beta;

alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Zthetas, &one);

if(offsets==NULL)	//usual case
{
for(i=0;i<ns;i++){probs[i]=pow(1+exp(-Zthetas[i]),-1);}
}
else	//allow for offsets
{
for(i=0;i<ns;i++){probs[i]=pow(1+exp(-Zthetas[i]-offsets[i]),-1);}
}

for(i=0;i<ns;i++)
{
if(probs[i]==0){probs[i]=1e-10;}
if(probs[i]==1){probs[i]=1-1e-10;}
}
}

void reg_covar_log(double *Y, double *Z, int ns, int num_covars, int num_tops, double *offsets, double *thetas, double *thetasds, double *thetapvas, double *Yadj, int type, double *covher, double *topher, double prev, double tol, int maxiter)
//regress Y on Z logistically, then maybe fill Yadj and get fixed effect heritabilities
//there may be padded values - in theory should adjust covher and topher, but currently rarely use padding
//type=0 - fill Yadj with function of liabilities, type=1 - fill Yadj with new thresholds
//type=2 - fill Yadj with probabilities, type=3 - fill Yadj with probabilities*(1-probabilities)
//type=4 - fill Yadj with Y - probabilities (dont think types 3 or 4 are used)
{
int i, j, count, count2, num_fixed, one=1;
double sum, sum2, sumsq, sumsq2, var, var2, value, value2, relax;
double ascer, like, like2, likeold, fixedvar, covarvar, alpha, beta;
double *Z2, *Zthetas, *probs, *probs2, *AI, *AI2, *AI3, *BI, *thetadiffs, *Ks, *thresh, *regs, *thetas2;


num_fixed=num_covars+num_tops;

//allocate variables

Z2=malloc(sizeof(double)*ns*num_fixed);
Zthetas=malloc(sizeof(double)*ns);
probs=malloc(sizeof(double)*ns);
probs2=malloc(sizeof(double)*ns);
AI=malloc(sizeof(double)*num_fixed*num_fixed);
AI2=malloc(sizeof(double)*num_fixed);
AI3=malloc(sizeof(double)*num_fixed*num_fixed);
BI=malloc(sizeof(double)*num_fixed);
thetadiffs=malloc(sizeof(double)*num_fixed);

if(prev!=-9999)
{
Ks=malloc(sizeof(double)*ns);
thresh=malloc(sizeof(double)*ns);
regs=malloc(sizeof(double)*ns);
}

//get ascer
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
ascer=sum/ns;

//solve for all fixed - first covariate starts at logit ascer, rest are zero
thetas[0]=log(ascer/(1-ascer));
for(j=1;j<num_fixed;j++){thetas[j]=0;}

//get starting probabilities
get_log_probs(probs, Zthetas, Z, thetas, ns, num_fixed, offsets);

//get likelihood
like=0;for(i=0;i<ns;i++){like+=Y[i]*log(probs[i])+(1-Y[i])*log(1-probs[i]);}

//update via newton-raphson
count=0;
relax=1;
while(1)
{
if(count==maxiter){printf("Error, logistic regression did not converge within %d iterations; consider reducing the number of fixed effects\n\n", maxiter);exit(1);}

//get Z2=Z*probs*(1-probs)
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*probs[i]*(1-probs[i]);}
}

//Bi = 1st deriv = sum(Yi-pi)Zi
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
BI[j]=0;for(i=0;i<ns;i++){BI[j]+=(Y[i]-probs[i])*Z[i+j*ns];}
}

//AI = -2nd deriv = sum(pi(1-pi)ZiZi)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_fixed);

//invert AI and get proposed move
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, AI, &num_fixed, BI, &one, &beta, thetadiffs, &one);

likeold=like;
relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs, update probabilities and likelihood
for(j=0;j<num_fixed;j++){thetas[j]+=relax*thetadiffs[j];}
get_log_probs(probs2, Zthetas, Z, thetas, ns, num_fixed, offsets);
like2=0;for(i=0;i<ns;i++){like2+=Y[i]*log(probs2[i])+(1-Y[i])*log(1-probs2[i]);}

if(like2>like-tol)	//accept move
{
like=like2;
for(i=0;i<ns;i++){probs[i]=probs2[i];}
break;
}
else	//move back and next turn try smaller move
{
for(j=0;j<num_fixed;j++){thetas[j]-=relax*thetadiffs[j];}
relax*=.5;
}
}

//see if breaking
if(fabs(like-likeold)<tol){break;}

count++;
}

//recompute inverse AI
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*probs[i]*(1-probs[i]);}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_fixed);
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);

//get sds and pvalues
for(j=0;j<num_fixed;j++)
{
if(AI[j+j*num_fixed]>0){thetasds[j]=pow(AI[j+j*num_fixed],.5);thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);}
else{thetasds[j]=-9999;thetapvas[j]=-9999;}
}

if(type==2)	//set Yadj to probs
{
for(i=0;i<ns;i++){Yadj[i]=probs[i];}
}
if(type==3)	//set Yadj to probs*(1-probs)
{
for(i=0;i<ns;i++){Yadj[i]=probs[i]*(1-probs[i]);}
}
if(type==4)	//set Yadj to Y-probs
{
for(i=0;i<ns;i++){Yadj[i]=Y[i]-probs[i];}
}

if(prev!=-9999)	//get thresholds, fill Yadj for types 0 and 1, maybe get covarvar and topher 
{
//get Ks (population probs of being case given covariates) and corresponding thresholds - probs was set above
value=prev*(1-ascer)/ascer/(1-prev);
for(i=0;i<ns;i++){Ks[i]=value*probs[i]/(1+value*probs[i]-probs[i]);thresh[i]=normal_inv(1-Ks[i]);}

if(type==0)	//get PCGC regressions, then set Yadj to (Y-probs)/root(probs*(1-probs)) divided by regressions
{
for(i=0;i<ns;i++)
{
value=pow(2*M_PI,-.5)*exp(-.5*pow(thresh[i],2))*(1-probs[i]*(ascer-prev)/ascer/(1-prev));
value2=pow(probs[i]*(1-probs[i]),.5)*(Ks[i]+(1-Ks[i])*prev*(1-ascer)/ascer/(1-prev));
regs[i]=value/value2;
Yadj[i]=(Y[i]-probs[i])*pow(probs[i]*(1-probs[i]),-.5)/regs[i];
}
}
if(type==1)	//set Yadj to revised thresholds
{
for(i=0;i<ns;i++){Yadj[i]=thresh[i];}
}

if(covher!=NULL)	//solve for just covariates, and get covher and topher
{
//get variance due to all fixed effects - I think only valid to consider non-padded values
sum=0;sum2=0;sumsq=0;sumsq2=0;count=0;count2=0;
for(i=0;i<ns;i++)
{
if(Y[i]==0){sum+=thresh[i];sumsq+=pow(thresh[i],2);count++;}
if(Y[i]==1){sum2+=thresh[i];sumsq2+=pow(thresh[i],2);count2++;}
}
var=sumsq/count-pow(sum/count,2);
var2=sumsq2/count2-pow(sum2/count2,2);
fixedvar=prev*(1-prev)*pow(sum2/count2-sum/count,2)+prev*var2+(1-prev)*var;

//get variance for just covariates
if(num_tops==0)	//already done
{covarvar=fixedvar;}
else	//must estimate thresholds for just covariates
{
thetas2=malloc(sizeof(double)*num_covars);
thetas2[0]=log(ascer/(1-ascer));
for(j=1;j<num_covars;j++){thetas2[j]=0;}

//get starting probabilities and likelihood
get_log_probs(probs, Zthetas, Z, thetas2, ns, num_covars, offsets);
like=0;for(i=0;i<ns;i++){like+=Y[i]*log(probs[i])+(1-Y[i])*log(1-probs[i]);}

count=0;
relax=1;
while(1)
{
if(count==maxiter){printf("Error, logistic regression did not converge within %d iterations; consider reducing the number of fixed effects\n\n", maxiter);exit(1);}

//get Z2=Z*probs*(1-probs)
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_covars;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*probs[i]*(1-probs[i]);}
}

//Bi = 1st deriv = sum(Yi-pi)Zi
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_covars;j++)
{
BI[j]=0;for(i=0;i<ns;i++){BI[j]+=(Y[i]-probs[i])*Z[i+j*ns];}
}

//AI = -2nd deriv = sum(pi(1-pi)ZiZi)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_covars, &num_covars, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_covars);

//invert AI and get proposed move
(void)eigen_invert(AI, num_covars, AI2, -1, AI3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_covars, &num_covars, &alpha, AI, &num_covars, BI, &one, &beta, thetadiffs, &one);

likeold=like;
relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs, update probabilities and likelihood
for(j=0;j<num_covars;j++){thetas2[j]+=relax*thetadiffs[j];}
get_log_probs(probs2, Zthetas, Z, thetas2, ns, num_covars, offsets);
like2=0;for(i=0;i<ns;i++){like2+=Y[i]*log(probs2[i])+(1-Y[i])*log(1-probs2[i]);}

if(like2>like-tol)	//accept move
{
like=like2;
for(i=0;i<ns;i++){probs[i]=probs2[i];}
break;
}
else	//move back and next turn try smaller move
{
for(j=0;j<num_covars;j++){thetas2[j]-=relax*thetadiffs[j];}
relax*=.5;
}
}

//see if breaking
if(fabs(like-likeold)<tol){break;}

count++;
}

//update thresholds
value=prev*(1-ascer)/ascer/(1-prev);
for(i=0;i<ns;i++){Ks[i]=value*probs[i]/(1+value*probs[i]-probs[i]);thresh[i]=normal_inv(1-Ks[i]);}

//get variance
sum=0;sum2=0;sumsq=0;sumsq2=0;count=0;count2=0;
for(i=0;i<ns;i++)
{
if(Y[i]==0){sum+=thresh[i];sumsq+=pow(thresh[i],2);count++;}
if(Y[i]==1){sum2+=thresh[i];sumsq2+=pow(thresh[i],2);count2++;}
}
var=sumsq/count-pow(sum/count,2);
var2=sumsq2/count2-pow(sum2/count2,2);
covarvar=prev*(1-prev)*pow(sum2/count2-sum/count,2)+prev*var2+(1-prev)*var;

free(thetas2);
}	//end of solving for just covariates

*covher=covarvar/(1+fixedvar);
*topher=(fixedvar-covarvar)/(1+fixedvar-covarvar);
}
}	//end of prev!=-9999

free(Z2);free(Zthetas);free(probs);free(probs2);free(AI);free(AI2);free(AI3);free(BI);free(thetadiffs);
if(prev!=-9999){free(Ks);free(thresh);free(regs);}
}	//end of reg_covar_log

////////

void reg_covar_log_missing(double *Y, double *Z, int ns, int num_fixed, double *offsets, double *thetas, double *thetasds, double *thetapvas, double *nullprobs, double *nullweights, double *Yadj, double tol, int maxiter, double missingvalue, int type)
//regress Y on Z logistically and fill nullprobs, nullweights, Yadj (type=0 - (Y-mu); type=1 - (Y-mu)/W, padding if necessary)
{
int i, j, count, indcount, one=1;
double sum, ascer, like, like2, likeold, relax, alpha, beta;

int *indexer;
double *Z2, *Zthetas, *probs, *probs2, *AI, *AI2, *AI3, *BI, *thetadiffs;


indexer=malloc(sizeof(int)*ns);
Z2=malloc(sizeof(double)*ns*num_fixed);
Zthetas=malloc(sizeof(double)*ns);
probs=malloc(sizeof(double)*ns);
probs2=malloc(sizeof(double)*ns);
AI=malloc(sizeof(double)*num_fixed*num_fixed);
AI2=malloc(sizeof(double)*num_fixed);
AI3=malloc(sizeof(double)*num_fixed*num_fixed);
BI=malloc(sizeof(double)*num_fixed);
thetadiffs=malloc(sizeof(double)*num_fixed);

//index non missing
indcount=0;
for(i=0;i<ns;i++)
{
if(Y[i]!=missingvalue){indexer[indcount]=i;indcount++;}
}

//get ascer
sum=0;for(i=0;i<indcount;i++){sum+=Y[indexer[i]];}
ascer=sum/indcount;

//first covariate starts at logit ascer, rest are zero
thetas[0]=log(ascer/(1-ascer));
for(j=1;j<num_fixed;j++){thetas[j]=0;}

//get starting probabilities - use all ind
get_log_probs(probs, Zthetas, Z, thetas, ns, num_fixed, offsets);

//get likelihood - use non-missing ind
like=0;for(i=0;i<indcount;i++){like+=Y[indexer[i]]*log(probs[indexer[i]])+(1-Y[indexer[i]])*log(1-probs[indexer[i]]);}

//update via newton-raphson
count=0;
relax=1;
while(1)
{
if(count==maxiter){printf("Error, logistic regression did not converge within %d iterations; consider reducing the number of fixed effects\n\n", maxiter);exit(1);}

//get Z2=Z*probs*(1-probs) - set to zero for missing ind (makes computing AI easier)
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=0;}
for(i=0;i<indcount;i++){Z2[indexer[i]+j*ns]=Z[indexer[i]+j*ns]*probs[indexer[i]]*(1-probs[indexer[i]]);}
}

//Bi = 1st deriv = sum(Yi-pi)Zi - use non-missing ind
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
BI[j]=0;for(i=0;i<indcount;i++){BI[j]+=(Y[indexer[i]]-probs[indexer[i]])*Z[indexer[i]+j*ns];}
}

//AI = -2nd deriv = sum(pi(1-pi)ZiZi)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_fixed);

//invert AI and get proposed move
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, AI, &num_fixed, BI, &one, &beta, thetadiffs, &one);

likeold=like;
relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs, update probabilities and likelihood
for(j=0;j<num_fixed;j++){thetas[j]+=relax*thetadiffs[j];}
get_log_probs(probs2, Zthetas, Z, thetas, ns, num_fixed, offsets);
like2=0;for(i=0;i<indcount;i++){like2+=Y[indexer[i]]*log(probs2[indexer[i]])+(1-Y[indexer[i]])*log(1-probs2[indexer[i]]);}

if(like2>like-tol)	//accept move
{
like=like2;
for(i=0;i<ns;i++){probs[i]=probs2[i];}
break;
}
else	//move back and next turn try smaller move
{
for(j=0;j<num_fixed;j++){thetas[j]-=relax*thetadiffs[j];}
relax*=.5;
}
}

//see if breaking
if(fabs(like-likeold)<tol){break;}

count++;
}

//recompute inverse AI
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=0;}
for(i=0;i<indcount;i++){Z2[indexer[i]+j*ns]=Z[indexer[i]+j*ns]*probs[indexer[i]]*(1-probs[indexer[i]]);}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_fixed);
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);

//get sds and pvalues
for(j=0;j<num_fixed;j++)
{
if(AI[j+j*num_fixed]>0){thetasds[j]=pow(AI[j+j*num_fixed],.5);thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);}
else{thetasds[j]=-9999;thetapvas[j]=-9999;}
}

//load nullprobs
for(i=0;i<ns;i++){nullprobs[i]=probs[i];}

//load nullprobs
for(i=0;i<ns;i++){nullweights[i]=probs[i]*(1-probs[i]);}

//load Yadj, setting missing values to zero
for(i=0;i<ns;i++){Yadj[i]=0;}
for(i=0;i<indcount;i++){Yadj[indexer[i]]=Y[indexer[i]]-probs[indexer[i]];}
if(type==1)	//scale by weights
{
for(i=0;i<ns;i++){Yadj[i]=Yadj[i]/nullweights[i];}
}

free(indexer);free(Z2);free(Zthetas);free(probs);free(probs2);free(AI);free(AI2);free(AI3);free(BI);free(thetadiffs);
}	//end of reg_covar_log_missing

////////

void reg_covar_log_lite(double *Y, double *Z, int ns, int num_fixed, double *offsets, double *stats, double *thetas2, double tol, int maxiter)
//regress Y on Z logistically - will be no missing
{
int i, j, count, one=1;
double sum, relax; 
double like, like2, likeold, alpha, beta;
double *thetas, *Zthetas, *probs, *probs2, *Z2, *AI, *AI2, *AI3, *BI, *thetadiffs;


thetas=malloc(sizeof(double)*num_fixed);
Zthetas=malloc(sizeof(double)*ns);
probs=malloc(sizeof(double)*ns);
probs2=malloc(sizeof(double)*ns);
Z2=malloc(sizeof(double)*ns*num_fixed);
AI=malloc(sizeof(double)*num_fixed*num_fixed);
AI2=malloc(sizeof(double)*num_fixed);
AI3=malloc(sizeof(double)*num_fixed*num_fixed);
BI=malloc(sizeof(double)*num_fixed);
thetadiffs=malloc(sizeof(double)*num_fixed);

if(thetas2==NULL)	//use default starting values
{
//first covariate starts at logit(ascer), rest are zero
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
thetas[0]=log(sum/(ns-sum));
for(j=1;j<num_fixed;j++){thetas[j]=0;}
}
else	//use provided values - must be single-predictor regression
{
//first num_fixed-1 values provided, last set to zero
for(j=0;j<num_fixed-1;j++){thetas[j]=thetas2[j];}
thetas[num_fixed-1]=0;
}

//get starting probabilities and likelihood
get_log_probs(probs, Zthetas, Z, thetas, ns, num_fixed, offsets);
like=0;for(i=0;i<ns;i++){like+=Y[i]*log(probs[i])+(1-Y[i])*log(1-probs[i]);}

count=0;
relax=1;
while(1)
{
if(count==maxiter){printf("Error, logistic regression did not converge within %d iterations; consider reducing the number of fixed effects\n\n", maxiter);exit(1);}

//get Z2=Z*probs*(1-probs)
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*probs[i]*(1-probs[i]);}
}

//Bi = 1st deriv = sum(Yi-pi)Zi
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
BI[j]=0;for(i=0;i<ns;i++){BI[j]+=(Y[i]-probs[i])*Z[i+j*ns];}
}

//AI = -2nd deriv = sum(pi(1-pi)ZiZi)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_fixed);

//invert AI and get proposed move
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, AI, &num_fixed, BI, &one, &beta, thetadiffs, &one);

likeold=like;
relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs, update probabilities and likelihood
for(j=0;j<num_fixed;j++){thetas[j]+=relax*thetadiffs[j];}
get_log_probs(probs2, Zthetas, Z, thetas, ns, num_fixed, offsets);
like2=0;for(i=0;i<ns;i++){like2+=Y[i]*log(probs2[i])+(1-Y[i])*log(1-probs2[i]);}

if(like2>like-tol)	//accept move
{
like=like2;
for(i=0;i<ns;i++){probs[i]=probs2[i];}
break;
}
else	//move back and next turn try smaller move
{
for(j=0;j<num_fixed;j++){thetas[j]-=relax*thetadiffs[j];}
relax*=.5;
}
}

//see if breaking
if(fabs(like-likeold)<tol){break;}

count++;
}

//recompute inverse AI
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*probs[i]*(1-probs[i]);}
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z2, &ns, &beta, AI, &num_fixed);
(void)eigen_invert(AI, num_fixed, AI2, -1, AI3, 1);

stats[0]=thetas[num_fixed-1];
if(AI[num_fixed*num_fixed-1]>0)
{
stats[1]=pow(AI[num_fixed*num_fixed-1],.5);
stats[2]=stats[0]/stats[1];
stats[3]=erfc(fabs(stats[0]/stats[1])*M_SQRT1_2);
}
else{stats[1]=-9999;stats[2]=-9999;stats[3]=-9999;}

free(thetas);free(Zthetas);free(probs);free(probs2);free(Z2);free(AI);free(AI2);free(AI3);free(BI);free(thetadiffs);
}	//end of reg_covar_log_lite

///////////////////////////

void get_Z2(double *Z2, double *Z, double *Z3, int ns, int num_fixed)
//set X = inv(ZTZ) Z
{
double alpha, beta;
double *ZTZ, *ZTZ2, *ZTZ3;


ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z3, &ns, &beta, ZTZ, &num_fixed);
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);

alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_fixed, &num_fixed, &alpha, Z, &ns, ZTZ, &num_fixed, &beta, Z2, &ns);

free(ZTZ);free(ZTZ2);free(ZTZ3);
}

////////

void reg_covar_matrix(double *X, double *Z, int ns, int length, int num_fixed)
//replace X with residual from regressing on Z
{
double alpha, beta;
double *ZTZ, *ZTZ2, *ZTX;


if(length>0)
{
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTX=malloc(sizeof(double)*num_fixed*length);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
dgemm_("T", "N", &num_fixed, &length, &ns, &alpha, Z, &ns, X, &ns, &beta, ZTX, &num_fixed);

(void)eigen_invert(ZTZ, num_fixed, ZTZ2, length, ZTX, 1);

alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &length, &num_fixed, &alpha, Z, &ns, ZTX, &num_fixed, &beta, X, &ns);

free(ZTZ);free(ZTZ2);free(ZTX);
}
}

////////

void reg_covar_weighted(double *X, double *Z, int ns, int length, int num_fixed, double *weights)
//replace X with residual from regressing on Z with weights
{
int i, j;
double alpha, beta;
double *Z2, *ZTZ, *ZTZ2, *ZTX;


if(length>0)
{
Z2=malloc(sizeof(double)*ns*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTX=malloc(sizeof(double)*num_fixed*length);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*weights[i];}
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z2, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
dgemm_("T", "N", &num_fixed, &length, &ns, &alpha, Z2, &ns, X, &ns, &beta, ZTX, &num_fixed);

(void)eigen_invert(ZTZ, num_fixed, ZTZ2, length, ZTX, 1);

alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &length, &num_fixed, &alpha, Z, &ns, ZTX, &num_fixed, &beta, X, &ns);

free(Z2);free(ZTZ);free(ZTZ2);free(ZTX);
}
}

////////

void reg_covar_transpose(double *X, double *Z, int ns, int length, int num_fixed, double *weights)
//replace X with residual from regressing on Z with weights "backwards" (get t(H) X instead of H X)
{
int i, j;
double alpha, beta;
double *Z2, *ZTZ, *ZTZ2, *ZTX;


if(length>0)
{
Z2=malloc(sizeof(double)*ns*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTX=malloc(sizeof(double)*num_fixed*length);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*weights[i];}
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z2, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
dgemm_("T", "N", &num_fixed, &length, &ns, &alpha, Z, &ns, X, &ns, &beta, ZTX, &num_fixed);

(void)eigen_invert(ZTZ, num_fixed, ZTZ2, length, ZTX, 1);

alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &length, &num_fixed, &alpha, Z2, &ns, ZTX, &num_fixed, &beta, X, &ns);

free(Z2);free(ZTZ);free(ZTZ2);free(ZTX);
}
}

///////////////////////////

void reg_covar_thetas(double *thetas, double *X, double *Z, int ns, int length, int num_fixed, double *weights, int type)
//divides reg_covar_matrix into two - type=0 - compute thetas, type=1 - use thetas, type=2 - both
{
int i, j;
double alpha, beta;
double *Z2, *ZTZ, *ZTZ2;


if(length>0)
{
if(type==0||type==2)	//compute thetas
{
if(weights==NULL)	//without weights
{
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
dgemm_("T", "N", &num_fixed, &length, &ns, &alpha, Z, &ns, X, &ns, &beta, thetas, &num_fixed);

(void)eigen_invert(ZTZ, num_fixed, ZTZ2, length, thetas, 1);

free(ZTZ);free(ZTZ2);
}
else	//with weights
{
Z2=malloc(sizeof(double)*ns*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*weights[i];}
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z2, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
dgemm_("T", "N", &num_fixed, &length, &ns, &alpha, Z2, &ns, X, &ns, &beta, thetas, &num_fixed);

(void)eigen_invert(ZTZ, num_fixed, ZTZ2, length, thetas, 1);

free(Z2);free(ZTZ);free(ZTZ2);
}
}

if(type==1||type==2)	//use thetas
{
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &length, &num_fixed, &alpha, Z, &ns, thetas, &num_fixed, &beta, X, &ns);
}
}
}

///////////////////////////

double comp_like(int ns, double *residuals, double resvar, double pen, int dichot, double *nullweights, int type)
//dichot=0 - linear, dichot=1 - quasi-logistic, dichot=2 - logistic (not currently used)
{
int i;
double sum, sumsq, value;


if(dichot==0)
{
sumsq=0;for(i=0;i<ns;i++){sumsq+=pow(residuals[i],2);}
value=-.5*ns*log(2*M_PI*resvar)-.5*sumsq/resvar-pen;
}
if(dichot==1)
{
sum=0;for(i=0;i<ns;i++){sum+=log(nullweights[i]);}
sumsq=0;for(i=0;i<ns;i++){sumsq+=pow(residuals[i],2)*nullweights[i];}
value=-.5*ns*log(2*M_PI*resvar)-.5*sumsq/resvar+.5*sum-pen;
}
if(dichot==2)
{
sum=0;for(i=0;i<ns;i++){sum+=log(nullweights[i]);}
sumsq=0;for(i=0;i<ns;i++){sumsq+=pow(residuals[i],2)*nullweights[i];}
value=-.5*ns*log(2*M_PI)-.5*sumsq+.5*sum-pen;
}

if(type!=-9999){printf("Model %d has like %f pen %f sumsq %f\n", type+1, value, pen, sumsq);}

if(value!=value&&type==-9999){printf("Partial problem pen %f sumsq %f\n", pen, sumsq);}

return(value);
}

///////////////////////////

void set_lambdas(int p, double *lambdas, double *lambdas2, double *lambdas3, double *lambdas4, int length, double *exps, double varphen, double her, double *trylams, double *tryps, double *tryp2s, double *tryp3s, double *tryp4s, double *tryf2s, int type)
//type=1 - lasso-sparse, type=2 - lasso-shrink, type=3 - ridge, type=4 - bolt, type=5 bayesr-sparse, type=6 - bayesr-shrink, type=7 - elastic
//have already set all lambdas to zero
{
int j;
double value, value2;


if(type==1)	//lasso-sparse - lambda is set to match the lassosum paper
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[(size_t)p*length+j]=trylams[p]*pow(length*exps[j],-.5);}
}
}

if(type==2)	//lasso-shrink - exp(betaj^2) = exps*varphen*her = 2/lambdas^2
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[(size_t)p*length+j]=pow(2.0/exps[j]/her/varphen,.5);}
}
}

if(type==3)	//ridge - exp(betaj^2) = exps*varphen*her = lambdas
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[(size_t)p*length+j]=exps[j]*her*varphen;}
}
}

if(type==4)	//bolt - exp(betaj^2) = exps*varphen*her = plambdas + (1-p)lambdas2 - might have f2=0, in which case lambdas2=0
{
value=(1-tryf2s[p])/tryps[p];
value2=tryf2s[p]/tryp2s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas[(size_t)p*length+j]=exps[j]*her*varphen*value;lambdas2[(size_t)p*length+j]=exps[j]*her*varphen*value2;}
}
}

if(type==5)	//bayesr-sparse - exp(betaj^2) = exps*varphen*her = p2lambdas2 + p3lambdas3 + p4lambdas4, where lambdas2=lambdas4/100, lambdas3=lambdas4/10
{
value=tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0)
{
value2=exps[j]*her*varphen/value;
lambdas2[(size_t)p*length+j]=value2/100;
lambdas3[(size_t)p*length+j]=value2/10;
lambdas4[(size_t)p*length+j]=value2;
}
}
}

if(type==6)	//bayesr-shrink - exp(betaj^2) = exps*varphen*her = plambdas + p2lambdas2 + p3lambdas3 + p4lambdas4
{
value=tryps[p]/1000+tryp2s[p]/100+tryp3s[p]/10+tryp4s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0)
{
value2=exps[j]*her*varphen/value;
lambdas[(size_t)p*length+j]=value2/1000;
lambdas2[(size_t)p*length+j]=value2/100;
lambdas3[(size_t)p*length+j]=value2/10;
lambdas4[(size_t)p*length+j]=value2;
}
}
}

if(type==7)	//elastic - exp(betaj^2) = exps*varphen*her = 2p/lambdas^2 + 2p/lambdas2^2 + (1-p)lambdas3
//lasso model has p3=0 (and f2=0), ridge model has p3=1 (and f2=1)
{
if(tryp3s[p]==0)	//lasso - use lambda and lambda2
{
for(j=0;j<length;j++)
{
if(exps[j]>0)
{
lambdas[(size_t)p*length+j]=pow(2.0/exps[j]/her/varphen,.5);lambdas2[(size_t)p*length+j]=pow(2.0/exps[j]/her/varphen,.5);
}
}
}
if(tryp3s[p]==1)	//ridge - just use lambda3
{
for(j=0;j<length;j++)
{
if(exps[j]>0){lambdas3[(size_t)p*length+j]=exps[j]*her*varphen;}
}
}
if(tryp3s[p]>0&&tryp3s[p]<1)	//elastic - use all three lambdas (note lambdas=lambdas2)
{
value=2*(1-tryp3s[p])/(1-tryf2s[p]);
value2=tryf2s[p]/tryp3s[p];
for(j=0;j<length;j++)
{
if(exps[j]>0)
{
lambdas[(size_t)p*length+j]=pow(value/exps[j]/her/varphen,.5);
lambdas2[(size_t)p*length+j]=pow(value/exps[j]/her/varphen,.5);
lambdas3[(size_t)p*length+j]=exps[j]*her*varphen*value2;
}
}
}
}
}	//end of set_lambdas

///////////////////////////

double get_postmean(double sum, double lam, double lam2, double lam3, double lam4, double dsq, double resvar, double pp, double pp2, double pp3, double pp4, double *pen, int type, double *prob)
//type=1 - lasso-sparse, type=2 - lasso-shrink, type=3 - ridge, type=4 - bolt, type=5 bayesr-sparse, type=6 - bayesr-shrink, type=7 - elastic
//when using gaussians, means are sum/A, vars are resvar/A, where A is (dsq+resvar/lam), mults are prop*sqrt(var/lam)*exp(mean^2/2var)
//previously reported posterior variance (now instead report post prob associated)
{
double mean, mean2, mean3, mean4, var, var2, var3, var4, frac, frac2, frac3, frac4, area, area2, trun, trun2, brun, brun2;
double max, value, value2, value3, value4, postmean, postvar;

if(type==1)	//lasso with posterior mode
{
if(fabs(sum)>resvar*lam)
{
if(sum>0){postmean=(sum-resvar*lam)/dsq;}
else{postmean=(sum+resvar*lam)/dsq;}
}
else{postmean=0;}

//will not be using penalty

if(prob!=NULL){*prob=(postmean!=0);}
}

if(type==2)	//lasso with posterior mean - divide prior as positive exp then negative exp
{
//these are the means and variances of the non-truncated distributions
mean=(sum-resvar*lam)/dsq;
mean2=(sum+resvar*lam)/dsq;
var=resvar/dsq;
var2=resvar/dsq;

//for positive normal, area (denominator) is 1-Phi(-mean/var^.5) = Phi(mean/var^.5) = .5*erfc(-mean/var^.5 root(1/2))
//for negative normal, area (denominator) is Phi(-mean2/var2^.5) = .5*erfc(mean2/var2^.5 root(1/2))
area=.5*erfc(-mean*pow(var,-.5)*M_SQRT1_2);
area2=.5*erfc(mean2*pow(var2,-.5)*M_SQRT1_2);

//for means and variances, require the following
value=exp(-.5*pow(-mean,2)/var)*pow(2*M_PI,-.5);	//on wiki page, this is phi(alpha)
value2=exp(-.5*pow(-mean2,2)/var2)*pow(2*M_PI,-.5);	//on wiki page, this is phi(beta)

//get mean and variance of positive truncated distribution
//if area < 1e-16, will replace with point mass at zero
if(area>1e-16)
{
trun=mean+value/area*pow(var,.5);
brun=var*(1-mean*pow(var,-.5)*value/area-pow(value/area,2));
}
else{trun=0;brun=0;}

//get mean and variance of negative truncated distribution
//if area2 < 1e-16, will replace with point mass at zero
if(area2>1e-16)
{
trun2=mean2-value2/area2*pow(var2,.5);
brun2=var2*(1+mean2*pow(var2,-.5)*value2/area2-pow(-value2/area2,2));
}
else{trun2=0;brun2=0;}

//get fractions
if(area>1e-16&&area2>1e-16)	//normal case
{
max=log(area)+pow(mean,2)/var;
if(log(area2)+pow(mean2,2)/var2>max){max=log(area2)+pow(mean2,2)/var2;}
//here ignore pow(var,.5)=pow(var2,.5) in respective lines
value=exp(log(area)+.5*pow(mean,2)/var-.5*max);
value2=exp(log(area2)+.5*pow(mean2,2)/var2-.5*max);
}
else	//one or both areas are small - will set to 1-0, 0-1 or 1-1
{
if(area>1e-16&&area2<=1e-16){value=1;value2=0;}
if(area<=1e-16&&area2>1e-16){value=0;value2=1;}
if(area<=1e-16&&area2<=1e-16){value=1;value2=1;}
}

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*trun+frac2*trun2;
postvar=frac*(brun+pow(trun,2))+frac2*(brun2+pow(trun2,2))-pow(postmean,2);

if(pen!=NULL)	
{
*pen+=.5*dsq*postvar/resvar;
if(area>1e-16)	//using truncated normal (although maybe frac can be zero?)
{
if(frac>0){*pen+=frac*(log(frac/0.5/lam/area)-.5*log(2*M_PI*var)-.5*(brun+pow(trun,2)-2*trun*mean+pow(mean,2))/var+lam*trun);}
}
//else will have trun=brun=0, and probably frac=0, so all terms disappear
if(area2>1e-16)	//using truncated normal (although maybe frac can be zero?)
{
if(frac2>0){*pen+=frac2*(log(frac2/0.5/lam2/area2)-.5*log(2*M_PI*var2)-.5*(brun2+pow(trun2,2)-2*trun2*mean2+pow(mean2,2))/var2-lam2*trun2);}
}
//else will have trun2=brun2=0, and probably frac2=0, so all terms disappear
}

if(prob!=NULL){*prob=1;}
}

if(type==3)	//ridge - just one gaussian - possible that lam is zero or negative
{
if(lam>0)
{
postmean=sum/(dsq+resvar/lam);
postvar=resvar/(dsq+resvar/lam);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
*pen-=.5*(1+log(postvar/lam)-(postvar+pow(postmean,2))/lam);

if(*pen!=*pen){printf("dsq %f postvar %f resvar %f lam %f\n", dsq, postvar, resvar, lam);}
}

if(prob!=NULL){*prob=1;}
}
else
{
postmean=0;
if(prob!=NULL){*prob=0;}
}
}

if(type==4)	//bolt - note that lam2 can be zero (but pp always positive)
{
if(lam2>0)	//have two gaussians
{
mean=sum/(dsq+resvar/lam);
mean2=sum/(dsq+resvar/lam2);
var=resvar/(dsq+resvar/lam);
var2=resvar/(dsq+resvar/lam2);

max=pow(mean,2)/var;
if(pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
value=pp*pow(var/lam,.5)*exp(.5*pow(mean,2)/var-.5*max);
value2=pp2*pow(var2/lam2,.5)*exp(.5*pow(mean2,2)/var2-.5*max);

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*mean+frac2*mean2;
postvar=frac*(var+pow(mean,2))+frac2*(var2+pow(mean2,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
*pen-=.5*frac2*(1+log(var2/lam2)-(var2+pow(mean2,2))/lam2);
}

if(prob!=NULL){*prob=frac;}
}
else	//have gaussian and point mass
{
mean=sum/(dsq+resvar/lam);
var=resvar/(dsq+resvar/lam);

max=pow(mean,2)/var;
if(0>max){max=0;}
value=pp*pow(var/lam,.5)*exp(.5*pow(mean,2)/var-.5*max);
value2=pp2*exp(-.5*max);

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*mean;
postvar=frac*(var+pow(mean,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
}

if(prob!=NULL){*prob=frac;}
}
}

if(type==8)	//ldpred - saving stats
{
mean=sum/(dsq+resvar/lam);
var=resvar/(dsq+resvar/lam);

max=pow(mean,2)/var;
if(0>max){max=0;}
value=pp*pow(var/lam,.5)*exp(.5*pow(mean,2)/var-.5*max);
value2=pp2*exp(-.5*max);

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*mean;
postvar=frac*(var+pow(mean,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
}

prob[0]=frac;
prob[1]=mean;
prob[2]=var;
}

if(type==5)	//bayesr with three gaussians and a point mass (pp, pp2 and pp3 can each be zero)
{
mean2=sum/(dsq+resvar/lam2);
mean3=sum/(dsq+resvar/lam3);
mean4=sum/(dsq+resvar/lam4);
var2=resvar/(dsq+resvar/lam2);
var3=resvar/(dsq+resvar/lam3);
var4=resvar/(dsq+resvar/lam4);

max=0;
if(pp2>0&&pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
if(pp3>0&&pow(mean3,2)/var3>max){max=pow(mean3,2)/var3;}
if(pp4>0&&pow(mean4,2)/var4>max){max=pow(mean4,2)/var4;}

value=0;value2=0;value3=0;value4=0;
if(pp>0){value=pp*exp(-.5*max);}
if(pp2>0){value2=pp2*pow(var2/lam2,.5)*exp(.5*pow(mean2,2)/var2-.5*max);}
if(pp3>0){value3=pp3*pow(var3/lam3,.5)*exp(.5*pow(mean3,2)/var3-.5*max);}
if(pp4>0){value4=pp4*pow(var4/lam4,.5)*exp(.5*pow(mean4,2)/var4-.5*max);}

frac=value/(value+value2+value3+value4);
frac2=value2/(value+value2+value3+value4);
frac3=value3/(value+value2+value3+value4);
frac4=value4/(value+value2+value3+value4);
postmean=frac2*mean2+frac3*mean3+frac4*mean4;
postvar=frac2*(var2+pow(mean2,2))+frac3*(var3+pow(mean3,2))+frac4*(var4+pow(mean4,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
if(frac3>0){*pen+=frac3*log(frac3/pp3);}
if(frac4>0){*pen+=frac4*log(frac4/pp4);}
*pen-=.5*frac2*(1+log(var2/lam2)-(var2+pow(mean2,2))/lam2);
*pen-=.5*frac3*(1+log(var3/lam3)-(var3+pow(mean3,2))/lam3);
*pen-=.5*frac4*(1+log(var4/lam4)-(var4+pow(mean4,2))/lam4);
}

if(prob!=NULL){*prob=frac4;}
}

if(type==6)	//bayesr with four gaussians
{
mean=sum/(dsq+resvar/lam);
mean2=sum/(dsq+resvar/lam2);
mean3=sum/(dsq+resvar/lam3);
mean4=sum/(dsq+resvar/lam4);
var=resvar/(dsq+resvar/lam);
var2=resvar/(dsq+resvar/lam2);
var3=resvar/(dsq+resvar/lam3);
var4=resvar/(dsq+resvar/lam4);

//pp will always be >0
max=pow(mean,2)/var;
if(pp2>0&&pow(mean2,2)/var2>max){max=pow(mean2,2)/var2;}
if(pp3>0&&pow(mean3,2)/var3>max){max=pow(mean3,2)/var3;}
if(pp4>0&&pow(mean4,2)/var4>max){max=pow(mean4,2)/var4;}

value=0;value2=0;value3=0;value4=0;
if(pp>0){value=pp*exp(-.5*max);}
if(pp2>0){value2=pp2*pow(var2/lam2,.5)*exp(.5*pow(mean2,2)/var2-.5*max);}
if(pp3>0){value3=pp3*pow(var3/lam3,.5)*exp(.5*pow(mean3,2)/var3-.5*max);}
if(pp4>0){value4=pp4*pow(var4/lam4,.5)*exp(.5*pow(mean4,2)/var4-.5*max);}

frac=value/(value+value2+value3+value4);
frac2=value2/(value+value2+value3+value4);
frac3=value3/(value+value2+value3+value4);
frac4=value4/(value+value2+value3+value4);
postmean=frac*mean+frac2*mean2+frac3*mean3+frac4*mean4;
postvar=frac*(var+pow(mean,2))+frac2*(var2+pow(mean2,2))+frac3*(var3+pow(mean3,2))+frac4*(var4+pow(mean4,2))-pow(postmean,2);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(frac>0){*pen+=frac*log(frac/pp);}
if(frac2>0){*pen+=frac2*log(frac2/pp2);}
if(frac3>0){*pen+=frac3*log(frac3/pp3);}
if(frac4>0){*pen+=frac4*log(frac4/pp4);}
*pen-=.5*frac*(1+log(var/lam)-(var+pow(mean,2))/lam);
*pen-=.5*frac2*(1+log(var2/lam2)-(var2+pow(mean2,2))/lam2);
*pen-=.5*frac3*(1+log(var3/lam3)-(var3+pow(mean3,2))/lam3);
*pen-=.5*frac4*(1+log(var4/lam4)-(var4+pow(mean4,2))/lam4);
}

if(prob!=NULL){*prob=frac4;}
}

if(type==7)	//elastic - remember using lam=lam2 and pp=pp2 for lasso, lam3 and pp3 for ridge
{
if(pp3==0)	//lasso - copied from above
{
//these are the means and variances of the non-truncated distributions
mean=(sum-resvar*lam)/dsq;
mean2=(sum+resvar*lam2)/dsq;
var=resvar/dsq;
var2=resvar/dsq;

//for positive normal, area (denominator) is 1-Phi(-mean/var^.5) = Phi(mean/var^.5) = .5*erfc(-mean/var^.5 root(1/2))
//for negative normal, area (denominator) is Phi(-mean2/var2^.5) = .5*erfc(mean2/var2^.5 root(1/2))
area=.5*erfc(-mean*pow(var,-.5)*M_SQRT1_2);
area2=.5*erfc(mean2*pow(var2,-.5)*M_SQRT1_2);

//for means and variances, require the following
value=exp(-.5*pow(-mean,2)/var)*pow(2*M_PI,-.5);	//on wiki page, this is phi(alpha)
value2=exp(-.5*pow(-mean2,2)/var2)*pow(2*M_PI,-.5);	//on wiki page, this is phi(beta)

//get mean and variance of positive truncated distribution
//if area < 1e-16, will replace with point mass at zero
if(area>1e-16)
{
trun=mean+value/area*pow(var,.5);
brun=var*(1-mean*pow(var,-.5)*value/area-pow(value/area,2));
}
else{trun=0;brun=0;}

//get mean and variance of negative truncated distribution
//if area2 < 1e-16, will replace with point mass at zero
if(area2>1e-16)
{
trun2=mean2-value2/area2*pow(var2,.5);
brun2=var2*(1+mean2*pow(var2,-.5)*value2/area2-pow(-value2/area2,2));
}
else{trun2=0;brun2=0;}

//get fractions
if(area>1e-16&&area2>1e-16)	//normal case
{
max=log(area)+pow(mean,2)/var;
if(log(area2)+pow(mean2,2)/var2>max){max=log(area2)+pow(mean2,2)/var2;}
//here ignore pow(var,.5)=pow(var2,.5) in respective lines
value=exp(log(area)+.5*pow(mean,2)/var-.5*max);
value2=exp(log(area2)+.5*pow(mean2,2)/var2-.5*max);
}
else	//one or both areas are small - will set to 1-0, 0-1 or 1-1
{
if(area>1e-16&&area2<=1e-16){value=1;value2=0;}
if(area<=1e-16&&area2>1e-16){value=0;value2=1;}
if(area<=1e-16&&area2<=1e-16){value=1;value2=1;}
}

frac=value/(value+value2);
frac2=value2/(value+value2);
postmean=frac*trun+frac2*trun2;
postvar=frac*(brun+pow(trun,2))+frac2*(brun2+pow(trun2,2))-pow(postmean,2);

if(postmean!=postmean||isinf(postmean)){printf("Error adww44, please tell doug values %f %f means %f %f vars %f %f\n\n", value, value2, trun, trun2, brun, brun2);}

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(area>1e-16)	//using truncated normal (although maybe frac can be zero?)
{
if(frac>0){*pen+=frac*(log(frac/pp/lam/area)-.5*log(2*M_PI*var)-.5*(brun+pow(trun,2)-2*trun*mean+pow(mean,2))/var+lam*trun);}
}
//else will have trun=brun=0, and probably frac=0, so all terms disappear
if(area2>1e-16)	//using truncated normal (although maybe frac can be zero?)
{
if(frac2>0){*pen+=frac2*(log(frac2/pp2/lam2/area2)-.5*log(2*M_PI*var2)-.5*(brun2+pow(trun2,2)-2*trun2*mean2+pow(mean2,2))/var2-lam2*trun2);}
}
//else will have trun2=brun2=0, and probably frac2=0, so all terms disappear
}

if(isinf(*pen)){printf("aracs %f %f means %f %f, vars %f %f brusn %f %f and %f %f area %e %e reg %f, lams %f %f logs %f %f\n", frac, frac2,mean, mean2, var, var2, trun, trun2, brun, brun2, area, area2, sum/dsq, lam, lam2, frac/pp/area/lam, log(lam));exit(1);}

if(prob!=NULL){*prob=1;}
}

if(pp3==1)	//ridge - copied from above
{
postmean=sum/(dsq+resvar/lam3);
postvar=resvar/(dsq+resvar/lam3);

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
*pen-=.5*(1+log(postvar/lam3)-(postvar+pow(postmean,2))/lam3);
}

if(prob!=NULL){*prob=1;}
}

if(pp3>0&&pp3<1)	//elastic
{
//these are the means and variances of the non-truncated distributions
mean=(sum-resvar*lam)/dsq;
mean2=(sum+resvar*lam2)/dsq;
var=resvar/dsq;
var2=resvar/dsq;

//for positive normal, area (denominator) is 1-Phi(-mean/var^.5) = Phi(mean/var^.5) = .5*erfc(-mean/var^.5 root(1/2))
//for negative normal, area (denominator) is Phi(-mean2/var2^.5) = .5*erfc(mean2/var2^.5 root(1/2))
area=.5*erfc(-mean*pow(var,-.5)*M_SQRT1_2);
area2=.5*erfc(mean2*pow(var2,-.5)*M_SQRT1_2);

//for means and variances, require the following
value=exp(-.5*pow(-mean,2)/var)*pow(2*M_PI,-.5);	//on wiki page, this is phi(alpha)
value2=exp(-.5*pow(-mean2,2)/var2)*pow(2*M_PI,-.5);	//on wiki page, this is phi(beta)

//get mean and variance of positive truncated distribution
//if area < 1e-16, will replace with point mass at zero
if(area>1e-16)
{
trun=mean+value/area*pow(var,.5);
brun=var*(1-mean*pow(var,-.5)*value/area-pow(value/area,2));
}
else{trun=0;brun=0;}

//get mean and variance of negative truncated distribution
//if area2 < 1e-16, will replace with point mass at zero
if(area2>1e-16)
{
trun2=mean2-value2/area2*pow(var2,.5);
brun2=var2*(1+mean2*pow(var2,-.5)*value2/area2-pow(-value2/area2,2));
}
else{trun2=0;brun2=0;}

//get mean and variance corresponding to ridge prior
mean3=sum/(dsq+resvar/lam3);
var3=resvar/(dsq+resvar/lam3);

//get fractions
max=log(area)+pow(mean,2)/var;
if(log(area2)+pow(mean2,2)/var2>max){max=log(area2)+pow(mean2,2)/var2;}
if(pow(mean3,2)/var3>max){max=pow(mean3,2)/var3;}

value=pp*lam*pow(2*M_PI*var,.5)*exp(log(area)+.5*pow(mean,2)/var-.5*max);
value2=pp2*lam2*pow(2*M_PI*var2,.5)*exp(log(area2)+.5*pow(mean2,2)/var2-.5*max);
value3=pp3*pow(var3/lam3,.5)*exp(.5*pow(mean3,2)/var3-.5*max);

frac=value/(value+value2+value3);
frac2=value2/(value+value2+value3);
frac3=value3/(value+value2+value3);

postmean=frac*trun+frac2*trun2+frac3*mean3;
postvar=frac*(brun+pow(trun,2))+frac2*(brun2+pow(trun2,2))+frac3*(var3+pow(mean3,2))-pow(postmean,2);

if(postmean!=postmean||isinf(postmean)){
printf("sum is %f dsq %f\n", sum, dsq);
printf("lams %e %e %e pps %e %e %e and %f\n", lam, lam2, lam3, pp3, var3, var3/lam3, value3);
printf("Error adww55, please tell doug values %f %f %f means %f %f %f vars %e %e %e areas %e %e reg %f diff %f and %f\n\n", value, value2, value3, mean, mean2, mean3, var, var2, var3, area, area2, sum/dsq, .5*pow(mean,2)/var-.5*max, exp(.5*pow(mean2,2)/var2-.5*max));}

if(pen!=NULL)
{
*pen+=.5*dsq*postvar/resvar;
if(area>1e-16)	//using truncated normal (although maybe frac can be zero?)
{
if(frac>0){*pen+=frac*(log(frac/pp/lam/area)-.5*log(2*M_PI*var)-.5*(brun+pow(trun,2)-2*trun*mean+pow(mean,2))/var+lam*trun);}
}
//else will have trun=brun=0, and probably frac=0, so all terms disappear
if(area2>1e-16)	//using truncated normal (although maybe frac can be zero?)
{
if(frac2>0){*pen+=frac2*(log(frac2/pp2/lam2/area2)-.5*log(2*M_PI*var2)-.5*(brun2+pow(trun2,2)-2*trun2*mean2+pow(mean2,2))/var2-lam2*trun2);}
}
//else will have trun2=brun2=0, and probably frac2=0, so all terms disappear
if(frac3>0){*pen+=frac3*log(frac3/pp3);}
*pen-=.5*frac3*(1+log(var3/lam3)-(var3+pow(mean3,2))/lam3);

if(isinf(*pen)){printf("bracs %f %f means %f %f, vars %f %f brusn %f %f and %f %f area %e %e reg %f, lams %f %f logs %f %f\n", frac, frac2,mean, mean2, var, var2, trun, trun2, brun, brun2, area, area2, sum/dsq, lam, lam2, frac/pp/area/lam, log(lam));exit(1);}
}

if(prob!=NULL){*prob=frac+frac2;}
}
}

if(pen!=NULL)
{
if(*pen!=*pen||isinf(*pen)){printf("Warning, pen is nan, please tell Doug - sum %f lam %e lam2 %e lam3 %e lam4 %e dsq %f resvar %f pp %f pp2 %f pp3 %f pp4 %f type %d\n", sum, lam, lam2, lam3, lam4, dsq, resvar, pp, pp2, pp3, pp4, type);}
}

return(postmean);
}	//end of get_postmean

///////////////////////////

double inter_hers(int nhers, double *polates, double *tryhers, int step)
{
int left, right;
double value;


//int j;for(j=0;j<nhers;j++){printf("j %d her %f ratio %f\n", j+1, tryhers[j], polates[j*step]);}


//work out closest knots
if(polates[0*step]<polates[(nhers-1)*step])	//positive gradient
{
if(polates[0*step]<1&&polates[(nhers-1)*step]>1)	//knot inside
{
left=0;
while(polates[(left+1)*step]<1){left++;}
right=left+1;
}
else	//knot outside
{
if(polates[0*step]>1){left=0;right=1;}
else{left=nhers-2;right=nhers-1;}
}
}
else	//negative gradient 
{
if(polates[0*step]>1&&polates[(nhers-1)*step]<1)	//knot inside
{
left=0;
while(polates[(left+1)*step]>1){left++;}
right=left+1;
}
else	//knot outside
{
if(polates[0*step]<1){left=0;right=1;}
else{left=nhers-2;right=nhers-1;}
}
}

value=(tryhers[right]*(1-polates[left*step])-tryhers[left]*(1-polates[right*step]))/(polates[right*step]-polates[left*step]);
if(value!=value||isinf(value)){value=-9999;}

return(value);
}

///////////////////////////

void elastic_net(double *effs, double *cors, double alpha, int length, double *XTX, double *XTY, double *XTY2, int num_lambdas, double neff, double tol, int maxiter)
//penalty is lambda1 sum (abs(betas)) + lambda2/2 sum (betas^2)
//lambda1 = lambda*alpha, lambda2 = lambda*(1-alpha) - alpha=0 is ridge
{
int j, j2, p, count, ecount;
double max, value, value2;
double lambdamax, *lambdas, lambda1, lambda2;
double *betas, sum, sumsq, pen, like, likeold;


lambdas=malloc(sizeof(double)*num_lambdas);
betas=malloc(sizeof(double)*length);

//get max (abs(XTY))
max=0;
for(j=0;j<length;j++)
{
if(fabs(XTY[j])>max){max=fabs(XTY[j]);}
}

//now lambdamax
if(alpha>0){lambdamax=max/alpha;}
else{lambdamax=max;}

//get sequence
value=log(.00001)/num_lambdas;
for(p=0;p<num_lambdas;p++){lambdas[p]=lambdamax*exp(value*(1+p));}

//set betas to zero
for(j=0;j<length;j++){betas[j]=0;}

ecount=0;
for(p=0;p<num_lambdas;p++)
{
lambda1=alpha*lambdas[p];
lambda2=(1-alpha)*lambdas[p];

//get likelihood
sumsq=0.5;
for(j=0;j<length;j++)
{
sumsq-=betas[j]*XTY[j];
for(j2=0;j2<length;j2++){sumsq+=0.5*betas[j]*XTX[j+j2*length]*betas[j2];}
}

pen=0;
for(j=0;j<length;j++){pen+=lambda1*fabs(betas[j])+.5*lambda2*pow(betas[j],2);}
likeold=neff*(sumsq+pen);

count=0;
while(1)
{
count++;

//update betas
for(j=0;j<length;j++)
{
sum=XTY[j];
for(j2=0;j2<length;j2++)
{
if(j2!=j){sum-=XTX[j+j2*length]*betas[j2];}
}

betas[j]=0;
if(sum>lambda1){betas[j]=(sum-lambda1)/(1+lambda2);}
if(sum<-lambda1){betas[j]=(sum+lambda1)/(1+lambda2);}
}

//get likelihood
sumsq=0.5;
for(j=0;j<length;j++)
{
sumsq-=betas[j]*XTY[j];
for(j2=0;j2<length;j2++){sumsq+=0.5*betas[j]*XTX[j+j2*length]*betas[j2];}
}

pen=0;
for(j=0;j<length;j++){pen+=lambda1*fabs(betas[j])+.5*lambda2*pow(betas[j],2);}
like=neff*(sumsq+pen);

if(fabs(like-likeold)<tol){break;}
likeold=like;

if(count==maxiter)
{
if(ecount<5){printf("Warning, Run %d did not converge within %d iterations\n", p+1, maxiter);}
ecount++;
break;
}
}

//load up betas
for(j=0;j<length;j++){effs[j+p*length]=betas[j];}

if(cors!=NULL)	//compute correlation
{
//get sumj (betaj XjTY)
value=0;
for(j=0;j<length;j++){value+=effs[j+p*length]*XTY2[j];}

//get sumjj2 (betaj betaj2 XjTXj2)
value2=0;
for(j=0;j<length;j++)
{
for(j2=0;j2<length;j2++){value2+=effs[j+p*length]*XTX[j+j2*length]*effs[j2+p*length];}
}

//get cor
if(value2>0){cors[p]=value*pow(value2,-.5);}
else{cors[p]=-9999;}
}
}	//end of p loop

if(ecount>0){printf("In total, %d runs did not converge\n\n", ecount);}

free(lambdas);free(betas);
}

///////////////////////////

double get_like_trun(double obs, double real, double sd, double thresh)
{
//returns truncated likelihood L(obs|N(real,sd) and obs>thresh)
double va, vb, pa, pb, value;

va=(-thresh-real)/sd;
pa=.5*erfc(-va*M_SQRT1_2);

vb=(thresh-real)/sd;
pb=.5*erfc(vb*M_SQRT1_2);

value=-.5*pow((obs-real)/sd,2)-log(pa+pb);

return(value);
}

////////

double get_first_trun(double obs, double real, double sd, double thresh)
{
//returns the first derivative of the truncated likelihood
double va, vb, pa, pb, value;

va=(-thresh-real)/sd;
pa=.5*erfc(-va*M_SQRT1_2);

vb=(thresh-real)/sd;
pb=.5*erfc(vb*M_SQRT1_2);

value=(obs-real)*pow(sd,-2)-pow(pa+pb,-1)*pow(2*M_PI,-.5)/sd*(-exp(-.5*pow(va,2))+exp(-.5*pow(vb,2)));

return(value);
}

////////

double get_second_trun(double obs, double real, double sd, double thresh)
{
//returns the second derivative of the truncated likelihood
double va, vb, pa, pb, value, value2;

va=(-thresh-real)/sd;
pa=.5*erfc(-va*M_SQRT1_2);

vb=(thresh-real)/sd;
pb=.5*erfc(vb*M_SQRT1_2);

value2=pow(2*M_PI,-.5)/sd*(-exp(-.5*pow(va,2))+exp(-.5*pow(vb,2)));
value=-real*pow(sd,-2)+pow(pa+pb,-2)*pow(value2,2)-pow(pa+pb,-1)*pow(2*M_PI,-.5)*pow(sd,-3)*(exp(-.5*pow(va,2))*(thresh+real)+exp(-.5*pow(vb,2))*(thresh-real));

return(value);
}

////////

double secant_trun(double obs, double sd, double thresh, double tol, int maxiter, int type)
{
//finds value that maximizes truncated likelihood
//type=0 - quiet, type=1 - warn
int count;
double x0, x1, x2, f0, f1;

x0=0;
x1=obs;

f0=get_first_trun(obs, x0, sd, thresh);
f1=get_first_trun(obs, x1, sd, thresh);

count=0;
while(fabs(f1)>tol)
{
x2=x1-f1*(x1-x0)/(f1-f0);
x0=x1;
f0=f1;
x1=x2;
f1=get_first_trun(obs, x1, sd, thresh);

count++;
if(count==maxiter)
{
if(type==1){printf("Failed to converge, values are %f %f and deriv is %f thre %f\n", obs, x1, f1, thresh);}
break;
}
}	//end of while loop

return(x1);
}

///////////////////////////


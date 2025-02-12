/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Solves generalized REML - variance matrix constructed from noise + kinships + regions (region can be a gene)
//shortcut=0 is standard way, shortcut=1 uses an eigendecomposition, shortcut=2 is fast-reml
//shortcut=3 is when using only the diagonals of kinship matrices

///////////////////////////

void multi_reml(int ns, int num_covars, int num_envs, int num_tops, int num_kins, int num_regs, double *Y, double *Z, double **Mkins, float **Mkins_single, double *kintraces, double *kinsums, double *sweights, double *X, int Xtotal, int *Xstarts, int *Xends, double *Xsums, double prev, int np, char *hersfile, int hestart, int shortcut, double *U, double *E, int num_vects, int ldlt, int discenv, char *oversfile, double **ssums, int constrain, double tol, int maxiter, int memsave, int maxthreads, char **kinstems, char **ids3, char *outfile, char **ids1, char **ids2, int num_preds, char **allpreds, char *allal1, char *allal2, int *rkeeppreds, int **regindex, double *rcentres, double *rmults, double *rweights, int *tkeeppreds, double *tcentres, double *vstarts, double missingvalue)
{
size_t scount, stotal;
int i, i2, j, j2, k, k2, k3, r, g, count, token, flag, one=1, info, *ipiv, lwork;
float wkopt_single, *work_single, alpha_single, beta_single;
double value, value2, sum, sum2, sumsq, max, alpha, beta;

int num_fixed, total, total2, cflag, rflag, *fixed, *fixedsave, nlost, nlost2;
double nfree, varnull, varmin, gam, gam2, gam3, relax;
double likenull, like, likeold, diff, lrtstat, lrtpva, covher, topher, factor;
double *scales, *vars, *vardiffs, *varsds, *hers, *hersds, *shares, *sharesds;

double *AI, *AI2, *AI3, *BI, *J, *JAI, *JAIJT, *cohers, *cohers2;
double *Z2, *ZTY, *ZTZ, *ZTZ2, *ZTZZTY, detZTZ;
double detV, *ZTVZ, *ZTVZ2, *ZTVZ3, detZTVZ, *PY, **KPY, **PKPY, *traces;
double *ZTVY, *thetas, *thetasds, *thetapvas, *Yadj;
double **mg, **mg2, **effects;

float *R_single, **KR_single, *V_full, *V_single, *RHS_single, *PY_single, *KPY_single, *RHS2_single;
double *V, *VZ, *VY, *VR, *VZZTVZ, *P, *ZTVR, *PR, *VKPY, *ZTVKPY;
double *kin_diags;
double *UTY, *UTZ, *D, detD, *BUTZ, *H, *HUTY, *HKPY;
double *PX, *XTPY;
double *UTX, *DUTX, detC, *XTVX, *XTVX2, *XTVX3, detXTVX, *F, *FUTZ, *FUTY, *FUTX, *HUTX, *FKPY;


char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500], filename10[500];
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7, *output8, *output9, *output10;

//variables for Andys software
int *Sindexer, *Snums;
double *Sscales, *SJ, *Shers, *Shersds;
char **wantids;


//set num_fixed, nfree, total and stotal
num_fixed=num_covars+num_envs+num_tops;
nfree=ns-num_fixed;
total=1+num_kins+num_regs;
total2=num_fixed+1+num_vects;
stotal=(size_t)ns*ns;

//allocate variables

fixed=malloc(sizeof(int)*total);
fixedsave=malloc(sizeof(int)*total);
vars=malloc(sizeof(double)*total);
vardiffs=malloc(sizeof(double)*total);
varsds=malloc(sizeof(double)*total);
scales=malloc(sizeof(double)*total);
hers=malloc(sizeof(double)*(total+1));
hersds=malloc(sizeof(double)*(total+1));
shares=malloc(sizeof(double)*total);
sharesds=malloc(sizeof(double)*total);

AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);
J=malloc(sizeof(double)*total*total);
JAI=malloc(sizeof(double)*total*total);
JAIJT=malloc(sizeof(double)*total*total);
cohers=malloc(sizeof(double)*total*total);
cohers2=malloc(sizeof(double)*total*total);

Z2=malloc(sizeof(double)*ns*num_fixed);
ZTY=malloc(sizeof(double)*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZZTY=malloc(sizeof(double)*num_fixed);

ZTVZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTVZ2=malloc(sizeof(double)*num_fixed);
ZTVZ3=malloc(sizeof(double)*num_fixed*num_fixed);
PY=malloc(sizeof(double)*ns);
KPY=malloc(sizeof(double*)*total);
for(k=0;k<total;k++){KPY[k]=malloc(sizeof(double)*ns);}
PKPY=malloc(sizeof(double*)*total);
for(k=0;k<total;k++){PKPY[k]=malloc(sizeof(double)*ns);}
traces=malloc(sizeof(double)*total);

ZTVY=malloc(sizeof(double)*num_fixed);
thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*ns);

if(num_kins+num_regs>0)
{
mg=malloc(sizeof(double*)*(num_kins+num_regs));
for(k=0;k<num_kins+num_regs;k++){mg[k]=malloc(sizeof(double)*ns);}
mg2=malloc(sizeof(double*)*(num_kins+num_regs));
for(k=0;k<num_kins+num_regs;k++){mg2[k]=malloc(sizeof(double)*ns);}
}

if(num_regs>0||num_tops>0)
{
effects=malloc(sizeof(double*)*(num_regs+4));
for(r=0;r<num_regs+4;r++){effects[r]=malloc(sizeof(double)*num_preds);}
}

if(shortcut==0)
{
V=malloc(sizeof(double)*ns*ns);
VZ=malloc(sizeof(double)*ns*num_fixed);
VZZTVZ=malloc(sizeof(double)*ns*num_fixed);
P=malloc(sizeof(double)*ns*ns);
}

if(shortcut==1)
{
UTY=malloc(sizeof(double)*ns);
UTZ=malloc(sizeof(double)*ns*num_fixed);
D=malloc(sizeof(double)*ns);
BUTZ=malloc(sizeof(double)*ns*num_fixed);
H=malloc(sizeof(double)*num_fixed*ns);
HUTY=malloc(sizeof(double)*num_fixed);
HKPY=malloc(sizeof(double)*num_fixed);
}

if(shortcut==2)
{
R_single=malloc(sizeof(float)*ns*num_vects);
KR_single=malloc(sizeof(float*)*num_kins);
for(k=0;k<num_kins;k++){KR_single[k]=malloc(sizeof(float)*ns*num_vects);}
V_full=malloc(sizeof(float)*ns*(ns+1));
V_single=V_full+ns;
RHS_single=malloc(sizeof(float)*ns*total2);
PY_single=malloc(sizeof(float)*ns);
KPY_single=malloc(sizeof(float)*ns);
RHS2_single=malloc(sizeof(float)*ns*total);

VZ=malloc(sizeof(double)*ns*num_fixed);
VY=malloc(sizeof(double)*ns);
VR=malloc(sizeof(double)*ns*num_vects);
VZZTVZ=malloc(sizeof(double)*ns*num_fixed);
ZTVR=malloc(sizeof(double)*num_fixed*num_vects);
PR=malloc(sizeof(double)*ns*num_vects);
VKPY=malloc(sizeof(double)*ns*total);
ZTVKPY=malloc(sizeof(double)*num_fixed);
}

if(shortcut==3)
{
kin_diags=malloc(sizeof(double)*num_kins*ns);
D=malloc(sizeof(double)*ns);
BUTZ=malloc(sizeof(double)*ns*num_fixed);
H=malloc(sizeof(double)*num_fixed*ns);
HUTY=malloc(sizeof(double)*num_fixed);
HKPY=malloc(sizeof(double)*num_fixed);
}

if(Xtotal>0)
{
PX=malloc(sizeof(double)*ns*Xtotal);
XTPY=malloc(sizeof(double)*Xtotal);

if(shortcut==1)
{
UTX=malloc(sizeof(double)*ns*Xtotal);
DUTX=malloc(sizeof(double)*ns*Xtotal);
XTVX=malloc(sizeof(double)*Xtotal*Xtotal);
XTVX2=malloc(sizeof(double)*Xtotal);
XTVX3=malloc(sizeof(double)*Xtotal*Xtotal);
F=malloc(sizeof(double)*Xtotal*ns);
FUTZ=malloc(sizeof(double)*Xtotal*num_fixed);
FUTY=malloc(sizeof(double)*Xtotal);
FUTX=malloc(sizeof(double)*Xtotal*Xtotal);
HUTX=malloc(sizeof(double)*num_fixed*Xtotal);
FKPY=malloc(sizeof(double)*Xtotal);
}
}

//fill some variables

if(shortcut==1)	//get UTY, UTZ and maybe UTX
{
if(num_kins==1)
{
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &ns, &alpha, U, &ns, Y, &one, &beta, UTY, &one);
dgemm_("T", "N", &ns, &num_fixed, &ns, &alpha, U, &ns, Z, &ns, &beta, UTZ, &ns);
if(Xtotal>0)
{dgemm_("T", "N", &ns, &Xtotal, &ns, &alpha, U, &ns, X, &ns, &beta, UTX, &ns);}
}
else
{
for(i=0;i<ns;i++)
{
UTY[i]=Y[i];
for(j=0;j<num_fixed;j++){UTZ[i+j*ns]=Z[i+j*ns];}
for(j=0;j<Xtotal;j++){UTX[i+j*ns]=X[i+j*ns];}
}
}
}

if(shortcut==2)	//load random vectors into R_single, and get KR_single for each kinship
{
for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){R_single[i+g*ns]=rnorm_safe();}
}

printf("Multiplying each kinship matrix by the random vectors\n\n");
for(k=0;k<num_kins;k++)
{
if(memsave==0)	//kins already stored
{
alpha_single=1.0;beta_single=0.0;
ssymm_("L", "U", &ns, &num_vects, &alpha_single, Mkins_single[k], &ns, R_single, &ns, &beta_single, KR_single[k], &ns);
}
else	//must read kins
{
read_kins(kinstems[k], NULL, V_full, 1.0, ns, ids3, 4, maxthreads);

alpha_single=1.0;beta_single=0.0;
ssymm_("L", "U", &ns, &num_vects, &alpha_single, V_full, &ns, R_single, &ns, &beta_single, KR_single[k], &ns);
}
}
}

if(shortcut==3)	//read kin diagonals
{
for(k=0;k<num_kins;k++)
{(void)read_kin_trace(kinstems[k], ns, ids3, kin_diags+k*ns, 1, 1, 1);}
}

////////

//solve model with just fixed effects to get varnull, varmin and null likelihood
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){Z2[i+j*ns]=Z[i+j*ns]*sweights[i];}
}

alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z2, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z2, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
for(j=0;j<num_fixed;j++){ZTZZTY[j]=ZTY[j];}
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, 1, ZTZZTY, 1);

sumsq=0;
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
for(j=0;j<num_fixed;j++){sumsq-=ZTY[j]*ZTZZTY[j];}
varnull=sumsq/nfree;
varmin=0.0001*varnull;

sum=0;for(i=0;i<ns;i++){sum+=log(sweights[i]);}
likenull=-.5*(nfree+nfree*log(2*M_PI*varnull)+detZTZ-sum);

//set scales
sum=0;for(i=0;i<ns;i++){sum+=pow(sweights[i],-1);}
scales[0]=sum/ns;
for(k=0;k<num_kins;k++){scales[1+k]=kintraces[k];}
for(r=0;r<num_regs;r++){scales[1+num_kins+r]=1;}

//set starting vars, hers, shares and fixed
if(total==1){hers[0]=1;}
else
{
if(strcmp(hersfile,"blank")!=0)	//read from file - have already checked correct size and their sum
{
read_values(hersfile, hers+1, num_kins+num_regs, NULL, 1, 0, 0);
sum=0;for(k=1;k<total;k++){sum+=hers[k];}
hers[0]=1-sum;
}
else
{
if(hestart==1)	//set based on he regression
{
printf("Computing starting heritabilities\n");
if(shortcut==0||shortcut==1){he_starts(hers, ns, num_covars, num_envs, num_tops, num_kins, num_regs, Y, Z, Mkins, NULL, kintraces, X, Xtotal, Xstarts, Xends, Xsums, memsave, maxthreads, kinstems, ids3, missingvalue);}
else{he_starts(hers, ns, num_covars, num_envs, num_tops, num_kins, num_regs, Y, Z, NULL, Mkins_single, kintraces, X, Xtotal, Xstarts, Xends, Xsums, memsave, maxthreads, kinstems, ids3, missingvalue);}

//make sure no genetic heritabilities are very small or negative (noise checked in next step)
for(k=1;k<total;k++)
{
if(hers[k]<0.01){hers[0]+=hers[k]-0.01;hers[k]=0.01;}
}

//ensure total her is not too large 0.95
sum=0;for(k=1;k<total;k++){sum+=hers[k];}
if(sum>.95)
{
hers[0]=.05;
for(k=1;k<total;k++){hers[k]=hers[k]*.95/sum;}
}
}
else	//set agnostically
{
hers[0]=.5;
for(k=1;k<total;k++){hers[k]=.5/(total-1);}
}
}
}	//end of total!=1

for(k=0;k<total;k++){vars[k]=hers[k]*varnull/scales[k];}
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}
for(k=0;k<total;k++){fixed[k]=0;}

//prepare to screen and file print
printf("Iter\t");
for(k=0;k<num_kins;k++){printf("Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){printf("Her_R%d\t", r+1);}
printf("Her_All\tLikelihood\tDifference\tTarget\tNum_Constrained\n");

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Iteration\t");
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d\t", r+1);}
fprintf(output, "Her_All\tTotal_Variance\tLikelihood\tDifference\tTarget\tNum_Constrained\n");
fclose(output);

////////

//now iterate
count=0;
cflag=1;	//flips to 0 if likelihood does not converge
rflag=0;	//0 - normal, 1 - reduced, 2 - transitioning, 3 - normal (after reduced), 4 - about to stop
while(1)
{
//compute detV, invZTVZ, detZTVZ, PY and gam (plus other variables depending on shortcut)
#include "remllike.c"

//get likelihood
like=-.5*gam-.5*detZTVZ-.5*detV-.5*nfree*log(2*M_PI);

if(count>0)	//see if move was good
{
if(like>likeold-tol)	//good move - decide next move type
{
switch(rflag)
{
case 0:rflag=0;break;
case 1:rflag=2;break;
case 2:rflag=3;break;
case 3:rflag=3;break;
}
}
else	//bad move - decide next move type, return to previous state, and maybe recompute likelihood terms
{
switch(rflag)
{
case 0:rflag=1;break;
case 1:rflag=1;break;
case 2:rflag=4;break;
case 3:rflag=4;break;
}

for(k=0;k<total;k++){fixed[k]=fixedsave[k];}
for(k=0;k<total;k++){vars[k]-=relax*vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}

if(count==maxiter||rflag==4)	//about to finish, so need to recompute likelihood terms
{
printf("Warning, the last move reduced the likelihood, so will return to previous state\n");
#include "remllike.c"
like=-.5*gam-.5*detZTVZ-.5*detV-.5*nfree*log(2*M_PI);
}
else	//will do a reduced move - so no need to recompute likelihood terms
{
printf("Warning, the last move reduced the likelihood, so will try a smaller move\n");
like=likeold;
}
}

diff=like-likeold;
}	//end of count>0

likeold=like;
for(k=0;k<total;k++){fixedsave[k]=fixed[k];}

if(count==maxiter||rflag==0||rflag==2||rflag==3||rflag==4)	//need to compute derivatives (not required if doing a reduced move)
{
//compute PX, XTPY, KPY, PKPY, (inverse) AI and BI (plus other variables, depending on shortcut)
#include "remlderiv.c"
}

//print update
sum=0;for(k=0;k<num_kins+num_regs;k++){sum+=hers[1+k];}
value=0;for(k=0;k<total;k++){value+=vars[k];}
nlost=0;for(k=0;k<total;k++){nlost+=(fixed[k]>=3);}
nlost2=0;for(k=0;k<total;k++){nlost2+=(fixed[k]==1||fixed[k]==2);}

if(count==0){printf("Start\t");}
else{printf("%d\t", count);}
for(k=0;k<num_kins+num_regs;k++){printf("%.4f\t", hers[1+k]);}
printf("%.4f\t", sum);
if(rflag==0||rflag==2||rflag==3){printf("%.2f\t", like);}
else{printf("%.2f*\t", like);}
if(count==0){printf("n/a\t\t%.6f\t", tol);}
else{printf("%.6f\t%.6f\t", diff, tol);}
if(nlost2==0){printf("%d\n", nlost);}
else{printf("%d*\n", nlost);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
fprintf(output, "%d\t", count);
for(k=0;k<num_kins+num_regs;k++){fprintf(output, "%.6f\t", hers[1+k]);}
if(count==0){fprintf(output, "%.6f\t%.6f\t%.6f\tNA\t%.6f\t%d\n", sum, value, like, tol, nlost);}
else{fprintf(output, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n", sum, value, like, diff, tol, nlost);}
fclose(output);

//see if breaking (normally can only break if rflag=0|3 and nlost2=0, unless at iter limit)

if(num_kins+num_regs==0){break;}	//null model
if(nlost>=num_kins+num_regs){printf("All heritabilities are constrained\n");break;}
if(rflag==4){cflag=0;break;}	//giving up

if(count>0)
{
if(fabs(diff)<tol&&(rflag==0||rflag==3)&&nlost2==0){break;}
}

if(count==maxiter)
{
if(shortcut==0||shortcut==1||shortcut==3){printf("\nWarning, REML failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)\n", maxiter, tol);}
if(shortcut==2){printf("\nWarning, REML failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e), or using \"--repetitions\" to increase the number of random samplings (currently %d)\n", maxiter, tol, num_vects);}
if(constrain==0)
{printf("Additionally, it might help to use \"--constrain YES\" to constrain heritabilities to [0,1]\n");}
else
{printf("Additionally, it might help to use \"--constrain NO\" to allow heritabilities outside [0,1]\n");}
cflag=0;
break;
}

////////

//update variances using NR

if(rflag==0||rflag==2||rflag==3)	//doing a normal move
{
relax=1;

//get proposed moves
alpha=1.0;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, vardiffs, &one);

if(constrain==1)	//new variances can not be negative
{
for(k=0;k<total;k++)
{
if(fixed[k]<3)	//free to update
{
if(vars[k]+vardiffs[k]<varmin){vardiffs[k]=varmin-vars[k];fixed[k]++;}
else{fixed[k]=0;}
if(fixed[k]==3){vardiffs[k]=-vars[k];}
}}
}

//make sure moves not too large
max=0;
for(k=0;k<total;k++)
{
if(fabs(vardiffs[k])>max){max=fabs(vardiffs[k]);}
}
if(max>.1*varnull)	//reduce (all) moves
{
for(k=0;k<total;k++){vardiffs[k]*=.1*varnull/max;}
}
}
else	//doing a reduced version of previous move
{
if(shortcut==2){relax*=.1;}	//will have large sample size
else{relax*=.5;}
}

//move
for(k=0;k<total;k++){vars[k]+=relax*vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}

count++;
}	//end of while loop
printf("\n");

//get some stats
lrtstat=2*(like-likenull);
lrtpva=.5*erfc(pow(lrtstat,.5)*M_SQRT1_2);
if(lrtstat<0){lrtpva=.75;}
//if(hers[0]<=0){lrtpva=1;}

if(vstarts!=NULL)	//must return genetic variances, scaled so they sum to one (will only have kinships)
{
sum=0;for(k=0;k<num_kins;k++){sum+=vars[1+k];}
for(k=0;k<num_kins;k++){vstarts[k]=vars[1+k]/sum;}
}

////////

//get fixed coefficients (ZTinvVZ)^-1 ZTinvVY with variance matrix (ZTinvVZ)^-1
alpha=1.0;beta=0.0;
if(shortcut==0){dgemv_("T", &ns, &num_fixed, &alpha, VZ, &ns, Y, &one, &beta, ZTVY, &one);}
if(shortcut==1){dgemv_("T", &ns, &num_fixed, &alpha, BUTZ, &ns, UTY, &one, &beta, ZTVY, &one);}
//for shortcut=2 already have ZTVY
if(shortcut==3){dgemv_("T", &ns, &num_fixed, &alpha, BUTZ, &ns, Y, &one, &beta, ZTVY, &one);}
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTVZ, &num_fixed, ZTVY, &one, &beta, thetas, &one);

for(j=0;j<num_fixed;j++)
{
if(ZTVZ[j+j*num_fixed]>=0)
{
thetasds[j]=pow(ZTVZ[j+j*num_fixed],.5);
thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);
}
else{thetasds[j]=-9999;thetapvas[j]=-9999;}
}

//set Yadj to contribution of fixed effects
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Yadj, &one);

for(k=0;k<num_kins+num_regs;k++)	//random effects are Kg invV Yadj = g KPY - get also g PY (redundant for regions)
{
if(shortcut==0||shortcut==2||shortcut==3)
{
for(i=0;i<ns;i++){mg[k][i]=vars[1+k]*KPY[1+k][i];mg2[k][i]=vars[1+k]*PY[i];}
}
if(shortcut==1)
{
alpha=vars[1+k];beta=0.0;
dgemv_("N", &ns, &ns, &alpha, U, &ns, KPY[1+k], &one, &beta, mg[k], &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, PY, &one, &beta, mg2[k], &one);
}

//add contribution of random effects to Yadj
for(i=0;i<ns;i++){Yadj[i]+=mg[k][i];}
}

////////

if(num_regs>0||num_tops>0)	//load up effects for regions (equal to XTPY var/sum) and tops (in thetas)
{
//effects contain effects for regions and tops, then sum, centres, predictor used?
for(r=0;r<num_regs+4;r++)
{
for(j=0;j<num_preds;j++){effects[r][j]=0;}
}

for(r=0;r<num_regs;r++)	//note, will have removed trivial regional predictors and those with zero weight
{
if(vars[1+num_kins+r]!=0)
{
for(j=Xstarts[r];j<Xends[r];j++)
{
j2=regindex[r][1+j-Xstarts[r]];
effects[r][rkeeppreds[j2]]+=XTPY[j]*vars[1+num_kins+r]/Xsums[r]*rmults[j2]*pow(rweights[j2],.5);
effects[num_regs+1][rkeeppreds[j2]]+=effects[r][rkeeppreds[j2]];
effects[num_regs+2][rkeeppreds[j2]]=rcentres[j2];
effects[num_regs+3][rkeeppreds[j2]]=1;
}
}}

for(j=0;j<num_tops;j++)
{
effects[num_regs][tkeeppreds[j]]=thetas[num_covars+num_envs+j];
effects[num_regs+1][tkeeppreds[j]]+=effects[num_regs][tkeeppreds[j]];
effects[num_regs+2][tkeeppreds[j]]=tcentres[j];
effects[num_regs+3][tkeeppreds[j]]=1;
}
}

////////

//get SEs - for transformed variances, must compute J invAI JT, where Jij=dnewi/dvarsj
//values of AI corresponding to fixed components will have been set to zero; do same for J

//load up SEs of vars direct from AI
for(k=0;k<total;k++)
{
if(AI[k+k*total]>=0){varsds[k]=pow(AI[k+k*total],.5);}
else{varsds[k]=-9999;}
}

//get SEs of hers - Jij=delta*scalei/sum-scalei*scalej*vari/sum^2 (where sum across all elements)
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}

for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++)
{
if(fixed[k]<3&&fixed[k2]<3){J[k+k2*total]=-scales[k]*vars[k]*scales[k2]*pow(sum,-2);}
else{J[k+k2*total]=0;}
}
if(fixed[k]<3){J[k+k*total]+=scales[k]/sum;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

//load up (not sure first heritability makes sense, but include anyway)
for(k=0;k<total;k++)
{
if(JAIJT[k+k*total]>=0){hersds[k]=pow(JAIJT[k+k*total],.5);}
else{hersds[k]=-9999;}
}

//save details for her_all as final element of hers / hersds
hers[total]=0;for(k=1;k<total;k++){hers[total]+=hers[k];}
sum=0;
for(k=1;k<total;k++)
{
for(k2=1;k2<total;k2++){sum+=JAIJT[k+k2*total];}
}
if(sum>=0){hersds[total]=pow(sum,.5);}
else{hersds[total]=-9999;}

//save JAIJT for printing
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){cohers[k+k2*total]=JAIJT[k+k2*total];}
}

//get SEs of shares - Jij=delta*scalei/sum-scalei*scalej*vari/sum^2 (where sum excludes noise term)
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}

//first row and column (easier to simply set new0=var0)
if(fixed[0]<3){J[0]=1;}
else{J[0]=0;}
for(k=1;k<total;k++){J[k]=0;J[k*total]=0;}

//rest
for(k=1;k<total;k++)
{
for(k2=1;k2<total;k2++)
{
if(fixed[k]<3&&fixed[k2]<3){J[k+k2*total]=-scales[k]*vars[k]*scales[k2]*pow(sum,-2);}
else{J[k+k2*total]=0;}
}
if(fixed[k]<3){J[k+k*total]+=scales[k]/sum;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

//load up
for(k=0;k<total;k++)
{
if(JAIJT[k+k*total]>=0){sharesds[k]=pow(JAIJT[k+k*total],.5);}
else{sharesds[k]=-9999;}
}

//save JAIJT for printing
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){cohers2[k+k2*total]=JAIJT[k+k2*total];}
}

///////////////////////////

//get variance due to covariates, environments and top predictors - will already have ZTY

sumsq=-ZTY[0]/ns*ZTY[0];
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
value=-ZTY[0]/ns*ZTY[0];for(j=0;j<num_covars+num_envs;j++){value+=ZTY[j]*thetas[j];}
value2=0;for(j=num_covars+num_envs;j<num_fixed;j++){value2+=ZTY[j]*thetas[j];}
covher=value/sumsq;
topher=value2/(sumsq-value);

if(num_covars>1){printf("Proportion of variance explained by the %d covariates: %.4f\n", num_covars, covher);}
if(num_tops==1){printf("Proportion of variance explained by the top predictor: %.4f\n", topher);}
if(num_tops>1){printf("Proportion of variance explained by the %d top predictors: %.4f\n", num_tops, topher);}
if(num_covars>1||num_tops>1){printf("\n");}

//adjust for tops
for(k=0;k<total+1;k++)
{
hers[k]*=(1-topher);
if(hersds[k]!=-9999){hersds[k]*=(1-topher);}
}

//save stuff

flag=0;
for(k=0;k<num_kins;k++){flag+=(kinsums[k]==-9999);}
if(flag==0)
{
sum=0;
for(k=0;k<num_kins;k++){sum+=kinsums[k];}
for(r=0;r<num_regs;r++){sum+=Xsums[r];}
}

sprintf(filename2,"%s.reml", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
if(num_kins+num_regs>0){fprintf(output2, "Blupfile %s.indi.blp\n", outfile);}
else{fprintf(output2, "Blupfile none\n");}
if(num_regs>0||num_tops>0){fprintf(output2, "Regfile %s.reg.blup\n", outfile);}
else{fprintf(output2, "Regfile none\n");}
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
if(cflag==1){fprintf(output2,"Converged YES\n");}
else{fprintf(output2,"Converged NO\n");}
fprintf(output2, "Null_Likelihood %.4f\nAlt_Likelihood %.4f\n", likenull, like);
if(num_kins+num_regs==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[1+k], hersds[1+k], kinsums[k], hers[1+k]/kinsums[k]*1000000, hersds[1+k]/kinsums[k]*1000000);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[1+num_kins+r], hersds[1+num_kins+r], Xsums[r], hers[1+num_kins+r]/Xsums[r]*1000000, hersds[1+num_kins+r]/Xsums[r]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
if(num_kins+num_regs==0){fprintf(output2, "Her_All %.6f NA NA NA NA\n", topher);}
else{fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output2, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[1+k], hersds[1+k]);}
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[1+num_kins+r], hersds[1+num_kins+r]);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f NA NA NA\n", hers[total]+topher, hersds[total]);
}
fclose(output2);

sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Component Effect SE P\n");
fprintf(output3, "Intercept %.4e %.4e %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.4e %.4e %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);} 
for(j=0;j<num_envs;j++){fprintf(output3, "Enviromental_%d %.4e %.4e %.4e\n",j, thetas[num_covars+j], thetasds[num_covars+j], thetapvas[num_covars+j]);}
fclose(output3);

sprintf(filename4,"%s.share", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SE Expected Enrichment SE\n");
if(flag==0)
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, shares[1+k], sharesds[1+k], kinsums[k]/sum, shares[1+k]/kinsums[k]*sum, sharesds[1+k]/kinsums[k]*sum);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f %.6f %.6f %.6f\n", r+1, shares[1+num_kins+r], sharesds[1+num_kins+r], Xsums[r]/sum, shares[1+num_kins+r]/Xsums[r]*sum, sharesds[1+num_kins+r]/Xsums[r]*sum);}
}
else
{
for(k=0;k<num_kins;k++){fprintf(output4, "Share_K%d %.6f %.6f NA NA NA\n", k+1, shares[1+k], sharesds[1+k]);}
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f NA NA NA\n", r+1, shares[1+num_kins+r], sharesds[1+num_kins+r]);}
}
fclose(output4);

sprintf(filename5,"%s.vars", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Component Variance SE\n");
for(k=0;k<num_kins;k++){fprintf(output5, "Var_K%d %.6f %.6f\n", k+1, vars[1+k], varsds[1+k]);}
for(r=0;r<num_regs;r++){fprintf(output5, "Var_R%d %.6f %.6f\n", r+1, vars[1+num_kins+r], varsds[1+num_kins+r]);}
fprintf(output5, "Var_E %.6f %.6f\n", vars[0], varsds[0]);
fclose(output5);

sprintf(filename6,"%s.indi.res", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "ID1\tID2\tPhenotype\tFitted\tResidual\n");
for(i=0;i<ns;i++){fprintf(output6, "%s\t%s\t%f\t%f\t%f\n", ids1[i], ids2[i], Y[i], Yadj[i], Y[i]-Yadj[i]);}
fclose(output6);

if(num_kins+num_regs>0)
{
sprintf(filename7,"%s.indi.blp", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(i=0;i<ns;i++)
{
fprintf(output7, "%s\t%s\t", ids1[i], ids2[i]);
for(k=0;k<num_kins;k++){fprintf(output7, "%.6f\t%.6f\t", mg2[k][i], mg[k][i]);}
for(r=0;r<num_regs;r++){fprintf(output7, "0\t%.6f\t", mg[num_kins+r][i]);}
fprintf(output7, "\n");
}
fclose(output7);
}

if(num_regs>0||num_tops>0)
{
sprintf(filename8,"%s.reg.blup", outfile);
if((output8=fopen(filename8,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename8);exit(1);}
fprintf(output8, "Predictor\tA1\tA2\tCentre\t");
for(r=0;r<num_regs;r++){fprintf(output8, "Region%d\t", r+1);}
if(num_tops>0){fprintf(output8, "Top_Preds\n");}
else{fprintf(output8, "\n");}
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{
fprintf(output8,"%s\t%c\t%c\t%.6f\t", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]);
for(r=0;r<num_regs;r++){fprintf(output8, "%.6f\t", effects[r][j]);}
if(num_tops>0){fprintf(output8, "%.6f\n", effects[num_regs][j]);}
else{fprintf(output8,"\n");}
}
}
fclose(output8);

sprintf(filename9,"%s.reg.score", outfile);
if((output9=fopen(filename9,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename9);exit(1);}
fprintf(output9, "Predictor\tA1\tA2\tCentre\tEffect\n");
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output9,"%s\t%c\t%c\t%.6f\t%.6f\n", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j], effects[num_regs+1][j]);}
}
fclose(output9);
}

sprintf(filename10,"%s.cross", outfile);
if((output10=fopen(filename10,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename10);exit(1);}
for(k=0;k<num_kins;k++){fprintf(output10, "Her_K%d\t", k+1);}
for(r=0;r<num_regs;r++){fprintf(output10, "Her_R%d\t", r+1);}
fprintf(output10, "\n");
for(k=1;k<total;k++)
{
for(k2=1;k2<total;k2++)
{fprintf(output10, "%.6f\t", cohers[k+k2*total]);}
fprintf(output10, "\n");
}
fclose(output10);

////////

if(prev!=-9999)
{
factor=get_factor(Y, ns, prev, -9999, outfile);

sprintf(filename,"%s.liab", filename2);	//.reml
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Num_Kinships %d\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_kins, num_regs, num_tops, num_covars, num_envs);
if(num_kins+num_regs>0){fprintf(output, "Blupfile %s.indi.blp.liab\n", outfile);}
else{fprintf(output, "Blupfile none\n");}
if(num_regs>0||num_tops>0){fprintf(output, "Regfile %s.reg.blup.liab\n", outfile);}
else{fprintf(output, "Regfile none\n");}
fprintf(output, "Coeffsfile %s.coeff.liab\nCovar_Heritability %.4f\n", outfile, covher*factor);
fprintf(output, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);
if(cflag==1){fprintf(output,"Converged YES\n");}
else{fprintf(output,"Converged NO\n");}
fprintf(output, "Null_Likelihood %.4f\nAlt_Likelihood %.4f\n", likenull, like);
if(num_kins+num_regs==1){fprintf(output, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output, "Component Heritability SE Size Mega_Intensity SE\n");
if(flag==0)	//might have null model
{
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d %.6f %.6f %.2f %.6f %.6f\n", k+1, hers[1+k]*factor, hersds[1+k]*factor, kinsums[k], hers[1+k]/kinsums[k]*1000000*factor, hersds[1+k]/kinsums[k]*1000000*factor);}
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[1+num_kins+r]*factor, hersds[1+num_kins+r]*factor, Xsums[r], hers[1+num_kins+r]/Xsums[r]*1000000*factor, hersds[1+num_kins+r]/Xsums[r]*1000000*factor);}
fprintf(output, "Her_Top %.6f NA NA NA NA\n", topher*factor);
if(num_kins+num_regs==0){fprintf(output, "Her_All %.6f NA NA NA NA\n", topher*factor);}
else{fprintf(output, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);}
}
else	//must be non-null
{
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d %.6f %.6f NA NA NA\n", k+1, hers[1+k]*factor, hersds[1+k]*factor);}
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d %.6f %.6f NA NA NA\n", r+1, hers[1+num_kins+r]*factor, hersds[1+num_kins+r]*factor);}
fprintf(output, "Her_Top %.6f NA NA NA NA\n", topher*factor);
fprintf(output, "Her_All %.6f %.6f NA NA NA\n", hers[total]*factor+topher*factor, hersds[total]*factor);
}
fclose(output);

sprintf(filename,"%s.liab", filename3);	//.coeff
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Effect SE P\n");
fprintf(output, "Intercept %.4e %.4e %.4e\n", thetas[0]*pow(factor,.5), thetasds[0]*pow(factor,.5), thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.4e %.4e %.4e\n",j, thetas[j]*pow(factor,.5), thetasds[j]*pow(factor,.5), thetapvas[j]);} 
for(j=0;j<num_envs;j++){fprintf(output, "Enviromental_%d %.4e %.4e %.4e\n",j, thetas[num_covars+j]*pow(factor,.5), thetasds[num_covars+j]*pow(factor,.5), thetapvas[num_covars+j]);}
fclose(output);

sprintf(filename,"%s.liab", filename7);	//.indi.blp
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(i=0;i<ns;i++)
{
fprintf(output, "%s\t%s\t", ids1[i], ids2[i]);
for(k=0;k<num_kins;k++){fprintf(output, "%.6f\t%.6f\t", mg2[k][i]*pow(factor,.5), mg[k][i]*pow(factor,.5));}
for(r=0;r<num_regs;r++){fprintf(output, "0\t%.6f\t", mg[num_kins+r][i]*pow(factor,.5));}
fprintf(output, "\n");
}
fclose(output);

if(num_regs>0||num_tops>0)
{
sprintf(filename,"%s.liab", filename8);	//.reg.blup
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\t");
for(r=0;r<num_regs;r++){fprintf(output, "Region%d\t", r+1);}
if(num_tops>0){fprintf(output, "Top_Preds\n");}
else{fprintf(output, "\n");}
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output,"%s\t%c\t%c\t%.6f\t", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]);}
for(r=0;r<num_regs;r++){fprintf(output, "%.6f\t", effects[r][j]*pow(factor,.5));}
if(num_tops>0){fprintf(output, "%.6f\n", effects[num_regs][j]*pow(factor,.5));}
else{fprintf(output,"\n");}
}
fclose(output);

sprintf(filename,"%s.liab", filename9);	//.reg.score
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\tEffect\n");
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output9,"%s\t%c\t%c\t%.6f\t%.6f\n", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]*pow(factor,.5), effects[num_regs+1][j]*pow(factor,.5));}
}
fclose(output);
}
}	//end of binary

if(shortcut==3)	//print inverse variances
{
sprintf(filename,"%s.weights", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "ID1\tID2\tWeight\n");
for(i=0;i<ns;i++)
{
if(D[i]>0){fprintf(output, "%s\t%s\t%f\t\n", ids1[i], ids2[i], pow(D[i],-1));}
}
fclose(output);
}

if(strcmp(oversfile,"blank")!=0)	//compute and save cats and enrichments - cant have regions
{
sprintf(filename,"%s.cats", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Component Heritability SE\n");
for(k=0;k<num_kins;k++)
{
sum=0;sum2=0;
for(k2=0;k2<num_kins;k2++)
{
sum+=hers[1+k2]*ssums[k][k2];
for(k3=0;k3<num_kins;k3++){sum2+=cohers[1+k2+(1+k3)*total]*ssums[k][k2]*ssums[k][k3];}
}
if(sum2>0){value=pow(sum2,.5);}
else{value=-9999;}
fprintf(output, "Cat_K%d %.6f %.6f\n", k+1, sum, value);
}
fclose(output);

sprintf(filename,"%s.enrich", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Component Share SE Expected Enrichment SE\n");
for(k=0;k<num_kins;k++)
{
sum=0;sum2=0;
for(k2=0;k2<num_kins;k2++)
{
sum+=shares[1+k2]*ssums[k][k2];
for(k3=0;k3<num_kins;k3++){sum2+=cohers2[1+k2+(1+k3)*total]*ssums[k][k2]*ssums[k][k3];}
}
if(sum2>0){value=pow(sum2,.5);}
else{value=-9999;}
fprintf(output, "Enrich_K%d %.6f %.6f %.6f %.6f %.6f\n", k+1, sum, value, ssums[k][num_kins+1], sum/ssums[k][num_kins+1], value/ssums[k][num_kins+1]);
}
fclose(output);
}

printf("Main results saved in %s", filename2);
if(prev!=-9999){printf(", with a liability version saved in %s.liab", filename2);}
printf("\n\n");

////////

if(discenv==1)	//there are 1+num_envs+num_regs genetic kinships - get proportion of variance these explain
//there is no need to compute the entire jacobian, as we only need the top line (gen sums / total sums)
{
printf("Computing heritabilities for %d subgroups\n", num_envs);

Sindexer=malloc(sizeof(int)*ns);
Snums=malloc(sizeof(int)*(1+num_envs));
Sscales=malloc(sizeof(double)*total);
SJ=malloc(sizeof(double)*total);
Shers=malloc(sizeof(double)*(1+num_envs));
Shersds=malloc(sizeof(double)*(1+num_envs));
if(memsave==1){wantids=malloc(sizeof(char*)*ns);}

//first do for all individuals - can use existing scales
sum=0;
for(k=0;k<1+num_envs;k++){sum+=scales[1+k]*vars[1+k];}
for(r=0;r<num_regs;r++){sum+=scales[1+num_kins+r]*vars[1+num_kins+r];}
sum2=0;for(k=0;k<total;k++){sum2+=scales[k]*vars[k];}

//new1 = sum/sum2 (equals sum scale*var for gens / sum scale*var for all)
for(k=0;k<total;k++){SJ[k]=-scales[k]*sum*pow(sum2,-2);}
for(k=0;k<1+num_envs;k++){SJ[1+k]+=scales[1+k]/sum2;}
for(r=0;r<num_regs;r++){SJ[1+num_kins+r]+=scales[1+num_kins+r]/sum2;}

Snums[num_envs]=ns;
Shers[num_envs]=sum/sum2;
value=0;for(k=0;k<total;k++){for(k2=0;k2<total;k2++){value+=AI[k+k2*total]*SJ[k]*SJ[k2];}}
if(value>=0){Shersds[num_envs]=pow(value,.5);}
else{Shersds[num_envs]=-9999;}

//now for each subgroup
for(j=0;j<num_envs;j++)
{
//work out individuals in subgroup
count=0;
for(i=0;i<ns;i++)
{
if(Z[i+(num_covars+num_tops+j)*ns]==1){Sindexer[count]=i;count++;}
}

//get traces across these individuals
Sscales[0]=1;
for(k=0;k<num_kins;k++)
{
if(memsave==0)	//can get directly
{
sum=0;for(i=0;i<count;i++){sum+=Mkins[k][(size_t)Sindexer[i]*ns+Sindexer[i]];}
Sscales[1+k]=sum/count;
}
else	//load individuals into wantids, then read from file
{
for(i=0;i<count;i++){wantids[i]=ids3[Sindexer[i]];}
Sscales[1+k]=read_kin_trace(kinstems[k], count, wantids, NULL, 1, 0, maxthreads);
}
}
for(r=0;r<num_regs;r++)
{
sum=0;
for(i=0;i<count;i++)
{
for(k=Xstarts[r];k<Xends[r];k++){sum+=pow(X[Sindexer[i]+k*ns],2);}
}
Sscales[1+num_kins+r]=sum/Xsums[r]/count;
}

sum=0;
for(k=0;k<1+num_envs;k++){sum+=Sscales[1+k]*vars[1+k];}
for(r=0;r<num_regs;r++){sum+=Sscales[1+num_kins+r]*vars[1+num_kins+r];}
sum2=0;for(k=0;k<total;k++){sum2+=Sscales[k]*vars[k];}

//new1 = sum/sum2 (equals sum Sscale*var for gens / sum Sscale*var for all)
for(k=0;k<total;k++){SJ[k]=-Sscales[k]*sum*pow(sum2,-2);}
for(k=0;k<1+num_envs;k++){SJ[1+k]+=Sscales[1+k]/sum2;}
for(r=0;r<num_regs;r++){SJ[1+num_kins+r]+=Sscales[1+num_kins+r]/sum2;}

Snums[j]=count;
Shers[j]=sum/sum2;
value=0;for(k=0;k<total;k++){for(k2=0;k2<total;k2++){value+=AI[k+k2*total]*SJ[k]*SJ[k2];}}
if(value>=0){Shersds[j]=pow(value,.5);}
else{Shersds[j]=-9999;}
}	//end of j loop

sprintf(filename,"%s.subgroups", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Num_Samples Heritability SE\n");
for(j=0;j<num_envs;j++){fprintf(output, "Her_Sub%d %d %.6f %.6f\n", j+1, Snums[j], Shers[j], Shersds[j]);}
fprintf(output, "Her_All %d %.6f %.6f\n", Snums[num_envs], Shers[num_envs], Shersds[num_envs]);
fclose(output);

printf("Estimates saved in %s\n\n", filename);
free(Sindexer);free(Snums);free(Sscales);free(Shers);free(Shersds);free(SJ);
if(memsave==1){free(wantids);}
}	//end of discenv==1

////////

free(fixed);free(fixedsave);free(vars);free(vardiffs);free(varsds);free(scales);free(hers);free(hersds);free(shares);free(sharesds);
free(AI);free(AI2);free(AI3);free(BI);free(J);free(JAI);free(JAIJT);free(cohers);free(cohers2);
free(Z2);free(ZTY);free(ZTZ);free(ZTZ2);free(ZTZZTY);
free(ZTVZ);free(ZTVZ2);free(ZTVZ3);free(PY);
for(k=0;k<total;k++){free(KPY[k]);free(PKPY[k]);}free(KPY);free(PKPY);free(traces);
free(ZTVY);free(thetas);free(thetasds);free(thetapvas);free(Yadj);
if(num_kins+num_regs>0)
{
for(k=0;k<num_kins+num_regs;k++){free(mg[k]);}free(mg);
for(k=0;k<num_kins+num_regs;k++){free(mg2[k]);}free(mg2);
}
if(num_regs>0||num_tops>0)
{
for(r=0;r<num_regs+4;r++){free(effects[r]);}free(effects);
}

if(shortcut==0){free(V);free(VZ);free(VZZTVZ);free(P);}
if(shortcut==1){free(UTY);free(UTZ);free(D);free(BUTZ);free(H);free(HUTY);free(HKPY);}
if(shortcut==2)
{
free(R_single);for(k=0;k<num_kins;k++){free(KR_single[k]);}free(KR_single);
free(V_full);free(RHS_single);free(PY_single);free(KPY_single);free(RHS2_single);
free(VZ);free(VY);free(VR);free(VZZTVZ);free(ZTVR);free(PR);free(VKPY);free(ZTVKPY);
}
if(shortcut==3)
{
free(kin_diags);free(D);free(BUTZ);free(H);free(HUTY);free(HKPY);
}

if(Xtotal>0)
{
free(PX);free(XTPY);
if(shortcut==1)
{
free(UTX);free(DUTX);free(XTVX);free(XTVX2);free(XTVX3);free(F);free(FUTZ);free(FUTY);free(FUTX);free(HUTX);free(FKPY);}
}

}	//end of multi_reml

///////////////////////////


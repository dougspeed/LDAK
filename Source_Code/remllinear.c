/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Mixed-model linear regression
//will have null model (stats==NULL) or normal (stats!=NULL)

///////////////////////////

void linear_reml(int ns, int num_fixed, double *Y, double *Z, double *U, double *E, double *UTY, double *UTZ, double *kintraces, double *stats, double *vstarts, double *ts, double *tsds, double *tpvas, int constrain, double tol, int maxiter, int exact)
{
int i, i2, j, j2, k, k2, r, g, count, token, one=1, info;
double sum, sumsq, max, alpha, beta;

int total, total2, rflag, *fixed, *fixedsave, nlost, nlost2;
double nfree, varnull, varmin, gam, gam2, gam3, relax;
double like, likeold, diff;
double *scales, *vars, *vardiffs, *varsds, *hers, *hersds;

double *AI, *AI2, *AI3, *BI, *J, *JAI, *JAIJT;
double *ZTY, *ZTZ, *ZTZ2, *ZTZZTY;
double detV, *ZTVZ, *ZTVZ2, *ZTVZ3, detZTVZ, *PY, **KPY, **PKPY, *traces;
double *ZTVY, *thetas, *thetasds, *thetapvas;

double *D, detD, *BUTZ, *H, *HUTY, *HKPY;

//variables to make compatible with remllikes.c and remlderiv.c
int num_kins=1, num_regs=0, Xtotal=0, *Xstarts=NULL, *Xends=NULL, shortcut=1, num_vects=-9999, ldlt=-9999, memsave=-9999, maxthreads=-9999;
float **Mkins_single=NULL;
double **Mkins=NULL, *sweights, *X=NULL, *Xsums=NULL;
char **kinstems=NULL, **ids3=NULL;
sweights=malloc(sizeof(double)*ns);
for(i=0;i<ns;i++){sweights[i]=1;}

size_t scount=0, stotal=0;
int *ipiv, lwork;
float wkopt_single, *work_single=NULL, alpha_single, beta_single;
float *R_single=NULL, **KR_single=NULL, *V_full=NULL, *V_single=NULL, *RHS_single=NULL;
float *PY_single=NULL, *KPY_single=NULL, *RHS2_single=NULL;
double *V=NULL, *VZ=NULL, *VY=NULL, *VR=NULL, *VZZTVZ=NULL, *P=NULL, *ZTVR=NULL, *PR=NULL, *VKPY=NULL, *ZTVKPY=NULL;
double *kin_diags=NULL;
double *PX=NULL, *XTPY=NULL;
double *UTX=NULL, *DUTX=NULL, detC, *XTVX=NULL, *XTVX2=NULL, *XTVX3=NULL, detXTVX;
double *F=NULL, *FUTZ=NULL, *FUTY=NULL, *FUTX=NULL, *HUTX=NULL, *FKPY=NULL;


//set nfree and total
nfree=ns-num_fixed;
total=1+num_kins+num_regs;

//allocate variables

fixed=malloc(sizeof(int)*total);
fixedsave=malloc(sizeof(int)*total);
vars=malloc(sizeof(double)*total);
vardiffs=malloc(sizeof(double)*total);
varsds=malloc(sizeof(double)*total);
scales=malloc(sizeof(double)*total);
hers=malloc(sizeof(double)*total);
hersds=malloc(sizeof(double)*total);

AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);
J=malloc(sizeof(double)*total*total);
JAI=malloc(sizeof(double)*total*total);
JAIJT=malloc(sizeof(double)*total*total);

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

D=malloc(sizeof(double)*ns);
BUTZ=malloc(sizeof(double)*ns*num_fixed);
H=malloc(sizeof(double)*num_fixed*ns);
HUTY=malloc(sizeof(double)*num_fixed);
HKPY=malloc(sizeof(double)*num_fixed);

////////

//solve model with just fixed effects to get varnull and varmin
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
for(j=0;j<num_fixed;j++){ZTZZTY[j]=ZTY[j];}
(void)eigen_invert(ZTZ, num_fixed, ZTZ2, 1, ZTZZTY, 1);

sumsq=0;
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
for(j=0;j<num_fixed;j++){sumsq-=ZTY[j]*ZTZZTY[j];}
varnull=sumsq/nfree;
varmin=0.0001*varnull;

//set scales (remember total=2)
scales[0]=1;scales[1]=kintraces[0];

//set starting vars, hers and fixed (remember total=2)
if(stats==NULL)	//set heritabilties agnostically
{hers[0]=0.5;hers[1]=.5;}
else	//set based on vstarts
{hers[0]=vstarts[0];hers[1]=vstarts[1];}

for(k=0;k<total;k++){vars[k]=hers[k]*varnull/scales[k];}
for(k=0;k<total;k++){fixed[k]=0;}

if(stats==NULL)	//will screen print
{printf("Iter\tHer_K1\tHer_All\tLikelihood\tDifference\tTarget\tNum_Constrained\n");}

////////

//now iterate
count=0;
rflag=0;	//0 if normal moves, 1 if reduced moves, 2 if transitioning from reduced to normal moves
while(1)
{
//compute invV, detV, invZTVZ, detZTVZ, P, PY and gam
#include "remllike.c"

if(exact==2){break;}	//only required BUTZ and ZTVZ

//get likelihood
like=-.5*gam-.5*detZTVZ-.5*detV-.5*nfree*log(2*M_PI);

if(count>0)	//set diff and decide what type of move to do
{
if(like>likeold-tol)	//move was fine, so next move will be normal or transitioning
{
diff=like-likeold;
if(rflag==1){rflag=2;}
else{rflag=0;}
}
else	//move was poor, so return to previous state and next move will be reduced
{
if(stats==NULL){printf("Warning, the last move reduced the likelihood, so have returned to the previous state\n");}

for(k=0;k<total;k++){fixed[k]=fixedsave[k];}
for(k=0;k<total;k++){vars[k]-=vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}
like=likeold;
diff=0;
rflag=1;
}
}
likeold=like;
for(k=0;k<total;k++){fixedsave[k]=fixed[k];}

//compute PX, XTPY, KPY, PKPY, (inverse) AI and BI
#include "remlderiv.c"

nlost=0;for(k=0;k<total;k++){nlost+=(fixed[k]>=3);}
nlost2=0;for(k=0;k<total;k++){nlost2+=(fixed[k]==1||fixed[k]==2);}

if(stats==NULL)	//print update
{
if(count==0){printf("Start\t");}
else{printf("%d\t", count);}
printf("%.4f\t%.4f\t", hers[1], hers[1]);
if(rflag==0||rflag==2){printf("%.2f\t", like);}
else{printf("%.2f*\t", like);}
if(count==0){printf("n/a\t\t%.6f\t", tol);}
else{printf("%.6f\t%.6f\t", diff, tol);}
if(nlost2==0){printf("%d\n", nlost);}
else{printf("%d*\n", nlost);}
}

//see if breaking (normally can only break if rflag=0 and nlost2=0, unless at iter limit)

if(nlost>=num_kins){break;}	//all heritabilities constrained
if(count>0)
{
if(fabs(diff)<tol&&rflag==0&&nlost2==0){break;}
}
if(count==maxiter)
{printf("Warning, REML failed to converge; if this happens many times, consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)\n", maxiter, tol);break;}

////////

//update variances using NR

//decide how far to move
if(rflag==0||rflag==2){relax=1;}
else{relax*=.5;}

//get proposed moves, ensuring not too large
alpha=relax;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, vardiffs, &one);

max=0;
for(k=0;k<total;k++)
{
if(fabs(vardiffs[k])>max){max=fabs(vardiffs[k]);}
}
if(max>.1*varnull)	//then reduce moves
{
relax*=.1*varnull/max;
for(k=0;k<total;k++){vardiffs[k]*=.1*varnull/max;}
}

if(constrain==1)	//variances can not be negative
{
for(k=0;k<total;k++)
{
if(fixed[k]<3)	//free to update
{
if(vars[k]+vardiffs[k]<varmin){vardiffs[k]=varmin-vars[k];fixed[k]++;}
else{fixed[k]=0;}
if(fixed[k]==3){vardiffs[k]=-vars[k];}
}
}}

//now move
for(k=0;k<total;k++){vars[k]+=vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}

count++;
}	//end of while loop

////////

//get fixed coefficients (ZTinvVZ)^-1 ZTinvVY with variance matrix (ZTinvVZ)^-1
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, BUTZ, &ns, UTY, &one, &beta, ZTVY, &one);
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

if(stats==NULL)	//save thetas and heritabilities / variances
{
for(j=0;j<num_fixed;j++){ts[j]=thetas[j];tsds[j]=thetasds[j];tpvas[j]=thetapvas[j];}
vstarts[0]=hers[0];vstarts[1]=hers[1];vstarts[2]=vars[0];vstarts[3]=vars[1];
}
else
{
//load up stats
stats[0]=thetas[num_fixed-1];
stats[1]=thetasds[num_fixed-1];
if(stats[1]!=-9999){stats[2]=stats[0]/stats[1];}
stats[3]=thetapvas[num_fixed-1];
}

free(sweights);
free(fixed);free(fixedsave);free(vars);free(vardiffs);free(varsds);free(scales);free(hers);free(hersds);
free(AI);free(AI2);free(AI3);free(BI);free(J);free(JAI);free(JAIJT);
free(ZTY);free(ZTZ);free(ZTZ2);free(ZTZZTY);
free(ZTVZ);free(ZTVZ2);free(ZTVZ3);free(PY);
for(k=0;k<total;k++){free(KPY[k]);free(PKPY[k]);}free(KPY);free(PKPY);free(traces);
free(ZTVY);free(thetas);free(thetasds);free(thetapvas);
free(D);free(BUTZ);free(H);free(HUTY);free(HKPY);
}	//end of adv_reml

///////////////////////////


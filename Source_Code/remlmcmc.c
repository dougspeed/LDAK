/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//MCMC-based REML - variance matrix constructed from noise + kinships

///////////////////////////

void mcmc_reml(int ns, int num_covars, int num_tops, int num_kins, int num_vects, double *Y, double *Z, double **Mkins, double *kintraces, double *kinsums, double *R, int np, int hestart, int compact, int constrain, double tol, int maxiter, int memsave, char **kinstems, char **ids3, char *outfile, char **ids1, char **ids2, int num_preds, char **allpreds, char *allal1, char *allal2, int *tkeeppreds, double *tcentres)
{
size_t scount, stotal;
int i, j, j2, k, k2, r, g, count, token, flag, one=1, info;
double value, value2, sum, sum2, sumsq, max, alpha, beta;

int num_fixed, total, total2, cflag, rflag, *fixed, *fixedsave, nlost, nlost2;
double nfree, varnull, varmin, gam, gam2, gam3, relax;
double likenull, like, likeold, diff, lrtstat, lrtpva, covher, topher, factor;
double *scales, *vars, *vardiffs, *varsds, *hers, *hersds, *shares, *sharesds;

double *AI, *AI2, *AI3, *AIsave, *BI, *BIsave, *J, *JAI, *JAIJT, *cohers;
double *ZTY, *ZTZ, *ZTZ2, *ZTZZTY, detZTZ;
double *vect, *vect2;
double *V, *ZTVZ, *ZTVZ2, *ZTVZ3, detZTVZ, *PY, **KPY, **PKPY, *traces;
double *ZTVY, *thetas, *thetasds, *thetapvas, *Yadj;
double **mg, **mg2, **effects;

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500], filename10[500];
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7, *output8, *output9, *output10;


//set num_fixed, nfree, total and stotal
num_fixed=num_covars+num_tops;
nfree=ns-num_fixed;
total=1+num_kins;
total2=num_vects+num_fixed+1;
stotal=(size_t)ns*ns;

//allocate variables

fixed=malloc(sizeof(int)*total);
fixedsave=malloc(sizeof(int)*total);
vars=malloc(sizeof(double)*total);
vardiffs=malloc(sizeof(double)*total);
varsds=malloc(sizeof(double)*total);
scales=malloc(sizeof(double)*total);
hers=malloc(sizeof(double)*total);
hersds=malloc(sizeof(double)*total);
shares=malloc(sizeof(double)*total);
sharesds=malloc(sizeof(double)*total);

AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
AIsave=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);
BIsave=malloc(sizeof(double)*total);
J=malloc(sizeof(double)*total*total);
JAI=malloc(sizeof(double)*total*total);
JAIJT=malloc(sizeof(double)*total*total);
cohers=malloc(sizeof(double)*total*total);

ZTY=malloc(sizeof(double)*num_fixed);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZZTY=malloc(sizeof(double)*num_fixed);

vect=malloc(sizeof(double)*ns*(num_vects+1+num_fixed));
vect2=malloc(sizeof(double)*ns*(num_vects+1+num_fixed));

V=malloc(sizeof(double)*ns*ns);
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

mg=malloc(sizeof(double*)*num_kins);
for(k=0;k<num_kins;k++){mg[k]=malloc(sizeof(double)*ns);}
mg2=malloc(sizeof(double*)*num_kins);
for(k=0;k<num_kins;k++){mg2[k]=malloc(sizeof(double)*ns);}

if(num_tops>0)
{
effects=malloc(sizeof(double*)*4);
for(r=0;r<4;r++){effects[r]=malloc(sizeof(double)*num_preds);}
}

////////

//solve model with just fixed effects to get varnull, varmin and null likelihood
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);
for(j=0;j<num_fixed;j++){ZTZZTY[j]=ZTY[j];}
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, 1, ZTZZTY, 1);

sumsq=0;
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
for(j=0;j<num_fixed;j++){sumsq-=ZTY[j]*ZTZZTY[j];}
varnull=sumsq/nfree;
varmin=0.0001*varnull;
likenull=-.5*(nfree+nfree*log(2*M_PI*varnull)+detZTZ);

//set scales
scales[0]=1;
for(k=0;k<num_kins;k++){scales[1+k]=kintraces[k];}

//set starting vars, hers, shares and fixed - total always >1
if(hestart==1)	//set based on he regression
{
printf("Computing starting heritabilities\n");
he_starts(hers, ns, num_covars, 0, num_tops, num_kins, 0, Y, Z, Mkins, NULL, kintraces, kinsums, NULL, -9999, NULL, NULL, NULL, tol, memsave, kinstems, ids3, missingvalue);
}
else	//set agnostically
{
hers[0]=.5;shares[0]=1;
for(k=1;k<total;k++){hers[k]=.5/(total-1);shares[k]=1.0/(total-1);}
}
for(k=0;k<total;k++){vars[k]=hers[k]*varnull/scales[k];}
for(k=0;k<total;k++){fixed[k]=0;}

//prepare to screen and file print
printf("Iter\t");
for(k=0;k<num_kins;k++){printf("Her_K%d\t", k+1);}
printf("Her_All\tDifference\tTarget\tNum_Constrained\n");

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Iteration\t");
for(k=0;k<num_kins;k++){fprintf(output, "Her_K%d\t", k+1);}
fprintf(output, "Her_All\tTotal_Variance\tDifference\tTarget\tNum_Constrained\n");
fclose(output);

////////

//now iterate
count=0;
cflag=1;
rflag=0;	//0 if normal moves, 1 if reduced moves, 2 if transitioning from reduced to normal moves
while(1)
{
//load up V - will change this to just upper part and single precision ???
for(scount=0;scount<stotal;scount++){V[scount]=0;}
for(i=0;i<ns;i++){V[(size_t)i*ns+i]=vars[0];}

for(k=0;k<num_kins;k++)	//add on varK for kinships
{
if(vars[1+k]!=0)
{
if(memsave==0)	//have kinships saved in Mkins
{
for(scount=0;scount<stotal;scount++){V[scount]+=vars[1+k]*Mkins[k][scount];}
}
else{read_kins(kinstems[k], V, NULL, vars[1+k], ns, ids3, 3);}
}
}

//get 

//construct vectors with E(zzt)=V
for(g=0;g<num_vects;g++)
{
//start with noise
for(i=0;i<ns;i++){vect[i+g*ns]=vars[0]*R[i+g*ns];}

for(k=0;k<num_kins;k++)	//add on for kinships
{
if(vars[1+k]!=0)
{
for(i=0;i<ns;i++){vect[i+g*ns]=vars[1+k]*R[i+(g+(1+k)*num_vects)*ns];}
}}
}

//put Z and Y into end of vects
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){vect[i+(num_vects+j)*ns]=Z[i+j*ns];}
}
for(i=0;i<ns;i++){vect[i+(num_vects+num_fixed)*ns]=Y[i];}

double *vect3=malloc(sizeof(double)*ns*total2);
double *vect4=malloc(sizeof(double)*ns*total2);

for(i=0;i<5;i++){printf("%f %f %f | %f %f %f\n", vect[i], vect[i+num_vects*ns], vect[i+(num_vects+1)*ns], vect2[i], vect2[i+num_vects*ns], vect[i+(num_vects+1)*ns]);}

//compute invV vect
for(g=0;g<num_vects+1+num_fixed;g++)
{
for(i=0;i<ns;i++){vect2[i+num_vects*ns]=0;}
}
cg_solve(V, ns, vect2, vect, num_vects+1+num_fixed, tol);

//construct PY = invVY - invVZ (ZTinvVZ) ZTinvVY, so require invVY and invVZ
for(i=0;i<ns;i++){PY[i]=vect2[i+num_vects*ns];}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &total2, &ns, &alpha, V, &ns, vect2, &ns, &beta, vect3, &ns);

for(i=0;i<5;i++){printf("%f %f %f | %f %f %f\n", vect2[i], vect2[i+num_vects*ns], vect2[i+(num_vects+1)*ns], vect3[i], vect3[i+num_vects*ns], vect3[i+(num_vects+1)*ns]);}

//repeat

//construct vectors with E(zzt)=V
for(g=0;g<num_vects;g++)
{
//start with noise
for(i=0;i<ns;i++){vect[i+g*ns]=vars[0]*R[i+g*ns];}

for(k=0;k<num_kins;k++)	//add on for kinships
{
if(vars[1+k]!=0)
{
for(i=0;i<ns;i++){vect[i+g*ns]=vars[1+k]*R[i+(g+(1+k)*num_vects)*ns];}
}}
}

//put Z and Y into end of vects
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){vect[i+(num_vects+j)*ns]=Z[i+j*ns];}
}
for(i=0;i<ns;i++){vect[i+(num_vects+num_fixed)*ns]=Y[i];}

ldlt_invert(V, ns, num_vects+num_fixed+1, vect, &info, 1);

alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &total2, &ns, &alpha, V, &ns, vect, &ns, &beta, vect3, &ns);

for(i=0;i<5;i++){printf("%f %f %f | %f %f %f\n", vect[i], vect[i+num_vects*ns], vect[i+(num_vects+1)*ns], vect3[i], vect3[i+num_vects*ns], vect3[i+(num_vects+1)*ns]);}
exit(1);

//while for traces require invVR

/*
, compute invV, PY and gam (plus other variables depending on shortcut)
#include "remllike.c"

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
printf("Warning, the last move reduced the likelihood, so have returned to the previous state\n");
for(k=0;k<total;k++){fixed[k]=fixedsave[k];}
for(k=0;k<total;k++){vars[k]-=vardiffs[k];}
sum=0;for(k=0;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){hers[k]=scales[k]*vars[k]/sum;}
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}
like=likeold;
diff=0;
rflag=1;
}
}
likeold=like;
for(k=0;k<total;k++){fixedsave[k]=fixed[k];}

if(rflag==0||rflag==2)	//compute derivatives
{
//compute PX, XTPY, KPY, PKPY, (inverse) AI and BI (plus other variables, depending on shortcut)
#include "remlderiv.c"

//save derivatives in case returning
for(k=0;k<total;k++)
{
BIsave[k]=BI[k];
for(k2=0;k2<total;k2++){AIsave[k+k2*total]=AI[k+k2*total];}
}
}
else	//recover saved values
{
for(k=0;k<total;k++)
{
BI[k]=BIsave[k];
for(k2=0;k2<total;k2++){AI[k+k2*total]=AIsave[k+k2*total];}
}
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
if(rflag==0||rflag==2){printf("%.2f\t", like);}
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

//see if breaking (normally can only break if rflag=0 and nlost2=0, unless at iter limit)

if(num_kins+num_regs==0){break;}	//null model
if(nlost>=num_kins+num_regs){printf("All heritabilities are constrained\n");break;}
if(count>0)
{
if(fabs(diff)<tol&&rflag==0&&nlost2==0){break;}
}
if(count==maxiter)
{
printf("\nWarning, REML failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)\n", maxiter, tol);
if(constrain==0)
{printf("; additionally, it might help to use \"--constrain YES\" to constrain heritabilities to [0,1]\n");}
else
{printf("; additionally, it might help to use \"--constrain NO\" to allow heritabilities outside [0,1]\n");}
cflag=0;
break;
}

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
sum=0;for(k=1;k<total;k++){sum+=scales[k]*vars[k];}
for(k=0;k<total;k++){shares[k]=scales[k]*vars[k]/sum;}

*/
count++;
}	//end of while loop
printf("\n");

/*

//get some stats
lrtstat=2*(like-likenull);
lrtpva=.5*erfc(pow(lrtstat,.5)*M_SQRT1_2);
if(lrtstat<0){lrtpva=.75;}
if(hers[0]<=0){lrtpva=1;}

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

//set Yadj to fixed effects
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Yadj, &one);

for(k=0;k<num_kins+num_regs;k++)	//random effects are Kg invV Yadj = g KPY - get also g PY (redundant for regions)
{
if(shortcut==0)
{
for(i=0;i<ns;i++){mg[k][i]=vars[1+k]*KPY[1+k][i];mg2[k][i]=vars[1+k]*PY[i];}
}
else
{
alpha=vars[1+k];beta=0.0;
dgemv_("N", &ns, &ns, &alpha, U, &ns, KPY[1+k], &one, &beta, mg[k], &one);
dgemv_("N", &ns, &ns, &alpha, U, &ns, PY, &one, &beta, mg2[k], &one);
}

//add random effects to Yadj
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
effects[num_regs][tkeeppreds[j]]=thetas[num_covars+j];
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

//load up
for(k=1;k<total;k++)
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

*/
}	//end of mcmc_reml

///////////////////////////


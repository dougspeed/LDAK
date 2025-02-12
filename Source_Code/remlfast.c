/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Solves special case of REML - when variance matrix constructed from just regions (at least one)

///////////////////////////

void fast_reml(int ns, int num_covars, int num_envs, int num_tops, int num_regs, double* Y, double* Z, double* X, int Xtotal, int *Xstarts, int *Xends, double *Xsums, double *Xnss, double *Xrhos, double prev, int np, int constrain, double shrink, double strip, double tol, int maxiter, char *outfile, char **ids1, char **ids2, int num_preds, char **allpreds, char *allal1, char *allal2, int *rkeeppreds, int **regindex, double *rcentres, double *rmults, double *rweights, int *tkeeppreds, double *tcentres)
{
int i, j, j2, r, r2, count, token, one=1, info;
double value, value2, sum, sumsq, max, alpha, beta;

int num_fixed, total, cflag, rflag, *fixed, *fixedsave, nlost, nlost2;
double nfree, scale, varnull, varmin, gam, gam2, gam3, relax;
double likenull, like, likeold, diff, lrtstat, lrtpva, covher, topher, factor;
double *vars, *vardiffs, *varsds, *hers, *hersds, *shares, *sharesds;

double *AI, *AI2, *AI3, *BI, *J, *JAI, *JAIJT, *cohers;
double *ZTY, *ZTX, *ZTZ, *ZTZ2, *ZTZ3, detZTZ, *ZTZZTY, *ZTZZTX;
double YTCY, *XTCY, *XTCX,  *T, *T2, *T3, detT, *TD, **ITD, **TITD, *traces;
double *B, *B2, *B3, *ZTXB, *XTY, *ZTVY, *ZTVZ, *ZTVZ2, *ZTVZ3, *thetas, *thetasds, *thetapvas, *Yadj;
double **mg, **mg2, **effects;

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], filename6[500], filename7[500], filename8[500], filename9[500], filename10[500];
FILE *output, *output2, *output3, *output4, *output5, *output6, *output7, *output8, *output9, *output10;



//set num_fixed, nfree, scale and total
num_fixed=num_covars+num_envs+num_tops;
if(Xnss==NULL)
{
nfree=ns-num_fixed;
scale=1;
}
else
{
sum=0;for(j=0;j<Xtotal;j++){sum+=Xnss[j];}
nfree=sum/Xtotal;
scale=nfree/ns;
}
total=1+num_regs;

//allocate variables - will always have Xtotal>0

fixed=malloc(sizeof(int)*total);
fixedsave=malloc(sizeof(int)*total);
vars=malloc(sizeof(double)*total);
vardiffs=malloc(sizeof(double)*total);
varsds=malloc(sizeof(double)*total);
hers=malloc(sizeof(double)*total);
hersds=malloc(sizeof(double)*total);
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

ZTY=malloc(sizeof(double)*num_fixed);
ZTX=malloc(sizeof(double)*num_fixed*Xtotal);
ZTZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZ2=malloc(sizeof(double)*num_fixed);
ZTZ3=malloc(sizeof(double)*num_fixed*num_fixed);
ZTZZTY=malloc(sizeof(double)*num_fixed);
ZTZZTX=malloc(sizeof(double)*num_fixed*Xtotal);

XTCY=malloc(sizeof(double)*Xtotal);
XTCX=malloc(sizeof(double)*Xtotal*Xtotal);
T=malloc(sizeof(double)*Xtotal*Xtotal);
T2=malloc(sizeof(double)*Xtotal);
T3=malloc(sizeof(double)*Xtotal*Xtotal);
TD=malloc(sizeof(double)*Xtotal);

ITD=malloc(sizeof(double*)*total);
for(r=0;r<total;r++){ITD[r]=malloc(sizeof(double)*Xtotal);}
TITD=malloc(sizeof(double*)*total);
for(r=0;r<total;r++){TITD[r]=malloc(sizeof(double)*Xtotal);}
traces=malloc(sizeof(double)*total);

B=malloc(sizeof(double)*Xtotal*Xtotal);
B2=malloc(sizeof(double)*Xtotal);
B3=malloc(sizeof(double)*Xtotal*Xtotal);
ZTXB=malloc(sizeof(double)*num_fixed*Xtotal);
XTY=malloc(sizeof(double)*Xtotal);
ZTVY=malloc(sizeof(double)*num_fixed);
ZTVZ=malloc(sizeof(double)*num_fixed*num_fixed);
ZTVZ2=malloc(sizeof(double)*num_fixed);
ZTVZ3=malloc(sizeof(double)*num_fixed*num_fixed);
thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*ns);

mg=malloc(sizeof(double*)*num_regs);
for(r=0;r<num_regs;r++){mg[r]=malloc(sizeof(double)*ns);}
mg2=malloc(sizeof(double*)*num_regs);
for(r=0;r<num_regs;r++){mg2[r]=malloc(sizeof(double)*ns);}

if(num_regs>0||num_tops>0)
{
effects=malloc(sizeof(double*)*(num_regs+4));
for(r=0;r<num_regs+4;r++){effects[r]=malloc(sizeof(double)*num_preds);}
}

//fill some variables

//get ZTY, ZTX, ZTZ
alpha=scale;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTY, &one);
dgemm_("T", "N", &num_fixed, &Xtotal, &ns, &alpha, Z, &ns, X, &ns, &beta, ZTX, &num_fixed);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTZ, &num_fixed);

//invert ZTZ and get ZTZZTY and ZTZZTX
detZTZ=eigen_invert(ZTZ, num_fixed, ZTZ2, -1, ZTZ3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTZ, &num_fixed, ZTY, &one, &beta, ZTZZTY, &one);
dgemm_("N", "N", &num_fixed, &Xtotal, &num_fixed, &alpha, ZTZ, &num_fixed, ZTX, &num_fixed, &beta, ZTZZTX, &num_fixed);

//get YTCY (when using summary statistics, this equals var(Y)*nfree), varnull, varmin and likenull
YTCY=0;
for(i=0;i<ns;i++){YTCY+=scale*pow(Y[i],2);}
for(j=0;j<num_fixed;j++){YTCY-=ZTY[j]*ZTZZTY[j];}
varnull=YTCY/nfree;
varmin=0.0001*varnull;
likenull=-.5*(nfree+nfree*log(2*M_PI*varnull)+detZTZ);

//get XTCX (when using summary statistics, this equals covar(X)*nfree)
alpha=scale;beta=0.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &ns, &alpha, X, &ns, X, &ns, &beta, XTCX, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &num_fixed, &alpha, ZTX, &num_fixed, ZTZZTX, &num_fixed, &beta, XTCX, &Xtotal);

//get XTCY
if(Xnss==NULL)
{
alpha=scale;beta=0.0;
dgemv_("T", &ns, &Xtotal, &alpha, X, &ns, Y, &one, &beta, XTCY, &one);
alpha=-1.0;beta=1.0;
dgemv_("T", &num_fixed, &Xtotal, &alpha, ZTX, &num_fixed, ZTZZTY, &one, &beta, XTCY, &one);
}
else
{
for(j=0;j<Xtotal;j++){XTCY[j]=Xrhos[j]*pow(XTCX[j+j*Xtotal]*YTCY,.5);}
}

if(shrink!=1)	//alter XTCX (must be using summary statistics, so XTCX=XTX)
{
if(shrink<1)	//deflate off diagonal terms
{
for(j=0;j<Xtotal;j++)
{
for(j2=j+1;j2<Xtotal;j2++){XTCX[j+j2*Xtotal]*=shrink;XTCX[j2+j*Xtotal]*=shrink;}
}
}
if(shrink>1)	//inflate diagonal terms
{
for(j=0;j<Xtotal;j++){XTCX[j+j*Xtotal]*=shrink;}
}
}

if(strip>0)	//decompose XTCX, then set evalues to zero and recompose
//for moment will strip whole matrix, but might be better to strip diagonal blocks
{eigen_strip(XTCX, 0, Xtotal, Xtotal, strip);}

////////

//set starting vars, shares and fixed
hers[0]=.5;shares[0]=1;
for(r=1;r<total;r++){hers[r]=.5/num_regs;shares[r]=1.0/num_regs;}
for(r=0;r<total;r++){vars[r]=hers[r]*varnull;}
for(r=0;r<total;r++){fixed[r]=0;}

//prepare to screen and file print
printf("Iter\t");
for(r=0;r<num_regs;r++){printf("Her_R%d\t", r+1);}
printf("Her_All\tLikelihood\tDifference\tTarget\tNum_Constrained\n");

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Iteration\t");
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d\t", r+1);}
fprintf(output, "Her_All\tTotal_Variance\tLikelihood\tDifference\tTarget\tNum_Constrained\n");
fclose(output);

////////

//now iterate
count=0;
cflag=1;
rflag=0;	//0 if normal moves, 1 if reduced moves, 2 if transitioning from reduced to normal moves
while(1)
{
//get likelihood
//get T = XTCX + We/g - will blank for zero variance terms
for(r=0;r<num_regs;r++)
{
for(j=Xstarts[r];j<Xends[r];j++)	//fill row j
{
if(fixed[1+r]<3)	//copy in XTCXj or blank, adding eWr/vr to diagonal
{
for(r2=0;r2<num_regs;r2++)
{
for(j2=Xstarts[r2];j2<Xends[r2];j2++)
{
if(fixed[1+r2]<3){T[j+j2*Xtotal]=XTCX[j+j2*Xtotal];}
else{T[j+j2*Xtotal]=0;}
}}
T[j+j*Xtotal]+=Xsums[r]*vars[0]/vars[1+r];
}
else	//blank row except for diagonal
{
for(j2=0;j2<Xtotal;j2++){T[j+j2*Xtotal]=0;}
T[j+j*Xtotal]=1;
}
}}

//invert T - save T in case ldlt fails
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){T3[j+j2*Xtotal]=T[j+j2*Xtotal];}
}
detT=ldlt_invert(T, Xtotal, -1, NULL, &info, 1);
if(info!=0)	//eigen will definitely work, so can use T3 as work space
{
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){T[j+j2*Xtotal]=T3[j+j2*Xtotal];}
}
detT=eigen_invert(T, Xtotal, T2, -1, T3, 1);
}

//remove ones from diagonal for blanked components
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]>=3)
{
for(j=Xstarts[r];j<Xends[r];j++){T[j+j*Xtotal]=0;}
}
}

//get invTD and gam = YTCY - DTinvTD, where D=XTCY (does not matter that XTCY not blanked)
alpha=1.0;beta=0.0;
dgemv_("N", &Xtotal, &Xtotal, &alpha, T, &Xtotal, XTCY, &one, &beta, TD, &one);
gam=YTCY;for(j=0;j<Xtotal;j++){gam-=TD[j]*XTCY[j];}

//get likelihood
value=0;
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]<3){value+=(Xends[r]-Xstarts[r])*log(Xsums[r]*fabs(vars[0]/vars[1+r]));}
}
like=-.5*(gam/vars[0]+nfree*log(2*M_PI*vars[0])+detT-value+detZTZ);

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
for(r=0;r<total;r++){fixed[r]=fixedsave[r];}
for(r=0;r<total;r++){vars[r]-=vardiffs[r];}
sum=0;for(r=0;r<total;r++){sum+=vars[r];}
for(r=0;r<total;r++){hers[r]=vars[r]/sum;}
sum=0;for(r=1;r<total;r++){sum+=vars[r];}
for(r=0;r<total;r++){shares[r]=vars[r]/sum;}
like=likeold;
diff=0;
rflag=1;
}
}
likeold=like;
for(r=0;r<total;r++){fixedsave[r]=fixed[r];}

//get ITD[0] = W/g TD - get even if noise constrained
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]<3)
{
for(j=Xstarts[r];j<Xends[r];j++){ITD[0][j]=TD[j]*Xsums[r]/vars[1+r];}
}
}

for(r=0;r<num_regs;r++)	//get I_k TD - similarly, get even if variance constrained
{
for(j=0;j<Xtotal;j++){ITD[1+r][j]=0;}
for(j=Xstarts[r];j<Xends[r];j++){ITD[1+r][j]=TD[j];}
}

for(r=0;r<total;r++)	//get TITD[1+r] = invT ITD
{
alpha=1.0;beta=0.0;
dgemv_("N", &Xtotal, &Xtotal, &alpha, T, &Xtotal, ITD[r], &one, &beta, TITD[r], &one);
}

//get traces[0] = tr(invT W/g)
traces[0]=0;
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]<3)
{
for(j=Xstarts[r];j<Xends[r];j++){traces[0]+=T[j+j*Xtotal]*Xsums[r]/vars[1+r];}
}
}

for(r=0;r<num_regs;r++)	//get traces[1+r] = tr(invT I_k)
{
traces[1+r]=0;
for(j=Xstarts[r];j<Xends[r];j++){traces[1+r]+=T[j+j*Xtotal];}
}

////////

//fill up AI = -2nd deriv and BI = 1st deriv - be careful about zero/-tive vars

//start assuming all terms fixed (AI diag, BI=0)
for(r=0;r<total;r++)
{
for(r2=0;r2<total;r2++){AI[r+r2*total]=0;}
AI[r+r*total]=1;BI[r]=0;
}

//now fill up
if(fixed[0]<3)
{
//first for noise x noise
gam2=0;for(j=0;j<Xtotal;j++){gam2+=ITD[0][j]*TD[j];}
gam3=0;for(j=0;j<Xtotal;j++){gam3+=TITD[0][j]*ITD[0][j];}
AI[0]=.5*gam*pow(vars[0],-3)-.5*gam2*pow(vars[0],-2)-.5*gam3/vars[0];

for(r=0;r<num_regs;r++)	//now noise x regions
{
if(fixed[1+r]<3)
{
gam3=0;for(j=0;j<Xtotal;j++){gam3+=TITD[0][j]*ITD[1+r][j];}
AI[1+r]=.5*gam3*Xsums[r]*pow(vars[1+r],-2);
AI[(1+r)*total]=AI[1+r];
}
}

BI[0]=.5*gam*pow(vars[0],-2)-.5*gam2/vars[0]-.5*traces[0]-.5*nfree/vars[0];
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]<3){BI[0]+=.5*(Xends[r]-Xstarts[r])/vars[0];}
}
}

//now for regions x regions
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]<3)
{
//diagonal
gam2=0;for(j=0;j<Xtotal;j++){gam2+=ITD[1+r][j]*TD[j];}
gam3=0;for(j=0;j<Xtotal;j++){gam3+=TITD[1+r][j]*ITD[1+r][j];}
AI[1+r+(1+r)*total]=.5*gam2*Xsums[r]*pow(vars[1+r],-3)-.5*gam3*pow(Xsums[r],2)*vars[0]*pow(vars[1+r],-4);

for(r2=0;r2<r;r2++)	//off diagonals
{
if(fixed[1+r2]<3)
{
gam3=0;for(j=0;j<Xtotal;j++){gam3+=TITD[1+r][j]*ITD[1+r2][j];}
AI[1+r+(1+r2)*total]=-.5*gam3*Xsums[r]*Xsums[r2]*vars[0]*pow(vars[1+r],-2)*pow(vars[1+r2],-2);
AI[1+r2+(1+r)*total]=AI[1+r+(1+r2)*total];
}}
BI[1+r]=.5*gam2*Xsums[r]*pow(vars[1+r],-2)+.5*traces[1+r]*Xsums[r]*vars[0]*pow(vars[1+r],-2)-.5*(Xends[r]-Xstarts[r])/vars[1+r];
}
}

//for stability, scale AI so has trace one
sum=0;for(r=0;r<total;r++){sum+=AI[r+r*total];}
for(r=0;r<total;r++)
{
for(r2=0;r2<total;r2++){AI[r+r2*total]*=pow(sum,-1);}
}

//invert then set fixed to zero
(void)eigen_invert(AI, total, AI2, -1, AI3, 1);
for(r=0;r<total;r++)
{
if(fixed[r]>=3){AI[r+r*total]=0;}
}

//undo scaling of AI
for(r=0;r<total;r++)
{
for(r2=0;r2<total;r2++){AI[r+r2*total]*=pow(sum,-1);}
}

//print update
sum=0;for(r=0;r<num_regs;r++){sum+=hers[1+r];}
value=0;for(r=0;r<total;r++){value+=vars[r];}
nlost=0;for(r=0;r<total;r++){nlost+=(fixed[r]>=3);}
nlost2=0;for(r=0;r<total;r++){nlost2+=(fixed[r]==1||fixed[r]==2);}

if(count==0){printf("Start\t");}
else{printf("%d\t", count);}
for(r=0;r<num_regs;r++){printf("%.4f\t", hers[1+r]);}
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
for(r=0;r<num_regs;r++){fprintf(output, "%.6f\t", hers[1+r]);}
if(count==0){fprintf(output, "%.6f\t%.6f\t%.6f\tNA\t%.6f\t%d\n", sum, value, like, tol, nlost);}
else{fprintf(output, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n", sum, value, like, diff, tol, nlost);}
fclose(output);

//see if breaking (normally can only break if rflag=0 and nlost2=0, unless at iter limit)

if(nlost>=num_regs){printf("All heritabilities are constrained\n");break;}
if(count>0)
{
if(fabs(diff)<tol&&rflag==0&&nlost2==0){break;}
}
if(count==maxiter)
{
printf("\nWarning, REML failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)", maxiter, tol);
if(constrain==0)
{printf(", and using \"--constrain YES\" to constrain heritabilities to [0,1]\n");}
else
{printf(", and using \"--constrain NO\" to allow heritabilities outside [0,1]\n");}
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
for(r=0;r<total;r++)
{
if(fabs(vardiffs[r])>max){max=fabs(vardiffs[r]);}
}
if(max>.1*varnull)	//then reduce moves
{
relax*=.1*varnull/max;
for(r=0;r<total;r++){vardiffs[r]*=.1*varnull/max;}
}

if(constrain==1)	//variances can not be negative
{
for(r=0;r<total;r++)
{
if(fixed[r]<3)	//free to update
{
if(vars[r]+vardiffs[r]<varmin){vardiffs[r]=varmin-vars[r];fixed[r]++;}
else{fixed[r]=0;}
if(fixed[r]==3){vardiffs[r]=-vars[r];}
}
}}

//now move
for(r=0;r<total;r++){vars[r]+=vardiffs[r];}
sum=0;for(r=0;r<total;r++){sum+=vars[r];}
for(r=0;r<total;r++){hers[r]=vars[r]/sum;}
sum=0;for(r=1;r<total;r++){sum+=vars[r];}
for(r=0;r<total;r++){shares[r]=vars[r]/sum;}

count++;
}	//end of while loop
printf("\n");

//get some stats
lrtstat=2*(like-likenull);
lrtpva=.5*erfc(pow(lrtstat,.5)*M_SQRT1_2);
if(lrtstat<0){lrtpva=.75;}
//if(hers[0]<=0){lrtpva=1;}

////////

if(Xnss==NULL)	//get fixed and random effects
{
//get fixed coefficients (ZTinvVZ)^-1 ZTinvVY with variance matrix (ZTinvVZ)^-1

//start with B = XTX + We/g, so that invV = (I - XBXT) / vars[0]
alpha=1.0;beta=0.0;	//not using summaries so no need for scale
dgemm_("T", "N", &Xtotal, &Xtotal, &ns, &alpha, X, &ns, X, &ns, &beta, B, &Xtotal);

//add on diagonal terms and blank (latter means no need to later blank ZTY, ZTX, ZTZ)
for(r=0;r<num_regs;r++)
{
for(j=Xstarts[r];j<Xends[r];j++)	//deal with row j
{
if(fixed[1+r]<3)	//leave or blank, adding eWr/vr to diagonal
{
for(r2=0;r2<num_regs;r2++)
{
for(j2=Xstarts[r2];j2<Xends[r2];j2++)
{
if(fixed[1+r2]>=3){B[j+j2*Xtotal]=0;}
}}
B[j+j*Xtotal]+=Xsums[r]*vars[0]/vars[1+r];
}
else	//blank row except for diagonal
{
for(j2=0;j2<Xtotal;j2++){B[j+j2*Xtotal]=0;}
B[j+j*Xtotal]=1;
}
}}

//invert B - save B in case ldlt fails
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){B3[j+j2*Xtotal]=B[j+j2*Xtotal];}
}
(void)ldlt_invert(B, Xtotal, -1, NULL, &info, 1);
if(info!=0)	//eigen will definitely work, so can use B3 as work space
{
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){B[j+j2*Xtotal]=B3[j+j2*Xtotal];}
}
(void)eigen_invert(B, Xtotal, B2, -1, B3, 1);
}

//remove ones from diagonal for blanked components
for(r=0;r<num_regs;r++)
{
if(fixed[1+r]>=3)
{
for(j=Xstarts[r];j<Xends[r];j++){B[j+j*Xtotal]=0;}
}
}

//get ZTXB and XTY (not using summaries, so will have X and Y)
alpha=1.0;beta=0.0;	//not using summaries so no need for scale
dgemm_("N", "N", &num_fixed, &Xtotal, &Xtotal, &alpha, ZTX, &num_fixed, B, &Xtotal, &beta, ZTXB, &num_fixed);
dgemv_("T", &ns, &Xtotal, &alpha, X, &ns, Y, &one, &beta, XTY, &one);

//now ZTVY = (ZTY - ZTXB XTY) and ZTVZ = (ZTZ - ZTXB XTZ) - these actually ZTVY and ZTVZ times vars[0]
alpha=1.0;beta=0.0;	//not using summaries so no need for scale
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, Y, &one, &beta, ZTVY, &one);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, Z, &ns, &beta, ZTVZ, &num_fixed);
alpha=-1.0;beta=1.0;
dgemv_("N", &num_fixed, &Xtotal, &alpha, ZTXB, &num_fixed, XTY, &one, &beta, ZTVY, &one);
dgemm_("N","T", &num_fixed, &num_fixed, &Xtotal, &alpha, ZTXB, &num_fixed, ZTX, &num_fixed, &beta, ZTVZ, &num_fixed);

//get theta = (ZTinvVZ)^-1 ZTinvVY (the vars[0] will cancel out)
(void)eigen_invert(ZTVZ, num_fixed, ZTVZ2, -1, ZTVZ3, 1);
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &num_fixed, &alpha, ZTVZ, &num_fixed, ZTVY, &one, &beta, thetas, &one);

for(j=0;j<num_fixed;j++)	//now have to multiply inverse variance matrix by vars[0]
{
if(ZTVZ[j+j*num_fixed]>=0)
{
thetasds[j]=pow(ZTVZ[j+j*num_fixed]*vars[0],.5);
thetapvas[j]=erfc(fabs(thetas[j]/thetasds[j])*M_SQRT1_2);
}
else{thetasds[j]=-9999;thetapvas[j]=-9999;}
}

//set Yadj to fixed effects
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &num_fixed, &alpha, Z, &ns, thetas, &one, &beta, Yadj, &one);

//random effects are Xr invT XrT CY
for(r=0;r<num_regs;r++)
{
token=Xends[r]-Xstarts[r];
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &token, &alpha, X+ns*Xstarts[r], &ns, TD+Xstarts[r], &one, &beta, mg[r], &one);

//add random effects to Yadj
for(i=0;i<ns;i++){Yadj[i]+=mg[r][i];}
}
}	//end of getting fixed and random effects

////////

//load up effects for regions (saved in TD) and tops (saved in thetas)
{
//effects contain effects for regions, tops, sum, centres, predictor used?
for(r=0;r<num_regs+4;r++)
{
for(j=0;j<num_preds;j++){effects[r][j]=0;}
}

for(r=0;r<num_regs;r++)	//note, will have removed trivial regional predictors and those with zero weight
{
if(vars[1+r]!=0)
{
for(j=Xstarts[r];j<Xends[r];j++)
{
j2=regindex[r][1+j-Xstarts[r]];
effects[r][rkeeppreds[j2]]+=TD[j]*rmults[j2]*pow(rweights[j2],.5);
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
for(r=0;r<total;r++)
{
if(AI[r+r*total]>=0){varsds[r]=pow(AI[r+r*total],.5);}
else{varsds[r]=-9999;}
}

//get SEs of hers - Jij=delta/sum-vari/sum^2 (where sum across all elements)
sum=0;for(r=0;r<total;r++){sum+=vars[r];}

for(r=0;r<total;r++)
{
for(r2=0;r2<total;r2++)
{
if(fixed[r]<3&&fixed[r2]<3){J[r+r2*total]=-vars[r]*pow(sum,-2);}
else{J[r+r2*total]=0;}
}
if(fixed[r]<3){J[r+r*total]+=1.0/sum;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

//load up
for(r=1;r<total;r++)
{
if(JAIJT[r+r*total]>=0){hersds[r]=pow(JAIJT[r+r*total],.5);}
else{hersds[r]=-9999;}
}

//save details for her_all as final element of hers / hersds
hers[total]=0;for(r=1;r<total;r++){hers[total]+=hers[r];}
sum=0;
for(r=1;r<total;r++)
{
for(r2=1;r2<total;r2++){sum+=JAIJT[r+r2*total];}
}
if(sum>=0){hersds[total]=pow(sum,.5);}
else{hersds[total]=-9999;}

//save JAIJT for printing
for(r=0;r<total;r++)
{
for(r2=0;r2<total;r2++){cohers[r+r2*total]=JAIJT[r+r2*total];}
}

//get SEs of shares - Jij=delta/sum-vari/sum^2 (where sum excludes noise term)
sum=0;for(r=1;r<total;r++){sum+=vars[r];}

//first row and column (easier to simply set new0=var0)
if(fixed[0]<3){J[0]=1;}
else{J[0]=0;}
for(r=1;r<total;r++){J[r]=0;J[r*total]=0;}

//rest
for(r=1;r<total;r++)
{
for(r2=1;r2<total;r2++)
{
if(fixed[r]<3&&fixed[r2]<3){J[r+r2*total]=-vars[r]*pow(sum,-2);}
else{J[r+r2*total]=0;}
}
if(fixed[r]<3){J[r+r*total]+=1.0/sum;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

//load up
for(r=0;r<total;r++)
{
if(JAIJT[r+r*total]>=0){sharesds[r]=pow(JAIJT[r+r*total],.5);}
else{sharesds[r]=-9999;}
}

///////////////////////////

if(Xnss==NULL)	//get variance due to covariates and top predictors - will already have ZTY
{
sumsq=-ZTY[0]/ns*ZTY[0];
for(i=0;i<ns;i++){sumsq+=pow(Y[i],2);}
value=-ZTY[0]/ns*ZTY[0];for(j=0;j<num_covars;j++){value+=ZTY[j]*thetas[j];}
value2=0;for(j=num_covars;j<num_fixed;j++){value2+=ZTY[j]*thetas[j];}
covher=value/sumsq;
topher=value2/(sumsq-value);

if(num_covars>1){printf("Proportion of variance explained by the %d covariates: %.4f\n", num_covars, covher);}
if(num_tops==1){printf("Proportion of variance explained by the top predictor: %.4f\n", topher);}
if(num_tops>1){printf("Proportion of variance explained by the %d top predictors: %.4f\n", num_tops, topher);}
if(num_covars>1||num_tops>1){printf("\n");}
}
else{covher=0;topher=0;}

//adjust for tops
for(r=0;r<total+1;r++)
{
hers[r]*=(1-topher);
if(hersds[r]!=-9999){hersds[r]*=(1-topher);}
}

//save stuff

sum=0;for(r=0;r<num_regs;r++){sum+=Xsums[r];}

sprintf(filename2,"%s.reml", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Num_Kinships 0\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_regs, num_tops, num_covars, num_envs);
if(Xnss==NULL){fprintf(output2, "Blupfile %s.indi.blp\n", outfile);}
else{fprintf(output2, "Blupfile none\n");}
fprintf(output2, "Regfile %s.reg.blup\n", outfile);
if(Xnss==NULL){fprintf(output2, "Coeffsfile %s.coeff\n", outfile);}
else{fprintf(output2, "Coeffsfile none\n");}
fprintf(output2, "Covar_Heritability %.4f\n", covher);
if(Xnss==NULL){fprintf(output2, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);}
else{fprintf(output2, "Num_Genotyped_Samples %d\nAverage_Sample_Size %.2f\n", ns, nfree);}
if(cflag==1){fprintf(output2,"Converged YES\n");}
else{fprintf(output2,"Converged NO\n");}
fprintf(output2, "Null_Likelihood %.4f\nAlt_Likelihood %.4f\n", likenull, like);
if(num_regs==1){fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output2, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output2, "Component Heritability SE Size Mega_Intensity SE\n");
for(r=0;r<num_regs;r++){fprintf(output2, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[1+r], hersds[1+r], Xsums[r], hers[1+r]/Xsums[r]*1000000, hersds[1+r]/Xsums[r]*1000000);}
fprintf(output2, "Her_Top %.6f NA NA NA NA\n", topher);
fprintf(output2, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]+topher, hersds[total], sum, hers[total]/sum*1000000, hersds[total]/sum*1000000);
fclose(output2);

if(Xnss==NULL)
{
sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Component Effect SE P\n");
fprintf(output3, "Intercept %.4e %.4e %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.4e %.4e %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);} 
for(j=0;j<num_envs;j++){fprintf(output3, "Enviromental_%d %.4e %.4e %.4e\n",j, thetas[num_covars+j], thetasds[num_covars+j], thetapvas[num_covars+j]);}
fclose(output3);
}

sprintf(filename4,"%s.share", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Share SE Expected Enrichment SE\n");
for(r=0;r<num_regs;r++){fprintf(output4, "Share_R%d %.6f %.6f %.6f %.6f %.6f\n", r+1, shares[1+r], sharesds[1+r], Xsums[r]/sum, shares[1+r]/Xsums[r]*sum, sharesds[1+r]/Xsums[r]*sum);}
fclose(output4);

sprintf(filename5,"%s.vars", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Component Variance SE\n");
for(r=0;r<num_regs;r++){fprintf(output5, "Var_R%d %.6f %.6f\n", r+1, vars[1+r], varsds[1+r]);}
fprintf(output5, "Var_E %.6f %.6f\n", vars[0], varsds[0]);
fclose(output5);

if(Xnss==NULL)
{
sprintf(filename6,"%s.indi.res", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "ID1\tID2\tPhenotype\tFitted\tResidual\n");
for(i=0;i<ns;i++){fprintf(output6, "%s\t%s\t%f\t%f\t%f\n", ids1[i], ids2[i], Y[i], Yadj[i], Y[i]-Yadj[i]);}
fclose(output6);

sprintf(filename7,"%s.indi.blp", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(i=0;i<ns;i++)
{
fprintf(output7, "%s\t%s\t", ids1[i], ids2[i]);
for(r=0;r<num_regs;r++){fprintf(output7, "0\t%.6f\t", mg[r][i]);}
fprintf(output7, "\n");
}
fclose(output7);
}

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

sprintf(filename10,"%s.cross", outfile);
if((output10=fopen(filename10,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename10);exit(1);}
for(r=0;r<num_regs;r++){fprintf(output10, "Her_R%d\t", r+1);}
fprintf(output10, "\n");
for(r=1;r<total;r++)
{
for(r2=1;r2<total;r2++)
{fprintf(output10, "%.6f\t", cohers[r+r2*total]);}
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
fprintf(output, "Num_Kinships 0\nNum_Regions %d\nNum_Top_Predictors %d\nNum_Covariates %d\nNum_Environments %d\n", num_regs, num_tops, num_covars, num_envs);
if(Xnss==NULL){fprintf(output, "Blupfile %s.indi.blp.liab\n", outfile);}
else{fprintf(output, "Blupfile none\n");}
fprintf(output, "Regfile %s.reg.blup.liab\n", outfile);
if(Xnss==NULL){fprintf(output, "Coeffsfile %s.coeff.liab\n", outfile);}
else{fprintf(output, "Coeffsfile none\n");}
fprintf(output, "Covar_Heritability %.4f\n", covher*factor);
if(Xnss==NULL){fprintf(output, "Total_Samples %d\nWith_Phenotypes %d\n", ns, np);}
else{fprintf(output, "Num_Genotyped_Samples %d\nAverage_Sample_Size %.2f\n", ns, nfree);}
if(cflag==1){fprintf(output,"Converged YES\n");}
else{fprintf(output,"Converged NO\n");}
fprintf(output, "Null_Likelihood %.4f\nAlt_Likelihood %.4f\n", likenull, like);
if(num_regs==1){fprintf(output, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);}
else{fprintf(output, "LRT_Stat %.4f\nLRT_P NA\n", lrtstat);}

fprintf(output, "Component Heritability SE Size Mega_Intensity SE\n");
for(r=0;r<num_regs;r++){fprintf(output, "Her_R%d %.6f %.6f %.2f %.6f %.6f\n", r+1, hers[1+r]*factor, hersds[1+r]*factor, Xsums[r], hers[1+r]/Xsums[r]*1000000*factor, hersds[1+r]/Xsums[r]*1000000*factor);}
fprintf(output, "Her_Top %.6f NA NA NA NA\n", topher*factor);
fprintf(output, "Her_All %.6f %.6f %.2f %.6f %.6f\n", hers[total]*factor+topher*factor, hersds[total]*factor, sum, hers[total]/sum*1000000*factor, hersds[total]/sum*1000000*factor);
fclose(output);

if(Xnss==NULL)
{
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
for(r=0;r<num_regs;r++){fprintf(output, "0\t%.6f\t", mg[r][i]*pow(factor,.5));}
fprintf(output, "\n");
}
fclose(output);
}

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
{
fprintf(output,"%s\t%c\t%c\t%.6f\t", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]);
for(r=0;r<num_regs;r++){fprintf(output, "%.6f\t", effects[r][j]*pow(factor,.5));}
if(num_tops>0){fprintf(output, "%.6f\n", effects[num_regs][j]*pow(factor,.5));}
else{fprintf(output,"\n");}
}
}
fclose(output);

sprintf(filename,"%s.liab", filename9);	//.reg.score
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\tEffect\n");
for(j=0;j<num_preds;j++)
{
if(effects[num_regs+3][j]==1)
{fprintf(output,"%s\t%c\t%c\t%.6f\t%.6f\n", allpreds[j], allal1[j], allal2[j], effects[num_regs+2][j]*pow(factor,.5), effects[num_regs+1][j]*pow(factor,.5));}
}
fclose(output);
}	//end of binary

printf("Main results saved in %s", filename2);
if(prev!=-9999){printf(", with a liability version saved in %s.liab", filename2);}
printf("\n\n");

////////

free(fixed);free(fixedsave);free(vars);free(vardiffs);free(varsds);free(hers);free(hersds);free(shares);free(sharesds);
free(AI);free(AI2);free(AI3);free(BI);free(J);free(JAI);free(JAIJT);free(cohers);
free(ZTY);free(ZTX);free(ZTZ);free(ZTZ2);free(ZTZ3);free(ZTZZTY);free(ZTZZTX);
free(XTCY);free(XTCX);free(T);free(T2);free(T3);free(TD);
for(r=0;r<total;r++){free(ITD[r]);free(TITD[r]);}free(ITD);free(TITD);free(traces);
free(B);free(B2);free(B3);free(ZTXB);free(XTY);
free(ZTVY);free(ZTVZ);free(ZTVZ2);free(ZTVZ3);free(thetas);free(thetasds);free(thetapvas);free(Yadj);

for(r=0;r<num_regs;r++){free(mg[r]);}free(mg);
for(r=0;r<num_regs;r++){free(mg2[r]);}free(mg2);
for(r=0;r<num_regs+4;r++){free(effects[r]);}free(effects);

}	//end of fast_reml

///////////////////////////


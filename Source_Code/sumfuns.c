/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the fGNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//summary functions

//////////////////////////

void solve_sums_pairs(double *stats, int length, double *ldscores, double *stags, double *schis, int num_blocks)
//lite version that only estimates ratio of intersects
{
int j, p, count, start, end, one=1;
double value, sum, inter, alpha, beta;

int total;
double tol=0.001, relax, like, like2, likeold, diff;
double thetas[2], thetadiffs[2], *exps, *exps2;
double *sT, *sT2, *AI, *AI2, *AI3, *BI;


total=2;

exps=malloc(sizeof(double)*length);
exps2=malloc(sizeof(double)*length);

sT=malloc(sizeof(double)*length*total);
sT2=malloc(sizeof(double)*length*total);
AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);

for(p=0;p<num_blocks+1;p++)
{
if(p<num_blocks)	//work out which predictors to exclude
{
start=(double)(p)/num_blocks*length;
end=(double)(p+1)/num_blocks*length;
}

//set sT (contain ldscores followed by 1)
for(j=0;j<length;j++){sT[j]=ldscores[j];sT[j+length]=1.0;}
if(p<num_blocks)	//blank block terms
{
for(j=start;j<end;j++){sT[j]=0;sT[j+length]=0;}
}

//solve for first column of schis

//set starting model (assumes no causal variation and intercept 1)
thetas[0]=0;
thetas[1]=1;
for(j=0;j<length;j++){exps[j]=1.0;}

//set starting likelihood
sum=0;value=0;
for(j=0;j<length;j++){sum+=schis[j]/exps[j]/stags[j];value+=log(schis[j]*exps[j])/stags[j];}
if(p<num_blocks)	//blank block terms
{
for(j=start;j<end;j++){sum-=schis[j]/exps[j]/stags[j];value-=log(schis[j]*exps[j])/stags[j];}
}
like=-.5*sum-.5*value;

//ready to iterate
count=0;
while(1)
{
//load up sT2=(stat-exps/2)/exps^3/stags*sT
for(j=0;j<length;j++)
{
value=(schis[j]-.5*exps[j])*pow(exps[j],-3)/stags[j];
sT2[j]=value*sT[j];
sT2[j+length]=value*sT[j+length];
}

//get AI (-2nd deriv) and BI (1st deriv)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sT, &length, sT2, &length, &beta, AI, &total);
BI[0]=0;for(j=0;j<length;j++){BI[0]+=.5*(schis[j]-exps[j])*pow(exps[j],-2)/stags[j]*sT[j];}
BI[1]=0;for(j=0;j<length;j++){BI[1]+=.5*(schis[j]-exps[j])*pow(exps[j],-2)/stags[j]*sT[j+length];}

(void)eigen_invert(AI, total, AI2, -1, AI3, 1);

//get proposed move
alpha=1.0;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, thetadiffs, &one);

relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs and get expectations
thetas[0]+=relax*thetadiffs[0];
thetas[1]+=relax*thetadiffs[1];

alpha=1.0;beta=0.0;
dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta, exps2, &one);

for(j=0;j<length;j++)
{
if(exps2[j]<=0){exps2[j]=1e-6;}
}

//get corresponding likelihood
sum=0;value=0;
for(j=0;j<length;j++){sum+=schis[j]/exps2[j]/stags[j];value+=log(schis[j]*exps2[j])/stags[j];}
if(p<num_blocks)	//blank block terms
{
for(j=start;j<end;j++){sum-=schis[j]/exps2[j]/stags[j];value-=log(schis[j]*exps2[j])/stags[j];}
}
like2=-.5*sum-.5*value;

if(like2>like-tol)	//accept move
{
like=like2;
for(j=0;j<length;j++){exps[j]=exps2[j];}
break;
}
else	//move back and next turn try smaller move
{
thetas[0]-=relax*thetadiffs[0];
thetas[1]-=relax*thetadiffs[1];
relax*=.5;
}
}

if(count>0)	//see if breaking
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;

count++;
}

inter=thetas[1];

//solve for second column of chis

//set starting model (assumes no causal variation and intercept 1)
thetas[0]=0;
thetas[1]=1;
for(j=0;j<length;j++){exps[j]=1.0;}

//set starting likelihood
sum=0;value=0;
for(j=0;j<length;j++){sum+=schis[j+length]/exps[j]/stags[j];value+=log(schis[j+length]*exps[j])/stags[j];}
if(p<num_blocks)	//blank block terms
{
for(j=start;j<end;j++){sum-=schis[j+length]/exps[j]/stags[j];value-=log(schis[j+length]*exps[j])/stags[j];}
}
like=-.5*sum-.5*value;

//ready to iterate
count=0;
while(1)
{
//load up sT2=(stat-exps/2)/exps^3/stags*sT
for(j=0;j<length;j++)
{
value=(schis[j+length]-.5*exps[j])*pow(exps[j],-3)/stags[j];
sT2[j]=value*sT[j];
sT2[j+length]=value*sT[j+length];
}

//get AI (-2nd deriv) and BI (1st deriv)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sT, &length, sT2, &length, &beta, AI, &total);
BI[0]=0;for(j=0;j<length;j++){BI[0]+=.5*(schis[j+length]-exps[j])*pow(exps[j],-2)/stags[j]*sT[j];}
BI[1]=0;for(j=0;j<length;j++){BI[1]+=.5*(schis[j+length]-exps[j])*pow(exps[j],-2)/stags[j]*sT[j+length];}

(void)eigen_invert(AI, total, AI2, -1, AI3, 1);

//get proposed move
alpha=1.0;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, thetadiffs, &one);

relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs and get expectations
thetas[0]+=relax*thetadiffs[0];
thetas[1]+=relax*thetadiffs[1];

alpha=1.0;beta=0.0;
dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta, exps2, &one);

for(j=0;j<length;j++)
{
if(exps2[j]<=0){exps2[j]=1e-6;}
}

//get corresponding likelihood
sum=0;value=0;
for(j=0;j<length;j++){sum+=schis[j+length]/exps2[j]/stags[j];value+=log(schis[j+length]*exps2[j])/stags[j];}
if(p<num_blocks)	//blank block terms
{
for(j=start;j<end;j++){sum-=schis[j+length]/exps2[j]/stags[j];value-=log(schis[j+length]*exps2[j])/stags[j];}
}
like2=-.5*sum-.5*value;

if(like2>like-tol)	//accept move
{
like=like2;
for(j=0;j<length;j++){exps[j]=exps2[j];}
break;
}
else	//move back and next turn try smaller move
{
thetas[0]-=relax*thetadiffs[0];
thetas[1]-=relax*thetadiffs[1];
relax*=.5;
}
}

if(count>0)	//see if breaking
{
diff=like-likeold;
if(fabs(diff)<tol){break;}
}
likeold=like;

count++;
}

//load ratio into stats
stats[p]=thetas[1]/inter;
}	//end of p loop

free(exps);free(exps2);
free(sT);free(sT2);free(AI);free(AI2);free(AI3);free(BI);
}

////////

void solve_sums(double *stats, double *likes, double *cohers, double *cohers2, double *influs, int num_parts, int gcon, int cept, int num_blocks, int length, int ncv, int *cvindex, double *cvexps, double *stags, double **svars, double **ssums, double *snss, double *schis, double tol, int maxiter, int chisol, int sflag, char *filename)
//sflag=0 - normal, sflag=1 - first pass, sflag=2 - second pass
//sflag=3 - just get expectations and likelihood, sflag=4 - LDSC, sflag=5 - divide+updating, sflag=6 - quiet
{
int j, j2, p, q, q2, q3, count, count2, start, end, mark, one=1;
double value, value2, sum, sum2, sumsq, mean, mean2, var, alpha, beta;

int total, cflag, rflag;
double scale, gc, sumhers, relax, likenull, like, like2, likeold, diff;
double *thetas, *thetadiffs, *exps, *exps2, *jacks;
double *sW, *sX, *sY, *sXTX, *sXTX2, *sXTY, *sXTXs, *sXTYs, *sT, *sTb, *sT2, *AI, *AI2, *AI3, *BI, *J, *JAI, *JAIJT;

FILE *output;


//assuming statj = c * (1 + nja + njvj b), where c is gc, a is intercept/n, b is h2SNP/ssums
//easier to write as (statj-1) = nj/nvj cbn + (c-1) + nj/n can = sT theta, where n is average nj

//set total and maybe num_blocks
total=num_parts+gcon+cept;
if(num_blocks==-1||num_blocks>length){num_blocks=length;}

//allocate variables

thetas=malloc(sizeof(double)*total);
thetadiffs=malloc(sizeof(double)*total);
exps=malloc(sizeof(double)*length);
exps2=malloc(sizeof(double)*length);
if(num_blocks!=-9999){jacks=malloc(sizeof(double)*(total+1+num_parts)*num_blocks);}

sW=malloc(sizeof(double)*length);
sX=malloc(sizeof(double)*total*length);
sY=malloc(sizeof(double)*length);
sXTX=malloc(sizeof(double)*total*total);
sXTX2=malloc(sizeof(double)*total*total);
sXTY=malloc(sizeof(double)*total);
sXTXs=malloc(sizeof(double)*total*total);
sXTYs=malloc(sizeof(double)*total);

sT=malloc(sizeof(double)*length*total);
sTb=malloc(sizeof(double)*length*total);
sT2=malloc(sizeof(double)*length*total);
AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);
J=malloc(sizeof(double)*total*total);
JAI=malloc(sizeof(double)*total*total);
JAIJT=malloc(sizeof(double)*total*total);

//set variables that stay the same throughout

//get scale - weighted average sample size
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=snss[j]/stags[j];sum2+=pow(stags[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=snss[j2]/stags[j2];sum2-=pow(stags[j2],-1);}
scale=sum/sum2;

//set sT and sTb (latter might not be used)
for(j=0;j<length;j++)
{
for(q=0;q<num_parts;q++){sT[j+q*length]=snss[j]/scale*svars[q][j];}
if(gcon==1){sT[j+num_parts*length]=1;}
if(cept==1){sT[j+(num_parts+gcon)*length]=snss[j]/scale;}
}

//sTb matches sT, except blanked out for cv predictors
for(j=0;j<length;j++)
{
for(q=0;q<total;q++){sTb[j+q*length]=sT[j+q*length];}
}
for(j=0;j<ncv;j++)
{
for(q=0;q<total;q++){sTb[cvindex[j]+q*length]=0;}
}

////////

//get null likelihood
if(cept+gcon>0)	//null model will be E[Sj]=constant
{
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=schis[j]/stags[j];sum2+=pow(stags[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/stags[j2];sum2-=pow(stags[j2],-1);}
value2=sum/sum2;
}
else{value2=1;}

sum=0;sum2=0;value=0;
for(j=0;j<length;j++){sum+=schis[j]/value2/stags[j];sum2+=pow(stags[j],-1);value+=log(schis[j]*value2)/stags[j];}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/value2/stags[j2];sum2-=pow(stags[j2],-1);value-=log(schis[j2]*value2)/stags[j2];}
likenull=-.5*sum-.5*value-.5*sum2*log(2*M_PI);

//set starting model

if(sflag==0||sflag==1)	//starting model assumes no causal variation (plus gcon=1 and cept=1)
{
for(q=0;q<total;q++){thetas[q]=0;}
}

if(sflag==2||sflag==3||sflag==5||sflag==6)	//starting model saved in stats
{
gc=1;if(gcon==1){gc=stats[num_parts];}
for(q=0;q<num_parts;q++){thetas[q]=stats[q]*gc/ssums[q][q]*scale;}
if(gcon==1){thetas[num_parts]=(gc-1);}
if(cept==1){thetas[num_parts+gcon]=stats[num_parts+gcon]*gc;}
}

if(sflag==4)	//starting model copies ldsc - theta = (sum_j stat -1) / sum_jq (n_j var_jq) 
{
value=0;value2=0;
for(j=0;j<length;j++)
{
value+=schis[j]-1;
for(q=0;q<num_parts;q++){value2+=snss[j]*svars[q][j];}
}
for(j=0;j<ncv;j++)
{
value-=schis[cvindex[j]]-1;
for(q=0;q<num_parts;q++){value2-=snss[cvindex[j]]*svars[q][cvindex[j]];}
}

for(q=0;q<num_parts;q++){thetas[q]=value/value2;}
if(gcon==1){thetas[num_parts]=0;}
if(cept==1){thetas[num_parts+gcon]=0;}
}

if(sflag!=6)	//prepare to screen and file print
{
printf("Iter\tHer_All\t");
if(gcon==1){printf("Scaling\t");}
if(cept==1){printf("Intercept\t");}
printf("Likelihood\tDifference\tTarget\n");

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
if(sflag==1){fprintf(output,"Solving for first %d categories\n", num_parts);}
if(sflag==2){fprintf(output,"Solving for all %d categories\n", num_parts);}
fprintf(output,"Iter\tHer_All\t");
if(gcon==1){fprintf(output, "Scaling\t");}
if(cept==1){fprintf(output, "Intercept\t");}
fprintf(output, "Likelihood\tDifference\tTarget\n");
fclose(output);
}

////////

//now iterate - rflag indicates type of move
count=0;
cflag=1;
rflag=0;	//0 if using multiNR (or least-squares), 1 if about to do singleNR, 2 if just done singleNR
while(1)
{
if(chisol==0||count==0)	//get exps and likelihood (already have for NR solver after first iteration)
{
for(j=0;j<length;j++){exps[j]=1;}
alpha=1.0;beta=1.0;
dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta, exps, &one);
for(j=0;j<length;j++)
{
if(exps[j]<=0){exps[j]=1e-6;}
}

sum=0;sum2=0;value=0;
for(j=0;j<length;j++){sum+=schis[j]/exps[j]/stags[j];sum2+=pow(stags[j],-1);value+=log(schis[j]*exps[j])/stags[j];}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/exps[j2]/stags[j2];sum2-=pow(stags[j2],-1);value-=log(schis[j2]*exps[j2])/stags[j2];}
like=-.5*sum-.5*value-.5*sum2*log(2*M_PI);
}

if(count>0)	//set diff
{diff=like-likeold;}
likeold=like;

if(sflag!=6)	//print update
{
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
sumhers=0;for(q=0;q<num_parts;q++){sumhers+=thetas[q]/gc*ssums[q][q]/scale;}

if(count==0){printf("Start\t");}
else{printf("%d\t", count);}
printf("%.4f\t", sumhers);
if(gcon==1){printf("%.4f\t", gc);}
if(cept==1){printf("%.4f\t", 1+thetas[num_parts+gcon]/gc);}
printf("%.2f\t", like);
if(count==0){printf("n/a\t\t%.6f\n", tol);}
else{printf("%.6f\t%.6f\n", diff, tol);}

if(sflag!=6)
{
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "%d\t%.6f\t", count, sumhers);
if(gcon==1){fprintf(output, "%.6f\t", gc);}
if(cept==1){fprintf(output, "%.6f\t", 1+thetas[num_parts+gcon]/gc);}
if(count==0){fprintf(output, "%.6f\tNA\t%.6f\n", like, tol);}
else{fprintf(output, "%.6f\t%.6f\t%.6f\n", like, diff, tol);}
fclose(output);
}
}

//see if breaking (normally can only break if rflag=0, unless at iter limit)

if(sflag==3){break;}	//only wanted expectations and likelihood
if(sflag==4&&num_parts>1&&count==1){break;}	//multi-tagging ldsc uses a single iteration
if(count>0)
{
if(fabs(diff)<tol&&(rflag==0||chisol==0)){break;}
}
if(count==maxiter){printf("\nWarning, the optimizer failed to converge within %d iterations\n", maxiter);cflag=0;break;}

////////

//update thetas

if(chisol==0)	//use least squares (will always accept proposed move)
{
//load up sW, sX and sY - regressing Sj-1 on nj*vjk, 1, nj
for(j=0;j<length;j++)
{
sW[j]=stags[j]*pow(exps[j],2);
value=pow(sW[j],-0.5);
for(q=0;q<num_parts;q++){sX[j+q*length]=snss[j]/scale*svars[q][j]*value;}
if(gcon==1){sX[j+num_parts*length]=value;}
if(cept==1){sX[j+(num_parts+gcon)*length]=snss[j]/scale*value;}
sY[j]=(schis[j]-1)*value;
}

for(j=0;j<ncv;j++)	//blank out sX for cv predictors
{
for(q=0;q<total;q++){sX[cvindex[j]+q*length]=0;}
}

//new theta is (sXTX)^-1 sXTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sX, &length, sX, &length, &beta, sXTX, &total);
dgemv_("T", &length, &total, &alpha, sX, &length, sY, &one, &beta, sXTY, &one);

for(q=0;q<total;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total, sXTX2, 1, thetas, 1);
}
else	//use either single nr or multi nr (only move if good for likelihood)
{
if(rflag==1)	//single nr
{
for(q=total-1;q>=0;q--)	//do backwards (might help if last category is base)
{
//get derivs for q - use sTb so as not to include contributions of cv predictors
value=0;value2=0;
for(j=0;j<length;j++)
{
value+=.5*(schis[j]-exps[j])/stags[j]*pow(exps[j],-2)*sTb[j+q*length];
value2+=(schis[j]-.5*exps[j])/stags[j]*pow(exps[j],-3)*pow(sTb[j+q*length],2);
}

//get proposed move
thetadiffs[q]=value/value2;

//move and get expectations
thetas[q]+=thetadiffs[q];

for(j=0;j<length;j++){exps2[j]=1;}
alpha=1.0;beta=1.0;
dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta, exps2, &one);

for(j=0;j<length;j++)
{
if(exps2[j]<=0){exps2[j]=1e-6;}
}

//get corresponding likelihood
sum=0;sum2=0;value=0;
for(j=0;j<length;j++){sum+=schis[j]/exps2[j]/stags[j];sum2+=pow(stags[j],-1);value+=log(schis[j]*exps2[j])/stags[j];}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/exps2[j2]/stags[j2];sum2-=pow(stags[j2],-1);value-=log(schis[j2]*exps2[j2])/stags[j2];}
like2=-.5*sum-.5*value-.5*sum2*log(2*M_PI);

if(like2>like-tol)	//accept move
{
like=like2;
for(j=0;j<length;j++){exps[j]=exps2[j];}
}
else	//move back
{thetas[q]-=thetadiffs[q];}
}	//end of q loop
rflag=2;
}
else	//multi nr - use BI and AI
{
//load up sT2=(stat-exps/2)*w/exps^3*sTb - will be blanked for cvpredictors
for(j=0;j<length;j++)
{
value=(schis[j]-.5*exps[j])/stags[j]*pow(exps[j],-3);
for(q=0;q<total;q++){sT2[j+q*length]=value*sTb[j+q*length];}
}

//get AI (-2nd deriv) and BI (1st deriv) - cvpredictors already taken care of (by blanking sTb and ST2)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sTb, &length, sT2, &length, &beta, AI, &total);
for(q=0;q<total;q++)
{
BI[q]=0;for(j=0;j<length;j++){BI[q]+=.5*(schis[j]-exps[j])/stags[j]*pow(exps[j],-2)*sTb[j+q*length];}
}

(void)eigen_invert(AI, total, AI2, -1, AI3, 1);

//get proposed move
alpha=1.0;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, thetadiffs, &one);

rflag=1;
relax=1;
while(relax>0.0001)
{
//move relax*thetadiffs and get expectations
for(q=0;q<total;q++){thetas[q]+=relax*thetadiffs[q];}

for(j=0;j<length;j++){exps2[j]=1;}
alpha=1.0;beta=1.0;
dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta, exps2, &one);

for(j=0;j<length;j++)
{
if(exps2[j]<=0){exps2[j]=1e-6;}
}

//get corresponding likelihood
sum=0;sum2=0;value=0;
for(j=0;j<length;j++){sum+=schis[j]/exps2[j]/stags[j];sum2+=pow(stags[j],-1);value+=log(schis[j]*exps2[j])/stags[j];}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/exps2[j2]/stags[j2];sum2-=pow(stags[j2],-1);value-=log(schis[j2]*exps2[j2])/stags[j2];}
like2=-.5*sum-.5*value-.5*sum2*log(2*M_PI);

if(like2>like-tol)	//accept move
{
like=like2;
for(j=0;j<length;j++){exps[j]=exps2[j];}
rflag=0;
break;
}
else	//move back and next turn try smaller move
{
for(q=0;q<total;q++){thetas[q]-=relax*thetadiffs[q];}
relax*=.5;
}
}
}
}	//end of NR

count++;
}	//end of while loop
if(sflag!=6){printf("\n");}

//load up first column of stats which contain Q hers, gc, na, sum of Q hers, Q cats
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
sumhers=0;for(q=0;q<num_parts;q++){sumhers+=thetas[q]/gc*ssums[q][q]/scale;}

for(q=0;q<num_parts;q++){stats[q]=thetas[q]/gc*ssums[q][q]/scale;}
if(gcon==1){stats[num_parts]=gc;}
if(cept==1){stats[num_parts+gcon]=thetas[num_parts+gcon]/gc;}
stats[total]=sumhers;
for(q=0;q<num_parts;q++)
{
stats[total+1+q]=0;
for(q2=0;q2<num_parts;q2++){stats[total+1+q]+=thetas[q2]/gc*ssums[q][q2]/scale;}
}

////////

if(likes!=NULL)	//save likelihoods - currently redundant if using cv predictors
{
//already have chisq likelihoods
likes[0]=likenull;
likes[1]=like;

//normal likelihood using starting weights
if(cept+gcon>0)	//for null fit E[Sj]=constant
{
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=schis[j]/stags[j];sum2+=pow(stags[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/stags[j2];sum2-=pow(stags[j2],-1);}
value2=sum/sum2;
}
else{value2=1;}

sum=0;sum2=0;
for(j=0;j<length;j++){sum+=pow(schis[j]-value2,2)/stags[j];sum2+=pow(stags[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=pow(schis[j2]-value2,2)/stags[j2];sum2-=pow(stags[j2],-1);}
likes[2]=-.5*sum2*(1+log(2*M_PI*sum/sum2));

sum=0;sum2=0;
for(j=0;j<length;j++){sum+=pow(schis[j]-exps[j],2)/stags[j];sum2+=pow(stags[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=pow(schis[j2]-exps[j2],2)/stags[j2];sum2-=pow(stags[j2],-1);}
likes[3]=-.5*sum2*(1+log(2*M_PI*sum/sum2));

//normal likelihood using final weights
for(j=0;j<length;j++){sW[j]=stags[j]*pow(exps[j],2);}

if(cept+gcon>0)	//for null fit E[Sj]=constant
{
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=schis[j]/sW[j];sum2+=pow(sW[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=schis[j2]/sW[j2];sum2-=pow(sW[j2],-1);}
value2=sum/sum2;
}
else{value2=1;}

sum=0;sum2=0;
for(j=0;j<length;j++){sum+=pow(schis[j]-value2,2)/sW[j];sum2+=pow(sW[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=pow(schis[j2]-value2,2)/sW[j2];sum2-=pow(sW[j2],-1);}
likes[4]=-.5*sum2*(1+log(2*M_PI*sum/sum2));

sum=0;sum2=0;
for(j=0;j<length;j++){sum+=pow(schis[j]-exps[j],2)/sW[j];sum2+=pow(sW[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=pow(schis[j2]-exps[j2],2)/sW[j2];sum2-=pow(sW[j2],-1);}
likes[5]=-.5*sum2*(1+log(2*M_PI*sum/sum2));

//save cflag
likes[6]=cflag;
}

//get negative counts
count=0;count2=0;
sum=0;sum2=0;
for(j=0;j<length;j++)
{
value=0;for(q=0;q<num_parts;q++){value+=svars[q][j]*thetas[q];}
if(value<0){count++;sum+=pow(stags[j],-1);}
value2=1+snss[j]/scale*value;
if(gcon==1){value2+=thetas[num_parts];}
if(cept==1){value2+=snss[j]*thetas[num_parts+gcon]/scale;}
if(value2<=0){count2++;sum2+=pow(stags[j],-1);}
}

if(count>10000||count2>10)
{printf("Warning, %d (%d) of the %d predictors have negative expected heritability (test statistic); this suggests an over-complicated heritability model\n\n", count, count2, length);}

if(likes!=NULL)	//save them
{
likes[7]=count;likes[8]=sum;
likes[9]=count2;likes[10]=sum2;
}

if(ncv>0)	//save expectations for cv predictors
{
for(j=0;j<ncv;j++){cvexps[j]=exps[cvindex[j]];}
}

////////

if(ncv==0&&(sflag==0||sflag==2||sflag==4||sflag==5||sflag==6))	//get SEs
{
if(num_blocks==-9999)	//get from second derivative of likelihood (still valid to use thetas)
{
//get inverse AI for current state
for(j=0;j<length;j++)
{
value=(schis[j]-.5*exps[j])/stags[j]*pow(exps[j],-3);
for(q=0;q<total;q++){sT2[j+q*length]=value*sTb[j+q*length];}
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sTb, &length, sT2, &length, &beta, AI, &total);
(void)eigen_invert(AI, total, AI2, -1, AI3, 1);

//AI provides variances for thetas = (cb*scale, (c-1), ca*scale) - for transformed variances, must compute J invAI JT, where Jij=dnewi/dthetaj

//make sure we have gc
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}

//first want variances of (hers, c, na) = (b ssums, c, a*scale) = (t1*ssums/(t2+1)/scale, t2+1, t3/(t2+1))
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){J[q+q2*total]=0;}
}
for(q=0;q<num_parts;q++){J[q+q*total]=ssums[q][q]/gc/scale;}
if(gcon==1)
{
for(q=0;q<num_parts;q++){J[q+num_parts*total]=-thetas[q]*ssums[q][q]/scale*pow(gc,-2);}
J[num_parts+num_parts*total]=1;
}
if(cept==1)
{
if(gcon==1){J[num_parts+gcon+num_parts*total]=-thetas[num_parts+gcon]*pow(gc,-2);}
J[num_parts+gcon+(num_parts+gcon)*total]=1.0/gc;
}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

for(q=0;q<total;q++)
{
if(JAIJT[q+q*total]>=0){stats[q+total+1+num_parts]=pow(JAIJT[q+q*total],.5);}
else{stats[q+total+1+num_parts]=-9999;}
}

if(cohers!=NULL)
{
for(q=0;q<num_parts;q++)
{
for(q2=0;q2<num_parts;q2++){cohers[q+q2*num_parts]=JAIJT[q+q2*total];}
}
}

//can also get variance of sumhers
sum=0;
for(q=0;q<num_parts;q++)
{
for(q2=0;q2<num_parts;q2++){sum+=JAIJT[q+q2*total];}
}
if(sum>=0){stats[total+total+1+num_parts]=pow(sum,.5);}
else{stats[total+total+1+num_parts]=-9999;}

//and variances of category heritabilities
for(q=0;q<num_parts;q++)
{
sum=0;
for(q2=0;q2<num_parts;q2++)
{
for(q3=0;q3<num_parts;q3++)
{sum+=JAIJT[q2+q3*total]*ssums[q][q2]/ssums[q2][q2]*ssums[q][q3]/ssums[q3][q3];}
}
if(sum>=0){stats[total+1+q+total+1+num_parts]=pow(sum,.5);}
else{stats[total+1+q+total+1+num_parts]=-9999;}
}

if(cohers2!=NULL)	//compute variance of (category - exp x total), where exp = ssums[q][np+1]
{
for(q=0;q<num_parts;q++)
{
//get cov(category, total) = cov( sum her_q2 s_q_q2/s_q2_q2, sum her_q3)
sum=0;
for(q2=0;q2<num_parts;q2++)
{
for(q3=0;q3<num_parts;q3++)
{sum+=JAIJT[q2+q3*total]*ssums[q][q2]/ssums[q2][q2];}
}

//can compute variance if variance of category, total and sum>=0
if(stats[total+1+q+total+1+num_parts]&&stats[total+total+1+num_parts]>0&&sum>=0)
{
value=pow(stats[total+1+q+total+1+num_parts],2)-2*ssums[q][total+1]*sum+pow(ssums[q][total+1]*stats[total+total+1+num_parts],2);
if(value>=0){cohers2[q]=pow(value,.5);}
else{cohers2[q]=-9999;}
}
else{cohers2[q]=-9999;}
}
}

//next want variances of (shares, c, ca) = (b ssums / sum b ssums, c/avsum, ca) - last two irrelevant
sum=0;for(q=0;q<num_parts;q++){sum+=thetas[q]*ssums[q][q];}
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){J[q+q2*total]=0;}
}
for(q=0;q<num_parts;q++)
{
for(q2=0;q2<num_parts;q2++){J[q+q2*total]=-thetas[q]*ssums[q][q]*ssums[q2][q2]*pow(sum,-2);}
J[q+q*total]+=ssums[q][q]/sum;
}
if(gcon==1){J[num_parts+num_parts*total]=1;}
if(cept==1){J[num_parts+gcon+(num_parts+gcon)*total]=1;}

alpha=1.0;beta=0.0;
dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total, &beta, JAI, &total);
dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total, &beta, JAIJT, &total);

for(q=0;q<total;q++)
{
if(JAIJT[q+q*total]>=0){stats[q+2*(total+1+num_parts)]=pow(JAIJT[q+q*total],.5);}
else{stats[q+2*(total+1+num_parts)]=-9999;}
}

//and variances of enrichments
for(q=0;q<num_parts;q++)
{
sum=0;
for(q2=0;q2<num_parts;q2++)
{
for(q3=0;q3<num_parts;q3++)
{sum+=JAIJT[q2+q3*total]*ssums[q][q2]/ssums[q2][q2]*ssums[q][q3]/ssums[q3][q3];}
}
if(sum>=0){stats[total+1+q+2*(total+1+num_parts)]=pow(sum,.5);}
else{stats[total+1+q+2*(total+1+num_parts)]=-9999;}
}
}
else	//jackknifing, using non-iterative least squares
{
//save sXTX and sXTY based on final expectations
for(j=0;j<length;j++)
{
sW[j]=stags[j]*pow(exps[j],2);
for(q=0;q<num_parts;q++){sX[j+q*length]=snss[j]/scale*svars[q][j]*pow(sW[j],-0.5);}
if(gcon==1){sX[j+num_parts*length]=pow(sW[j],-0.5);}
if(cept==1){sX[j+(num_parts+gcon)*length]=snss[j]/scale*pow(sW[j],-0.5);}
sY[j]=(schis[j]-1)*pow(sW[j],-0.5);
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sX, &length, sX, &length, &beta, sXTX, &total);
dgemv_("T", &length, &total, &alpha, sX, &length, sY, &one, &beta, sXTY, &one);

for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){sXTXs[q+q2*total]=sXTX[q+q2*total];}
sXTYs[q]=sXTY[q];
}

for(p=0;p<num_blocks;p++)
{
if(p%100000==0&&sflag!=6){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
start=(double)(p)/num_blocks*length;
end=(double)(p+1)/num_blocks*length;

//reset sXTY and sXTY, then subtract and solve
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){sXTX[q+q2*total]=sXTXs[q+q2*total];}
sXTY[q]=sXTYs[q];
}

alpha=-1.0;beta=1.0;count=end-start;
dgemm_("T", "N", &total, &total, &count, &alpha, sX+start, &length, sX+start, &length, &beta, sXTX, &total);
dgemv_("T", &count, &total, &alpha, sX+start, &length, sY+start, &one, &beta, sXTY, &one);
for(q=0;q<total;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total, sXTX2, 1, thetas, 1);

//save
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
sumhers=0;for(q=0;q<num_parts;q++){sumhers+=thetas[q]/gc*ssums[q][q]/scale;}

mark=p*(total+1+num_parts);
for(q=0;q<num_parts;q++){jacks[q+mark]=thetas[q]/gc*ssums[q][q]/scale;}
if(gcon==1){jacks[num_parts+mark]=gc;}
if(cept==1){jacks[num_parts+gcon+mark]=thetas[num_parts+gcon]/gc;}
jacks[total+mark]=sumhers;
for(q=0;q<num_parts;q++)
{
jacks[total+1+q+mark]=0;
for(q2=0;q2<num_parts;q2++){jacks[total+1+q+mark]+=thetas[q2]/gc*ssums[q][q2]/scale;}
}
}	//end of p loop
if(sflag!=6){printf("\n");}

//get sds for all stats (as is)
for(q=0;q<total+1+num_parts;q++)
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+1+num_parts);
sum+=jacks[q+mark];sumsq+=pow(jacks[q+mark],2);
}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
stats[q+total+1+num_parts]=pow(var,.5);
}
 
//now sds after dividing by sum of hers (for annotation shares, considered instead dividing by base)
for(q=0;q<total+1+num_parts;q++)
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+1+num_parts);
sum+=jacks[q+mark]/jacks[total+mark];sumsq+=pow(jacks[q+mark]/jacks[total+mark],2);
}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
stats[q+2*(total+1+num_parts)]=pow(var,.5);
}

if(cohers!=NULL)	//now coheritabilities - will be duplication of variances
{
for(q=0;q<num_parts;q++)
{
for(q2=0;q2<num_parts;q2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+1+num_parts);
sum+=jacks[q+mark];sum2+=jacks[q2+mark];sumsq+=jacks[q+mark]*jacks[q2+mark];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
cohers[q+q2*num_parts]=var;
}}
}

if(cohers2!=NULL)	//compute variance of (category - exp x total), where exp = ssums[q][np+1]
{
//first need top corner of JAIJT (covars of hers) - probably stored in cohers, but cant be sure
for(q=0;q<num_parts;q++)
{
for(q2=0;q2<num_parts;q2++)
{
sum=0;sum2=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+1+num_parts);
sum+=jacks[q+mark];sum2+=jacks[q2+mark];sumsq+=jacks[q+mark]*jacks[q2+mark];
}
mean=sum/num_blocks;mean2=sum2/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-mean*mean2);
JAIJT[q+q2*total]=var;
}}

//now get cohers2
for(q=0;q<num_parts;q++)
{
//get cov(category, total) = cov( sum her_q2 s_q_q2/s_q2_q2, sum her_q3)
sum=0;
for(q2=0;q2<num_parts;q2++)
{
for(q3=0;q3<num_parts;q3++)
{sum+=JAIJT[q2+q3*total]*ssums[q][q2]/ssums[q2][q2];}
}

//can compute variance if variance of category, total and sum>=0
if(stats[total+1+q+total+1+num_parts]&&stats[total+total+1+num_parts]>0&&sum>=0)
{
value=pow(stats[total+1+q+total+1+num_parts],2)-2*ssums[q][total+1]*sum+pow(ssums[q][total+1]*stats[total+total+1+num_parts],2);
if(value>=0){cohers2[q]=pow(value,.5);}
else{cohers2[q]=-9999;}
}
else{cohers2[q]=-9999;}
}
}
}	//end of jackknifing
}	//end of getting SEs

if(influs!=NULL)	//get factors for scaling hers to influences (can no longer use thetas)
{
//first get sum (stat/c-1-nja)^2/stags, weighted sum sq of test stats under null (allowing for gc/cept)
gc=1;if(gcon==1){gc=stats[num_parts];}
value=0;if(cept==1){value=stats[num_parts+gcon]/scale;}
sumsq=0;for(j=0;j<length;j++){sumsq+=pow(schis[j]/gc-1-snss[j]*value,2)/stags[j];}

for(q=0;q<num_parts;q++)	//get sum (stat/c-1-nja)(njvj)/stags /sumsq /ssums[q][q]
{
sum=0;for(j=0;j<length;j++){sum+=(schis[j]/gc-1-snss[j]*value)*snss[j]*svars[q][j]/stags[j];}
influs[q]=sum/sumsq/ssums[q][q];
}
}

free(thetas);free(thetadiffs);free(exps);free(exps2);
if(num_blocks!=-9999){free(jacks);}
free(sW);free(sX);free(sY);free(sXTX);free(sXTX2);free(sXTY);free(sXTXs);free(sXTYs);
free(sT);free(sTb);free(sT2);free(AI);free(AI2);free(AI3);free(BI);free(J);free(JAI);free(JAIJT);
}	//end of solve_sums

//////////////////////////

void solve_cors(double *stats, int num_parts, int gcon, int cept, int num_blocks, int length, double *stags, double **svars, double **ssums, double *snss, double *schis, double *srhos, double *snss2, double *schis2, double *srhos2, double tol, int maxiter, char *filename)
{
int j, q, q2, p, count, start, end, mark, one=1;
double value, sum, sum2, sum3, sumsq, mean, var, alpha, beta;

int total, total2, total3, total4;
double scale, scale2, scale3, gc, gc2, gc3, sumold, diff, sumhers, sumhers2, sumhers3, *snss3, *schis3, *exps, *exps2, *exps3;
double *sW, *sX, *sY, *sXTX, *sXTX2, *sXTY, *sXTXs, *sXTYs, *thetas, *jacks;

FILE *output;


//set totals and maybe num_blocks
total=2*(num_parts+gcon+cept)+num_parts+1;
total2=num_parts+gcon+cept;
total3=num_parts+1;
total4=num_parts+2;
if(num_blocks==-1||num_blocks>length){num_blocks=length;}

//allocate variables

snss3=malloc(sizeof(double)*length);
schis3=malloc(sizeof(double)*length);
exps=malloc(sizeof(double)*length);
exps2=malloc(sizeof(double)*length);
exps3=malloc(sizeof(double)*length);

sW=malloc(sizeof(double)*length);
sX=malloc(sizeof(double)*length*total4);
sY=malloc(sizeof(double)*length);
sXTX=malloc(sizeof(double)*total4*total4);
sXTX2=malloc(sizeof(double)*total4);
sXTY=malloc(sizeof(double)*total4);
sXTXs=malloc(sizeof(double)*total4*total4);
sXTYs=malloc(sizeof(double)*total4);
thetas=malloc(sizeof(double)*total4);
jacks=malloc(sizeof(double)*(total+4+num_parts)*num_blocks);

//set variables that stay the same throughout

//snss3 = root (snss*snss2)
for(j=0;j<length;j++){snss3[j]=pow(snss[j],.5)*pow(snss2[j],.5);}

//product of z-statistics
for(j=0;j<length;j++)
{
if(srhos[j]*srhos2[j]>0){schis3[j]=pow(schis[j],.5)*pow(schis2[j],.5);}
else{schis3[j]=-pow(schis[j],.5)*pow(schis2[j],.5);}
}

//get scale, scale2 and scale3 - weighted average sample sizes
sum=0;sum2=0;sum3=0;
for(j=0;j<length;j++){sum+=snss[j]/stags[j];sum2+=snss2[j]/stags[j];sum3+=pow(stags[j],-1);}
scale=sum/sum3;
scale2=sum2/sum3;
scale3=pow(scale,.5)*pow(scale2,.5);

//deal with progress file
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

////////

//estimate first total2 parameter
printf("Estimating heritabilies for Trait 1\n");
printf("Iter\tHer_All\t");
if(gcon==1){printf("Scaling\t");}
if(cept==1){printf("Intercept\t");}
printf("Difference\tTarget\n");
printf("Start\t0.0000\t");
if(gcon==1){printf("1.0000\t");}
if(cept==1){printf("1.0000\t");}
printf("n/a\t\t%.6f\n", tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Estimating heritabilies for Trait 1\n");
fprintf(output,"Iter\tHer_All\t");
if(gcon==1){fprintf(output, "Scaling\t");}
if(cept==1){fprintf(output, "Intercept\t");}
fprintf(output, "Difference\tTarget\n");
fclose(output);

//null model blank
sumold=0;
for(j=0;j<length;j++){exps[j]=1;}

count=0;
while(1)
{
//load up sW, sX and sY
for(j=0;j<length;j++)
{
sW[j]=stags[j]*pow(exps[j],2);
value=pow(sW[j],-.5);
for(q=0;q<num_parts;q++){sX[j+q*length]=snss[j]/scale*svars[q][j]*value;}
if(gcon==1){sX[j+num_parts*length]=value;}
if(cept==1){sX[j+(num_parts+gcon)*length]=snss[j]/scale*value;}
sY[j]=(schis[j]-1)*value;
}

//get sXTX and sXTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX, &length, &beta, sXTX, &total2);
dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTY, &one);

//solve for all parameters
for(q=0;q<total2;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total2, sXTX2, 1, thetas, 1);

//get diff
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
sumhers=0;for(q=0;q<num_parts;q++){sumhers+=thetas[q]/gc*ssums[q][q]/scale;}
diff=sumhers-sumold;
sumold=sumhers;

//print update - have just got gc and sumhers
printf("%d\t%.4f\t", count+1, sumhers);
if(gcon==1){printf("%.4f\t", gc);}
if(cept==1){printf("%.4f\t", 1+thetas[num_parts+gcon]/gc);}
printf("%.6f\t%.6f\n", diff, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "%d\t%.6f\t", count+1, sumhers);
if(gcon==1){fprintf(output, "%.6f\t", gc);}
if(cept==1){fprintf(output, "%.6f\t", 1+thetas[num_parts+gcon]/gc);}
fprintf(output, "%.6f\t%.6f\n", diff, tol);
fclose(output);

//update exps - can use sX, but remember to multiply by sW[j] and add on one
alpha=1.0;beta=0.0;
dgemv_("N", &length, &total2, &alpha, sX, &length, thetas, &one, &beta, exps, &one);
for(j=0;j<length;j++){exps[j]=1+exps[j]*pow(sW[j],.5);}
for(j=0;j<length;j++)
{
if(exps[j]<=0){exps[j]=1e-6;}
}

//see if breaking
if(fabs(diff)<tol){break;}
if(count==maxiter){printf("Warning, the optimizer failed to converge within %d iterations\n", maxiter);break;}

count++;
}

//put cat hers, sum of hers, gc and cept into first column of stats (sum goes right at end)
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
sumhers=0;for(q=0;q<num_parts;q++){sumhers+=thetas[q]/gc*ssums[q][q]/scale;}
//for(q=0;q<num_parts;q++){stats[q]=thetas[q]/gc*ssums[q][q]/scale;}
for(q=0;q<num_parts;q++)
{
stats[q]=0;
for(q2=0;q2<num_parts;q2++){stats[q]+=thetas[q2]/gc*ssums[q][q2]/scale;}
}
if(gcon==1){stats[num_parts]=gc;}
if(cept==1){stats[num_parts+gcon]=1+thetas[num_parts+gcon]/gc;}
stats[total]=sumhers;

//jackknife
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX, &length, &beta, sXTXs, &total2);
dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTYs, &one);

for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*length;
end=(double)(p+1)/num_blocks*length;
if(p%500000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
count=end-start;

//reset sXTX and sXTY
for(q=0;q<total2;q++)
{
for(q2=0;q2<total2;q2++){sXTX[q+q2*total2]=sXTXs[q+q2*total2];}
sXTY[q]=sXTYs[q];
}

//then subtract and solve
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &total2, &total2, &count, &alpha, sX+start, &length, sX+start, &length, &beta, sXTX, &total2);
dgemv_("T", &count, &total2, &alpha, sX+start, &length, sY+start, &one, &beta, sXTY, &one);
for(q=0;q<total2;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total2, sXTX2, 1, thetas, 1);

//put cat hers, sum of hers, gc and cept into jacks
mark=p*(total+4+num_parts);
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
sumhers=0;for(q=0;q<num_parts;q++){sumhers+=thetas[q]/gc*ssums[q][q]/scale;}
//for(q=0;q<num_parts;q++){jacks[q+mark]=thetas[q]/gc*ssums[q][q]/scale;}
for(q=0;q<num_parts;q++)
{
jacks[q+mark]=0;
for(q2=0;q2<num_parts;q2++){jacks[q+mark]+=thetas[q2]/gc*ssums[q][q2]/scale;}
}
if(gcon==1){jacks[num_parts+mark]=gc;}
if(cept==1){jacks[num_parts+gcon+mark]=1+thetas[num_parts+gcon]/gc;}
jacks[total+mark]=sumhers;
}
printf("\n");

////////

//estimate second total2 parameters
printf("Estimating heritabilies for Trait 2\n");
printf("Iter\tHer_All\t");
if(gcon==1){printf("Scaling\t");}
if(cept==1){printf("Intercept\t");}
printf("Difference\tTarget\n");
printf("Start\t0.0000\t");
if(gcon==1){printf("1.0000\t");}
if(cept==1){printf("1.0000\t");}
printf("n/a\t\t%.6f\n", tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Estimating heritabilies for Trait 2\n");
fprintf(output,"Iter\tHer_All\t");
if(gcon==1){fprintf(output, "Scaling\t");}
if(cept==1){fprintf(output, "Intercept\t");}
fprintf(output, "Difference\tTarget\n");
fclose(output);

//null model blank
sumold=0;
for(j=0;j<length;j++){exps2[j]=1;}

count=0;
while(1)
{
//load up sW, sX and sY
for(j=0;j<length;j++)
{
sW[j]=stags[j]*pow(exps2[j],2);
value=pow(sW[j],-.5);
for(q=0;q<num_parts;q++){sX[j+q*length]=snss2[j]/scale2*svars[q][j]*value;}
if(gcon==1){sX[j+num_parts*length]=value;}
if(cept==1){sX[j+(num_parts+gcon)*length]=snss2[j]/scale2*value;}
sY[j]=(schis2[j]-1)*value;
}

//get sXTX and sXTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX, &length, &beta, sXTX, &total2);
dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTY, &one);

//solve for all parameters
for(q=0;q<total2;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total2, sXTX2, 1, thetas, 1);

//get diff
gc2=1;if(gcon==1){gc2=thetas[num_parts]+1;}
sumhers2=0;for(q=0;q<num_parts;q++){sumhers2+=thetas[q]/gc2*ssums[q][q]/scale2;}
diff=sumhers2-sumold;
sumold=sumhers2;

//print update - have just got gc2 and sumhers2
printf("%d\t%.4f\t", count+1, sumhers2);
if(gcon==1){printf("%.4f\t", gc2);}
if(cept==1){printf("%.4f\t", 1+thetas[num_parts+gcon]/gc2);}
printf("%.6f\t%.6f\n", diff, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "%d\t%.6f\t", count+1, sumhers2);
if(gcon==1){fprintf(output, "%.6f\t", gc2);}
if(cept==1){fprintf(output, "%.6f\t", 1+thetas[num_parts+gcon]/gc2);}
fprintf(output, "%.6f\t%.6f\n", diff, tol);
fclose(output);

//update exps2 - can use sX, but remember to multiply by sW[j] and add on one
alpha=1.0;beta=0.0;
dgemv_("N", &length, &total2, &alpha, sX, &length, thetas, &one, &beta, exps2, &one);
for(j=0;j<length;j++){exps2[j]=1+exps2[j]*pow(sW[j],.5);}
for(j=0;j<length;j++)
{
if(exps2[j]<=0){exps2[j]=1e-6;}
}

//see if breaking
if(fabs(diff)<tol){break;}
if(count==maxiter){printf("Warning, the optimizer failed to converge within %d iterations\n", maxiter);break;}

count++;
}

//put cat hers, sum of hers, gc and cept into first column of stats (sum goes right at end)
gc2=1;if(gcon==1){gc2=thetas[num_parts]+1;}
sumhers2=0;for(q=0;q<num_parts;q++){sumhers2+=thetas[q]/gc2*ssums[q][q]/scale2;}
//for(q=0;q<num_parts;q++){stats[total2+q]=thetas[q]/gc2*ssums[q][q]/scale2;}
for(q=0;q<num_parts;q++)
{
stats[total2+q]=0;
for(q2=0;q2<num_parts;q2++){stats[total2+q]+=thetas[q2]/gc2*ssums[q][q2]/scale2;}
}
if(gcon==1){stats[total2+num_parts]=gc2;}
if(cept==1){stats[total2+num_parts+gcon]=1+thetas[num_parts+gcon]/gc2;}
stats[total+1]=sumhers2;

//jackknife
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX, &length, &beta, sXTXs, &total2);
dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTYs, &one);

for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*length;
end=(double)(p+1)/num_blocks*length;
if(p%500000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
count=end-start;

//reset sXTX and sXTY
for(q=0;q<total2;q++)
{
for(q2=0;q2<total2;q2++){sXTX[q+q2*total2]=sXTXs[q+q2*total2];}
sXTY[q]=sXTYs[q];
}

//then subtract and solve
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &total2, &total2, &count, &alpha, sX+start, &length, sX+start, &length, &beta, sXTX, &total2);
dgemv_("T", &count, &total2, &alpha, sX+start, &length, sY+start, &one, &beta, sXTY, &one);
for(q=0;q<total2;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total2, sXTX2, 1, thetas, 1);

//put cat hers, sum of hers, gc and cept into jacks
mark=p*(total+4+num_parts);
gc2=1;if(gcon==1){gc2=thetas[num_parts]+1;}
sumhers2=0;for(q=0;q<num_parts;q++){sumhers2+=thetas[q]/gc2*ssums[q][q]/scale2;}
//for(q=0;q<num_parts;q++){jacks[total2+q+mark]=thetas[q]/gc2*ssums[q][q]/scale2;}
for(q=0;q<num_parts;q++)
{
jacks[total2+q+mark]=0;
for(q2=0;q2<num_parts;q2++){jacks[total2+q+mark]+=thetas[q2]/gc2*ssums[q][q2]/scale2;}
}
if(gcon==1){jacks[total2+num_parts+mark]=gc2;}
if(cept==1){jacks[total2+num_parts+gcon+mark]=1+thetas[num_parts+gcon]/gc2;}
jacks[total+1+mark]=sumhers2;
}
printf("\n");

////////

//estimate final total3 parameters
printf("Estimating coheritabilies\n");
printf("Iter\tCoh_All\t");
printf("Difference\tTarget\n");
printf("Start\t0.0000\t");
printf("n/a\t\t%.6f\n", tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output,"Estimating coheritabilies\n");
fprintf(output,"Iter\tCoher_All\t");
fprintf(output, "Difference\tTarget\n");
fclose(output);

//null model blank
sumold=0;
for(j=0;j<length;j++){exps3[j]=0;}

count=0;
while(1)
{
//load up sW, sX and sY
for(j=0;j<length;j++)
{
sW[j]=stags[j]*(exps[j]*exps2[j]+pow(exps3[j],2));
if(sW[j]<=0){sW[j]=1e-6;}
for(q=0;q<num_parts;q++){sX[j+q*length]=snss3[j]/scale3*svars[q][j]*pow(sW[j],-.5);}
sX[j+num_parts*length]=pow(sW[j],-.5);
sY[j]=schis3[j]*pow(sW[j],-.5);
}

//get sXTX and sXTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total3, &total3, &length, &alpha, sX, &length, sX, &length, &beta, sXTX, &total3);
dgemv_("T", &length, &total3, &alpha, sX, &length, sY, &one, &beta, sXTY, &one);

//solve for all parameters
for(q=0;q<total3;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total3, sXTX2, 1, thetas, 1);

//get diff
gc=1;if(gcon==1){gc=stats[num_parts];}
gc2=1;if(gcon==1){gc2=stats[total2+num_parts];}
sumhers3=0;for(q=0;q<num_parts;q++){sumhers3+=thetas[q]*pow(gc*gc2,-.5)*ssums[q][q]/scale3;}
diff=sumhers3-sumold;
sumold=sumhers3;

//print update - have just got sumhers3
printf("%d\t%.4f\t", count+1, sumhers3);
printf("%.6f\t%.6f\n", diff, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "%d\t%.6f\t", count+1, sumhers3);
fprintf(output, "%.6f\t%.6f\n", diff, tol);
fclose(output);

//update exps3 - can use sX, but remember to multiply by sW[j] (but no need to add on one)
alpha=1.0;beta=0.0;
dgemv_("N", &length, &total3, &alpha, sX, &length, thetas, &one, &beta, exps3, &one);
for(j=0;j<length;j++){exps3[j]=exps3[j]*pow(sW[j],.5);}

//see if breaking
if(fabs(diff)<tol){break;}
if(count==maxiter){printf("Warning, the optimizer failed to converge within %d iterations\n", maxiter);break;}

count++;
}

//put cat cohers, sum of cohers and the weird term into first column of stats (sum goes right at end)
gc=1;if(gcon==1){gc=stats[num_parts];}
gc2=1;if(gcon==1){gc2=stats[total2+num_parts];}
gc3=pow(gc*gc2,.5);
sumhers3=0;for(q=0;q<num_parts;q++){sumhers3+=thetas[q]/gc3*ssums[q][q]/scale3;}
//for(q=0;q<num_parts;q++){stats[2*total2+q]=thetas[q]/gc3*ssums[q][q]/scale3;}
for(q=0;q<num_parts;q++)
{
stats[2*total2+q]=0;
for(q2=0;q2<num_parts;q2++){stats[2*total2+q]+=thetas[q2]/gc3*ssums[q][q2]/scale3;}
}
stats[2*total2+num_parts]=thetas[num_parts]/gc3;
stats[total+2]=sumhers3;

//and now also correlations
stats[total+3]=stats[total+2]*pow(stats[total]*stats[total+1],-.5);
for(q=0;q<num_parts;q++){stats[total+4+q]=stats[2*total2+q]*pow(stats[q]*stats[total2+q],-.5);}

//jackknife
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total3, &total3, &length, &alpha, sX, &length, sX, &length, &beta, sXTXs, &total3);
dgemv_("T", &length, &total3, &alpha, sX, &length, sY, &one, &beta, sXTYs, &one);

for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*length;
end=(double)(p+1)/num_blocks*length;
if(p%500000==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}
count=end-start;

//reset sXTX and sXTY
for(q=0;q<total3;q++)
{
for(q2=0;q2<total3;q2++){sXTX[q+q2*total3]=sXTXs[q+q2*total3];}
sXTY[q]=sXTYs[q];
}

//then subtract and solve
alpha=-1.0;beta=1.0;
dgemm_("T", "N", &total3, &total3, &count, &alpha, sX+start, &length, sX+start, &length, &beta, sXTX, &total3);
dgemv_("T", &count, &total3, &alpha, sX+start, &length, sY+start, &one, &beta, sXTY, &one);
for(q=0;q<total3;q++){thetas[q]=sXTY[q];}
(void)eigen_invert(sXTX, total3, sXTX2, 1, thetas, 1);

//put cat cohers, sum of cohers and the weird term into jacks
mark=p*(total+4+num_parts);
gc=1;if(gcon==1){gc=jacks[num_parts+mark];}
gc2=1;if(gcon==1){gc2=jacks[total2+num_parts+mark];}
gc3=pow(gc*gc2,.5);
sumhers3=0;for(q=0;q<num_parts;q++){sumhers3+=thetas[q]/gc3*ssums[q][q]/scale3;}
//for(q=0;q<num_parts;q++){jacks[2*total2+q+mark]=thetas[q]/gc3*ssums[q][q]/scale3;}
for(q=0;q<num_parts;q++)
{
jacks[2*total2+q+mark]=0;
for(q2=0;q2<num_parts;q2++){jacks[2*total2+q+mark]+=thetas[q2]/gc3*ssums[q][q2]/scale3;}
}
jacks[2*total2+num_parts+mark]=thetas[num_parts]/gc3;
jacks[total+2+mark]=sumhers3;

//and now also correlations
mark=p*(total+4+num_parts);
jacks[total+3+mark]=jacks[total+2+mark]*pow(jacks[total+mark]*jacks[total+1+mark],-.5);
for(q=0;q<num_parts;q++){jacks[total+4+q+mark]=jacks[2*total2+q+mark]*pow(jacks[q+mark]*jacks[total2+q+mark],-.5);}
}
printf("\n");

////////

//now get all sds
for(q=0;q<total+4+num_parts;q++)
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+4+num_parts);
sum+=jacks[q+mark];
sumsq+=pow(jacks[q+mark],2);
}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
stats[q+total+4+num_parts]=pow(var,.5);
}

printf("Trait 1 heritability: %.4f (%.4f)\n", stats[total], stats[total+total+4+num_parts]);
printf("Trait 2 heritability: %.4f (%.4f)\n", stats[total+1], stats[total+1+total+4+num_parts]);
printf("Coheritability: %.4f (%.4f)\n", stats[total+2], stats[total+2+total+4+num_parts]);
printf("Correlation: %.4f (%.4f)\n\n", stats[total+3], stats[total+3+total+4+num_parts]);

if(cept==1)	//get stats for ratio of intercepts
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+4+num_parts);
sum+=jacks[num_parts+gcon+mark]/jacks[total2+num_parts+gcon+mark];
sumsq+=pow(jacks[num_parts+gcon+mark]/jacks[total2+num_parts+gcon+mark],2);
}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Ratio\t%.6f\t%.6f\n", stats[num_parts+gcon]/stats[total2+num_parts+gcon], pow(var,.5));
fclose(output);
}

free(snss3);free(schis3);free(exps);free(exps2);free(exps3);
free(sW);free(sX);free(sY);free(sXTX);free(sXTX2);free(sXTY);free(sXTXs);free(sXTYs);free(thetas);free(jacks);
}	//end of solve_cors

//////////////////////////

/*
//code for log least squares

//first divide non-gc thetas by gc / replace gc theta with log
gc=1;if(gcon==1){gc=thetas[num_parts]+1;}
for(q=0;q<num_parts;q++){thetas[q]=thetas[q]/gc;}
if(gcon==1){thetas[num_parts]=log(gc);}
if(cept==1){thetas[num_parts+1]=thetas[num_parts+1]/gc;}

//load up sW, sX and sY - regressing logSj-dg-log(2)-log(expj) on nj*vjk/expj, 1, nj/expsj
for(j=0;j<length;j++)
{
sW[j]=stags[j];
for(q=0;q<num_parts;q++){sX[j+q*length]=snss[j]*svars[q][j]/exps[j]*pow(sW[j],-0.5);}
if(gcon==1){sX[j+num_parts*length]=pow(sW[j],-0.5);}
if(cept==1){sX[j+(num_parts+gcon)*length]=snss[j]/exps[j]*pow(sW[j],-0.5);}
sY[j]=(log(schis[j])-digamma(.5)-log(2)-log(exps[j]))*pow(sW[j],-0.5);
}

if(ncv>0)	//must blank out sX for cv predictors
{
for(j=0;j<ncv;j++)
{
for(q=0;q<token;q++){sX[cvindex[j]+q*length]=0;}
}
}

//solve to update theta - except for gc theta, these should be added
alpha=1.0;beta=0.0;
dgemm_("T", "N", &token, &token, &length, &alpha, sX, &length, sX, &length, &beta, sXTX, &token);
dgemv_("T", &length, &token, &alpha, sX, &length, sY, &one, &beta, sXTY, &one);
for(q=0;q<token;q++){thetas2[q]=sXTY[q];}
(void)eigen_invert(sXTX, token, sXTX2, 1, thetas2, 1);

//except for gc theta, we have found differences
for(q=0;q<num_parts;q++){thetas[q]+=thetas2[q];}
if(gcon==1){thetas[num_parts]=thetas[num_parts];}
if(cept==1){thetas[num_parts+gcon]+=thetas2[num_parts+gcon];}

//now multiply non-gc thetas by gc / replace gc theta with gc-1
gc=1;if(gcon==1){gc=exp(thetas[num_parts]);}
for(q=0;q<num_parts;q++){thetas[q]=thetas[q]*gc;}
if(gcon==1){thetas[num_parts]=gc-1;}
if(cept==1){thetas[num_parts+1]=thetas[num_parts+1]*gc;}
*/

//////////////////////////

/*
//tried to make mle solver for correlations - but did not work

void compute_exps(double *exps, double *exps2, double *exps3, double *sa, double *sb, double *sc, double *sd, double *se, double *sf, double *sZ, double *sZ2, double *sT, double *sT2, double *sT3, double *thetas, int num_parts, int length, int gcon, int cept)
{
int j, q, total2, total3, one=1;
double value, alpha, beta;

total2=num_parts+gcon+cept;
total3=num_parts+1;


//exps
for(j=0;j<length;j++){exps[j]=1;}
alpha=1.0;beta=1.0;
dgemv_("N", &length, &total2, &alpha, sT, &length, thetas, &one, &beta, exps, &one);
for(j=0;j<length;j++)
{
if(exps[j]<=0){exps[j]=1e-6;}
}

//exps2
for(j=0;j<length;j++){exps2[j]=1;}
alpha=1.0;beta=1.0;
dgemv_("N", &length, &total2, &alpha, sT2, &length, thetas+total2, &one, &beta, exps2, &one);
for(j=0;j<length;j++)
{
if(exps2[j]<=0){exps2[j]=1e-6;}
}

//exps3
alpha=1.0;beta=0.0;
dgemv_("N", &length, &total3, &alpha, sT3, &length, thetas+2*total2, &one, &beta, exps3, &one);

//as, bs, cs, ds, es, fs
sa[j]=exps2[j]/value;
sb[j]=-exps3[j]/value;
sc[j]=exps[j]/value;
sd[j]=value;
se[j]=sa[j]*sZ[j]+sb[j]*sZ2[j];
sf[j]=sb[j]*sZ[j]+sc[j]*sZ2[j];
}

}	//end of compute_exps

////////

void compute_derivs(double *AI, double *AI2, double *AI3, double *BI, double *sa, double *sb, double *sc, double *se, double *sf, double *sT, double *sT2, double *sT3, double *sT4, double *thetas, double *stags, int num_parts, int length, int gcon, int cept, int cflag)
{
int j, j2, q, q2, total, total2, total3, one=1;
double value, sum, alpha, beta;


total=num_parts+gcon+cept+num_parts+gcon+cept+num_parts+1;
total2=num_parts+gcon+cept;
total3=num_parts+1;

//get first derivatives
for(q=0;q<total;q++){BI[q]=0;}

if(cflag==0||cflag==2)
{
//first total2 thetas - multiply Ts by .5 (e^2 - a) / stags
for(q=0;q<total2;q++)
{
BI[q]=0;for(j=0;j<length;j++){BI[q]+=.5*(pow(se[j],2)-sa[j])/stags[j]*sT[j+q*length];}
}

//next total2 thetas - multiply Ts by .5 (f^2 - c) / stags
for(q=0;q<total2;q++)
{
BI[total2+q]=0;for(j=0;j<length;j++){BI[total2+q]+=.5*(pow(sf[j],2)-sc[j])/stags[j]*sT2[j+q*length];}
}
}

if(cflag==1||cflag==2)
{
//last total3 thetas - multiply Ts by (ef - b) / stags
for(q=0;q<total3;q++)
{
BI[2*total2+q]=0;for(j=0;j<length;j++){BI[2*total2+q]+=(se[j]*sf[j]-sb[j])/stags[j]*sT3[j+q*length];}
}
}

//now negative second derivatives
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){AI[q+q2*total]=0;}
AI[q+q*total]=1;
}

if(cflag==0||cflag==2)
{
//block 1,1 - multiply by (ae^2 - a^2/2)
for(j=0;j<length;j++)
{
value=(sa[j]*pow(se[j],2)-.5*pow(sa[j],2))/stags[j];
for(q=0;q<total2;q++){sT4[j+q*length]=value*sT[j+q*length];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sT, &length, sT4, &length, &beta, AI, &total);

//block 2,2 - multiply by (cf^2 - c^2/2)
for(j=0;j<length;j++)
{
value=(sc[j]*pow(sf[j],2)-.5*pow(sc[j],2))/stags[j];
for(q=0;q<total2;q++){sT4[j+q*length]=value*sT2[j+q*length];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sT2, &length, sT4, &length, &beta, AI+total2+total2*total, &total);

//blocks 1,2 and 2,1 - multiply by (bef - b^2/2)
for(j=0;j<length;j++)
{
value=(sb[j]*se[j]*sf[j]-.5*pow(sb[j],2))/stags[j];
for(q=0;q<total2;q++){sT4[j+q*length]=value*sT2[j+q*length];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total2, &length, &alpha, sT, &length, sT4, &length, &beta, AI+total2*total, &total);
for(j=0;j<total2;j++)
{
for(j2=total2;j2<2*total2;j2++){AI[j2+j*total]=AI[j+j2*total];}
}
}

if(cflag==1||cflag==2)
{
//block 3,3 - multiply by (af^2+2bef+ce^2 - b^2-ac)
for(j=0;j<length;j++)
{
value=(sa[j]*pow(sf[j],2)+2*sb[j]*se[j]*sf[j]+sc[j]*pow(se[j],2)-pow(sb[j],2)-sa[j]*sc[j])/stags[j];
for(q=0;q<total3;q++){sT4[j+q*length]=value*sT3[j+q*length];}
}
alpha=-1.0;beta=0.0;
dgemm_("T", "N", &total3, &total3, &length, &alpha, sT3, &length, sT4, &length, &beta, AI+2*total2+2*total2*total, &total);
}

if(cflag==2)
{
//blocks 1,3 and 3,1 - multiply by (aef+be2 - ab)
for(j=0;j<length;j++)
{
value=(sa[j]*se[j]*sf[j]+sb[j]*pow(se[j],2)-sa[j]*sb[j])/stags[j];
for(q=0;q<total3;q++){sT4[j+q*length]=value*sT3[j+q*length];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total3, &length, &alpha, sT, &length, sT4, &length, &beta, AI+2*total2*total, &total);
for(j=0;j<total2;j++)
{
for(j2=2*total2;j2<total;j2++){AI[j2+j*total]=AI[j+j2*total];}
}

//blocks 2,3 and 3,2 - multiply by (bf^2+cef - bc)
for(j=0;j<length;j++)
{
value=(sb[j]*pow(sf[j],2)+sc[j]*se[j]*sf[j]-sb[j]*sc[j])/stags[j];
for(q=0;q<total3;q++){sT4[j+q*length]=value*sT3[j+q*length];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total2, &total3, &length, &alpha, sT2, &length, sT4, &length, &beta, AI+total2+2*total2*total, &total);
for(j=total2;j<2*total2;j++)
{
for(j2=2*total2;j2<total;j2++){AI[j2+j*total]=AI[j+j2*total];}
}
}

//for stability, scale AI so has trace one
sum=0;for(j=0;j<total;j++){sum+=AI[j+j*total];}
for(j=0;j<total;j++)
{
for(j2=0;j2<total;j2++){AI[j+j2*total]*=pow(sum,-1);}
}

(void)eigen_invert(AI, total, AI2, -1, AI3, 1);

//undo scaling of AI
for(j=0;j<total;j++)
{
for(j2=0;j2<total;j2++){AI[j+j2*total]*=pow(sum,-1);}
}

for(q=0;q<total;q++)
{
//for(q2=0;q2<total;q2++){printf("%.3e ", AI[q+q2*total]);}
//printf(" >s %e\n", BI[q]);
}
}	//end of compute_derivs

////////

//(z,z2) is N(0,sigma), where sigma has components exps=exp(z^2), exps2=exp(z^2) and exps3=exp(zz2), 
//assume exp(z^2-1) = cnvjb + (c-1) + cna = [nvj 1 n] * (cb, c-1, ca) = sT theta
//assume exp(z2^2-1) = c2n2vjb2 + (c2-1) + c2n2a2 = [nvj 1 n] * (c2b2, c2-1, c2a2) = sT2 theta
//and exp(zz2) = (cc2)^.5 n3vjb3 + (cc2)^.5 frho = sT3 theta
//b, b2, b3 are the taus, n3 = (n n2)^.5, f is overlap fraction, rho is phenotypic correlation

//write inv sigma = (a b), let d be determinant of sigma, and inv sigma (z ) = (e)
//                  (b c)                                               (z2)   (f)

//set total, total2 and total3
total=num_parts+gcon+cept+num_parts+gcon+cept+num_parts+1;
total2=num_parts+gcon+cept;
total3=num_parts+1;

//allocate variables

snss3=malloc(sizeof(double)*length);
sZ=malloc(sizeof(double)*length);
sZ2=malloc(sizeof(double)*length);
sa=malloc(sizeof(double)*length);
sb=malloc(sizeof(double)*length);
sc=malloc(sizeof(double)*length);
sd=malloc(sizeof(double)*length);
se=malloc(sizeof(double)*length);
sf=malloc(sizeof(double)*length);

thetas=malloc(sizeof(double)*total);
thetadiffs=malloc(sizeof(double)*total);
exps=malloc(sizeof(double)*length);
exps2=malloc(sizeof(double)*length);
exps3=malloc(sizeof(double)*length);

sT=malloc(sizeof(double)*length*(num_parts+gcon+cept));
sT2=malloc(sizeof(double)*length*(num_parts+gcon+cept));
sT3=malloc(sizeof(double)*length*(num_parts+1));
sT4=malloc(sizeof(double)*length*(num_parts+2));
AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);

//snss3 = root (snss*snss2)
for(j=0;j<length;j++){snss3[j]=pow(snss[j],.5)*pow(snss2[j],.5);}

//z-statistics are signed roots of chis
for(j=0;j<length;j++)
{
if(srhos[j]>0){sZ[j]=pow(schis[j],.5);}
else{sZ[j]=-pow(schis[j],.5);}
if(srhos2[j]>0){sZ2[j]=pow(schis2[j],.5);}
else{sZ2[j]=-pow(schis2[j],.5);}
}

//get scale and scale2 - weighted average sample sizes
sum=0;sum2=0;sum3=0;
for(j=0;j<length;j++){sum+=snss[j]/stags[j];sum2+=snss2[j]/stags[j];sum3+=pow(stags[j],-1);}
scale=sum/sum3;
scale2=sum2/sum3;

//and avsum - average ssum (for scaling gcs and last term)
sum=0;for(q=0;q<num_parts;q++){sum+=ssums[q][num_parts];}
avsum=sum/num_parts;

//set sT, sT2 and sT3
for(j=0;j<length;j++)
{
for(q=0;q<num_parts;q++){sT[j+q*length]=snss[j]*svars[q][j];}
if(gcon==1){sT[j+num_parts*length]=avsum;}
if(cept==1){sT[j+(num_parts+gcon)*length]=snss[j];}

for(q=0;q<num_parts;q++){sT2[j+q*length]=snss2[j]*svars[q][j];}
if(gcon==1){sT2[j+num_parts*length]=1;}
if(cept==1){sT2[j+(num_parts+gcon)*length]=snss2[j];}

for(q=0;q<num_parts;q++){sT3[j+q*length]=snss3[j]*svars[q][j];}
sT3[j+num_parts*length]=avsum;
}

////////

//set starting model (remember, theta = b, c-1, ca, b2, c2-1, c2a2, (cc2)^.5b3, (cc2)^.5 frho)
for(q=0;q<total;q++){thetas[q]=0;}

//now iterate - rflag indicates type of move
count=0;
cflag=0;	//0 when doing individual traits, 1 when doing cross traits, 2 when doing everything
rflag=1;	//0 if using multiNR, 1 if about to do singleNR, 2 if just done singleNR
while(1)
{
//update expectations and likelihood
compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=(sZ[j]*se[j]+sZ2[j]*sf[j])/stags[j];sum2+=log(fabs(sd[j]))/stags[j];}
like=-.5*sum-.5*sum2;

if(count>0)	//set diff
{diff=like-likeold;}
likeold=like;

//see if switching cflag

if(count>0)
{
if(fabs(diff)<tol&&(rflag==0||cflag>0)){rflag=1;cflag++;}
if(cflag==3){break;}
}
if(count==maxiter){printf("\nWarning, the optimizer failed to converge within %d iterations\n", maxiter);break;}

////////

//update thetas using either single nr or multi nr (only move if good for likelihood)

if(rflag==1||cflag>0)	//single nr - for first iterations or if multi fails
{
if(cflag==0||cflag==2)
{
//first total2 thetas
for(q=0;q<total2;q++)
{
value=0;value2=0;
for(j=0;j<length;j++)
{
value+=.5*(pow(se[j],2)-sa[j])/stags[j]*sT[j+q*length];
value2+=(sa[j]*pow(se[j],2)-.5*pow(sa[j],2))/stags[j]*pow(sT[j+q*length],2);
}

//move, but return if reduces likelihood
thetadiffs[q]=value/value2;
thetas[q]+=thetadiffs[q];

compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=(sZ[j]*se[j]+sZ2[j]*sf[j])/stags[j];sum2+=log(fabs(sd[j]))/stags[j];}
like2=-.5*sum-.5*sum2;

if(like2>like-tol){like=like2;}
else
{
thetas[q]-=thetadiffs[q];
compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
}
}	//end of q loop

//second total2 thetas
for(q=0;q<total2;q++)
{
value=0;value2=0;
for(j=0;j<length;j++)
{
value+=.5*(pow(sf[j],2)-sc[j])/stags[j]*sT2[j+q*length];
value2+=(sc[j]*pow(sf[j],2)-.5*pow(sc[j],2))/stags[j]*pow(sT2[j+q*length],2);
}

//move, but return if reduces likelihood
thetadiffs[total2+q]=value/value2;
thetas[total2+q]+=thetadiffs[total2+q];

compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=(sZ[j]*se[j]+sZ2[j]*sf[j])/stags[j];sum2+=log(fabs(sd[j]))/stags[j];}
like2=-.5*sum-.5*sum2;

if(like2>like-tol){like=like2;}
else
{
thetas[total2+q]-=thetadiffs[total2+q];
compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
}
}	//end of q loop
}

if(cflag==1||cflag==2)
{
//last total3 thetas
//for(q=0;q<total3;q++)
for(q=total3-1;q>=0;q--)
{
value=0;value2=0;
for(j=0;j<length;j++)
{
value+=(se[j]*sf[j]-sb[j])/stags[j]*sT3[j+q*length];
value2+=(sa[j]*pow(sf[j],2)+2*sb[j]*se[j]*sf[j]+sc[j]*pow(se[j],2)-pow(sb[j],2)-sa[j]*sc[j])/stags[j]*pow(sT3[j+q*length],2);
}

//move, but return if reduces likelihood
thetadiffs[2*total2+q]=value/value2;
thetas[2*total2+q]+=thetadiffs[2*total2+q];

compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=(sZ[j]*se[j]+sZ2[j]*sf[j])/stags[j];sum2+=log(fabs(sd[j]))/stags[j];}
like2=-.5*sum-.5*sum2;

if(like2>like-tol){like=like2;}
else
{
thetas[2*total2+q]-=thetadiffs[2*total2+q];
compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
}
}	//end of q loop
}

rflag=2;
}
else	//multi nr
{
compute_derivs(AI, AI2, AI3, BI, sa, sb, sc, se, sf, sT, sT2, sT3, sT4, thetas, stags, num_parts, length, gcon, cept, cflag);

//get proposed move
alpha=1.0;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, thetadiffs, &one);

rflag=1;
relax=1;
while(relax>0.01)
{
//move relax*thetadiffs, but return if reduces likelihood
for(q=0;q<total;q++){thetas[q]+=relax*thetadiffs[q];}
compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=(sZ[j]*se[j]+sZ2[j]*sf[j])/stags[j];sum2+=log(fabs(sd[j]))/stags[j];}
like2=-.5*sum-.5*sum2;

if(like2>like-tol)	//accept move
{
like=like2;
rflag=0;
break;
}
else	//return, then next turn try smaller move
{
compute_exps(exps, exps2, exps3, sa, sb, sc, sd, se, sf, sZ, sZ2, sT, sT2, sT3, thetas, num_parts, length, gcon, cept);
for(q=0;q<total;q++){thetas[q]-=relax*thetadiffs[q];}
relax*=.5;
}
}
}

count++;
}	//end of while loop
printf("\n");

compute_derivs(AI, AI2, AI3, BI, sa, sb, sc, se, sf, sT, sT2, sT3, sT4, thetas, stags, num_parts, length, gcon, cept, 2);
*/

//////////////////////////


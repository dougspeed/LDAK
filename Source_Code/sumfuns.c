/*
Copyright 2026 Doug Speed.

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
start=(double)(p)*length/num_blocks;
end=(double)(p+1)*length/num_blocks;
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
while(relax>0.001)
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
while(relax>0.001)
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
//sflag=3 - just get expectations and likelihood, sflag=4 - LDSC, sflag=5 - divide+updating, sflag=6 - quiet (with starts)
//sflag=7,8 - quiet versions used by megaprs (8 uses starts)
//stats should have total+1+num_parts rows and three columns, where total=num_parts+gcon+cept
{
int j, j2, p, q, q2, q3, count, count2, start, end, mark, one=1;
double value, value2, sum, sum2, sumsq, mean, mean2, var, alpha, beta;

int total, cflag, rflag;
double scale, gc, sumhers, relax, likenull, like, like2, likeold, diff;
double *thetas, *thetadiffs, *exps, *exps2, *jacks;
double *sW, *sX, *sY, *sXTX, *sXTX2, *sXTY, *sXTXs, *sXTYs, *sT, *sT2, *AI, *AI2, *AI3, *BI, *J, *JAI, *JAIJT;

FILE *output;

//assuming statj = c * (1 + nja + njvj b), where c is gc, a is intercept/n, b is h2SNP/ssums
//easier to write as (statj-1) = nj/n vj cbn + (c-1) + nj/n can = sT theta, where n is average nj


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
sT=malloc(sizeof(double)*length*total);
sT2=malloc(sizeof(double)*length*total);
AI=malloc(sizeof(double)*total*total);
AI2=malloc(sizeof(double)*total);
AI3=malloc(sizeof(double)*total*total);
BI=malloc(sizeof(double)*total);
J=malloc(sizeof(double)*total*total);
JAI=malloc(sizeof(double)*total*total);
JAIJT=malloc(sizeof(double)*total*total);

if(chisol==0||num_blocks!=-9999)
{
sX=malloc(sizeof(double)*total*length);
sY=malloc(sizeof(double)*length);
sXTX=malloc(sizeof(double)*total*total);
sXTX2=malloc(sizeof(double)*total*total);
sXTY=malloc(sizeof(double)*total);
sXTXs=malloc(sizeof(double)*total*total);
sXTYs=malloc(sizeof(double)*total);
}

//set variables that stay the same throughout

//get scale - weighted average sample size
sum=0;sum2=0;
for(j=0;j<length;j++){sum+=snss[j]/stags[j];sum2+=pow(stags[j],-1);}
for(j=0;j<ncv;j++){j2=cvindex[j];sum-=snss[j2]/stags[j2];sum2-=pow(stags[j2],-1);}
scale=sum/sum2;

//set sT
for(j=0;j<length;j++)
{
for(q=0;q<num_parts;q++){sT[j+q*length]=snss[j]/scale*svars[q][j];}
if(gcon==1){sT[j+num_parts*length]=1;}
if(cept==1){sT[j+(num_parts+gcon)*length]=snss[j]/scale;}
}

//blank out sT for cv predictors
for(j=0;j<ncv;j++)
{
for(q=0;q<total;q++){sT[cvindex[j]+q*length]=0;}
}

////////

if(likes!=NULL) //get null likelihood
{
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
}

//set starting model

if(sflag==0||sflag==1||sflag==7)	//starting model assumes no causal variation (plus gcon=1 and cept=1)
{
for(q=0;q<total;q++){thetas[q]=0;}
}

if(sflag==2||sflag==3||sflag==5||sflag==6||sflag==8)	//use heritabilities (etc) saved in stats
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

if(sflag!=6&&sflag!=7&&sflag!=8)	//prepare to screen and file print
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

if(sflag!=6&&sflag!=7&&sflag!=8)	//print update
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

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "%d\t%.6f\t", count, sumhers);
if(gcon==1){fprintf(output, "%.6f\t", gc);}
if(cept==1){fprintf(output, "%.6f\t", 1+thetas[num_parts+gcon]/gc);}
if(count==0){fprintf(output, "%.6f\tNA\t%.6f\n", like, tol);}
else{fprintf(output, "%.6f\t%.6f\t%.6f\n", like, diff, tol);}
fclose(output);
}

//see if breaking (normally can only break if rflag=0, unless at iter limit or rflag=3)

if(sflag==3){break;}	//only wanted expectations and likelihood
if(sflag==4&&num_parts>1&&count==1){break;}	//multi-tagging ldsc uses a single iteration
if(rflag==3){cflag=0;break;}    //single nr did not update all components
if(count>0)
{
if(fabs(diff)<tol&&(rflag==0||chisol==0)){break;}
}
if(count==maxiter){printf("Warning, the optimizer failed to converge within %d iterations\n", maxiter);cflag=0;break;}

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
count2=0;
for(q=total-1;q>=0;q--)	//do backwards (might help if last category is base)
{
//get derivs for q
value=0;value2=0;
for(j=0;j<length;j++)
{
value+=.5*(schis[j]-exps[j])/stags[j]*pow(exps[j],-2)*sT[j+q*length];
value2+=(schis[j]-.5*exps[j])/stags[j]*pow(exps[j],-3)*pow(sT[j+q*length],2);
}

//get proposed move
thetadiffs[q]=value/value2;

relax=1;
while(relax>0.001)
{
//move relax*thetadiffs and get expectations
thetas[q]+=relax*thetadiffs[q];

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
count2++;
break;
}
else	//move back and next turn try smaller move
{
thetas[q]-=relax*thetadiffs[q];
relax*=.5;
}
}
}	//end of q loop

rflag=2;

//for sflag=6, will stop iterations if any components failed to update
if(sflag==6&&count2<total){rflag=3;break;}
}
else	//multi nr - use BI and AI
{
//load up sT2=(stat-exps/2)*w/exps^3*sT - will be blanked for cvpredictors
for(j=0;j<length;j++)
{
value=(schis[j]-.5*exps[j])/stags[j]*pow(exps[j],-3);
for(q=0;q<total;q++){sT2[j+q*length]=value*sT[j+q*length];}
}

//get AI (-2nd deriv) and BI (1st deriv) - cvpredictors already taken care of (by blanking sT and ST2)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sT, &length, sT2, &length, &beta, AI, &total);
for(q=0;q<total;q++)
{
BI[q]=0;for(j=0;j<length;j++){BI[q]+=.5*(schis[j]-exps[j])/stags[j]*pow(exps[j],-2)*sT[j+q*length];}
}

(void)eigen_invert(AI, total, AI2, -1, AI3, 1);

//get proposed move
alpha=1.0;beta=0.0;
dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta, thetadiffs, &one);

rflag=1;
relax=1;
while(relax>0.001)
{
//move relax*thetadiffs and get expectations
for(q=0;q<total;q++){
thetas[q]+=relax*thetadiffs[q];}

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
if(sflag!=6&&sflag!=7&&sflag!=8){printf("\n");}

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
{
if(sflag!=6&&sflag!=7&&sflag!=8)
{
if(num_parts==1){printf("Warning, %d (%d) of the %d predictors have negative expected heritability (test statistic)\n\n", count, count2, length);}
else{printf("Warning, %d (%d) of the %d predictors have negative expected heritability (test statistic); this suggests an over-complicated heritability model\n\n", count, count2, length);}
}
}

if(likes!=NULL)	//save them
{
likes[7]=count;likes[8]=sum;
likes[9]=count2;likes[10]=sum2;
}

if(ncv>0)	//compute expectations for cv predictors
{
//reset sT
for(j=0;j<length;j++)
{
for(q=0;q<num_parts;q++){sT[j+q*length]=snss[j]/scale*svars[q][j];}
if(gcon==1){sT[j+num_parts*length]=1;}
if(cept==1){sT[j+(num_parts+gcon)*length]=snss[j]/scale;}
}

//recompute exps
for(j=0;j<length;j++){exps[j]=1;}
alpha=1.0;beta=1.0;
dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta, exps, &one);
for(j=0;j<length;j++)
{
if(exps[j]<=0){exps[j]=1e-6;}
}

//extract exps for the cv predictors
for(j=0;j<ncv;j++){cvexps[j]=exps[cvindex[j]];}
}

////////

if(ncv==0&&(sflag==0||sflag==2||sflag==4||sflag==5||sflag==6))	//get SEs
{
if(num_blocks==-9999)	//get from second derivative of likelihood (still valid to use thetas)
{
if(cflag==0)    //failed to converge, so all SDs are NA
{
for(q=0;q<total+1+num_parts;q++){stats[q+total+1+num_parts]=-9999;stats[q+2*(total+1+num_parts)]=-9999;}
if(cohers2!=NULL)
{
for(q=0;q<num_parts;q++){cohers2[q]=-9999;}
}
}
else
{
//get inverse AI for current state
for(j=0;j<length;j++)
{
value=(schis[j]-.5*exps[j])/stags[j]*pow(exps[j],-3);
for(q=0;q<total;q++){sT2[j+q*length]=value*sT[j+q*length];}
}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &total, &total, &length, &alpha, sT, &length, sT2, &length, &beta, AI, &total);
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
//get cov(category, total) = cov( sum her_q3, sum her_q2 s_q_q2/s_q2_q2)
sum=0;
for(q2=0;q2<num_parts;q2++)
{
for(q3=0;q3<num_parts;q3++)
{sum+=JAIJT[q2+q3*total]*ssums[q][q2]/ssums[q2][q2];}
}

//can compute variance if variance of category and total >=0 (used to also require sum>=0)
if(stats[total+1+q+total+1+num_parts]&&stats[total+total+1+num_parts]>0)
{
value=pow(stats[total+1+q+total+1+num_parts],2)-2*ssums[q][num_parts+1]*sum+pow(ssums[q][num_parts+1]*stats[total+total+1+num_parts],2);
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
start=(double)(p)*length/num_blocks;
end=(double)(p+1)*length/num_blocks;
count=end-start;

//reset sXTY and sXTY, then subtract and solve
for(q=0;q<total;q++)
{
for(q2=0;q2<total;q2++){sXTX[q+q2*total]=sXTXs[q+q2*total];}
sXTY[q]=sXTYs[q];
}

alpha=-1.0;beta=1.0;
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

//can compute variance if variance of category and total >=0 (used to also require sum>=0)
if(stats[total+1+q+total+1+num_parts]&&stats[total+total+1+num_parts]>0&&sum>=0)
{
value=pow(stats[total+1+q+total+1+num_parts],2)-2*ssums[q][num_parts+1]*sum+pow(ssums[q][num_parts+1]*stats[total+total+1+num_parts],2);
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
free(sW);free(sT);free(sT2);free(AI);free(AI2);free(AI3);free(BI);free(J);free(JAI);free(JAIJT);
if(chisol==0||num_blocks!=-9999){free(sX);free(sY);free(sXTX);free(sXTX2);free(sXTY);free(sXTXs);free(sXTYs);}
}	//end of solve_sums

//////////////////////////

void solve_cors(double *stats, int num_parts, int gcon, int cept, int oversamp, int num_blocks, int length, double *stags, double **svars, double **ssums, double *snss, double *schis, double *srhos, double *snss2, double *schis2, double *srhos2, double tol, int maxiter, char *filename, int type, double *stats2)
//type=0 - regular tagging, type=1 - cross tagging (will move up svars2 for second and third parts)
//type=2 - quiet normal, type=3 - quiet cross tagging
//stats will contain total2(her1)+total2(her2)+total3(coher)+4(her1,her2,coher,cor)+np(cor cats) = total+4+np
{
int j, q, q2, p, count, start, end, mark, one=1;
double value, value2, value3, value4, sum, sum2, sum3, sumsq, mean, var, alpha, beta;
double prev, prev2, ascer, ascer2;

int total, total2, total3, total4;
double scale, scale2, scale3, gc, gc2, gc3, sumold, diff, sumhers, sumhers2, sumhers3, *snss3, *schis3, *exps, *exps2, *exps3;
double *sW, *sX, *sY, *sXTX, *sXTX2, *sXTY, *sXTXs, *sXTYs, *thetas, *jacks;

FILE *output;


count=0;for(j=0;j<length;j++){count+=(stags[j]==0);}
if(count>0){printf("Error AXXE4 %d, please tell Doug\n\n", count);exit(1);}

//set totals and maybe num_blocks
total2=num_parts+gcon+cept;	//number of components used when estimating hers
total3=num_parts+oversamp;	//number of components used when estimating coher
total=2*total2+total3;
total4=total2;
if(total3>total4){total4=total3;}
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

if(type!=2&&type!=3)	//deal with progress file
{
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);
}

////////

//estimate first total2 parameter
if(type!=2&&type!=3)
{
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
}

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

if(type!=2&&type!=3)	//print update - have just got gc and sumhers
{
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
}

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
start=(double)(p)*length/num_blocks;
end=(double)(p+1)*length/num_blocks;
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
if(type!=2&&type!=3){printf("\n");}

////////

if(type!=2&&type!=3)	//estimate second total2 parameters
{
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
}

//null model blank
sumold=0;
for(j=0;j<length;j++){exps2[j]=1;}

if(type==1||type==3)	//move second column of svars into first
{
for(j=0;j<length;j++){svars[0][j]=svars[1][j];}
}

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

if(type!=2&&type!=3)	//print update - have just got gc2 and sumhers2
{
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
}

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
start=(double)(p)*length/num_blocks;
end=(double)(p+1)*length/num_blocks;
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
if(type!=2&&type!=3){printf("\n");}

////////

//estimate final total3 parameters
if(type!=2&&type!=3)
{
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
}

//null model blank
sumold=0;
for(j=0;j<length;j++){exps3[j]=0;}

if(type==1||type==3)	//move third column of svars into first
{
for(j=0;j<length;j++){svars[0][j]=svars[2][j];}
}

count=0;
while(1)
{
//load up sW, sX and sY
for(j=0;j<length;j++)
{
sW[j]=stags[j]*(exps[j]*exps2[j]+pow(exps3[j],2));
if(sW[j]<=0){sW[j]=1e-6;}
for(q=0;q<num_parts;q++){sX[j+q*length]=snss3[j]/scale3*svars[q][j]*pow(sW[j],-.5);}
if(oversamp==1){sX[j+num_parts*length]=pow(sW[j],-.5);}
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

if(type!=2&&type!=3)	//print update - have just got sumhers3
{
printf("%d\t%.4f\t", count+1, sumhers3);
printf("%.6f\t%.6f\n", diff, tol);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "%d\t%.6f\t", count+1, sumhers3);
fprintf(output, "%.6f\t%.6f\n", diff, tol);
fclose(output);
}

//update exps3 - can use sX, but remember to multiply by sW[j] (but no need to add on one)
alpha=1.0;beta=0.0;
dgemv_("N", &length, &total3, &alpha, sX, &length, thetas, &one, &beta, exps3, &one);
for(j=0;j<length;j++){exps3[j]=exps3[j]*pow(sW[j],.5);}

//see if breaking
if(fabs(diff)<tol){break;}
if(count==maxiter){printf("Warning, the optimizer failed to converge within %d iterations\n", maxiter);break;}

count++;
}

//put cat cohers, sum of cohers and (maybe) the weird term into first column of stats (sum goes right at end)
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
if(oversamp==1){stats[2*total2+num_parts]=thetas[num_parts]/gc3;}
stats[total+2]=sumhers3;

//and now also correlations (overall and categories)
stats[total+3]=stats[total+2]*pow(stats[total]*stats[total+1],-.5);
for(q=0;q<num_parts;q++){stats[total+4+q]=stats[2*total2+q]*pow(stats[q]*stats[total2+q],-.5);}

//jackknife
alpha=1.0;beta=0.0;
dgemm_("T", "N", &total3, &total3, &length, &alpha, sX, &length, sX, &length, &beta, sXTXs, &total3);
dgemv_("T", &length, &total3, &alpha, sX, &length, sY, &one, &beta, sXTYs, &one);

for(p=0;p<num_blocks;p++)
{
start=(double)(p)*length/num_blocks;
end=(double)(p+1)*length/num_blocks;
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

//put cat cohers, sum of cohers and (maybe) the weird term into jacks
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
if(oversamp==1){jacks[2*total2+num_parts+mark]=thetas[num_parts]/gc3;}
jacks[total+2+mark]=sumhers3;

//and now also correlations
mark=p*(total+4+num_parts);
jacks[total+3+mark]=jacks[total+2+mark]*pow(jacks[total+mark]*jacks[total+1+mark],-.5);
for(q=0;q<num_parts;q++){jacks[total+4+q+mark]=jacks[2*total2+q+mark]*pow(jacks[q+mark]*jacks[total2+q+mark],-.5);}
}
if(type!=2&&type!=3){printf("\n");}

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

if(type!=2&&type!=3)
{
printf("Trait 1 heritability: %.4f (%.4f)\n", stats[total], stats[total+total+4+num_parts]);
printf("Trait 2 heritability: %.4f (%.4f)\n", stats[total+1], stats[total+1+total+4+num_parts]);
printf("Coheritability: %.4f (%.4f)\n", stats[total+2], stats[total+2+total+4+num_parts]);
printf("Correlation: %.4f (%.4f)\n\n", stats[total+3], stats[total+3+total+4+num_parts]);
}

if(cept==1)	//get stats for ratio of intercepts - cant remember why!
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

if(type!=2&&type!=3)
{
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n",filename);exit(1);}
fprintf(output, "Ratio\t%.6f\t%.6f\n", stats[num_parts+gcon]/stats[total2+num_parts+gcon], pow(var,.5));
fclose(output);
}
}

if(stats2!=NULL) //get stats for genetic separation 
{
prev=stats2[0], prev2=stats2[1], ascer=stats2[2], ascer2=stats2[3];

value=pow(1-prev,2)/ascer/(1-ascer);
value2=pow(1-prev2,2)/ascer2/(1-ascer2);
value3=(1-prev)*(1-prev2)*pow(ascer*(1-ascer)*ascer2*(1-ascer2),-.5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+4+num_parts);
value4=value*jacks[total+mark]+value2*jacks[total+1+mark]-2*value3*jacks[total+2+mark];
sum+=value4;
sumsq+=pow(value4,2);
}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));

stats2[0]=value*stats[total]+value2*stats[total+1]-2*value3*stats[total+2];
stats2[1]=pow(var,.5);
printf("Genetic_Separation: %.4f (%.4f)\n", stats2[0], stats2[1]);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{
mark=p*(total+4+num_parts);
value4=value*jacks[total+mark]+value2*jacks[total+1+mark]-2*value3*jacks[total+2+mark];
if(value4>0){sum+=erfc(-pow(value4/2,.5)*M_SQRT1_2)/2;sumsq+=pow(erfc(-pow(value4/2,.5)*M_SQRT1_2)/2,2);}
else{sum+=0.5;sum+=.25;}
}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));

value4=value*stats[total]+value2*stats[total+1]-2*value3*stats[total+2];
if(value4>0){stats2[2]=erfc(-pow(value4/2,.5)*M_SQRT1_2)/2;}
else{stats2[2]=0.5;}
stats2[3]=pow(var,.5);
printf("Max_AUC: %.4f (%.4f)\n\n", stats2[2], stats2[3]);
}

free(snss3);free(schis3);free(exps);free(exps2);free(exps3);
free(sW);free(sX);free(sY);free(sXTX);free(sXTX2);free(sXTY);free(sXTXs);free(sXTYs);free(thetas);free(jacks);
}	//end of solve_cors

//////////////////////////

void solve_sums_lite(double *sigma, double *omega, double *rjksums2, double *rjksums3, double *nss, double *chis, int length, int type, int q)
//type=0 - fix intercept to one; type=1 - estimate intercept; type=2 - fix intercept to current value
{
int j, count;
double value, sum, sum2, sumsq, mean, var, scale;

double *stags, *svars, *snss, *schis;
double sXTX[4], sXTX2[4], sXTY[2], est, est2, sd, sd2;


stags=malloc(sizeof(double)*length);
svars=malloc(sizeof(double)*length);
snss=malloc(sizeof(double)*length);
schis=malloc(sizeof(double)*length);

count=0;
for(j=0;j<length;j++)
{
if(nss[j]>0)	//have summaries
{
if(rjksums2[j]<=0){printf("Error ADD32, please tell doug %d %f\n", j+1, rjksums2[j]);exit(1);}

stags[count]=rjksums2[j];
svars[count]=rjksums3[j];
snss[count]=nss[j];
schis[count]=chis[j];
if(schis[count]>100){schis[count]=100;}
count++;
}
}

//get weighted average sample size
sum=0;sum2=0;
for(j=0;j<count;j++){sum+=snss[j]/stags[j];sum2+=pow(stags[j],-1);}
scale=sum/sum2;

if(type==1)    //estimate intercept
{
sXTX[0]=0;sXTX[1]=0;sXTX[3]=0;sXTY[0]=0;sXTY[1]=0;
sum=0;sumsq=0;
for(j=0;j<count;j++)
{
sXTX[0]+=pow(snss[j]/scale,2)/stags[j];
sXTX[1]+=pow(snss[j]/scale,2)*svars[j]/stags[j];
sXTX[3]+=pow(snss[j]/scale*svars[j],2)/stags[j];
sXTY[0]+=snss[j]/scale*(schis[j]-1)/stags[j];
sXTY[1]+=snss[j]/scale*svars[j]*(schis[j]-1)/stags[j];
sum+=(schis[j]-1)*pow(stags[j],-.5);
sumsq+=pow(schis[j]-1,2)/stags[j];
}
mean=sum/count;
var=sumsq/count-pow(mean,2);

//invert sXTX and solve
value=sXTX[0]*sXTX[3]-pow(sXTX[1],2);
sXTX2[0]=sXTX[3]/value;
sXTX2[1]=-sXTX[1]/value;
sXTX2[3]=sXTX[0]/value;
est=sXTX2[0]*sXTY[0]+sXTX2[1]*sXTY[1];
est2=(sXTX2[1]*sXTY[0]+sXTX2[3]*sXTY[1])/scale;

if(sXTX2[0]>0){sd=pow(var*sXTX2[0],.5);}
else{sd=100;}
if(sXTX2[3]>0){sd2=pow(var*sXTX2[3],.5)/scale;}
else{sd2=100;}

//printf("trait %d int is %f (%f), slope is %f (%f)\n", q+1, est, sd, est2, sd2);

if(sd<.2&&sd2<.2)   //estimates seem ok, will use
{
*sigma=1+est;
*omega=est2;
}
else    //estimates unreliable, will switch to type=0
{
printf("Warning, could not reliably estimate both the intercept and slope (suggesting there is limited LD or genetic signal); will fix the former to zero\n");
type=0;
}
}

if(type==0)    //fix intercept to one
{
sXTX[0]=0;sXTY[0]=0;
sum=0;sumsq=0;
for(j=0;j<count;j++)
{
sXTX[0]+=pow(snss[j]/scale*svars[j],2)/stags[j];
sXTY[0]+=snss[j]/scale*svars[j]*(schis[j]-1)/stags[j];
sum+=(schis[j]-1)*pow(stags[j],-.5);
sumsq+=pow(schis[j]-1,2)/stags[j];
}
mean=sum/count;
var=sumsq/count-pow(mean,2);

//invert sXTX and solve
sXTX2[0]=pow(sXTX[0],-1);
est=sXTX2[0]*sXTY[0]/scale;

if(sXTX2[0]>0){sd=pow(var*sXTX2[0],.5)/scale;}
else{sd=100;}

//printf("trait %d slope is %f (%f)\n", q+1, est, sd);

if(sd>.2){printf("Warning, could not reliably estimate the slope (suggesting there is limited LD or genetic signal),\n");}

*sigma=1;
*omega=est;
}

if(type==2)    //fix intercept to current value
{
value=(*sigma-1)/scale;
sXTX[0]=0;sXTY[0]=0;
sum=0;sumsq=0;
for(j=0;j<count;j++)
{
sXTX[0]+=pow(snss[j]/scale*svars[j],2)/stags[j];
sXTY[0]+=snss[j]/scale*svars[j]*(schis[j]-1-value*snss[j])/stags[j];
sum+=(schis[j]-1-value*snss[j])*pow(stags[j],-.5);
sumsq+=pow(schis[j]-1-value*snss[j],2)/stags[j];
}
mean=sum/count;
var=sumsq/count-pow(mean,2);

//invert sXTX and solve
sXTX2[0]=pow(sXTX[0],-1);
est=sXTX2[0]*sXTY[0]/scale;

if(sXTX2[0]>0){sd=pow(var*sXTX2[0],.5)/scale;}
else{sd=100;}

//if(sd>.2){printf("Warning, could not reliably estimate the slope (suggesting there is limited LD or genetic signal),\n");}

*omega=est;
}

free(stags);free(svars);free(snss);free(schis);
}

////////

void solve_cors_lite(double *sigma, double *omega, double *rjksums2a, double *rjksums2b, double *rjksums3, double *nss, double *nss2, double *chis, double *chis2, double *rhos, double *rhos2, int length, int type, int q, int q2)
//type=0 - fix intercept to one; type=1 - estimate intercept; type=2 - fix intercept to current value
{
int j, count;
double value, sum, sum2, sumsq, mean, var, scale;

double *stags, *svars, *snss, *sprods;
double sXTX[4], sXTX2[4], sXTY[2], est, est2, sd, sd2;


stags=malloc(sizeof(double)*length);
svars=malloc(sizeof(double)*length);
snss=malloc(sizeof(double)*length);
sprods=malloc(sizeof(double)*length);

count=0;
for(j=0;j<length;j++)
{
if(nss[j]>0&&nss2[j]>0)	//predictor in both correlations and have summaries
{
if(rjksums2a[j]<=0||rjksums2b[j]<=0){printf("Error BDD32, please tell doug %d %f %f\n", j+1, rjksums2a[j], rjksums2b[j]);exit(1);}

stags[count]=pow(rjksums2a[j]*rjksums2b[j],.5);
svars[count]=rjksums3[j];
snss[count]=pow(nss[j]*nss2[j],.5);
if(rhos[j]*rhos2[j]>0){sprods[count]=pow(chis[j]*chis2[j],.5);}
else{sprods[count]=-pow(chis[j]*chis2[j],.5);}
count++;
}
}

//get weighted average sample size
sum=0;sum2=0;
for(j=0;j<count;j++){sum+=snss[j]/stags[j];sum2+=pow(stags[j],-1);}
scale=sum/sum2;

if(type==1)    //estimate intercept
{
sXTX[0]=0;sXTX[1]=0;sXTX[3]=0;sXTY[0]=0;sXTY[1]=0;
sum=0;sumsq=0;
for(j=0;j<count;j++)
{
sXTX[0]+=1.0/stags[j];
sXTX[1]+=snss[j]/scale*svars[j]/stags[j];
sXTX[3]+=pow(snss[j]/scale*svars[j],2)/stags[j];
sXTY[0]+=sprods[j]/stags[j];
sXTY[1]+=snss[j]/scale*svars[j]*sprods[j]/stags[j];
sum+=sprods[j]*pow(stags[j],-.5);
sumsq+=pow(sprods[j],2)/stags[j];
}
mean=sum/count;
var=sumsq/count-pow(mean,2);

//invert sXTX and solve
value=sXTX[0]*sXTX[3]-pow(sXTX[1],2);
sXTX2[0]=sXTX[3]/value;
sXTX2[1]=-sXTX[1]/value;
sXTX2[3]=sXTX[0]/value;
est=sXTX2[0]*sXTY[0]+sXTX2[1]*sXTY[1];
est2=(sXTX2[1]*sXTY[0]+sXTX2[3]*sXTY[1])/scale;

if(sXTX2[0]>0){sd=pow(var*sXTX2[0],.5);}
else{sd=100;}
if(sXTX2[3]>0){sd2=pow(var*sXTX2[3],.5)/scale;}
else{sd2=100;}

//printf("traits %d and %d int is %f (%f), slope is %f (%f)\n", q+1, q2+1, est, sd, est2, sd2);

if(sd<.2&&sd2<.2)   //estimates seem ok, will use
{
*sigma=est;
*omega=est2;
}
else    //estimates unreliable, will switch to type=0
{
printf("Warning, could not reliably estimate both the intercept and slope (suggesting there is limited LD or genetic signal); will fix the former to zero\n");
type=0;
}
}

if(type==0)    //fix intercept to zero
{
sXTX[0]=0;sXTY[0]=0;
sum=0;sumsq=0;
for(j=0;j<count;j++)
{
sXTX[0]+=pow(snss[j]/scale*svars[j],2)/stags[j];
sXTY[0]+=snss[j]/scale*svars[j]*sprods[j]/stags[j];
sum+=sprods[j]*pow(stags[j],-.5);
sumsq+=pow(sprods[j],2)/stags[j];
}
mean=sum/count;
var=sumsq/count-pow(mean,2);

//invert sXTX and solve
sXTX2[0]=pow(sXTX[0],-1);
est=sXTX2[0]*sXTY[0]/scale;

if(sXTX2[0]>0){sd=pow(var*sXTX2[0],.5)/scale;}
else{sd=100;}

//printf("traits %d and %d slope is %f (%f)\n", q+1, q2+1, est, sd);

if(sd>.2){printf("Warning, could not reliably estimate the slope (suggesting there is limited LD or genetic signal),\n");}

*sigma=0;
*omega=est;
}

if(type==2)    //fix intercept to current value
{
value=*sigma;
sXTX[0]=0;sXTY[0]=0;
sum=0;sumsq=0;
for(j=0;j<count;j++)
{
sXTX[0]+=pow(snss[j]/scale*svars[j],2)/stags[j];
sXTY[0]+=snss[j]/scale*svars[j]*(sprods[j]-value)/stags[j];
sum+=(sprods[j]-value)*pow(stags[j],-.5);
sumsq+=pow(sprods[j]-value,2)/stags[j];
}
mean=sum/count;
var=sumsq/count-pow(mean,2);

//invert sXTX and solve
sXTX2[0]=pow(sXTX[0],-1);
est=sXTX2[0]*sXTY[0]/scale;

if(sXTX2[0]>0){sd=pow(var*sXTX2[0],.5)/scale;}
else{sd=100;}

//if(sd>.2){printf("Warning, could not reliably estimate the slope (suggesting there is limited LD or genetic signal),\n");}

*omega=est;
}

free(stags);free(svars);free(snss);free(sprods);
}

//////////////////////////


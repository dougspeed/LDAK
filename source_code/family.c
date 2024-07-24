/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Estimation of heritability for related pairs

///////////////////////////

void get_log_density(double *phen1, double *phen2, double *relats, double *relats2, int total, double *hers, double *logdens)
{
int i;
double cor, cor2, value;


value=log(2*M_PI);

for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;}
cor2=pow(cor,2);

logdens[i]=-value-.5*log(1-cor2)-.5*(pow(phen1[i],2)+pow(phen2[i],2)-2*cor*phen1[i]*phen2[i])*pow(1-cor2,-1);
}
}

////////

void get_probs(double *phen1, double *phen2, double *thresh1, double *thresh2, double *relats, double *relats2, int total, double *hers, double *probs, double *tops, double *sums, double prev1, double prev2, double ascer1, double ascer2)
//probs = tops / sums (when no ascertainment, sums=1 and probs=tops) - can not allow zero or one probabilities
{
int i;
double cor, value, value2, value3, p0, p1, p2, p3;


if(prev1==ascer1&&prev2==ascer2)
{
for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;}

value=bivnor(thresh1[i], thresh2[i], cor);
if(value==0){value=1e-10;}
if(value==1){value=1-1e-10;}

switch ((int)phen1[i]+2*(int)phen2[i])
{
case 0: probs[i]=gauss(thresh1[i])+gauss(thresh2[i])+value-1; break;
case 1: probs[i]=1-gauss(thresh1[i])-value; break;
case 2: probs[i]=1-gauss(thresh2[i])-value; break;
case 3: probs[i]=value; break;
}

tops[i]=probs[i];
sums[i]=1;
}
}
else
{
value2=prev1*(1-ascer1)/ascer1/(1-prev1);
value3=prev2*(1-ascer2)/ascer2/(1-prev2);

for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;}

value=bivnor(thresh1[i], thresh2[i], cor);
if(value==0){value=1e-10;}
if(value==1){value=1-1e-10;}

p0=gauss(thresh1[i])+gauss(thresh2[i])+value-1;
p1=1-gauss(thresh1[i])-value;
p2=1-gauss(thresh2[i])-value;
p3=value;

sums[i]=value2*value3*p0+value3*p1+value2*p2+p3;

switch ((int)phen1[i]+2*(int)phen2[i])
{
case 0: probs[i]=value2*value3*p0/sums[i];tops[i]=p0; break;
case 1: probs[i]=value3*p1/sums[i];tops[i]=p1; break;
case 2: probs[i]=value2*p2/sums[i];tops[i]=p2; break;
case 3: probs[i]=p3/sums[i];tops[i]=p3; break;
}
}
}
}

////////

void get_firsts_quant(double *phen1, double *phen2, double *relats, double *relats2, int total, double *hers, double *firsts)
{
int i, count;
double cor, cor2;


count=0;
for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;count++;}
cor2=pow(cor,2);

firsts[i]=cor*pow(1-cor2,-1)+(phen1[i]-cor*phen2[i])*(phen2[i]-cor*phen1[i])*pow(1-cor2,-2);
}

if(count>0){printf("Warning, %d correlations were greater than 0.9999\n",count);}
}

////////

void get_firsts_binary(double *phen1, double *phen2, double *thresh1, double *thresh2, double *relats, double *relats2, int total, double *hers, double *firsts, double *tops, double *sums, double prev1, double prev2, double ascer1, double ascer2)
//without ascertainment, first deriv is sign phi / top, where phi=pdf(thresh1,thresh2,cor)
//with ascertainment, subtract (s-1)^2 phi / sum, where s = prev*(1-ascer)/ascer/(1-prev)
{
int i, count;
double cor, cor2, value, scale;


if(prev1==ascer1&&prev2==ascer2){scale=0;}
else{scale=(prev1*(1-ascer1)/ascer1/(1-prev1)-1)*(prev2*(1-ascer2)/ascer2/(1-prev2)-1);}

count=0;
for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;count++;}
cor2=pow(cor,2);

value=.5*M_1_PI*pow(1-cor2,-.5)*exp(-.5*(pow(thresh1[i],2)+pow(thresh2[i],2)-2*cor*thresh1[i]*thresh2[i])*pow(1-cor2,-1));

if(phen1[i]+phen2[i]!=1)	//00 or 11
{firsts[i]=value/tops[i]-scale*value/sums[i];}
else	//10 or 01
{firsts[i]=-value/tops[i]-scale*value/sums[i];}
}

if(count>0){printf("Warning, %d correlations were greater than 0.9999\n",count);}
}

////////

void get_seconds_quant(double *phen1, double *phen2, double *relats, double *relats2, int total, double *hers, double *seconds)
{
int i, count;
double cor, cor2, value;


count=0;
for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;count++;}
cor2=pow(cor,2);

value=pow(phen1[i],2)+pow(phen2[i],2)-2*cor*phen1[i]*phen2[i];

seconds[i]=(1+cor2+4*cor*phen1[i]*phen2[i]-value)*pow(1-cor2,-2)-4*cor2*value*pow(1-cor2,-3);
}
}

////////

void get_seconds_binary(double *phen1, double *phen2, double *thresh1, double *thresh2, double *relats, double *relats2, int total, double *hers, double *seconds, double *tops, double *sums, double prev1, double prev2, double ascer1, double ascer2)
//without ascertainment, second deriv is -sign phi / top phi' -(phi/top)^2
//with ascertainment, add (phi (s-1)^2 / sum)^2 -(s-1)^2 phi phi' / sum
{
int i;
double cor, cor2, value, value2, scale;


if(prev1==ascer1&&prev2==ascer2){scale=0;}
else{scale=(prev1*(1-ascer1)/ascer1/(1-prev1)-1)*(prev2*(1-ascer2)/ascer2/(1-prev2)-1);}

for(i=0;i<total;i++)
{
cor=relats[i]*hers[0]+relats2[i]*hers[1];
if(cor>0.9999){cor=0.9999;}
cor2=pow(cor,2);

value=.5*M_1_PI*pow(1-cor2,-.5)*exp(-.5*(pow(thresh1[i],2)+pow(thresh2[i],2)-2*cor*thresh1[i]*thresh2[i])*pow(1-cor2,-1));

value2=cor*pow(1-cor2,-1)+(thresh1[i]-cor*thresh2[i])*(thresh2[i]-cor*thresh1[i])*pow(1-cor2,-2);

if(phen1[i]+phen2[i]!=1)	//00 or 11
{seconds[i]=value/tops[i]*value2-pow(value/tops[i],2)+pow(scale*value/sums[i],2)-scale*value*value2/sums[i];}
else	//10 or 01
{seconds[i]=-value/tops[i]*value2-pow(value/tops[i],2)+pow(scale*value/sums[i],2)-scale*value*value2/sums[i];}
}
}

///////////////////////////

int solve_fam(double *stats, int total, int total2, int start, int end, double *relats, double *relats2, double *phen1, double *phen2, double *weights, double *thresh1, double *thresh2, double prev1, double prev2, double ascer1, double ascer2, double tol, int maxiter, char *outfile, int constrain, int type, int flag)
//type=0 - quant, type=1 - binary
//flag=0, univariate main; flag=1 bivariate; flag=2 jackknife
{
int i, count, cflag, rflag, one=1, two=2;
double sum, sum2, value, value2, value3, alpha, beta;

double *firsts, *seconds, *logdens, *probs, *tops, *sums;

double relax, hers[3], herdiffs[2], hersds[3], likenull, like, likeold, like2, diff;
double AI[4], AI2[2], AI3[4], BI[2], lrtstat, lrtpva;

char filename[500];
FILE *output;


//allocate

firsts=malloc(sizeof(double)*total);
seconds=malloc(sizeof(double)*total);

if(type==0)
{
logdens=malloc(sizeof(double)*total);
}
else
{
probs=malloc(sizeof(double)*total);
tops=malloc(sizeof(double)*total);
sums=malloc(sizeof(double)*total);
}

if(flag==0||flag==1)
{
//prepare to screen and file print
printf("Iter\tGenetic\tEnvironment\tLikelihood\tDifference\tTarget\n");

sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output,"Iter\tGenetic\tEnvironment\tLikelihood\tDifference\tTarget\n");
fclose(output);
}

//get null likelihood
hers[0]=0;hers[1]=0;
if(type==0)
{
get_log_density(phen1, phen2, relats, relats2, total, hers, logdens);
likenull=0;for(i=0;i<total;i++){likenull+=logdens[i]*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){likenull-=logdens[i]*weights[i];}
}
}
else
{
get_probs(phen1, phen2, thresh1, thresh2, relats, relats2, total, hers, probs, tops, sums, prev1, ascer1, prev2, ascer2);
likenull=0;for(i=0;i<total;i++){likenull+=log(probs[i])*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){likenull-=log(probs[i])*weights[i];}
}
}

//will start values at zero
hers[0]=0.;hers[1]=0;
if(type==0)
{
get_log_density(phen1, phen2, relats, relats2, total, hers, logdens);
like=0;for(i=0;i<total;i++){like+=logdens[i]*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){like-=logdens[i]*weights[i];}
}
}
else
{
get_probs(phen1, phen2, thresh1, thresh2, relats, relats2, total, hers, probs, tops, sums, prev1, ascer1, prev2, ascer2);
like=0;for(i=0;i<total;i++){like+=log(probs[i])*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){like-=log(probs[i])*weights[i];}
}
}

//will first consider only genetic component, then add in environmental (if used)
count=0;
rflag=0;
while(rflag<total2)
{
cflag=1;
while(1)
{
likeold=like;

if(flag==0||flag==1)	//print update
{
if(count==0){printf("Start\t%.4f\t%.4f\t%.2f\tNA\t%.6f\n", hers[0], hers[1], like, tol);}
else{printf("%d\t%.4f\t%.4f\t%.2f\t%.6f\t%.6f\n", count, hers[0], hers[1], like, diff, tol);}

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n", filename);exit(1);}
if(count==0){fprintf(output,"%d\t%.4f\t%.4f\t%.2f\tNA\t%.6f\n", count, hers[0], hers[1], like, tol);}
else{fprintf(output,"%d\t%.4f\t%.4f\t%.2f\t%.6f\t%.6f\n", count, hers[0], hers[1], like, diff, tol);}
fclose(output);
}

//get firsts and seconds
if(type==0)
{
get_firsts_quant(phen1, phen2, relats, relats2, total, hers, firsts);
get_seconds_quant(phen1, phen2, relats, relats2, total, hers, seconds);
}
else
{
get_firsts_binary(phen1, phen2, thresh1, thresh2, relats, relats2, total, hers, firsts, tops, sums, prev1, ascer1, prev2, ascer2);
get_seconds_binary(phen1, phen2, thresh1, thresh2, relats, relats2, total, hers, seconds, tops, sums, prev1, ascer1, prev2, ascer2);
}

//fill derivatives for first heritability
value=0;for(i=0;i<total;i++){value+=relats[i]*firsts[i]*weights[i];}
value2=0;for(i=0;i<total;i++){value2+=pow(relats[i],2)*seconds[i]*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){value-=relats[i]*firsts[i]*weights[i];}
for(i=start;i<end;i++){value2-=pow(relats[i],2)*seconds[i]*weights[i];}
}
BI[0]=value;
AI[0]=-value2;

if(total2==1||rflag==0)	//only updating one heritability - will set AI[3]=1 so can invert
{
BI[1]=0;
AI[1]=0;
AI[2]=0;
AI[3]=1;
}
else	//updating two heritabilities
{
value=0;for(i=0;i<total;i++){value+=relats2[i]*firsts[i]*weights[i];}
value2=0;for(i=0;i<total;i++){value2+=relats[i]*relats2[i]*seconds[i]*weights[i];}
value3=0;for(i=0;i<total;i++){value3+=pow(relats2[i],2)*seconds[i]*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){value-=relats2[i]*firsts[i]*weights[i];}
for(i=start;i<end;i++){value2-=relats[i]*relats2[i]*seconds[i]*weights[i];}
for(i=start;i<end;i++){value3-=pow(relats2[i],2)*seconds[i]*weights[i];}
}
BI[1]=value;
AI[1]=-value2;
AI[2]=-value2;
AI[3]=-value3;
}

//invert AI
(void)eigen_invert(AI, 2, AI2, -1, AI3, 1);
if(total2==1){AI[3]=0;}

if(count>0)	//test for convergence
{
if(fabs(diff)<tol){break;}
}

if(count==maxiter)
{
if(rflag==total2-1){printf("Warning, the optimizer failed to converge within %d iterations; consider using \"--max-iter\" and/or \"--tolerance\" to change the iteration limit and tolerance", maxiter);}
cflag=0;
break;
}

//get proposed move
alpha=1.0;beta=0.0;
dgemv_("N", &two, &two, &alpha, AI, &two, BI, &one, &beta, herdiffs, &one);

//ensure moves not too large
if(herdiffs[0]<-0.2){herdiffs[0]=-.2;}
if(herdiffs[0]>0.2){herdiffs[0]=.2;}
if(herdiffs[1]<-0.2){herdiffs[1]=-.2;}
if(herdiffs[1]>0.2){herdiffs[1]=.2;}

if(constrain==1)	//new variances can not be negative
{
if(hers[0]+herdiffs[0]<0){herdiffs[0]=-hers[0];}
if(hers[1]+herdiffs[1]<0){herdiffs[1]=-hers[1];}
}

//ensure heritabilities and their sum remain within -1 and 1
if(hers[0]+herdiffs[0]<-0.99){herdiffs[0]=-0.99-hers[0];}
if(hers[0]+herdiffs[0]>0.99){herdiffs[0]=0.99-hers[0];}

if(hers[1]+herdiffs[1]<-0.99){herdiffs[1]=-0.99-hers[1];}
if(hers[1]+herdiffs[1]>0.99){herdiffs[1]=0.99-hers[1];}

sum=hers[0]+hers[1];
sum2=hers[0]+herdiffs[0]+hers[1]+herdiffs[1];

if(sum2>.99)	//too high
{
value=(.99-sum)/(sum2-sum);
herdiffs[0]*=value;
herdiffs[1]*=value;
}
if(sum2<-.99)	//too low
{
value=(sum+.99)/(sum-sum2);
herdiffs[0]*=value;
herdiffs[1]*=value;
}

relax=1;
while(relax>0.0001)
{
//try moving value*herdiff
hers[0]+=relax*herdiffs[0];
hers[1]+=relax*herdiffs[1];

//get likelihood
if(type==0)
{
get_log_density(phen1, phen2, relats, relats2, total, hers, logdens);
like2=0;for(i=0;i<total;i++){like2+=logdens[i]*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){like2-=logdens[i]*weights[i];}
}
}
else
{
get_probs(phen1, phen2, thresh1, thresh2, relats, relats2, total, hers, probs, tops, sums, prev1, ascer1, prev2, ascer2);
like2=0;for(i=0;i<total;i++){like2+=log(probs[i])*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){like2-=log(probs[i])*weights[i];}
}
}

//see whether to accept move, or switch to a smaller move
if(like2>like-tol)	//accept move
{
like=like2;
break;
}
else	//move back and next turn try smaller move
{
hers[0]-=relax*herdiffs[0];
hers[1]-=relax*herdiffs[1];
relax*=.5;
}
}

//refresh likelihood
if(type==0)
{
get_log_density(phen1, phen2, relats, relats2, total, hers, logdens);
like=0;for(i=0;i<total;i++){like+=logdens[i]*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){like-=logdens[i]*weights[i];}
}
}
else
{
get_probs(phen1, phen2, thresh1, thresh2, relats, relats2, total, hers, probs, tops, sums, prev1, ascer1, prev2, ascer2);
like=0;for(i=0;i<total;i++){like+=log(probs[i])*weights[i];}
if(flag==2)	//remove excluded block contributions
{
for(i=start;i<end;i++){like-=log(probs[i])*weights[i];}
}
}

//get difference
diff=like-likeold;

count++;
}	//end of inner while loop

rflag++;
}	//end of outer while loop
if(flag==0||flag==1){printf("\n");}

//get SDs 
if(AI[0]>=0){hersds[0]=pow(AI[0],.5);}
else{hersds[0]=pow(AI[0],.5);}
if(AI[3]>=0){hersds[1]=pow(AI[3],.5);}
else{hersds[1]=pow(AI[3],.5);}

//save details for her_all as final element of hers / hersds
hers[2]=hers[0]+hers[1];
sum=AI[0]+AI[1]+AI[2]+AI[3];
if(sum>=0){hersds[2]=pow(sum,.5);}
else{hersds[2]=-9999;}

//get test stats
lrtstat=2*(like-likenull);
lrtpva=erfc(pow(lrtstat,.5)*M_SQRT1_2);

if(flag==0)	//need to return hers, hersds, likes and lrts
{
stats[0]=hers[0];
stats[1]=hers[1];
stats[2]=hers[2];
stats[3]=hersds[0];
stats[4]=hersds[1];
stats[5]=hersds[2];
stats[6]=likenull;
stats[7]=like;
stats[8]=lrtstat;
stats[9]=lrtpva;
}
else	//need to just return hers
{
stats[0]=hers[0];
stats[1]=hers[1];
stats[2]=hers[2];
}

free(firsts);free(seconds);
if(type==0){free(logdens);}
else{free(probs);free(tops);free(sums);}

return(cflag);
}	//end of solve_fam

///////////////////////////

void mle_fam(int ns, int num_covars, double *Y, double *Z, char *relfile, double prev, double prev2, double tol, int maxiter, char **ids3, char *outfile, int constrain, int type, int cordups, int num_blocks, int flag)
//type=0 - quant, type=1 - binary
//flag=0, univariate; flag=1, bivariate
{
int i, j, p, start, end, count, total, total2, cflag;
double sum, sumsq, wsum, mean, mean2, mean3, var, var2, var3, value, value2, value3, value4, ascer, ascer2, factor;

int *indexer, *indexer2, *usedids, *order;
char **wantids, **wantids2, **wantids3;

double *thetas, *thetas2, *thetasds, *thetasds2, *thetapvas, *thetapvas2, *Yadj, *Yadj2, covher, covher2, topher, topher2;

int *tallys;
double *weights, *relats, *relats2, *relats3, *phen1, *phen2, *thresh1, *thresh2;

double *xweights, *xrelats, *xrelats2, *xphen1, *xphen2, *xthresh1, *xthresh2;

double *stats, *stats2;

double hers[3],hersds[3], likenull, like, lrtstat, lrtpva;

char filename2[500], filename3[500], filename4[500], filename5[500];
FILE *output2, *output3, *output4, *output5;


//deal with covariates (and maybe ascertainment)

thetas=malloc(sizeof(double)*num_covars);
thetasds=malloc(sizeof(double)*num_covars);
thetapvas=malloc(sizeof(double)*num_covars);
Yadj=malloc(sizeof(double)*ns);

if(flag==1)
{
thetas2=malloc(sizeof(double)*num_covars);
thetasds2=malloc(sizeof(double)*num_covars);
thetapvas2=malloc(sizeof(double)*num_covars);
Yadj2=malloc(sizeof(double)*ns);
}

if(type==1)	//deal with binary aspects
{
if(flag==0)	//one trait - get ascertainment and maybe set prevalence
{
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
ascer=sum/ns;

if(prev==-9999)
{
prev=ascer;
printf("Will assume the population prevalence matches the observed prevalence (%.6f); if you wish to specify a different value, use \"--prevalence\"\n\n", prev);
}
}
else	//two traits - get ascertainments and maybe set prevalences
{
sum=0;for(i=0;i<ns;i++){sum+=Y[i];}
ascer=sum/ns;
sum=0;for(i=0;i<ns;i++){sum+=Y[i+ns];}
ascer2=sum/ns;

if(prev==-9999)
{
prev=ascer;
printf("Will assume the population prevalence for Trait 1 matches its observed prevalence (%.6f); if you wish to specify a different value, use \"--prevalence\"\n\n", prev);
}
if(prev2==-9999)
{
prev2=ascer2;
printf("Will assume the population prevalence for Trait 2 matches its observed prevalence (%.6f); if you wish to specify a different value, use \"--prevalence2\"\n\n", prev2);
}
}
}

//regress trait 1 on covariates
if(type==0)	//linear model - Yadj contains standardized residuals
{reg_covar_lin(Y, Z, ns, num_covars, 0, thetas, thetasds, thetapvas, Yadj, 1, &covher, &topher);}
else	//logistic model - fill Yadj with new thresholds
{reg_covar_log(Y, Z, ns, num_covars, 0, NULL, thetas, thetasds, thetapvas, Yadj, 1, &covher, &topher, prev, 0.001, 100);}

if(flag==1)	//regress trait 2 on covariates
{
if(type==0){reg_covar_lin(Y+ns, Z, ns, num_covars, 0, thetas2, thetasds2, thetapvas2, Yadj2, 1, &covher2, &topher2);}
else{reg_covar_log(Y+ns, Z, ns, num_covars, 0, NULL, thetas2, thetasds2, thetapvas2, Yadj2, 1, &covher2, &topher2, prev2, 0.001, 100);}
}

if(num_covars>1)
{
if(flag==0){printf("Proportion of variance explained by the %d covariates: %.4f\n\n", num_covars, covher);}
else{printf("Proportions of variance explained by the %d covariates: %.4f and %.4f\n\n", num_covars, covher, covher2);}
}

////////

//deal with relfile

total=countrows(relfile);
total2=countcols(relfile)-4;

indexer=malloc(sizeof(int)*total);
indexer2=malloc(sizeof(int)*total);
usedids=malloc(sizeof(int)*total);

wantids=malloc(sizeof(char*)*total);
wantids2=malloc(sizeof(char*)*total);
relats=malloc(sizeof(double)*total);
relats2=malloc(sizeof(double)*total);

//read pairs and relatedness
printf("Reading details for %d related pairs from %s\n\n", total, relfile);
read_ids(relfile, NULL, NULL, wantids, total, NULL, 0, 0);
read_ids(relfile, NULL, NULL, wantids2, total, NULL, 0, 2);
read_values(relfile, relats, total, NULL, 5, 0, 0);
if(total2==1)
{
for(i=0;i<total;i++){relats2[i]=0;}
}
else{read_values(relfile, relats2, total, NULL, 6, 0, 0);}

for(i=0;i<total;i++)
{
if(strcmp(wantids[i], wantids2[i])==0)
{
read_ids(relfile, wantids, wantids2, NULL, total, NULL, 0, 0);
printf("Error reading Row %d of %s; the sample %s %s is listed twice IDs\n\n", i+1, relfile, wantids[i], wantids2[i]);exit(1);}
}

count=0;
for(i=0;i<total;i++){count+=(relats[i]>1||relats2[i]>1);}
if(count>0){printf("Warning, %d pairs of samples have similarity greater than one\n\n", count);}

//work out which pairs are both present

for(i=0;i<total;i++){usedids[i]=0;}

count=find_strings(wantids, total, ids3, ns, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count==0)
{
if(flag==0){printf("Error, none of the %d samples in Columns 1 and 2 of %s have phenotypes\n\n", total, relfile);exit(1);}
else{printf("Error, none of the %d samples in Columns 1 and 2 of %s have both phenotypes\n\n", total, relfile);exit(1);}
}
for(i=0;i<count;i++){usedids[indexer[i]]++;}

count=find_strings(wantids2, total, ids3, ns, indexer2, NULL, NULL, NULL, NULL, NULL, 3);
if(count==0)
{
if(flag==0){printf("Error, none of the %d samples in Columns 3 and 4 of %s have phenotypes\n\n", total, relfile);exit(1);}
else{printf("Error, none of the %d samples in Columns 3 and 4 of %s have both phenotypes\n\n", total, relfile);exit(1);}
}
for(i=0;i<count;i++){usedids[indexer2[i]]++;}

count=0;for(i=0;i<total;i++){count+=(usedids[i]==2);}
if(count==0){printf("Error, none of the %d pairs of samples in %s both have phenotypes\n\n", total, relfile);exit(1);}
if(count<total){printf("Warning, only %d of the %d pairs of samples in %s both have phenotypes\n\n", count, total, relfile);}

if(count<total)	//squeeze down ids and relats
{
count=0;
for(i=0;i<total;i++)
{
if(usedids[i]==2)
{
if(count!=i)
{
free(wantids[count]);copy_string(wantids,count,wantids[i]);
free(wantids2[count]);copy_string(wantids2,count,wantids2[i]);
relats[count]=relats[i];
relats2[count]=relats2[i];
}
count++;
}}
for(i=count;i<total;i++){free(wantids[i]);free(wantids2[i]);}
total=count;
}

if(num_blocks!=-9999)	//shuffle order of pairs
{
if(num_blocks==-1||num_blocks>total){num_blocks=total;}

order=malloc(sizeof(int)*total);
wantids3=malloc(sizeof(char*)*total);
relats3=malloc(sizeof(double)*total);

for(i=0;i<total;i++){order[i]=i;}
//permute_int(order,total);

for(i=0;i<total;i++){copy_string(wantids3,i,wantids[i]);free(wantids[i]);}
for(i=0;i<total;i++){copy_string(wantids,i,wantids3[order[i]]);free(wantids3[order[i]]);}
for(i=0;i<total;i++){copy_string(wantids3,i,wantids2[i]);free(wantids2[i]);}
for(i=0;i<total;i++){copy_string(wantids2,i,wantids3[order[i]]);free(wantids3[order[i]]);}

for(i=0;i<total;i++){relats3[i]=relats[i];}
for(i=0;i<total;i++){relats[i]=relats3[order[i]];}
for(i=0;i<total;i++){relats3[i]=relats2[i];}
for(i=0;i<total;i++){relats2[i]=relats3[order[i]];}

free(order);free(wantids3);free(relats3);
}

////////

//allocate

tallys=malloc(sizeof(double)*ns);
weights=malloc(sizeof(double)*total);

if(flag==0)	//one trait
{
phen1=malloc(sizeof(double)*total);
phen2=malloc(sizeof(double)*total);

if(type==1)
{
thresh1=malloc(sizeof(double)*total);
thresh2=malloc(sizeof(double)*total);
}

stats=malloc(sizeof(double)*10);
if(num_blocks!=-9999){stats2=malloc(sizeof(double)*3*num_blocks);}
}
else	//two traits
{
phen1=malloc(sizeof(double)*total*2);
phen2=malloc(sizeof(double)*total*2);
xweights=malloc(sizeof(double)*2*total);
xrelats=malloc(sizeof(double)*2*total);
xrelats2=malloc(sizeof(double)*2*total);
xphen1=malloc(sizeof(double)*2*total);
xphen2=malloc(sizeof(double)*2*total);

if(type==1)
{
thresh1=malloc(sizeof(double)*total*2);
thresh2=malloc(sizeof(double)*total*2);
xthresh1=malloc(sizeof(double)*2*total);
xthresh2=malloc(sizeof(double)*2*total);
}

stats=malloc(sizeof(double)*9*2);
stats2=malloc(sizeof(double)*9*num_blocks);
}

//load up responses (and maybe thresholds)

count=find_strings(wantids, total, ids3, ns, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count!=total){printf("Error 933E; please tell Doug %d %d %d\n\n", count, total, countrows(relfile));exit(1);}

count=find_strings(wantids2, total, ids3, ns, NULL, indexer2, NULL, NULL, NULL, NULL, 3);
if(count!=total){printf("Error 933F; please tell Doug %d %d %d\n\n", count, total, countrows(relfile));exit(1);}

//load trait 1
for(i=0;i<total;i++)
{
if(type==0)
{
phen1[i]=Yadj[indexer[i]];
phen2[i]=Yadj[indexer2[i]];
}
else
{
phen1[i]=Y[indexer[i]];
phen2[i]=Y[indexer2[i]];
thresh1[i]=Yadj[indexer[i]];
thresh2[i]=Yadj[indexer2[i]];
}
}

if(flag==1)	//load trait 2
{
for(i=0;i<total;i++)
{
if(type==0)
{
phen1[i+total]=Yadj2[indexer[i]];
phen2[i+total]=Yadj2[indexer2[i]];
}
else
{
phen1[i+total]=Y[indexer[i]+ns];
phen2[i+total]=Y[indexer2[i]+ns];
thresh1[i+total]=Yadj2[indexer[i]];
thresh2[i+total]=Yadj2[indexer2[i]];
}
}
}

//tally pairs
for(i=0;i<ns;i++){tallys[i]=0;}
for(i=0;i<total;i++){tallys[indexer[i]]++;tallys[indexer2[i]]++;}

if(cordups==1)	//work out weightings for each pair
{
//for(i=0;i<total;i++){weights[i]=pow(tallys[indexer[i]]*tallys[indexer2[i]],-.5);}
for(i=0;i<total;i++){weights[i]=.5*(pow(tallys[indexer[i]],-1)+pow(tallys[indexer2[i]],-1));}
}
else	//just give weighting one
{
for(i=0;i<total;i++){weights[i]=1.0;}
}
wsum=0;for(i=0;i<total;i++){wsum+=weights[i];}

if(flag==1)	//load cross stuff
{
for(i=0;i<total;i++)
{
//first phen of person 1 and second phen of person 2
xweights[0+i*2]=weights[i];
xrelats[0+i*2]=relats[i];
xrelats2[0+i*2]=relats2[i];
xphen1[0+i*2]=phen1[i];
xphen2[0+i*2]=phen2[i+total];
if(type==1)
{
xthresh1[0+i*2]=thresh1[i];
xthresh2[0+i*2]=thresh2[i+total];
}

//second phen of person 1 and first phen of person 2
xweights[1+i*2]=weights[i];
xrelats[1+i*2]=relats[i];
xrelats2[1+i*2]=relats2[i];
xphen1[1+i*2]=phen1[i+total];
xphen2[1+i*2]=phen2[i];
if(type==1)
{
xthresh1[1+i*2]=thresh1[i+total];
xthresh2[1+i*2]=thresh2[i];
}
}
}

////////

if(flag==0)	//normal solve (possibly with jackknifing)
{
cflag=solve_fam(stats, total, total2, -9999, -9999, relats, relats2, phen1, phen2, weights, thresh1, thresh2, prev, prev, ascer, ascer, tol, maxiter, outfile, constrain, type, 0);

hers[0]=stats[0];
hers[1]=stats[1];
hers[2]=stats[2];
hersds[0]=stats[3];
hersds[1]=stats[4];
hersds[2]=stats[5];
likenull=stats[6];
like=stats[7];
lrtstat=stats[8];
lrtpva=stats[9];

if(num_blocks!=-9999)	//jackknife
{
for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*total;
end=(double)(p+1)/num_blocks*total;
if(p%10==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}

(void)solve_fam(stats2+p*3, total, total2, start, end, relats, relats2, phen1, phen2, weights, thresh1, thresh2, prev, prev, ascer, ascer, tol, maxiter, outfile, constrain, type, 2);
}
printf("\n");

//replace hersds with jackknife estimates

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats2[p*3];sumsq+=pow(stats2[p*3],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
hersds[0]=pow(var,.5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats2[1+p*3];sumsq+=pow(stats2[1+p*3],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
hersds[1]=pow(var,.5);

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats2[2+p*3];sumsq+=pow(stats2[2+p*3],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
hersds[2]=pow(var,.5);
}

sprintf(filename2,"%s.mle", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output2, "Total_Pairs %d\n", total);
fprintf(output2, "Effective_Size %.2f\n", wsum);
if(prev!=-9999){fprintf(output2, "Prevalence %.6f\n", prev);}
else{fprintf(output2, "Prevalence NA\n");}
if(cflag==1){fprintf(output2,"Converged YES\n");}
else{fprintf(output2,"Converged NO\n");}
fprintf(output2, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
fprintf(output2, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);

fprintf(output2, "Component Heritability SD\n");
fprintf(output2, "Genetic %.4f %.4f\n", hers[0], hersds[0]);
if(total2==2){fprintf(output2, "Environmental %.4f %.4f\n", hers[1], hersds[1]);}
else{fprintf(output2, "Environmental 0 0\n");}
fprintf(output2, "Total %.4f %.4f\n", hers[2], hersds[2]);
fclose(output2);

sprintf(filename3,"%s.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Component Log_Odds SD P\n");
fprintf(output3, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output3);

if(type==0&&prev!=-9999)	//save liability version
{
factor=get_factor(Y, ns, prev, -9999, outfile);

sprintf(filename4,"%s.mle.liab", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher*factor);
fprintf(output4, "Total_Pairs %d\n", total);
fprintf(output4, "Prevalence %.6f\n", prev);
if(cflag==1){fprintf(output4,"Converged YES\n");}
else{fprintf(output4,"Converged NO\n");}
fprintf(output4, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
fprintf(output4, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);

fprintf(output4, "Component Heritability SD\n");
fprintf(output4, "Genetic %.4f %.4f\n", hers[0]*factor, hersds[0]*factor);
if(total2==2){fprintf(output4, "Environmental %.4f %.4f\n", hers[1]*factor, hersds[1]*factor);}
else{fprintf(output4, "Environmental 0 0\n");}
fprintf(output4, "Total %.4f %.4f\n", hers[2]*factor, hersds[2]*factor);
fclose(output4);
}

sprintf(filename5,"%s.mle.marginal", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Coeffsfile %s.coeff\nCovar_Heritability %.4f\n", outfile, covher);
fprintf(output5, "Total_Pairs %d\n", total);
if(prev!=-9999){fprintf(output5, "Prevalence %.6f\n", prev);}
else{fprintf(output5, "Prevalence NA\n");}
if(cflag==1){fprintf(output5,"Converged YES\n");}
else{fprintf(output5,"Converged NO\n");}
fprintf(output5, "Null_Likelihood %.6f\nAlt_Likelihood %.6f\n", likenull, like);
fprintf(output5, "LRT_Stat %.4f\nLRT_P %.4e\n", lrtstat, lrtpva);

fprintf(output5, "Component Heritability SD\n");
fprintf(output5, "Genetic %.4f %.4f\n", hers[0]*(1-covher), hersds[0]*(1-covher));
if(total2==2){fprintf(output5, "Environmental %.4f %.4f\n", hers[1]*(1-covher), hersds[1]*(1-covher));}
else{fprintf(output5, "Environmental 0 0\n");}
fprintf(output5, "Total %.4f %.4f\n", hers[2]*(1-covher), hersds[2]*(1-covher));
fclose(output5);

printf("Heritability estimates saved in %s", filename2);
if(type==0&&prev!=-9999){printf(", with a liability version saved in %s", filename4);}
printf("\n\n");
}
else	//bivariate solve (with jackknifing)
{
cflag=0;

printf("Estimating heritability for Trait 1\n\n");
cflag+=solve_fam(stats, total, total2, -9999, -9999, relats, relats2, phen1, phen2, weights, thresh1, thresh2, prev, prev, ascer, ascer, tol, maxiter, outfile, constrain, type, 1);

printf("Estimating heritability for Trait 2\n\n");
cflag+=solve_fam(stats+3, total, total2, -9999, -9999, relats, relats2, phen1+total, phen2+total, weights, thresh1+total, thresh2+total, prev2, prev2, ascer2, ascer2, tol, maxiter, outfile, constrain, type, 1);

printf("Estimating coheritability between Traits 1 and 2\n\n");
cflag+=solve_fam(stats+6, 2*total, total2, -9999, -9999, xrelats, xrelats2, xphen1, xphen2, xweights, xthresh1, xthresh2, prev, prev2, ascer, ascer2, tol, maxiter, outfile, constrain, type, 1);

for(p=0;p<num_blocks;p++)
{
start=(double)p/num_blocks*total;
end=(double)(p+1)/num_blocks*total;
if(p%10==0){printf("Performing Jackknife %d out of %d\n", p+1, num_blocks);}

//first trait
(void)solve_fam(stats2+p*9, total, total2, start, end, relats, relats2, phen1, phen2, weights, thresh1, thresh2, prev, prev, ascer, ascer, tol, maxiter, outfile, constrain, type, 2);

//second trait
(void)solve_fam(stats2+3+p*9, total, total2, start, end, relats, relats2, phen1+total, phen2+total, weights, thresh1+total, thresh2+total, prev2, prev2, ascer2, ascer2, tol, maxiter, outfile, constrain, type, 2);

//cross trait
(void)solve_fam(stats2+6+p*9, 2*total, total2, 2*start, 2*end, xrelats, xrelats2, xphen1, xphen2, xweights, xthresh1, xthresh2, prev, prev2, ascer, ascer2, tol, maxiter, outfile, constrain, type, 2);
}
printf("\n");

//fill second column of stats with jackknife estimates
for(j=0;j<9;j++)
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats2[j+p*9];sumsq+=pow(stats2[j+p*9],2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
stats[j+9]=pow(var,.5);
}

//genetic correlation
if(stats[0]>0&&stats[3]>0)
{
value=stats[6]*pow(stats[0]*stats[3],-.5);
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){value4=stats2[6+p*9]*pow(stats2[0+p*9]*stats2[3+p*9],-.5);sum+=value4;sumsq+=pow(value4,2);}
mean=sum/num_blocks;
var=(num_blocks-1)*(sumsq/num_blocks-pow(mean,2));
}
else
{value=-9999;var=0;}

//environmental correlation
if(stats[1]>0&&stats[4]>0)
{
value2=stats[7]*pow(stats[1]*stats[4],-.5);
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){value4=stats2[7+p*9]*pow(stats2[1+p*9]*stats2[4+p*9],-.5);sum+=value4;sumsq+=pow(value4,2);}
mean2=sum/num_blocks;
var2=(num_blocks-1)*(sumsq/num_blocks-pow(mean2,2));
}
else
{value2=-9999;var2=0;}

//total correlation
if(stats[2]>0&&stats[5]>0)
{
value3=stats[8]*pow(stats[2]*stats[5],-.5);
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){value4=stats2[8+p*9]*pow(stats2[2+p*9]*stats2[5+p*9],-.5);sum+=value4;sumsq+=pow(value4,2);}
mean3=sum/num_blocks;
var3=(num_blocks-1)*(sumsq/num_blocks-pow(mean3,2));
}
else
{value3=-9999;var3=0;}

sprintf(filename2,"%s.cross", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Coeffsfile %s.coeff\nTrait1_Covar_Heritability %.4f\nTrait1_Covar_Heritability %.4f\n", outfile, covher, covher2);
fprintf(output2, "Total_Pairs %d\n", total);
fprintf(output2, "Effective_Size %.2f\n", wsum);
if(prev!=-9999){fprintf(output2, "Trait1_Prevalence %.6f\n", prev);}
else{fprintf(output2, "Trait1_Prevalence NA\n");}
if(prev2!=-9999){fprintf(output2, "Trait2_Prevalence %.6f\n", prev2);}
else{fprintf(output2, "Trait2_Prevalence NA\n");}
if(cflag==3){fprintf(output2,"All_Converged YES\n");}
else{fprintf(output2,"Converged NO\n");}

fprintf(output2, "Component Value SD\n");

fprintf(output2, "Trait1_Genetic_Heritability %.4f %.4f\n", stats[0], stats[0+9]);
if(total2==2){fprintf(output2, "Trait1_Environmental_Contribution %.4f %.4f\n", stats[1], stats[1+9]);}
else{fprintf(output2, "Trait1_Environmental_Contribution 0 0\n");}
fprintf(output2, "Trait1_Total %.4f %.4f\n", stats[2], stats[2+9]);

fprintf(output2, "Trait2_Genetic_Heritability %.4f %.4f\n", stats[3], stats[3+9]);
if(total2==2){fprintf(output2, "Trait2_Environmental_Contribution %.4f %.4f\n", stats[4], stats[4+9]);}
else{fprintf(output2, "Trait2_Environmental_Contribution 0 0\n");}
fprintf(output2, "Trait2_Total %.4f %.4f\n", stats[5], stats[5+9]);

fprintf(output2, "Cross_Trait_Genetic_Coheritability %.4f %.4f\n", stats[6], stats[6+9]);
if(total2==2){fprintf(output2, "Cross_Trait_Environmental_Covariation %.4f %.4f\n", stats[7], stats[7+9]);}
else{fprintf(output2, "Cross_Trait_Environmental_Covariation 0 0\n");}
fprintf(output2, "Cross_Trait_Total %.4f %.4f\n", stats[8], stats[8+9]);

if(value!=-9999){fprintf(output2, "Cross_Trait_Genetic_Correlation %.4f %.4f\n", value, pow(var,.5));}
else{fprintf(output2, "Cross_Trait_Genetic_Correlation NA NA\n");}
if(value2!=-9999){fprintf(output2, "Cross_Trait_Environmental_Correlation %.4f %.4f\n", value2, pow(var2,.5));}
else{fprintf(output2, "Cross_Trait_Environmental_Correlation NA NA\n");}
if(value3!=-9999){fprintf(output2, "Cross_Trait_Total_Correlation %.4f %.4f\n", value3, pow(var3,.5));}
else{fprintf(output2, "Cross_Trait_Total_Correlation NA NA\n");}

fclose(output2);

sprintf(filename3,"%s.pheno1.coeff", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Component Log_Odds SD P\n");
fprintf(output3, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output3, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output3);

sprintf(filename4,"%s.pheno2.coeff", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Component Log_Odds SD P\n");
fprintf(output4, "Intercept %.6f %.6f %.4e\n", thetas2[0], thetasds2[0], thetapvas2[0]);
for(j=1;j<num_covars;j++){fprintf(output4, "Covariate_%d %.6f %.6f %.4e\n",j, thetas2[j], thetasds2[j], thetapvas2[j]);}
fclose(output4);

printf("Cross-trait estimates saved in %s\n\n", filename2);
}

////////

free(thetas);free(thetasds);free(thetapvas);free(Yadj);
if(flag==1){free(thetas2);free(thetasds2);free(thetapvas2);free(Yadj2);}
free(indexer);free(indexer2);free(usedids);
for(i=0;i<total;i++){free(wantids[i]);free(wantids2[i]);}free(wantids);free(wantids2);
free(relats);free(relats2);
free(tallys);free(weights);
if(flag==0)
{
free(phen1);free(phen2);
if(type==1){free(thresh1);free(thresh2);}
free(stats);
if(num_blocks!=-9999){free(stats2);}
}
else
{
free(phen1);free(phen2);free(xweights);free(xrelats);free(xrelats2);free(xphen1);free(xphen2);
if(type==1){free(thresh1);free(thresh2);free(xthresh1);free(xthresh2);}
free(stats);free(stats2);
}

}	//end of mle_fam

///////////////////////////


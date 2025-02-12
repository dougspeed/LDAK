/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Get per-predictor heritabilities using fast-he (data must be in binary format - also get cors and datasqs

///////////////////////////

printf("Estimating per-predictor heritabilities using Randomized Haseman-Elston Regression with %d random vectors\n\n", nmcmc);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating per-predictor heritabilities using Randomized Haseman-Elston Regression with %d random vectors\n", nmcmc);
fclose(output);

printf("Will divide the predictors into %d partitions (change this using \"--num-divides\")\n", ndivs);

if(ndivs>bittotal)
{
printf("Warning, the number of partitions has been reduced to %d, the total number of chunks\n", bittotal);
ndivs=bittotal;
}

if(maxcor==-9999)
{
maxcor=0.02;
if(maxcor<100.0/num_samples_use){maxcor=100.0/num_samples_use;}
}

printf("Will exclude chunks containing a predictor with estimated variance explained greater than %.4f (change this using \"--max-cor\")\n\n", maxcor);

////////

//allocate variables

total=nmcmc+num_resps_use;

kinsums=malloc(sizeof(double)*ndivs);

R=malloc(sizeof(double)*num_samples_use*total);
RTdata=malloc(sizeof(double)*bitsize*total);
MkinsD=malloc(sizeof(double)*num_samples_use*ndivs*num_pows);
MkinsR=malloc(sizeof(double)*num_samples_use*total*ndivs*num_pows);

KKtraces=malloc(sizeof(double)*ndivs*ndivs);
KKtraces2=malloc(sizeof(double)*ndivs);
KKtraces3=malloc(sizeof(double)*ndivs*ndivs);
KYtraces=malloc(sizeof(double)*ndivs);
KYtraces2=malloc(sizeof(double)*ndivs);

exps2=malloc(sizeof(double)*data_length*num_pows);

varexp=malloc(sizeof(double)*num_resps_use*num_pows);

cors2=malloc(sizeof(double)*num_resps_use*num_resps_use*num_pows);

ycounts=malloc(sizeof(int)*num_resps_use);

if(multi==1)
{
Wmat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat2=malloc(sizeof(double)*num_resps_use);
Wmat3=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat4=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat5=malloc(sizeof(double)*num_resps_use*num_resps_use);
Vmat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Vmat2=malloc(sizeof(double)*num_resps_use*num_resps_use);
Ymat=malloc(sizeof(double)*num_samples_use*num_resps_use);
}

////////

//have already filled Z and Yadj (Yadj will be scaled)

//fill start of R with random values, then end with adjusted phenotypes
for(g=0;g<nmcmc;g++)
{
for(i=0;i<num_samples_use;i++){R[(size_t)g*num_samples_use+i]=rnorm_safe();}
}
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){R[(size_t)(nmcmc+m)*num_samples_use+i]=Yadj[i+m*num_samples_use];}
}

//set MkinsD and MkinsR to zero
for(q=0;q<ndivs*num_pows;q++)
{
for(i=0;i<num_samples_use;i++){MkinsD[(size_t)q*num_samples_use+i]=0;}
}
for(k=0;k<total*ndivs*num_pows;k++)
{
for(i=0;i<num_samples_use;i++){MkinsR[(size_t)k*num_samples_use+i]=0;}
}

//will record number of chunks removed (both total, and per-phenotype)
wcount=0;
for(m=0;m<num_resps_use;m++){ycounts[m]=0;}

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%200==0)
{
printf("Calculating traces for Chunk %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating traces for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

//mark indicates the divide
mark=bit*ndivs/bittotal;

//read data, compute statistics, standardize and set missing to zero
if(dtype==1)	//fast way
{(void)read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, bedzeros, bedones, bedtwos, 2, maxthreads);}
else	//slow way
{
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1);
}

if(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999)	//perform qc
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
maf=centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]);
value=sqdevs[bitstart+j]/centres[bitstart+j]/(1-centres[bitstart+j]/2);}
if(minmaf!=-9999&&maf<minmaf){mults[bitstart+j]=-9999;}
if(maxmaf!=-9999&&maf>maxmaf){mults[bitstart+j]=-9999;}
if(minvar!=-9999&&sqdevs[bitstart+j]<minvar){mults[bitstart+j]=-9999;}
if(minobs!=-9999&&rates[bitstart+j]<minobs){mults[bitstart+j]=-9999;}
if(mininfo!=-9999&&value<mininfo){mults[bitstart+j]=-9999;}
}
}

//fill exps2 - one column for each power
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<num_pows;k++)
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
if(hwestand==1){exps2[bitstart+j+k*data_length]=weights[bitstart+j]*pow(centres[bitstart+j]*(1-centres[bitstart+j]/2),1+powers[k]);}
else{exps2[bitstart+j+k*data_length]=weights[bitstart+j]*pow(sqdevs[bitstart+j],1+powers[k]);}
}
else{exps2[bitstart+j+k*data_length]=0;}
}
}

if(num_fixed>1)	//adjust for covariates - will use non-weighted regression, even if binary
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

//compute t(data) R
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, R, &num_samples_use, &beta, RTdata, &bitlength);

for(m=0;m<num_resps_use;m++)	//look for big effect predictors
{
for(j=0;j<bitlength;j++)	//correlation is approx (XTY/n)^2/sd(Y) - XTY can be obtained from RTdata
{
value=pow(RTdata[(size_t)(nmcmc+m)*bitlength+j]/num_samples_use,2);
if(j==0){value2=value;best=0;}
if(value>value2){value2=value;best=j;}
}

if(value2>maxcor)	//blank the corresponding columns of RTdata and regress top predictor out of Yadj2
{
for(j=0;j<bitlength;j++){RTdata[(size_t)(nmcmc+m)*bitlength+j]=0;}

if(dichot==0){reg_covar_matrix(Yadj2+m*num_samples_use, data+(size_t)best*num_samples_use, num_samples_use, 1, 1);}
else{reg_covar_weighted(Yadj2+m*num_samples_use, data+(size_t)best*num_samples_use, num_samples_use, 1, 1, nullweights+m*num_samples_use);}

if(wcount<5)
{
if(mpheno!=-1){printf("Warning, excluding Chunk %d because it contains a predictor with estimated variance explained %.4f\n", bit+1, value2);}
else{printf("Warning, excluding Chunk %d for Phenotype %d because it contains a predictor with estimated variance explained %.4f\n", bit+1, m+1, value2);}
}
wcount++;

ycounts[m]++;
}
}

//add on contribution to diagonal for kinship mark for each power
#pragma omp parallel for private(k,j,i) schedule(static)
for(k=0;k<num_pows;k++)	
{
for(j=0;j<bitlength;j++)	//no need to worry about trivial predictors - will have exps2 zero
{
for(i=0;i<num_samples_use;i++)
{MkinsD[(size_t)(mark+k*ndivs)*num_samples_use+i]+=pow(data[(size_t)j*num_samples_use+i],2)*exps2[bitstart+j+k*data_length];}
}
}

//add on contributions to trace for kinship mark for each power - MkinsR contains K1R, then K2R, etc
for(k=0;k<num_pows;k++)	
{
//first load data x exps2 into data2
#pragma omp parallel for private(j,i) schedule(static)
for(j=0;j<bitlength;j++)
{
for(i=0;i<num_samples_use;i++)
{data2[(size_t)j*num_samples_use+i]=data[(size_t)j*num_samples_use+i]*exps2[bitstart+j+k*data_length];}
}

//now multiply data2 by t(data)R
alpha=1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total, &bitlength, &alpha, data2, &num_samples_use, RTdata, &bitlength, &beta, MkinsR+(size_t)(mark+k*ndivs)*total*num_samples_use, &num_samples_use);
}

//get cors and datasqs
if(dichot==0)	//using non-weighted regression, so only require one set
{
//have already adjusted for covariates

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+j-bitstart];}
}
else	//using weighted regression, so require one set per phenotype
{
for(m=0;m<num_resps_use;m++)
{
if(num_fixed>1)	//adjust for covariates - its ok that have already adjusted data for covariates
{reg_covar_weighted(data, Z, num_samples_use, bitlength, num_fixed, nullweights+m*num_samples_use);}

//put weighted version of data into data2
copy_matrix(num_samples_use, bitlength, data, data2, 1, nullweights+m*num_samples_use);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data2, &num_samples_use, &beta, cors+(size_t)(bitstart+m*data_length)*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j+m*data_length]=cors[(size_t)(j+m*data_length)*bitsize+j-bitstart];}
}
}
}	//end of bit loop

if(wcount>5){printf("In total, %d chunks were excluded\n", wcount);}
printf("\n");

count=0;for(j=0;j<num_preds_use;j++){count+=(mults[j]==-9999);}
if(count==num_preds_use)
{
if(minmaf==-9999&&maxmaf==-9999&&minvar==-9999&&minobs==-9999&&mininfo==-9999)
{printf("Error, all predictors are trivial (showed no variation)\n\n");}
else
{printf("Error, all predictors failed quality control\n\n");}
exit(1);
}
if(count>0)
{
if(minmaf==-9999&&maxmaf==-9999&&minvar==-9999&&minobs==-9999&&mininfo==-9999)
{printf("Warning, %d predictors were excluded because they were trivial (showed no variation)\n\n", count);}
else
{printf("Warning, %d predictors were excluded because they failed quality control\n\n", count);}
}

////////

//make sure exps2 sum to one
for(k=0;k<num_pows;k++)
{
sum=0;for(j=0;j<data_length;j++){sum+=exps2[j+k*data_length];}
for(j=0;j<data_length;j++){exps2[j+k*data_length]=exps2[j+k*data_length]/sum;}
}

////////

for(k=0;k<num_pows;k++)	//get solutions for kth power
{
//get the average diagonals of the kinship matrices
for(q=0;q<ndivs;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)(q+k*ndivs)*num_samples_use+i];}
kinsums[q]=sum/num_samples_use;
if(kinsums[q]==0){printf("Error, all the predictors in Block %d were excluded, so it is not possible to continue\n\n", q+1);exit(1);}
}

//set KKtraces to minus the cross-product of diagonals
alpha=-1.0;beta=0.0;
dgemm_("T", "N", &ndivs, &ndivs, &num_samples_use, &alpha, MkinsD+(size_t)k*ndivs*num_samples_use, &num_samples_use, MkinsD+(size_t)k*ndivs*num_samples_use, &num_samples_use, &beta, KKtraces, &ndivs);

//add average contribution from vectors
for(g=0;g<nmcmc;g++)
{
token=num_samples_use*total;
alpha=1.0/nmcmc;beta=1.0;
dgemm_("T", "N", &ndivs, &ndivs, &num_samples_use, &alpha, MkinsR+(size_t)(g+k*ndivs*total)*num_samples_use, &token, MkinsR+(size_t)(g+k*ndivs*total)*num_samples_use, &token, &beta, KKtraces, &ndivs);
}

//divide by traces and by two
for(q=0;q<ndivs;q++)
{
for(q2=0;q2<ndivs;q2++){KKtraces[q+q2*ndivs]=KKtraces[q+q2*ndivs]/kinsums[q]/kinsums[q2]/2;}
}

for(m=0;m<num_resps_use;m++)	//get estimate for phenotype m
{
//KYtraces is phenotype terms from MkinsR x Y, minus MkinsD x Y^2
token=num_samples_use*total;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &ndivs, &alpha, MkinsR+(size_t)(nmcmc+m+k*ndivs*total)*num_samples_use, &token, Yadj+m*num_samples_use, &one, &beta, KYtraces, &one);
for(q=0;q<ndivs;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)(q+k*ndivs)*num_samples_use+i]*pow(Yadj[i+m*num_samples_use],2);}
KYtraces[q]-=sum;
}

//divide by traces and by two
for(q=0;q<ndivs;q++){KYtraces[q]=KYtraces[q]/kinsums[q]/2;}

//copy traces and solve
for(q=0;q<ndivs;q++)
{
for(q2=0;q2<ndivs;q2++){KKtraces3[q+q2*ndivs]=KKtraces[q+q2*ndivs];}
KYtraces2[q]=KYtraces[q];
}
(void)eigen_invert(KKtraces3, ndivs, KKtraces2, 1, KYtraces2, 1);

//get heritability estimate
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces2[q];}
cors2[m+m*num_resps_use+k*num_resps_use*num_resps_use]=sum;

//compute variance explained (or something proportional to this)
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces[q]*KYtraces2[q];}
varexp[m+k*num_resps_use]=sum;

for(m2=m+1;m2<num_resps_use;m2++)	//get correlation with phenotype m2
{
//KYtraces is phenotype terms from MkinsR x Y2, minus MkinsD x Y*Y2
token=num_samples_use*total;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &ndivs, &alpha, MkinsR+(size_t)(nmcmc+m+k*ndivs*total)*num_samples_use, &token, Yadj+m2*num_samples_use, &one, &beta, KYtraces, &one);
for(q=0;q<ndivs;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)(q+k*ndivs)*num_samples_use+i]*Yadj[i+m*num_samples_use]*Yadj[i+m2*num_samples_use];}
KYtraces[q]-=sum;
}

//divide by traces and by two
for(q=0;q<ndivs;q++){KYtraces[q]=KYtraces[q]/kinsums[q]/2;}

//copy traces and solve
for(q=0;q<ndivs;q++)
{
for(q2=0;q2<ndivs;q2++){KKtraces3[q+q2*ndivs]=KKtraces[q+q2*ndivs];}
KYtraces2[q]=KYtraces[q];
}
(void)eigen_invert(KKtraces3, ndivs, KKtraces2, 1, KYtraces2, 1);

//get co-heritability estimate and save in upper triangle
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces2[q];}
cors2[m+m2*num_resps_use+k*num_resps_use*num_resps_use]=sum;

//get co-noise and save in lower triangle  - Exp(Y1i Y2i) = K[i,i] * sum + co-noise
sum=0;
for(i=0;i<num_samples_use;i++)
{
sum+=Yadj[i+m*num_samples_use]*Yadj[i+m2*num_samples_use];
for(q=0;q<ndivs;q++){sum-=MkinsD[(size_t)(q+k*ndivs)*num_samples_use+i]/kinsums[q]*KYtraces2[q];}
}
cors2[m2+m*num_resps_use+k*num_resps_use*num_resps_use]=sum/num_samples_use;
}
}
}	//end of k loop

if(multi==0)	//get top power for each phenotype separately
{
for(m=0;m<num_resps_use;m++)
{
for(k=0;k<num_pows;k++)
{
if(k==0){value=varexp[m+k*num_resps_use];Mtops[m]=0;}
if(varexp[m+k*num_resps_use]>value){value=varexp[m+k*num_resps_use];Mtops[m]=k;}
}
}
}
else	//get best power across all phenotypes
{
for(k=0;k<num_pows;k++)
{
sum=0;for(m=0;m<num_resps_use;m++){sum+=varexp[m+k*num_resps_use];}
if(k==0){value=sum;best=0;}
if(sum>value){value=sum;best=k;}
}
for(m=0;m<num_resps_use;m++){Mtops[m]=best;}
}

//get corresponding heritabilities
for(m=0;m<num_resps_use;m++)
{
value=cors2[m+m*num_resps_use+Mtops[m]*num_resps_use*num_resps_use];

if(power!=-9999)
{
if(mpheno!=-1){printf("Estimated heritability is %.4f\n", value);}
else{printf("Phenotype %d: estimated heritability %.4f\n", m+1, value);}

if(value<0.01){printf("Warning, this is very low, so has been increased to 0.01\n");value=0.01;}
if(value>maxher){printf("Warning, this is very high, so has been reduced to %.4f\n", maxher);value=maxher;}
}
else
{
if(mpheno!=-1){printf("Best power is %.4f, estimated heritability is %.4f\n", powers[Mtops[m]], value);}
else{printf("Phenotype %d: best power is %.4f, estimated heritability is %.4f\n", m+1, powers[Mtops[m]], value);}

if(value<0.01){printf("Warning, the latter is very low, so has been increased to 0.01\n");value=0.01;}
if(value>maxher){printf("Warning, the latter is very high, so has been reduced to %.4f\n", maxher);value=maxher;}
}

hers[m]=value;
}
printf("\n");

if(multi==1)
{
//compute Gmat and Emat (using best power)

for(m=0;m<num_resps_use;m++)
{
Gmat[m+m*num_resps_use]=hers[m];
Emat[m+m*num_resps_use]=1-hers[m];
for(m2=m+1;m2<num_resps_use;m2++)	//use mean from corresponding powers, checking correlations are between -1 and 1
{
sum=.5*(cors2[m+m2*num_resps_use+best*num_resps_use*num_resps_use]+cors2[m+m2*num_resps_use+best*num_resps_use*num_resps_use]);
value=sum*pow(hers[m]*hers[m2],-0.5);
if(value>1){value=0.99;}
if(value<-1){value=-0.99;}
Gmat[m+m2*num_resps_use]=value*pow(hers[m]*hers[m2],0.5);
Gmat[m2+m*num_resps_use]=value*pow(hers[m]*hers[m2],0.5);

sum=.5*(cors2[m2+m*num_resps_use+best*num_resps_use*num_resps_use]+cors2[m2+m*num_resps_use+best*num_resps_use*num_resps_use]);
value=sum*pow((1-hers[m])*(1-hers[m2]),-0.5);
if(value>1){value=0.99;}
if(value<-1){value=-0.99;}
Emat[m+m2*num_resps_use]=value*pow((1-hers[m])*(1-hers[m2]),0.5);
Emat[m2+m*num_resps_use]=value*pow((1-hers[m])*(1-hers[m2]),0.5);
}
}

//print and save Gmat and Emat

printf("Here is the estimated genetic covariance matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Gmat[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

printf("Here is the estimated noise covariance matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Emat[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

sprintf(filename3,"%s.transformed.genetic", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Gmat[m+m2*num_resps_use]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.transformed.noise", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Emat[m+m2*num_resps_use]);}
fprintf(output3, "\n");
}
fclose(output3);

////////

//decompose Emat = UEUT
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){Wmat[m+m2*num_resps_use]=Emat[m+m2*num_resps_use];}
}

lwork=-1;
dsyev_("V", "U", &num_resps_use, Wmat, &num_resps_use, Wmat2, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_resps_use, Wmat, &num_resps_use, Wmat2, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);

//get Wmat4=UE^-.5UT
for(m2=0;m2<num_resps_use;m2++)
{
if(Wmat2[m2]>0){value=pow(Wmat2[m2],-0.5);}
else{value=0;}
for(m=0;m<num_resps_use;m++){Wmat3[m+m2*num_resps_use]=Wmat[m+m2*num_resps_use]*value;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat3, &num_resps_use, Wmat, &num_resps_use, &beta, Wmat4, &num_resps_use);

//get Wmat5=UE^.5UT
for(m2=0;m2<num_resps_use;m2++)
{
if(Wmat2[m2]>0){value=pow(Wmat2[m2],0.5);}
else{value=0;}
for(m=0;m<num_resps_use;m++){Wmat3[m+m2*num_resps_use]=Wmat[m+m2*num_resps_use]*value;}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat3, &num_resps_use, Wmat, &num_resps_use, &beta, Wmat5, &num_resps_use);

//decompose Wmat4 Gmat Wmat4
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat4, &num_resps_use, Gmat, &num_resps_use, &beta, Vmat2, &num_resps_use);
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Vmat2, &num_resps_use, Wmat4, &num_resps_use, &beta, Vmat, &num_resps_use);

lwork=-1;
dsyev_("V", "U", &num_resps_use, Vmat, &num_resps_use, Vmat2, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_resps_use, Vmat, &num_resps_use, Vmat2, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);

//Umat is t(eigenvector) from second decomp times Wmat4
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Vmat, &num_resps_use, Wmat4, &num_resps_use, &beta, Umat, &num_resps_use);

//Umat2 (inverse of Umat) is Wmat5 times eigenvector from second decomp
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat5, &num_resps_use, Vmat, &num_resps_use, &beta, Umat2, &num_resps_use);

if(dougvar==1)
{
Umat[0]=1;Umat[1]=1;Umat[2]=1;Umat[3]=-1;
Umat2[0]=.5;Umat2[1]=.5;Umat2[2]=.5;Umat2[3]=-.5;
for(m=0;m<num_resps_use;m++)
{
//for(m2=0;m2<num_resps_use;m2++){Umat[m+m2*num_resps_use]=(m==m2);Umat2[m+m2*num_resps_use]=(m==m2);}
}
}

//compute transformed (adjusted) phenotypes, scaling so have variance one

for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Ymat[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Ymat, &num_samples_use, Umat, &num_resps_use, &beta, Yadj, &num_samples_use);

for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

value=pow(var,-0.5);
value2=pow(var,0.5);
for(m2=0;m2<num_resps_use;m2++){Umat[m+m2*num_resps_use]*=value;}
for(m2=0;m2<num_resps_use;m2++){Umat2[m2+m*num_resps_use]*=value2;}
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]*=value;}
}

printf("Here is the transformation matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Umat[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

sprintf(filename3,"%s.transformed.pheno", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output3,"%s %s", ids1[i], ids2[i]);
for(m=0;m<num_resps_use;m++){fprintf(output3, " %.6f", Yadj[i+m*num_samples_use]);}
fprintf(output3,"\n");
}
fclose(output3);

//compute transformed Yadj2

for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Ymat[i+m*num_samples_use]=Yadj2[i+m*num_samples_use];}
}

alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Ymat, &num_samples_use, Umat, &num_resps_use, &beta, Yadj2, &num_samples_use);

//get transformed Gmat and load up heritabilities

alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Umat, &num_resps_use, Gmat, &num_resps_use, &beta, Wmat4, &num_resps_use);
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat4, &num_resps_use, Umat, &num_resps_use, &beta, Wmat5, &num_resps_use);

for(m=0;m<num_resps_use;m++)
{
value=Wmat5[m+m*num_resps_use];
printf("Transformed phenotype %d: estimated heritability %.4f\n", m+1, value);
if(value<0.01){printf("Warning, this is very low, so has been increased to 0.01\n");value=0.01;}
if(value>maxher){printf("Warning, this is very high, so has been reduced to %.4f\n", maxher);value=maxher;}
hers[m]=value;
}
printf("\n");

printf("Here is the transformed genetic matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Wmat5[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

//add on transformed Emat to get variance matrix

alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Umat, &num_resps_use, Emat, &num_resps_use, &beta, Wmat4, &num_resps_use);
alpha=1.0;beta=1.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat4, &num_resps_use, Umat, &num_resps_use, &beta, Wmat5, &num_resps_use);

printf("Here is the transformed variance matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Wmat5[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

//save transformation and inverse, remembering scaling

sprintf(filename3,"%s.transformed.matrix", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Umat[m+m2*num_resps_use]*Mscales[m2]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.transformed.inverse", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Umat2[m+m2*num_resps_use]/Mscales[m]);}
fprintf(output3, "\n");
}
fclose(output3);
}

//load top exps2 into exps
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<data_length;j++){exps[j+m*data_length]=exps2[j+Mtops[m]*data_length];}
}

if(verbose==1)	//save
{
sprintf(filename3,"%s.hers", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Algorithm Phenotype Power Heritability Chunks_Excluded\n");
for(m=0;m<num_resps_use;m++){fprintf(output3,"RHE %d %.4f %.4f %d\n", m+1, powers[Mtops[m]], hers[m], ycounts[m]);}
fclose(output3);
}

////////

free(kinsums);
free(R);free(RTdata);
free(MkinsD);free(MkinsR);
free(KKtraces);free(KKtraces2);free(KKtraces3);free(KYtraces);free(KYtraces2);
free(exps2);
free(varexp);
free(cors2);
free(ycounts);
if(multi==1){free(Wmat);free(Wmat2);free(Wmat3);free(Wmat4);free(Wmat5);free(Vmat);free(Vmat2);free(Ymat);}

///////////////////////////


/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Estimate hers and power, and get exps and correlations, plus multivariate terms if multi=1

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

//allocate variables

total=nmcmc+num_resps_use;

R=malloc(sizeof(double)*num_samples_use*total);
RTdata=malloc(sizeof(double)*bitsize*total);

kinsums=malloc(sizeof(double)*ndivs);
MkinsD=malloc(sizeof(double)*num_samples_use*ndivs*num_pows);
MkinsR=malloc(sizeof(double)*num_samples_use*total*ndivs*num_pows);

KKtraces=malloc(sizeof(double)*ndivs*ndivs);
KKtraces2=malloc(sizeof(double)*ndivs);
KKtraces3=malloc(sizeof(double)*ndivs*ndivs);
KYtraces=malloc(sizeof(double)*ndivs);
KYtraces2=malloc(sizeof(double)*ndivs);

exps2=malloc(sizeof(double)*data_length*num_pows);

varexp=malloc(sizeof(double)*num_resps_use*num_pows);

cors2=malloc(sizeof(double)*num_resps_use*num_resps_use);
cors3=malloc(sizeof(double)*num_resps_use*num_resps_use);

cohers=malloc(sizeof(double)*num_resps*num_pows);
cohers2=malloc(sizeof(double)*num_resps*num_pows);

ycounts=malloc(sizeof(int)*num_resps_use);
Pscales=malloc(sizeof(double)*num_resps_use);

if(multi==1)
{
Luse=malloc(sizeof(int)*num_resps_use);
Lmat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Lmat2=malloc(sizeof(double)*num_resps_use);
Lmat3=malloc(sizeof(double)*num_resps_use);
Lmat4=malloc(sizeof(double)*num_resps_use);
Gmat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Emat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat2=malloc(sizeof(double)*num_resps_use);
Wmat3=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat4=malloc(sizeof(double)*num_resps_use*num_resps_use);
Wmat5=malloc(sizeof(double)*num_resps_use*num_resps_use);
Ymat=malloc(sizeof(double)*num_samples_use*num_resps_use);
}

//have already filled Z and Yadj (latter has been scaled to have variance one)

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

for(m=0;m<num_resps_use;m++){Pscales[m]=1;}
/*
if(multi==1) //could scale heritabilities by 1/callrate
{
for(m=0;m<num_resps_use;m++){Pscales[m]=pow((double)respcounts[m]/num_samples_use,-1);}
for(m=0;m<num_resps_use;m++){printf("phenotyep %d has callrate %f, scale her by %f\n", m+1, (double)respcounts[m]/num_samples_use, Pscales[m]);}
}
*/

if(multi==1)	//check for linear dependencies
{
//set Luse=1 (will update if phenotype redundant)
for(m=0;m<num_resps_use;m++){Luse[m]=1;}

//set Lmat=0 (will update columns corresponding to redundant phenotypes)
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){Lmat[m+m2*num_resps_use]=0;}
}

flag=0;
for(m=0;m<num_resps_use;m++)
{
//compute variance explained by first count phenotypes
if(m==0){value=0;}
else
{
//get XTX/n (divide by n so diagonal is one)
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &m, &m, &num_samples_use, &alpha, Yadj, &num_samples_use, Yadj, &num_samples_use, &beta, cors2, &m);

//get XTY/n (divide by n so consistent with above)
alpha=1.0/num_samples_use;beta=0.0;
dgemv_("T", &num_samples_use, &m, &alpha, Yadj, &num_samples_use, Yadj+m*num_samples_use, &one, &beta, Lmat3, &one);

//deal with already excluded phenotypes
for(m2=0;m2<m;m2++)
{
if(Luse[m2]==0)
{
for(m3=0;m3<m;m3++){cors2[m3+m2*m]=0;cors2[m2+m3*m]=0;}
cors2[m2+m2*m]=1;
Lmat3[m2]=0;
}
}

for(m2=0;m2<m;m2++){Lmat4[m2]=Lmat3[m2];}
(void)eigen_invert(cors2, m, Lmat2, 1, Lmat3, 1);

value=0;for(m2=0;m2<m;m2++){value+=Lmat3[m2]*Lmat4[m2];}
}

if(value>0.9999)	//dependent, so update Luse and Lmat
{
printf("Warning, Phenotype %d is (almost) perfectly explained by the previous %d phenotypes\n", m+1, m);
Luse[m]=0;
for(m2=0;m2<m;m2++){Lmat[m2+m*num_resps_use]=Lmat3[m2];}
flag=1;
}
}
if(flag==1){printf("\n");}
}

////////

//ready for bit loop - will record number of chunks removed (both total, and per-phenotype)
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

if(strcmp(sampwfile,"blank")!=0)	//multiply by W^1/2
{copy_matrix(num_samples_use, bitlength, data, data, 1, sweights);}

if(num_fixed>1||strcmp(sampwfile,"blank")!=0)	//adjust for covariates - will use non-weighted regression, even if binary
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

//compute t(data) R
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, R, &num_samples_use, &beta, RTdata, &bitlength);

for(m=0;m<num_resps_use;m++)	//look for big effect predictors
{
for(j=0;j<bitlength;j++)	//variance explained is (XTY/n)^2/var(X) - XTY can be obtained from RTdata
{
if(strcmp(sampwfile,"blank")==0)	//var(X) is approximately one
{var=1;}
else	//must compute var(X)
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=data[(size_t)j*num_samples_use+i];sumsq+=pow(data[(size_t)j*num_samples_use+i],2);}
mean=sum/num_samples_use;var=sumsq/num_samples_use-pow(mean,2);
}

value=pow(RTdata[(size_t)(nmcmc+m)*bitlength+j]/num_samples_use,2)/var;
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

count=0;for(j=0;j<data_length;j++){count+=(mults[j]==-9999);}
if(count==data_length)
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

for(k=0;k<num_pows;k++)	//get estimates of heritability (and maybe correlations) for kth power - remember to scale
{
//get the average diagonals of the (unscaled) kinship matrices
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

//invert KKtraces
(void)eigen_invert(KKtraces, ndivs, KKtraces2, -1, KKtraces3, 1);

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

//compute inv KKtraces x KYtraces
alpha=1.0;beta=0.0;
dgemv_("N", &ndivs, &ndivs, &alpha, KKtraces, &ndivs, KYtraces, &one, &beta, KYtraces2, &one);

//get heritability estimates (scaled)
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces2[q];}
cohers[m+k*num_resps_use]=sum*Pscales[m];

//get residual sumsq - note that total sumsq equals ([sum (Yi^2)]^2 - sum(Yi^4)])/2
sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(Yadj[i+m*num_samples_use],2);}
sum2=0;for(i=0;i<num_samples_use;i++){sum2+=pow(Yadj[i+m*num_samples_use],4);}
sumsq=(pow(sum,2)-sum2)/2;
for(q=0;q<ndivs;q++){sumsq-=KYtraces[q]*KYtraces2[q];}

//and variance (scaled)
sum=0;for(q=0;q<ndivs*ndivs;q++){sum+=KKtraces[q];}
var=sum*sumsq/num_samples_use*2/(num_samples_use-1);
cohers2[m+k*num_resps_use]=var*pow(Pscales[m],2);

//compute variance explained (or something proportional to this) - bit redundant, as could use values just above
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces[q]*KYtraces2[q];}
varexp[m+k*num_resps_use]=sum;
}	//end of m loop
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
value=cohers[m+Mtops[m]*num_resps_use];
var=cohers2[m+Mtops[m]*num_resps_use];

if(power!=-9999)
{
if(mpheno!=-1){printf("Estimated heritability is %.4f (SD %.4f)\n", value, pow(var,.5));}
else{printf("Phenotype %d has estimated heritability %.4f (SD %.4f)\n", m+1, value, pow(var,.5));}
}
else
{
if(mpheno!=-1){printf("Best power is %.4f and estimated heritability is %.4f (SD %.4f)\n", powers[Mtops[m]], value, pow(var,.5));}
else{printf("Phenotype %d has best power %.4f and estimated heritability %.4f (SD %.4f)\n", m+1, powers[Mtops[m]], value, pow(var,.5));}
}

if(value<0.01){printf("Warning, the heritability is very low, so has been increased to 0.01\n");value=0.01;}
if(value>maxher){printf("Warning, the heritability is very high, so has been reduced to %.4f\n", maxher);value=maxher;}

hers[m]=value;
}
printf("\n");

//save copy of hers in hers2
for(m=0;m<num_resps_use;m++){hers2[m]=hers[m];}

//load top exps2 into exps, making sure sum to one
for(m=0;m<num_resps_use;m++)
{
sum=0;for(j=0;j<data_length;j++){sum+=exps2[j+Mtops[m]*data_length];}
for(j=0;j<data_length;j++){exps[j+m*data_length]=exps2[j+Mtops[m]*data_length]/sum;}
}

if(verbose==1)	//save heritabilities
{
sprintf(filename3,"%s.hers", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Algorithm Phenotype Power Heritability Chunks_Excluded\n");
for(m=0;m<num_resps_use;m++){fprintf(output3,"RHE %d %.4f %.4f %d\n", m+1, powers[Mtops[m]], hers[m], ycounts[m]);}
fclose(output3);
}

////////

if(multi==1)	//get correlations (using best power)
{
printf("Computing cross-trait correlations\n");

//already have diagonals
for(m=0;m<num_resps_use;m++){cors2[m+m*num_resps_use]=cohers[m+best*num_resps_use];cors3[m+m*num_resps_use]=cohers2[m+best*num_resps_use];}

//recompute kinsums
for(q=0;q<ndivs;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)(q+best*ndivs)*num_samples_use+i];}
kinsums[q]=sum/num_samples_use;
}

//recompute KKtraces, but this time do not divide by two
alpha=-1.0;beta=0.0;
dgemm_("T", "N", &ndivs, &ndivs, &num_samples_use, &alpha, MkinsD+(size_t)best*ndivs*num_samples_use, &num_samples_use, MkinsD+(size_t)best*ndivs*num_samples_use, &num_samples_use, &beta, KKtraces, &ndivs);
for(g=0;g<nmcmc;g++)
{
token=num_samples_use*total;
alpha=1.0/nmcmc;beta=1.0;
dgemm_("T", "N", &ndivs, &ndivs, &num_samples_use, &alpha, MkinsR+(size_t)(g+best*ndivs*total)*num_samples_use, &token, MkinsR+(size_t)(g+best*ndivs*total)*num_samples_use, &token, &beta, KKtraces, &ndivs);
}
for(q=0;q<ndivs;q++)
{
for(q2=0;q2<ndivs;q2++){KKtraces[q+q2*ndivs]=KKtraces[q+q2*ndivs]/kinsums[q]/kinsums[q2];}
}

//invert KKtraces
(void)eigen_invert(KKtraces, ndivs, KKtraces2, -1, KKtraces3, 1);

for(m=0;m<num_resps_use;m++)
{
printf("Computing values for Trait %d\n", m+1);

for(m2=m+1;m2<num_resps_use;m2++)	//get correlation with phenotype m2
{
//KYtraces is phenotype terms from MkinsR x Y2, minus MkinsD x Y*Y2
token=num_samples_use*total;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &ndivs, &alpha, MkinsR+(size_t)(nmcmc+m+best*ndivs*total)*num_samples_use, &token, Yadj+m2*num_samples_use, &one, &beta, KYtraces, &one);
for(q=0;q<ndivs;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)(q+best*ndivs)*num_samples_use+i]*Yadj[i+m*num_samples_use]*Yadj[i+m2*num_samples_use];}
KYtraces[q]-=sum;
}

//divide by traces (but not by two)
for(q=0;q<ndivs;q++){KYtraces[q]=KYtraces[q]/kinsums[q];}

//compute inv KKtraces x KYtraces
alpha=1.0;beta=0.0;
dgemv_("N", &ndivs, &ndivs, &alpha, KKtraces, &ndivs, KYtraces, &one, &beta, KYtraces2, &one);

//get co-heritability estimate (scaled) and save in upper triangle of cor2
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces2[q];}
cors2[m+m2*num_resps_use]=sum*pow(Pscales[m]*Pscales[m2],.5);

//get residual sumsq (note that total sumsq equals [[sum (Y_i^2)]+[sum(Y'_i^2)] - [sum (Y_i^2)]^2 - [sum (Y'_i^2)]^2 - sum(Y_i^2 Y'_i^2)])
sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(Yadj[i+m*num_samples_use],2);}
sum2=0;for(i=0;i<num_samples_use;i++){sum2+=pow(Yadj[i+m2*num_samples_use],2);}
sum3=0;for(i=0;i<num_samples_use;i++){sum3+=pow(Yadj[i+m*num_samples_use]*Yadj[i+m2*num_samples_use],2);}
sumsq=sum*sum2-sum3;
for(q=0;q<ndivs;q++){sumsq-=KYtraces[q]*KYtraces2[q];}

//and variance (scaled)
sum=0;for(q=0;q<ndivs*ndivs;q++){sum+=KKtraces[q];}
var=sum*sumsq/num_samples_use/(num_samples_use-1);
cors3[m+m2*num_resps_use]=var*Pscales[m]*Pscales[m2];

//get co-noise and save in lower triangle  - Exp(Y1i Y2i) = K[i,i] * sum cohers + co-noise - will scale Exp(Y1i Y2i)
sum=ddot_(&num_samples_use, Yadj+m*num_samples_use, &one, Yadj+m2*num_samples_use, &one)/num_samples_use*pow(Pscales[m]*Pscales[m2],.5);
for(q=0;q<ndivs;q++){sum-=KYtraces2[q];}
cors2[m2+m*num_resps_use]=sum*pow(Pscales[m]*Pscales[m2],.5);

//set variance to same as that of co-heritability (assumes trait-trait covariance estimated perfectly)
cors3[m2+m*num_resps_use]=cors3[m+m2*num_resps_use];
}
}

//compute Gmat, Emat, Umat and Umat2
for(m=0;m<num_resps_use;m++)
{
Gmat[m+m*num_resps_use]=hers[m];
Emat[m+m*num_resps_use]=1-hers[m];
for(m2=m+1;m2<num_resps_use;m2++)	//fill values (will ensure correlations are between -1 and 1)
{
value=cors2[m+m2*num_resps_use]*pow(hers[m]*hers[m2],-0.5);
//if(value>1){value=1;}
//if(value<-1){value=-1;}
Gmat[m+m2*num_resps_use]=value*pow(hers[m]*hers[m2],0.5);
Gmat[m2+m*num_resps_use]=value*pow(hers[m]*hers[m2],0.5);

value=cors2[m2+m*num_resps_use]*pow((1-hers[m])*(1-hers[m2]),-0.5);
//if(value>1){value=1;}
//if(value<-1){value=-1;}
Emat[m+m2*num_resps_use]=value*pow((1-hers[m])*(1-hers[m2]),0.5);
Emat[m2+m*num_resps_use]=value*pow((1-hers[m])*(1-hers[m2]),0.5);

if(dougvar==1){Emat[m+m2*num_resps_use]=0;Emat[m2+m*num_resps_use]=0;}
}
}

printf("Here is the estimated correlation matrix (lower-triangle genetic, upper-triangle environmental):\n");
for(m=0;m<num_resps_use;m++){printf("\tP%d", m+1);}
printf("\n");
for(m=0;m<num_resps_use;m++)
{
printf("P_%d", m+1);
for(m2=0;m2<num_resps_use;m2++)
{
if(m==m2){printf("\t1");}
if(m>m2){printf("\t%.4f", Gmat[m+m2*num_resps_use]*pow(Gmat[m+m*num_resps_use]*Gmat[m2+m2*num_resps_use],-.5));}
if(m<m2){printf("\t%.4f", Emat[m+m2*num_resps_use]*pow(Emat[m+m*num_resps_use]*Emat[m2+m2*num_resps_use],-.5));}
}
printf("\n");
}
printf("\n");

sprintf(filename3,"%s.multivariate.correlations", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++)
{
if(m==m2){fprintf(output3, "1 ");}
if(m>m2){fprintf(output3, "%.4f ", Gmat[m+m2*num_resps_use]*pow(Gmat[m+m*num_resps_use]*Gmat[m2+m2*num_resps_use],-.5));}
if(m<m2){fprintf(output3, "%.4f ", Emat[m+m2*num_resps_use]*pow(Emat[m+m*num_resps_use]*Emat[m2+m2*num_resps_use],-.5));}
}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.multivariate.gencov", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++)
{fprintf(output3, "%.4f ", Gmat[m+m2*num_resps_use]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.multivariate.envcov", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++)
{fprintf(output3, "%.4f ", Emat[m+m2*num_resps_use]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.multivariate.covar", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m2=0;m2<num_resps_use;m2++)
{
for(m=0;m<num_resps_use;m++)
{fprintf(output3, "%.4f ", cors2[m+m2*num_resps_use]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.multivariate.se", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m2=0;m2<num_resps_use;m2++)
{
for(m=0;m<num_resps_use;m++)
{fprintf(output3, "%.4f ", pow(cors3[m+m2*num_resps_use],.5));}
fprintf(output3, "\n");
}
fclose(output3);

//deal with excluded phenotypes
for(m=0;m<num_resps_use;m++)
{
if(Luse[m]==0)	//blank rows and columns of Gmat and Emat
{
for(m2=0;m2<num_resps_use;m2++){Gmat[m+m2*num_resps_use]=0;Gmat[m2+m*num_resps_use]=0;}
Gmat[m+m*num_resps_use]=1;
for(m2=0;m2<num_resps_use;m2++){Emat[m+m2*num_resps_use]=0;Emat[m2+m*num_resps_use]=0;}
Emat[m+m*num_resps_use]=1;
}
}

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

//make sure all eigenvalues non-trivial (unlikely to be any trivial)
for(m=0;m<num_resps_use;m++)
{
if(Wmat2[m]<1e-10){Wmat2[m]=1e-10;}
}

//get Wmat4=UE^-.5UT
for(m2=0;m2<num_resps_use;m2++)
{
value=pow(Wmat2[m2],-0.5);
for(m=0;m<num_resps_use;m++){Wmat3[m+m2*num_resps_use]=Wmat[m+m2*num_resps_use]*value;}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat3, &num_resps_use, Wmat, &num_resps_use, &beta, Wmat4, &num_resps_use);

//and also Wmat5=UE^.5UT
for(m2=0;m2<num_resps_use;m2++)
{
value=pow(Wmat2[m2],0.5);
for(m=0;m<num_resps_use;m++){Wmat3[m+m2*num_resps_use]=Wmat[m+m2*num_resps_use]*value;}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat3, &num_resps_use, Wmat, &num_resps_use, &beta, Wmat5, &num_resps_use);

//decompose Wmat4 Gmat Wmat4
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat4, &num_resps_use, Gmat, &num_resps_use, &beta, Wmat3, &num_resps_use);
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat3, &num_resps_use, Wmat4, &num_resps_use, &beta, Wmat, &num_resps_use);

lwork=-1;
dsyev_("V", "U", &num_resps_use, Wmat, &num_resps_use, Wmat2, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;
work=malloc(sizeof(double)*lwork);
dsyev_("V", "U", &num_resps_use, Wmat, &num_resps_use, Wmat2, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug info %d\n\n", info);exit(1);}
free(work);

//turn around Wmat, so that start with most important eigenvectors
copy_matrix(num_resps_use, num_resps_use, Wmat, Wmat3, 0, NULL);
for(m2=0;m2<num_resps_use;m2++)
{
for(m=0;m<num_resps_use;m++){Wmat[m+m2*num_resps_use]=Wmat3[m+(num_resps_use-m2-1)*num_resps_use];}
}

//also turn around Wmat2
for(m=0;m<num_resps_use;m++){Wmat3[m]=Wmat2[m];}
for(m=0;m<num_resps_use;m++){Wmat2[m]=Wmat3[num_resps_use-m-1];}

if(mincor!=-9999)	//delete eigenvectors with small eigenvalues
{
sum=0;for(m=0;m<num_resps_use;m++){sum+=Wmat2[m];}
sum2=0;
for(m=0;m<num_resps_use;m++)
{
sum2+=Wmat2[m]/sum;
if(sum2<mincor)	//blank
{
printf("blank %d\n", m+1);
for(m2=0;m2<num_resps_use;m2++){Wmat[m2+m*num_resps_use]=0;}
}
}
}

//Umat is t(eigenvector) from second decomp times Wmat4
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat, &num_resps_use, Wmat4, &num_resps_use, &beta, Umat, &num_resps_use);

//Umat2 is Wmat5 times eigenvector from second decomp
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat5, &num_resps_use, Wmat, &num_resps_use, &beta, Umat2, &num_resps_use);

//deal with excluded phenotypes

for(m=0;m<num_resps_use;m++)
{
if(Luse[m]==0)	//blank contribution of phenotype to Umat
{
for(m2=0;m2<num_resps_use;m2++){Umat[m2+m*num_resps_use]=0;}
}
}

for(m=0;m<num_resps_use;m++)
{
if(Luse[m]==0)	//column m of Lmat says how phenotype depends on those that remain
{
//blank inverse row
for(m2=0;m2<num_resps_use;m2++){Umat2[m+m2*num_resps_use]=0;}

//add back on weighted sum of rows
for(m3=0;m3<num_resps_use;m3++)
{
value=Lmat[m3+m*num_resps_use];
for(m2=0;m2<num_resps_use;m2++){Umat2[m+m2*num_resps_use]+=value*Umat2[m3+m2*num_resps_use];}
}
}
}

//compute (provisional) transformed Yadj
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Yadj, &num_samples_use, Umat, &num_resps_use, &beta, Ymat, &num_samples_use);

//work out scalings that ensure transformed Ymat has variance one - use hers to indicate whether a transformation is trivial
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=Ymat[i+m*num_samples_use];sumsq+=pow(Ymat[i+m*num_samples_use],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

if(var>1e-6)
{
value=pow(var,-0.5);
value2=pow(var,0.5);
for(m2=0;m2<num_resps_use;m2++){Umat[m+m2*num_resps_use]*=value;}
for(m2=0;m2<num_resps_use;m2++){Umat2[m2+m*num_resps_use]*=value2;}
hers[m]=1;
}
else{hers[m]=0;}
}

//ensure redundant transformations are at the end
for(m=0;m<num_resps_use;m++)
{
Wmat2[m]=hers[m];
for(m2=0;m2<num_resps_use;m2++){Wmat3[m+m2*num_resps_use]=Umat[m+m2*num_resps_use];}
for(m2=0;m2<num_resps_use;m2++){Wmat4[m2+m*num_resps_use]=Umat2[m2+m*num_resps_use];}
}

//first load up non-redundant transforms, then redundant transforms
count=0;
for(m=0;m<num_resps_use;m++)
{
if(Wmat2[m]==1)
{
hers[count]=1;
for(m2=0;m2<num_resps_use;m2++){Umat[count+m2*num_resps_use]=Wmat3[m+m2*num_resps_use];}
for(m2=0;m2<num_resps_use;m2++){Umat2[m2+count*num_resps_use]=Wmat4[m2+m*num_resps_use];}
count++;
}
}
for(m=0;m<num_resps_use;m++)
{
if(Wmat2[m]==0)
{
hers[count]=0;
for(m2=0;m2<num_resps_use;m2++){Umat[count+m2*num_resps_use]=Wmat3[m+m2*num_resps_use];}
for(m2=0;m2<num_resps_use;m2++){Umat2[m2+count*num_resps_use]=Wmat4[m2+m*num_resps_use];}
count++;
}
}

//jasper would like adjusted Y
sprintf(filename3,"%s.multivariate.residuals", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
fprintf(output3, "FID IID");
for(m=0;m<num_resps_use;m++){fprintf(output3, " Residual_%d", m+1);}
fprintf(output3,"\n");
for(i=0;i<num_samples_use;i++)
{
fprintf(output3,"%s %s", ids1[i], ids2[i]);
for(m=0;m<num_resps_use;m++){fprintf(output3, " %.6f", Yadj[i+m*num_samples_use]);}
fprintf(output3,"\n");
}
fclose(output3);

//compute (final) transformed Yadj
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Ymat[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Ymat, &num_samples_use, Umat, &num_resps_use, &beta, Yadj, &num_samples_use);

//compute transformed Yadj2
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Ymat[i+m*num_samples_use]=Yadj2[i+m*num_samples_use];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Ymat, &num_samples_use, Umat, &num_resps_use, &beta, Yadj2, &num_samples_use);

printf("Here is the transformation matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Umat[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

sprintf(filename3,"%s.multivariate.transformation", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Umat[m+m2*num_resps_use]*Mscales[m2]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.multivariate.inverse", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Umat2[m+m2*num_resps_use]/Mscales[m]);}
fprintf(output3, "\n");
}
fclose(output3);

sprintf(filename3,"%s.multivariate.pheno", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
fprintf(output3, "FID IID");
for(m=0;m<num_resps_use;m++){fprintf(output3, " Transform_%d", m+1);}
fprintf(output3,"\n");
for(i=0;i<num_samples_use;i++)
{
fprintf(output3,"%s %s", ids1[i], ids2[i]);
for(m=0;m<num_resps_use;m++){fprintf(output3, " %.6f", Yadj[i+m*num_samples_use]);}
fprintf(output3,"\n");
}
fclose(output3);

//get transformed Gmat and load up heritabilities

alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Umat, &num_resps_use, Gmat, &num_resps_use, &beta, Wmat4, &num_resps_use);
dgemm_("N", "T", &num_resps_use, &num_resps_use, &num_resps_use, &alpha, Wmat4, &num_resps_use, Umat, &num_resps_use, &beta, Wmat5, &num_resps_use);

for(m=0;m<num_resps_use;m++)
{
if(hers[m]==1)
{
value=Wmat5[m+m*num_resps_use];
printf("Transformed phenotype %d has estimated heritability %.4f\n", m+1, value);
if(value<0.01){printf("Warning, this is very low, so has been increased to 0.01\n");value=0.01;}
if(value>maxher){printf("Warning, this is very high, so has been reduced to %.4f\n", maxher);value=maxher;}
hers[m]=value;
}
else{printf("Transformed phenotype %d is redundant (due to linear dependency)\n", m+1);}
}
printf("\n");

if(verbose==1)	//append transformed heritabilities
{
sprintf(filename3,"%s.hers", outfile);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++){fprintf(output3,"Transformation %d %.4f %.4f %d\n", m+1, powers[Mtops[m]], hers[m], ycounts[m]);}
fclose(output3);
}
}	//end of multi=1

////////

free(R);free(RTdata);
free(kinsums);free(MkinsD);free(MkinsR);
free(KKtraces);free(KKtraces2);free(KKtraces3);free(KYtraces);free(KYtraces2);
free(exps2);
free(varexp);
free(cors2);free(cors3);
free(cohers);free(cohers2);
free(ycounts);free(Pscales);
if(multi==1){free(Luse);free(Lmat);free(Lmat2);free(Lmat3);free(Lmat4);free(Gmat);free(Emat);free(Wmat);free(Wmat2);free(Wmat3);free(Wmat4);free(Wmat5);free(Ymat);}

///////////////////////////


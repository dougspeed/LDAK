/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Light version of getindhersb.c that only decides power

///////////////////////////

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

data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

total=nmcmc+num_resps_use;

kinsums=malloc(sizeof(double)*ndivs);

R=malloc(sizeof(double)*num_samples_use*total);
RTdata=malloc(sizeof(double)*bitsize*total);
MkinsD=malloc(sizeof(double)*num_samples_use*ndivs*num_pows);
MkinsR=malloc(sizeof(double)*num_samples_use*total*ndivs*num_pows);

KKtraces=malloc(sizeof(double)*ndivs*ndivs);
KKtraces2=malloc(sizeof(double)*ndivs);
KKtraces3=malloc(sizeof(double)*ndivs*ndivs);
KYtraces=malloc(sizeof(double)*ndivs*num_resps_use*num_pows);
KYtraces2=malloc(sizeof(double)*ndivs*num_resps_use*num_pows);

exps2=malloc(sizeof(double)*data_length*num_pows);

varexp=malloc(sizeof(double)*num_resps_use*num_pows);

ycounts=malloc(sizeof(int)*num_resps_use);

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
if(dtype==1&&dougvar==0)	//fast way
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

if(value2>maxcor)	//blank the column of RTdata
{
for(j=0;j<bitlength;j++){RTdata[(size_t)(nmcmc+m)*bitlength+j]=0;}

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
if(kinsums[q]==0){printf("Error, all the predictors in Block %d are trivial, so it is not possible to continue\n\n", q+1);exit(1);}
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

//KYtraces is phenotype terms from MkinsR x Y, minus MkinsD x Y^2
for(m=0;m<num_resps_use;m++)
{
token=num_samples_use*total;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &ndivs, &alpha, MkinsR+(size_t)(nmcmc+m+k*ndivs*total)*num_samples_use, &token, Yadj+m*num_samples_use, &one, &beta, KYtraces+(m+k*num_resps_use)*ndivs, &one);
for(q=0;q<ndivs;q++)
{
sum=0;for(i=0;i<num_samples_use;i++){sum+=MkinsD[(size_t)(q+k*ndivs)*num_samples_use+i]*pow(Yadj[i+m*num_samples_use],2);}
KYtraces[q+(m+k*num_resps_use)*ndivs]-=sum;
}
}

//divide by traces and by two
for(q=0;q<ndivs;q++)
{
for(q2=0;q2<ndivs;q2++){KKtraces[q+q2*ndivs]=KKtraces[q+q2*ndivs]/kinsums[q]/kinsums[q2]/2;}
for(m=0;m<num_resps_use;m++){KYtraces[q+(m+k*num_resps_use)*ndivs]=KYtraces[q+(m+k*num_resps_use)*ndivs]/kinsums[q]/2;}
}

for(m=0;m<num_resps_use;m++)	//get estimate for phenotype m
{
//copy traces
for(q=0;q<ndivs;q++)
{
for(q2=0;q2<ndivs;q2++){KKtraces3[q+q2*ndivs]=KKtraces[q+q2*ndivs];}
KYtraces2[q+(m+k*num_resps_use)*ndivs]=KYtraces[q+(m+k*num_resps_use)*ndivs];
}

//get estimates, and compute variance explained (or something proportional to this)
(void)eigen_invert(KKtraces3, ndivs, KKtraces2, 1, KYtraces2+(m+k*num_resps_use)*ndivs, 1);
sum=0;for(q=0;q<ndivs;q++){sum+=KYtraces[q+(m+k*num_resps_use)*ndivs]*KYtraces2[q+(m+k*num_resps_use)*ndivs];}
varexp[m+k*num_resps_use]=sum;
}
}	//end of k loop

for(m=0;m<num_resps_use;m++)	//get top power
{
for(k=0;k<num_pows;k++)
{
if(k==0){value=varexp[m];pedtops[m]=0;;}
if(varexp[m+k*num_resps_use]>value){value=varexp[m+k*num_resps_use];pedtops[m]=k;}
}

if(mpheno!=-1){printf("Best power is %.4f\n", powers[pedtops[m]]);}
else{printf("Phenotype %d: best power is %.4f\n", m+1, powers[pedtops[m]]);}
}
printf("\n");

free(data);free(data2);
free(kinsums);
free(R);free(RTdata);
free(MkinsD);free(MkinsR);
free(KKtraces);free(KKtraces2);free(KKtraces3);free(KYtraces);free(KYtraces2);
free(exps2);
free(varexp);
free(ycounts);

///////////////////////////


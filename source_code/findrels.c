/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Search data for relatives, then make pedigree PRS

///////////////////////////

//work out max number of vectors required

total=(2+ncal)*num_resps_use;
if(dichot==0)
{
if(nmcmc+num_resps_use>total){total=nmcmc+num_resps_use;}
}
else
{
if((1+nmcmc)*num_resps_use>total){total=(1+nmcmc)*num_resps_use;}
}

printf("Testing for structure using %d pedigree predictors\n", nped);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Testing for structure using %d pedigree predictors\n", nped);
fclose(output);

//work out memory

value=(double)num_samples_use/1024/1024/1024*8*(nped+total)+nped/1024/1024/1024*8*nped;
if(value>.001){printf("This usually takes at most a few minutes, and requires approximately %.1f Gb (if this is very high, consider using \"--num-pedigree-predictors\" to reduce the number of pedigree predictors\n\n", value);}
else{printf("This usually takes at most a few minutes, and requires less than 1 Gb\n\n");}

//set heritabilities

tryhers=malloc(sizeof(double)*num_hers1);
value=pow(maxher/0.001,1.0/(num_hers1-1));
for(j=0;j<num_hers1;j++){tryhers[j]=0.001*pow(value,j);}

//allocate (have already set total)

pedtops=malloc(sizeof(int)*num_resps_use);
pedscales=malloc(sizeof(double)*num_resps_use);
pedhers=malloc(sizeof(double*)*num_resps_use);
pedprs=malloc(sizeof(double)*num_samples_use*num_resps_use);
pedgammas=malloc(sizeof(double*)*num_resps_use);
pedsds=malloc(sizeof(double*)*num_resps_use);
pedmse=malloc(sizeof(double*)*num_resps_use);

dindex=malloc(sizeof(int)*data_length);
dindex2=malloc(sizeof(int)*data_length);
ddata=malloc(sizeof(double)*num_samples_use*nped);
dcentres=malloc(sizeof(double)*data_length);
dmults=malloc(sizeof(double)*data_length);
dsqdevs=malloc(sizeof(double)*data_length);
dpreds=malloc(sizeof(char*)*data_length);

XTCX=malloc(sizeof(double)*nped*nped);
XTCX2=malloc(sizeof(double)*num_samples_use);
kins=malloc(sizeof(double)*num_samples_use*bitsize);

firsts=malloc(sizeof(int)*maxpairs);
seconds=malloc(sizeof(int)*maxpairs);
kinships=malloc(sizeof(double)*maxpairs);

cX=malloc(sizeof(double)*num_samples_use*total);
cR=malloc(sizeof(double)*num_samples_use*total);

polates=malloc(sizeof(double)*num_hers1*num_resps_use);
polates2=malloc(sizeof(double)*num_hers1);

////////

//first compute stats for 20 x nped predictors

count=20*nped;
if(count<5000){count=5000;}
if(count>data_length){count=data_length;}

for(j=0;j<data_length;j++){dindex[j]=j;}
permute_int(dindex,data_length);
qsort(dindex,count,sizeof(int), compare_int);

for(j=0;j<count;j++){dindex2[j]=keeppreds_use[dindex[j]];}
for(j=0;j<count;j++){copy_string(dpreds,j,preds[dindex[j]]);}

printf("Calculating statistics for %d randomly-picked predictors (in order to find %d with relatively high variance)\n\n", count, nped);

bittotal2=(count-1)/nped+1;

for(bit=0;bit<bittotal2;bit++)
{
bitstart=bit*nped;
bitend=(bit+1)*nped;
if(bitend>count){bitend=count;}
bitlength=bitend-bitstart;

(void)read_data_fly(datafile, dtype, ddata, NULL, num_samples_use, keepsamps, bitstart, bitend, dindex2, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(ddata, dcentres+bitstart, dmults+bitstart, dsqdevs+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 0, dpreds+bitstart);
}

//sort based on variance, then get top nped
dptrs=malloc(sizeof(struct sorting_double)*count);

for(j=0;j<count;j++){dptrs[j].value=dsqdevs[j];dptrs[j].index=dindex[j];}
qsort(dptrs, count, sizeof(struct sorting_double), compare_sorting_double_rev);
for(j=0;j<nped;j++){dindex[j]=dptrs[j].index;}
qsort(dindex,nped,sizeof(int), compare_int);

for(j=0;j<count;j++){free(dpreds[j]);}
free(dptrs);

if(chr[dindex[nped-1]]==chr[dindex[0]]){printf("Error, all pedigree predictors are on Chromosome %d (this indicates that very few predictors are not on this chromosome)\n\n", chr[dindex[0]]);exit(1);}

for(j=0;j<nped;j++){dindex2[j]=keeppreds_use[dindex[j]];}
for(j=0;j<nped;j++){copy_string(dpreds,j,preds[dindex[j]]);}

//read pedigree predictors and standardize (and replace NAs with zero)
(void)read_data_fly(datafile, dtype, ddata, NULL, num_samples_use, keepsamps, 0, nped, dindex2, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(ddata, dcentres, dmults, dsqdevs, num_samples_use, nped, missingvalue, -1, 0, 0, NULL, 1, dpreds);

count=0;
for(j=0;j<nped;j++){count+=(dmults[j]==-9999);}
if(count==nped){printf("Error, all %d pedigree predictors are trivial\n\n", nped);}
if(count>0){printf("Warning, %d of the %d pedigree predictors are trivial\n\n", count, nped);}

if(num_fixed>1)	//adjust predictors for covariates (always use non-weighted residuals) then restandardize
{
reg_covar_matrix(ddata, Z, num_samples_use, nped, num_fixed);
stand_matrix_nomiss(ddata, num_samples_use, num_samples_use, nped);
}

//divide predictors by root nped-count (so that kinship matrix is XXT)
value=pow(nped-count,-.5);
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<nped;j++)
{
for(i=0;i<num_samples_use;i++){ddata[i+j*num_samples_use]*=value;}
}

//get t(X) X (which is num_samples_use / (nped-count) times correlation matrix)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &nped, &nped, &num_samples_use, &alpha, ddata, &num_samples_use, ddata, &num_samples_use, &beta, XTCX, &nped);

//test if significant structure (using cross-chromosome correlations)
sumsq=0;
count=0;
value=(double)(nped-count)/num_samples_use;
for(j2=0;j2<nped;j2++)
{
if(dmults[j2]!=-9999)
{
for(j=j2+1;j<nped;j++)
{
if(dmults[j]!=-9999)
{
if(chr[dindex[j]]!=chr[dindex[j2]]){sumsq+=pow(XTCX[j+j2*nped]*value,2);count++;}
}}
}}

if(count==0){printf("Error, it was not possible to test for significant structure, because there are no non-trivial cross-chromosome pairs of pedigree predictors\n\n");exit(1);}

//work out average squared correlation, then compare to mean under the null distribution
value=sumsq/count;
mean=1.0/(num_samples_use-1);
var=2*(num_samples_use-2)*pow(num_samples_use-1,-2)/(num_samples_use+1)/count;
value2=pow(value-mean,2)/var;
value3=erfc(pow(value2,.5)*M_SQRT1_2);
value4=(value-mean)*num_samples_use;

printf("Average cross-chromosome squared correlation is %.4e (chi-squared test statistic %.4f; p-value %.4e); estimated maximum average inflation of test statistics is %.4f\n\n", value, value2, value3, value4);

//see whether using pedigree PRS (or forcing ridge) - note that useped is currently set to zero
if(fastgwa==0)
{
if(value3>1e-4||value4<0.1){printf("This is relatively weak structure, so will continue as normal\n\n");}
else
{
if(dichot==0){printf("This is relatively strong structure, so will switch to ridge regression (you can avoid this by adding \"--check-pedigree NO\")\n\n");forceridge=1;}
else{printf("This is relatively strong structure, so will use a pedigree-based PRS\n\n");useped=1;}
}
}
else
{
if(value3>1e-4||value4<0.1){printf("This is relatively weak structure\n\n");}
else{printf("This is relatively strong structure\n\n");}
useped=1;
}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(useped==0){fprintf(output,"Will not construct a pedigree PRS\n");}
else{fprintf(output,"Will construct a pedigree PRS\n");}
fclose(output);

////////

if(useped==1)	//construct the pedigree prs - if dichot=1, will use alt logistic
{
//get diagonal of kinship matrix
for(i=0;i<num_samples_use;i++){XTCX2[i]=0;}
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<nped;j++)
{
for(i=0;i<num_samples_use;i++){XTCX2[i]+=pow(ddata[i+j*num_samples_use],2);}
}

//get sum of off-diagonal kinships = minus sum of diagonal kinships (sum of all kinships is zero)
sum=0;for(i=0;i<num_samples_use;i++){sum-=XTCX2[i];}

//get sum of squares of off-diagonal kinships = trace (XTX XTX) - sum of square of diagonal kinships
sumsq=0;
for(j2=0;j2<nped;j2++)
{
for(j=0;j<nped;j++){sumsq+=pow(XTCX[j+j2*nped],2);}
}
for(i=0;i<num_samples_use;i++){sumsq-=pow(XTCX2[i],2);}

//get mean and variance
mean=sum/num_samples_use/(num_samples_use-1);
var=sumsq/num_samples_use/(num_samples_use-1)-pow(mean,2);

//find cutoff so that will obtain about 0.1 n pairs by chance, or 0.05 (whichever higher)
value=mean-normal_inv(0.1/num_samples_use)*pow(var,.5);
if(value<0.05){value=0.05;}

printf("The off-diagonal kinships have mean %.4f and SD %.4f; will record sample pairs with estimated kinship above %.4f\n\n", mean, pow(var,.5), value);

bittotal2=(num_samples_use-1)/bitsize+1;

num_rels=0;
count=0;
for(bit=0;bit<bittotal2;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_samples_use){bitend=num_samples_use;}
bitlength=bitend-bitstart;

if(bit%200==0)
{
printf("Finding relatives for Chunk %d of %d\n", bit+1, bittotal2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Finding relatives for Chunk %d of %d\n", bit+1, bittotal2);
fclose(output);
}

//compute bitlength columns of kinship matrix
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &bitlength, &nped, &alpha, ddata, &num_samples_use, ddata+bitstart, &num_samples_use, &beta, kins, &num_samples_use);

for(i2=bitstart;i2<bitend;i2++)
{
firsts[num_rels]=i2;seconds[num_rels]=i2;
kinships[num_rels]=kins[(size_t)i2+(i2-bitstart)*num_samples_use];
num_rels++;

for(i=i2+1;i<num_samples_use;i++)
{
if(kins[(size_t)i+(i2-bitstart)*num_samples_use]>value)
{
firsts[num_rels]=i;seconds[num_rels]=i2;
kinships[num_rels]=kins[(size_t)i+(i2-bitstart)*num_samples_use];
num_rels++;
count++;
}

if(num_rels==maxpairs){break;}
}}
}	//end of bit loop

if(num_rels==maxpairs)
{
printf("Error, there are over %d significantly-related pairs; it will not be possible to continue (this has never happened before, so please can you tell Doug)\n\n", maxpairs);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Error, there are over %d significantly-related pairs; it will not be possible to continue (this has never happened before, so please tell Doug)\n\n", maxpairs);
fclose(output);

exit(1);
}

////////

if(count==0)	//no related pairs - set pedigree prs to zero
{
printf("Warning, there are no significantly-related pairs, the pedigree-based PRS will be set to zero\n\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Warning, there are no significantly-related pairs, the pedigree-based PRS will be set to zero\n");
fclose(output);

for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){pedprs[i+m*num_samples_use]=0;}
pedgammas[m]=1;
pedsds[m]=0;
}
}
else	//will make pedigree prs
{
printf("In total, there are %d related pairs (each sample has on average %.2f relatives)\n\n", count, 2.0*count/num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"In total, there are %d related pairs (each sample has on average %.2f relatives)\n", count, 2.0*count/num_samples_use);
fclose(output);

/*
sprintf(filename3,"%s.pairs.grm.id", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(i=0;i<num_samples_use;i++){fprintf(output3,"%s %s %s %s\n", ids1[i], ids2[i], ids2[i], ids2[i]);}
fclose(output3);

sprintf(filename3,"%s.pairs.grm.sp", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(i=0;i<num_rels;i++){fprintf(output3,"%d %d %f\n", firsts[i], seconds[i], kinships[i]);}
fclose(output3);

sprintf(filename3,"%s.pairs.grm", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
for(i=0;i<num_rels;i++){fprintf(output3,"%d %d %d %f\n", firsts[i]+1, seconds[i]+1, num_samples_use, kinships[i]);}
fclose(output3);

sprintf(cmd, "gzip -f %s.pairs.grm", outfile);
system(cmd);
*/

//ready to make the PRS (note that after this, we exit, so ok to modify Yadj and nullweights)

printf("Constructing the pedigree PRS\n\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Constructing the pedigree PRS\n");
fclose(output);

if(dichot==0)	//convenient to use scaled phenotypes
{
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);

pedscales[m]=pow(var,-0.5);
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=(Yadj[i+m*num_samples_use]-mean)*pedscales[m];}
}
}
else	//use unscaled phenotypes
{
for(m=0;m<num_resps_use;m++){pedscales[m]=1.0;}
}

//estimate heritability (have already set tryhers)

if(dichot==0){total=nmcmc+num_resps_use;}
else{total=(1+nmcmc)*num_resps_use;}

for(j=0;j<num_hers1;j++)	//evaluate ratio for jth heritability
{
if(dichot==0)	//only need one set of random vectors
{
//load up mcmc random vectors
for(g=0;g<nmcmc;g++)
{
for(i=0;i<num_samples_use;i++){cX[(size_t)g*num_samples_use+i]=0;}
for(i=0;i<num_samples_use;i++){cR[(size_t)g*num_samples_use+i]=gaussian[i+g*num_samples_use];}
}

//then the responses
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){cX[(size_t)(nmcmc+m)*num_samples_use+i]=0;}
for(i=0;i<num_samples_use;i++){cR[(size_t)(nmcmc+m)*num_samples_use+i]=Yadj[i+m*num_samples_use];}
}
}
else	//need one set of random vectors for each phenotype
{
for(m=0;m<num_resps_use;m++)
{
//load up mcmc random vectors
for(g=0;g<nmcmc;g++)
{
p=g+(nmcmc+1)*m;
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]=0;}
for(i=0;i<num_samples_use;i++){cR[(size_t)p*num_samples_use+i]=gaussian[i+g*num_samples_use];}
}

//then the response
p=nmcmc+(nmcmc+1)*m;
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]=0;}
for(i=0;i<num_samples_use;i++){cR[(size_t)p*num_samples_use+i]=Yadj[i+m*num_samples_use];}
}
}

sparse_cgd(num_samples_use, total, num_rels, firsts, seconds, kinships, cX, cR, num_resps_use, NULL, Yadj, dichot, nullweights, num_fixed, Z, 1, 0.001, tryhers[j], polates+j*num_resps_use, nmcmc, gaussian, -9999, NULL, NULL, NULL, NULL);
}

for(m=0;m<num_resps_use;m++)	//solve for mth phenotype
{
for(j=0;j<num_hers1;j++){polates2[j]=polates[m+j*num_resps_use];}

if(polates2[0]<1)	//will exclude the PRS
{
pedhers[m]=0;

if(mpheno!=-1){printf("The estimated pedigree heritability is below %.4f, will exclude the pedigree PRS\n", tryhers[0]);}
else{printf("Phenotype %d: the estimated pedigree heritability is below %.4f, will exclude the pedigree PRS\n", m+1, tryhers[0]);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output, "The estimated pedigree heritability is below %.4f, will exclude the pedigree PRS\n", tryhers[0]);}
else{fprintf(output, "Phenotype %d: the estimated pedigree heritability is below %.4f, will exclude the pedigree PRS\n", m+1, tryhers[0]);}
fclose(output);
}
else
{
pedhers[m]=inter_hers(num_hers1, polates2, tryhers, 1);

if(mpheno!=-1){printf("The estimated pedigree heritability is %.4f\n", pedhers[m]);}
else{printf("Phenotype %d: the estimated pedigree heritability is %.4f\n", m+1, pedhers[m]);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output, "The estimated pedigree heritability is %.4f\n", pedhers[m]);}
else{fprintf(output, "Phenotype %d: the estimated pedigree heritability is %.4f\n", m+1, pedhers[m]);}
fclose(output);
}
}
printf("\n");

//make pedigree PRS (for both full data and cv) and lambdas

total=(2+ncal)*num_resps_use;

for(m=0;m<num_resps_use;m++)	//have normal phenotypes then cv phenotypes, then calibration predictors
{
p=m*(2+ncal);
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]=0;}
for(i=0;i<num_samples_use;i++){cR[(size_t)p*num_samples_use+i]=Yadj[i+m*num_samples_use];}

p=1+m*(2+ncal);
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]=0;}
for(i=0;i<num_train;i++){cR[(size_t)p*num_samples_use+keeptrain[i]]=Yadj[keeptrain[i]+m*num_samples_use];}
for(i=0;i<num_test;i++){cR[(size_t)p*num_samples_use+keeptest[i]]=0;}

for(j=0;j<ncal;j++)
{
p=2+j+m*(2+ncal);
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]=0;}
for(i=0;i<num_samples_use;i++){cR[(size_t)p*num_samples_use+i]=cdata[i+(j+m*ncal)*num_samples_use];}
}
}

sparse_cgd(num_samples_use, total, num_rels, firsts, seconds, kinships, cX, cR, num_resps_use, pedhers, Yadj, dichot, nullweights, num_fixed, Z, 2, 0.0001, -9999, NULL, -9999, NULL, ncal, cdata, cmults, pedgammas, pedsds);

//extract prs (and if excluding the PRS, ensure scaling is 1)

for(m=0;m<num_resps_use;m++)
{
if(pedhers[m]>0)
{
p=m*(2+ncal);

if(dichot==0)	//prs is Y - invV Y (1-h2) - remember to scale
{
for(i=0;i<num_samples_use;i++)
{pedprs[i+m*num_samples_use]=(Yadj[i+m*num_samples_use]-cX[(size_t)p*num_samples_use+i]*(1-pedhers[m]))/pedscales[m];}
}
else	//prs is Y - invW invV Y - no need to scale
{
for(i=0;i<num_samples_use;i++)
{pedprs[i+m*num_samples_use]=Yadj[i+m*num_samples_use]-cX[(size_t)p*num_samples_use+i]/nullweights[i+m*num_samples_use];}
}
}
else
{
for(i=0;i<num_samples_use;i++){pedprs[i+m*num_samples_use]=0;}
pedgammas[m]=1;
pedsds[m]=0;
}
}

//compute MSE for cv PRS (remember we set phenotypes for test samples to zero) - for dichot=1 compute weighted MSE

for(m=0;m<num_resps_use;m++)
{
if(pedhers[m]>0)
{
p=1+m*(2+ncal);

sumsq=0;sumsq2=0;
for(i=0;i<num_test;i++)
{
if(respinds[keeptest[i]+m*num_samples_use]==1)
{
if(dichot==0)	//prs is Ycv - invV Y (1-h2)
{
sumsq+=pow(Yadj[keeptest[i]+m*num_samples_use]+cX[(size_t)p*num_samples_use+keeptest[i]]*(1-pedhers[m]),2);
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2);
}
else	//prs is Ycv - invW invV Y
{
sumsq+=pow(Yadj[keeptest[i]+m*num_samples_use]+cX[(size_t)p*num_samples_use+keeptest[i]]/nullweights[i+m*num_samples_use],2)*nullweights[i+m*num_samples_use];
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2)*nullweights[i+m*num_samples_use];
}
}
}
pedmse[m]=sumsq/sumsq2;
}
else
{pedmse[m]=1;}

if(mpheno!=-1){printf("The pedigree PRS has MSE %.4f and the scaling factor is %.4f\n", pedmse[m], pedgammas[m]);}
else{printf("Phenotype %d: the pedigree PRS has MSE %.4f and the scaling factor is %.4f\n", m+1, pedmse[m], pedgammas[m]);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output, "The pedigree PRS has MSE %.4f and the scaling factor is %.4f\n", pedmse[m], pedgammas[m]);}
else{fprintf(output, "Phenotype %d: the pedigree PRS has MSE %.4f and the scaling factor is %.4f\n", m+1, pedmse[m], pedgammas[m]);}
fclose(output);
}
printf("\n");

}	//end of count>0
}	//end of constructing pedigree PRS

////////

if(useped==1&&fastgwa==0)	//get best power
{
if(power!=-9999)
{
for(m=0;m<num_resps_use;m++){pedtops[m]=0;}
}
else	//do a light version of rhe
{
printf("Determining best power using Randomized Haseman-Elston Regression with %d random vectors\n\n", nmcmc);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Determining best power using Randomized Haseman-Elston Regression with %d random vectors\n", nmcmc);
fclose(output);

data_warn2(bitsize,2*num_samples_use);

value=(double)num_samples_use/1024/1024/1024*8*(nmcmc+num_resps_use)*ndivs*num_pows;
if(value>1){printf("Warning, to perform the analysis requires approximately %.1f Gb; sorry, this can not be reduced\n\n", value);}

#include "getindhersc.c"
}
}

if(useped==1||fastgwa==1)	//save results
{
//save root
sprintf(filename3,"%s.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Phenotypes %s\n", respfile);
if(fastgwa==0){fprintf(output3,"Analysis KVIK_Pedigree\n");}
else{fprintf(output3,"Analysis fastGWA\n");}
if(dichot==0){fprintf(output3,"Regression Linear\n");}
else{fprintf(output3,"Regression Logistic\n");}
if(mpheno!=-1){fprintf(output3,"One Phenotype\n");}
else{fprintf(output3,"Multiple Phenotypes\n");}
if(strcmp(covarfile,"blank")!=0){fprintf(output3,"Covariates %s\n", covarfile);}
else{fprintf(output3,"Covariates none\n");}
if(strcmp(topfile,"blank")!=0){fprintf(output3,"Top_Predictors %s\n", topfile);}
else{fprintf(output3,"Top_Predictors none\n");}
fprintf(output3,"Num_Samples_Used %d\n", num_samples_use);
fprintf(output3,"Num_Predictors_Used %d\n", data_length);
fprintf(output3,"Num_Chromosomes %d\n", num_chr);
fclose(output3);

for(m=0;m<num_resps_use;m++)	//save LOCO PRS and details (and maybe copy root)
{
if(mpheno!=-1){sprintf(filename3,"%s.loco.prs", outfile);}
else{sprintf(filename3,"%s.pheno%d.loco.prs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"FID IID Pedigree\n");

for(i=0;i<num_samples_use;i++)
{
if(respinds[i+m*num_samples_use]==1)
{
if(useped==1){fprintf(output3, "%s %s %.4f\n", ids1[i], ids2[i], pedprs[i+m*num_samples_use]);}
else{fprintf(output3, "%s %s 0\n", ids1[i], ids2[i]);}
}
}
fclose(output3);

if(mpheno!=-1){sprintf(filename3,"%s.loco.details", outfile);}
else{sprintf(filename3,"%s.pheno%d.loco.details", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(mpheno!=-1){fprintf(output3,"Phenotype_Number %d\n", mpheno+1);}
else{fprintf(output3,"Phenotype_Number %d\n", m+1);}
fprintf(output3,"Actual_Sample_Size %d\n", respcounts[m]);
fprintf(output3,"Approx_Effective_Sample_Size %d\n", respcounts[m]);
if(useped==1)
{
fprintf(output3,"Scaling_Estimate %.4f\n", pedgammas[m]);
fprintf(output3,"Scaling_SD %.4f\n", pedsds[m]);
}
else
{
fprintf(output3,"Scaling_Estimate 1\n");
fprintf(output3,"Scaling_SD 0\n");
}
if(fastgwa==0){fprintf(output3,"Power %.4f\n", powers[pedtops[m]]);}
else{fprintf(output3,"Power NA\n");}
fprintf(output3,"Heritability %.4f\n", pedhers[m]);
fclose(output3);

if(mpheno==-1)	//copy root
{
sprintf(cmd, "cp %s.root %s.pheno%d.root", outfile, outfile, m+1);
system(cmd);
}
}

if(mpheno!=-1){printf("LOCO results saved in %s.loco.prs and %s.loco.details\n\n", outfile, outfile);}
else{printf("LOCO results saved in %s.phenoX.loco.prs and %s.phenoX.loco.prs, where X is the phenotype number\n\n", outfile, outfile);}
}

free(tryhers);
free(pedtops);free(pedscales);free(pedhers);free(pedprs);free(pedgammas);free(pedsds);free(pedmse);
free(dindex);free(dindex2);free(ddata);free(dcentres);free(dmults);free(dsqdevs);
for(j=0;j<nped;j++){free(dpreds[j]);}free(dpreds);
free(XTCX);free(XTCX2);free(kins);
free(firsts);free(seconds);free(kinships);
free(cX);free(cR);
free(polates);free(polates2);

///////////////////////////


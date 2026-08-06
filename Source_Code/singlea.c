/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor linear regression - might permute, have top-preds, spa, sample weights (but no enviro)
//might transform (in which case, do not have top-preds)
//if using sample weights, will have adjpreds=2 (can use sample weights with spa, but not sure whether results will be sensible)

///////////////////////////

if(strcmp(transfile,"blank")!=0)	//read transformation and inverse
{
Umat=malloc(sizeof(double)*num_resps_use*num_resps_use);
Umat2=malloc(sizeof(double)*num_resps_use*num_resps_use);

sprintf(filename,"%s.transformation", transfile);
for(m=0;m<num_resps_use;m++){read_values(filename, Umat+m*num_resps_use, num_resps_use, NULL, m+1, 0, 0);}

sprintf(filename,"%s.inverse", transfile);
for(m=0;m<num_resps_use;m++){read_values(filename, Umat2+m*num_resps_use, num_resps_use, NULL, m+1, 0, 0);}

printf("transform matrix\n");
print_top(Umat, num_resps_use, num_resps_use, num_resps_use);
printf("transform inverse\n");
print_top(Umat2, num_resps_use, num_resps_use, num_resps_use);
printf("\n");
}

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//threshold is the p-value threshold used with adjpreds
threshold=0.05;
if(spatest==1&&spathresh>threshold){threshold=spathresh;}

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
datasqs=malloc(sizeof(double)*bitsize);

order=malloc(sizeof(int)*num_samples_use);
sweights=malloc(sizeof(double)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use*num_resps_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

Pscales=malloc(sizeof(double)*num_resps_use);
Pscales2=malloc(sizeof(double)*num_resps_use);
Pscales3=malloc(sizeof(double)*num_resps_use);
Pvarphens=malloc(sizeof(double)*num_resps_use);

thetas=malloc(sizeof(double)*num_fixed*num_resps_use);
thetasds=malloc(sizeof(double)*num_fixed*num_resps_use);
thetapvas=malloc(sizeof(double)*num_fixed*num_resps_use);
Yadj=malloc(sizeof(double)*num_samples_use*num_resps_use);

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4*bitsize*num_resps_use);
spastatus=malloc(sizeof(int)*bitsize*num_resps_use);

YTdata=malloc(sizeof(double)*bitsize*num_resps_use);
XTCX=malloc(sizeof(double)*bitsize);

if(adjpreds==1)
{
Z2=malloc(sizeof(double)*num_samples_use*num_fixed);
ZTdata=malloc(sizeof(double)*num_fixed);
}

if(spatest==1)
{
usedpreds=malloc(sizeof(int)*bitsize);

Pspamax=malloc(sizeof(double)*num_resps_use);
knots=malloc(sizeof(double)*num_knots*num_resps_use);
bins=malloc(sizeof(double)*num_bins);
CGF0=malloc(sizeof(double*)*num_knots*num_resps_use);
CGF1=malloc(sizeof(double*)*num_knots*num_resps_use);
CGF2=malloc(sizeof(double*)*num_knots*num_resps_use);
CGF3=malloc(sizeof(double*)*num_knots*num_resps_use);

for(j=0;j<num_knots*num_resps_use;j++)
{
CGF0[j]=malloc(sizeof(double)*num_bins);
CGF1[j]=malloc(sizeof(double)*num_bins);
CGF2[j]=malloc(sizeof(double)*num_bins);
CGF3[j]=malloc(sizeof(double)*num_bins);
}
}

if(strcmp(transfile,"blank")!=0)
{
Yadj2=malloc(sizeof(double)*num_samples_use*num_resps_use);
stats2=malloc(sizeof(double)*4*bitsize*num_resps_use);
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

if(strcmp(sampwfile,"blank")!=0)	//get square roots of sample weights (W^1/2)
{
readdoubles=malloc(sizeof(double)*num_samples_use);
read_sampwfile(sampwfile, readdoubles, num_samples_use, ids3, ids1, ids2);

//scale roots so the mean inverse is one (might help stability)
sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(readdoubles[i],-1);}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){readdoubles[i]*=mean;}

for(i=0;i<num_samples_use;i++){sweights[i]=pow(readdoubles[order[i]],0.5);}
free(readdoubles);
}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

//fill Y (maybe permuted and multiplying by W^1/2)
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Y[i+m*num_samples_use]=resp[order[i]+m*num_samples_use]*sweights[i];}
}

//fill Z (maybe permuted and multiplying by W^1/2)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use]*sweights[i];}
}

//decide how much to scale effect sizes, sds and spa test statistics
for(m=0;m<num_resps_use;m++)
{
//will scale by 1/p, 1/p and 1, where p is prop non-missing
Pscales[m]=pow((double)respcounts[m]/num_samples_use,-1);
Pscales2[m]=pow((double)respcounts[m]/num_samples_use,-1);
Pscales3[m]=1;
}

for(m=0;m<num_resps_use;m++)	//solve null models (and pad any missing values)
{
//get thetas, adjusted response, pad missing and get variance
reg_covar_lin_missing(Y+m*num_samples_use, Z, num_samples_use, num_fixed, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, Yadj+m*num_samples_use, missingvalue);
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);

//save coefficients
if(mpheno!=-1){sprintf(filename,"%s.coeff", outfile);}
else{sprintf(filename,"%s.pheno%d.coeff", outfile, m+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component\tEffect\tSE\tP\n");
fprintf(output, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0+m*num_fixed], thetasds[0+m*num_fixed], thetapvas[0+m*num_fixed]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[j+m*num_fixed], thetasds[j+m*num_fixed], thetapvas[j+m*num_fixed]);}
fclose(output);
}	//end of m loop

if(strcmp(transfile,"blank")!=0)	//replace Yadj with transformed phenotypes (and regress out covariates to be safe)
{
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Yadj2[i+m*num_samples_use]=Yadj[i+m*num_samples_use];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Yadj2, &num_samples_use, Umat, &num_resps_use, &beta, Yadj, &num_samples_use);

reg_covar_matrix(Yadj, Z, num_samples_use, num_resps_use, num_fixed);

//recompute Pvarphens, and get number of non trivial phenotypes
total=0;
for(m=0;m<num_resps_use;m++)
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;Pvarphens[m]=sumsq/num_samples_use-pow(mean,2);
total+=(Pvarphens[m]>0);
}
printf("%d out of %d phenotypes contribute\n\n", total, num_resps_use);
}

//set tindex - indicates whether each predictor is a top preds (and if so, which)
for(j=0;j<data_length;j++){tindex[j]=-9999;}
if(num_tops>0)	//find the tops
{
mark=0;
for(j=0;j<num_tops;j++)
{
while(tkeeppreds[j]>keeppreds_use[mark]){mark++;}
tindex[mark]=j;
}
}

if(spatest==0)	//convenient to set spastatus to zero
{
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<bitsize;j++){spastatus[j+m*bitsize]=0;}
}
}
else	//using empirical spa - can set bins and spamax
{
//set bins - evenly spaced from -2 to 2
for(k=0;k<num_bins;k++){bins[k]=-2+(double)k/(num_bins-1)*4;}

for(m=0;m<num_resps_use;m++)
{
if(spamax==-9999)	//set spamax based on variance of Yadj
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i+m*num_samples_use];sumsq+=pow(Yadj[i+m*num_samples_use],2);}
mean=sum/num_samples_use;var=sumsq/num_samples_use-pow(mean,2);
Pspamax[m]=2000*pow(num_samples_use*var,-.5);
}
else{Pspamax[m]=spamax;}
}
}

if(adjpreds==1)	//get Z2 = Z inv(ZTZ)
{get_Z2(Z2, Z, Z, num_samples_use, num_fixed);}

if(spatest==1)	//using empirical spa
{
for(m=0;m<num_resps_use;m++)
{
if(mpheno!=-1){printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -Pspamax[m], Pspamax[m], num_bins);}
else{printf("Phenotype %d: computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", m+1, num_knots, -Pspamax[m], Pspamax[m], num_bins);}

empirical_cumulants(Yadj+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
printf("\n");
}

////////

//deal with progress and on-the-fly files

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

for(m=0;m<num_resps_use;m++)
{
if(mpheno!=-1){sprintf(filename2,"%s.assoc", outfile);}
else{sprintf(filename2,"%s.pheno%d.assoc", outfile, m+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tMAF\tCallRate\tMachR2\tSPA_Status\n");
fclose(output2);

if(mpheno!=-1){sprintf(filename3,"%s.summaries", outfile);}
else{sprintf(filename3,"%s.pheno%d.summaries", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tA1\tA2\tZ\tn\tA1Freq\n");
fclose(output3);

if(mpheno!=-1){sprintf(filename4,"%s.pvalues", outfile);}
else{sprintf(filename4,"%s.pheno%d.pvalues", outfile, m+1);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor\tP\n");
fclose(output4);

if(mpheno!=-1){sprintf(filename5,"%s.score", outfile);}
else{sprintf(filename5,"%s.pheno%d.score", outfile, m+1);}
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor\tA1\tA2\tCentre\tALL\tP<0.1\tP<0.01\tP<0.001\tP<0.0001\tP<0.00001\tP<5e-8\n");
fclose(output5);
}

if(strcmp(transfile,"blank")!=0)
{
sprintf(filename6,"%s.factors", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "Predictor");
for(m=0;m<num_resps_use;m++){fprintf(output6,"\tStats_%d", m+1);}
fprintf(output6,"\n");
fclose(output6);

sprintf(filename7,"%s.combined", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
fprintf(output7, "Predictor\tChiSq_%d_Stat\tP\tChiSq_1_Stat\n", total);
fclose(output7);
}

//work out how many bits
bittotal=0;
bitend=0;
while(bitend<data_length)
{
bitstart=bitend;
bitend=bitstart+bitsize;
if(bitend>data_length){bitend=data_length;}
while(chr[bitend-1]>chr[bitstart]){bitend--;}

bittotal++;
}

//ready for bit loop
bit=0;
bitend=0;
while(bitend<data_length)
{
bitstart=bitend;
bitend=bitstart+bitsize;
if(bitend>data_length){bitend=data_length;}
while(chr[bitend-1]>chr[bitstart]){bitend--;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
printf("Performing linear regression for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing linear regression for Chunk %d of %d\n", bit+1, bittotal);
}

//read data, compute statistics, centre and set missing to zero
if(dtype==1)	//fast way
{
read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, bedzeros, bedones, bedtwos, 1, maxthreads);
}
else	//slow way
{
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1);
}

if(passqc==1)	//perform qc (will already have mults=-9999 for trivial predictors)
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
maf=centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]);
if(minmaf!=-9999&&maf<minmaf){mults[bitstart+j]=-9999;}
if(maxmaf!=-9999&&maf>maxmaf){mults[bitstart+j]=-9999;}
if(minvar!=-9999&&sqdevs[bitstart+j]<minvar){mults[bitstart+j]=-9999;}
if(minobs!=-9999&&rates[bitstart+j]<minobs){mults[bitstart+j]=-9999;}
if(mininfo!=-9999&&infos[bitstart+j]<mininfo){mults[bitstart+j]=-9999;}
}
}
}

if(strcmp(sampwfile,"blank")!=0)	//multiply by W^1/2
{copy_matrix(num_samples_use, bitlength, data, data, 1, sweights);}

//get XTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_resps_use, &num_samples_use, &alpha, data, &num_samples_use, Yadj, &num_samples_use, &beta, YTdata, &bitlength);

//ready to test
for(m=0;m<num_resps_use;m++)
{
if(m==0)	//only need to get XTCX for the first phenotype
{
if(adjpreds==2)	//regress out covariates
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

//get XTCX
if(adjpreds!=2)	//already have
{
for(j=0;j<bitlength;j++){XTCX[j]=num_samples_use*sqdevs[bitstart+j];}
}
else	//must compute
{
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}
}
}
}

if(Pvarphens[m]>0)
{
#pragma omp parallel for private(j, mark, value, value2, value3, i) schedule(static)
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//will be testing - remember Pscales and Pscales2
{
if(tindex[bitstart+j]!=-9999)	//a top predictor, so have already tested (did not pad, so no need to scale)
{
stats[0+mark]=thetas[num_covars+tindex[bitstart+j]+m*num_fixed];
stats[1+mark]=thetasds[num_covars+tindex[bitstart+j]+m*num_fixed];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=thetapvas[num_covars+tindex[bitstart+j]+m*num_fixed];
}
else	//not a top, so must test
{
if(sandwich==0)	//use standard estimates
{
//compute YTCX/XTCX and YCT(YC-XCvalue)/XTCX(n-nf-1)
value=YTdata[j+m*bitlength]/XTCX[j];
value2=(Pvarphens[m]*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);
}
else	//use sandwich estimates
{
//effect size remains YTCX/XTCX
value=YTdata[j+m*bitlength]/XTCX[j];

//variance is sum(Xi^2 ri^2)/XTCX^2, where ri are residuals
value3=0;
for(i=0;i<num_samples_use;i++)
{value3+=pow(data[(size_t)j*num_samples_use+i]*(Yadj[i+m*num_samples_use]-value*data[(size_t)j*num_samples_use+i]),2);}
value2=value3*pow(XTCX[j],-2);
}

if(strcmp(transfile,"blank")==0)	//will correct for padding now
{
stats[0+mark]=value*Pscales[m];
stats[1+mark]=pow(value2,.5)*Pscales2[m];
}
else	//will correct for padding after transforming
{
stats[0+mark]=value;
stats[1+mark]=pow(value2,.5);
}
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}	//end of testing
}	//end of j loop

if(adjpreds==1)	//revisit the most significant predictors
{
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//might consider
{
if(tindex[bitstart+j]==-9999&&stats[3+mark]<threshold)	//recalculate XTCX, then test statistics
{
//save original value
value4=XTCX[j];

//need to subtract Z2 x ZTX from data (its ok if have already done this for previous phenotypes)
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_fixed, &alpha, Z, &num_samples_use, data+(size_t)j*num_samples_use, &one, &beta, ZTdata, &one);

alpha=-1.0;beta=1.0;
dgemv_("N", &num_samples_use, &num_fixed, &alpha, Z2, &num_samples_use, ZTdata, &one, &beta, data+(size_t)j*num_samples_use, &one);

XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}

if(sandwich==0)	//use standard estimates
{
value=YTdata[j+m*bitlength]/XTCX[j];
value2=(Pvarphens[m]*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);
}
else	//use sandwich estimates
{
value=YTdata[j+m*bitlength]/XTCX[j];
value3=0;
for(i=0;i<num_samples_use;i++)
{value3+=pow(data[(size_t)j*num_samples_use+i]*(Yadj[i+m*num_samples_use]-value*data[(size_t)j*num_samples_use+i]),2);}
value2=value3*pow(XTCX[j],-2);
}

if(strcmp(transfile,"blank")==0)	//will correct for padding now
{
stats[0+mark]=value*Pscales[m];
stats[1+mark]=pow(value2,.5)*Pscales2[m];
}
else	//will correct for padding after transforming
{
stats[0+mark]=value;
stats[1+mark]=pow(value2,.5);
}
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

//restore original value
XTCX[j]=value4;
}}	//end of revisiting
}	//end of j loop
}	//end of adjpreds=1

if(spatest==1)	//compute spa test statistic - remember Pscales3 (bit silly, because Pscale3=1)
{
//work out which predictors to test
count=0;
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//might test
{
if(tindex[bitstart+j]==-9999&&stats[3+mark]<spathresh)	//will compute SPA
{usedpreds[count]=j;count++;}
else	//will report non-SPA results
{spastatus[j+m*bitlength]=0;}
}
else	//trivial, so will not test, but helps to set status to missing
{spastatus[j+m*bitlength]=-9999;}
}

while(1)
{
#pragma omp parallel for private(j2, j, mark) schedule(static)
for(j2=0;j2<count;j2++)
{
j=usedpreds[j2];
mark=4*j+m*4*bitlength;

spastatus[j+m*bitlength]=spa_test(YTdata[j+m*bitlength], data+(size_t)j*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, stats+mark, Pscales3[m]);
}

//see if necessary to increase range for any predictors (and recompute SPA)
count=0;
for(j=0;j<bitlength;j++)
{
if(spastatus[j+m*bitlength]==-2){usedpreds[count]=j;count++;}
}
if(count==0){break;}

Pspamax[m]*=5;
if(mpheno!=-1){printf("Increasing SPA range to %f and recomputing empirical CGF\n", Pspamax[m]);}
else{printf("Phenotype %d: increasing SPA range to %f and recomputing empirical CGF\n", m+1, Pspamax[m]);}
empirical_cumulants(Yadj+m*num_samples_use, num_samples_use, num_knots, knots+m*num_knots, num_bins, bins, CGF0+m*num_knots, CGF1+m*num_knots, CGF2+m*num_knots, CGF3+m*num_knots, Pspamax[m]);
}
}
}
else	//trivial (transformed) phenotype
{
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;
stats[0+mark]=0;stats[1+mark]=0;stats[2+mark]=0;stats[3+mark]=1;
}
}
}	//end of m loop

if(strcmp(transfile,"blank")!=0)	//get combined pvalues then transform statistics (will not have top predictors)
{
//save current values
sprintf(filename6,"%s.factors", outfile);
if((output6=fopen(filename6,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename6);exit(1);}

for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
fprintf(output6,"%s", preds[bitstart+j]);
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
fprintf(output6, "\t%.6f", stats[2+mark]);
}
fprintf(output6,"\n");
}
else
{
fprintf(output6,"%s", preds[bitstart+j]);
for(m=0;m<num_resps_use;m++){fprintf(output6, "\tNA");}
fprintf(output6,"\n");
}
}
fclose(output6);

sprintf(filename7,"%s.combined", outfile);
if((output7=fopen(filename7,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename7);exit(1);}

for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
sum=0;
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
sum+=pow(stats[2+mark],2);
}
value=pochisq(sum, total);
if(value>0)	//have a valid p-value, so convert to chisq(1) statistic
{value2=pow(normal_inv(value/2),2);}
else	//not valid, so use approximation from wiki
//https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Cumulative_distribution_function
{
value3=pow(sum/total,1.0/3)-1+pow(4.5*total,-1);
value2=pow(value3,2)*4.5*total;
}
fprintf(output7, "%s\t%.4f\t%.4e\t%.4f\n", preds[bitstart+j], sum, value, value2);
}
else{fprintf(output7, "%s\tNA\tNA\tNA\n", preds[bitstart+j]);}
}
fclose(output7);

#pragma omp parallel for private(j, m, m2, mark, mark3) schedule(static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
//first compute estimate and variance (save in stats2)
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
stats2[0+mark]=0;
stats2[1+mark]=0;
for(m2=0;m2<num_resps_use;m2++)
{
mark3=4*j+m2*4*bitlength;
stats2[0+mark]+=stats[0+mark3]*Umat2[m+m2*num_resps_use];
stats2[1+mark]+=pow(stats[1+mark3]*Umat2[m+m2*num_resps_use],2);
}
}

//now copy back into stats, remembering to scale
for(m=0;m<num_resps_use;m++)
{
mark=4*j+m*4*bitlength;
stats[0+mark]=stats2[0+mark]*Pscales[m];
stats[1+mark]=pow(stats2[1+mark],.5)*Pscales2[m];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}}
}

for(m=0;m<num_resps_use;m++)	//save
{
//reopen output files
if(mpheno!=-1){sprintf(filename2,"%s.assoc", outfile);}
else{sprintf(filename2,"%s.pheno%d.assoc", outfile, m+1);}
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
if(mpheno!=-1){sprintf(filename3,"%s.summaries", outfile);}
else{sprintf(filename3,"%s.pheno%d.summaries", outfile, m+1);}
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
if(mpheno!=-1){sprintf(filename4,"%s.pvalues", outfile);}
else{sprintf(filename4,"%s.pheno%d.pvalues", outfile, m+1);}
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
if(mpheno!=-1){sprintf(filename5,"%s.score", outfile);}
else{sprintf(filename5,"%s.pheno%d.score", outfile, m+1);}
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//save results
for(j=0;j<bitlength;j++)
{
mark=4*j+m*4*bitlength;

if(mults[bitstart+j]!=-9999)	//tested predictor
{
//print assoc
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4f\t%.4e\t%.4e\t%.4e\t", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
if(spastatus[j+m*bitlength]==0){fprintf(output2, "NOT_USED\n");}
if(spastatus[j+m*bitlength]==1){fprintf(output2, "SUCCESS\n");}
if(spastatus[j+m*bitlength]==2){fprintf(output2, "APPROX\n");}
if(spastatus[j+m*bitlength]==-1){fprintf(output2, "FAILED\n");}

//print summaries, pvalues and scores
fprintf(output3, "%s\t%s\t%s\t%.4f\t%.0f\t%.6f\n", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j], stats[2+mark], rates[bitstart+j]*respcounts[m], centres[bitstart+j]/2);
fprintf(output4, "%s\t%.4e\n", preds[bitstart+j], stats[3+mark]);
fprintf(output5, "%s\t%s\t%s\t%.6f\t%.6f", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j],  centres[bitstart+j], stats[0+mark]);
for(k=0;k<6;k++)
{
if(stats[3+mark]<cuts[k]){fprintf(output5, "\t%.6f", stats[0+mark]);}
else{fprintf(output5, "\t0");}
}
fprintf(output5, "\n");
}
else	//did not test, but will include in assoc if not doing qc
{
if(passqc==0)
{
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "NA\tNA\tNA\tNA\t");
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
fprintf(output2, "NOT_USED\n");
}
}
}	//end of j loop

fclose(output2);
fclose(output3);
fclose(output4);
fclose(output5);
}	//end of m loop

bit++;
}	//end of while loop
printf("\n");

fclose(output);

count=0;for(j=0;j<data_length;j++){count+=(mults[j]==-9999);}
if(count==data_length)
{
if(passqc==0)
{printf("Warning, all %d predictors were excluded because they were trivial (showed no variation)\n\n", count);}
else
{printf("Warning, all %d predictors were excluded because they were trivial or failed quality control\n\n", count);}
}
else
{
if(count>0)
{
if(passqc==0)
{printf("Warning, %d predictors were excluded because they were trivial (showed no variation)\n\n", count);}
else
{printf("Warning, %d predictors were excluded because they were trivial or failed quality control\n\n", count);}
}
}

if(mpheno!=-1){printf("Main results saved in %s.assoc, with a summary version in %s.summaries, p-values in %s.pvalues and score file in %s.score\n\n", outfile, outfile, outfile, outfile);}
else{printf("Main results saved in %s.phenoX.assoc, with a summary version in %s.phenoX.summaries, p-values in %s.phenoX.pvalues and score file in %s.phenoX.score, where X is the phenotype number\n\n", outfile, outfile, outfile, outfile);}

if(strcmp(transfile,"blank")!=0){free(Umat);free(Umat2);}
free(data);free(datasqs);
free(order);free(sweights);free(Y);free(Z);
free(Pscales);free(Pscales2);free(Pscales3);free(Pvarphens);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(tindex);free(stats);free(spastatus);
free(YTdata);free(XTCX);
if(adjpreds==1){free(Z2);free(ZTdata);}
if(spatest==1)
{
free(usedpreds);
for(j=0;j<num_knots*num_resps_use;j++){free(CGF0[j]);free(CGF1[j]);free(CGF2[j]);free(CGF3[j]);}
free(Pspamax);free(knots);free(bins);free(CGF0);free(CGF1);free(CGF2);free(CGF3);
}
if(strcmp(transfile,"blank")!=0){free(Yadj2);free(stats2);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


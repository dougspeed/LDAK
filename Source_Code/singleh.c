/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor logistic regression - might permute, have top-preds, spa and offsets (no enviro or sample weights)
//will not have (non-zero) offsets with mpheno=-1, nor missing phenotypic values / mpheno=-1 / spatest with scoretest=0

///////////////////////////

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//threshold is the p-value threshold used with adjpreds
threshold=0.05;
if(spatest==1&&spathresh<threshold){threshold=spathresh;}

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use*num_resps_use);
Z=malloc(sizeof(double)*num_samples_use*(num_fixed+1));

Pscales=malloc(sizeof(double)*num_resps_use);

thetas=malloc(sizeof(double)*num_fixed*num_resps_use);
thetasds=malloc(sizeof(double)*num_fixed*num_resps_use);
thetapvas=malloc(sizeof(double)*num_fixed*num_resps_use);

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4*bitsize);
spastatus=malloc(sizeof(int)*bitsize);

if(scoretest==1)
{
nullprobs=malloc(sizeof(double)*num_samples_use*num_resps_use);
nullweights=malloc(sizeof(double)*num_samples_use*num_resps_use);
Yadj=malloc(sizeof(double)*num_samples_use*num_resps_use);
YTdata=malloc(sizeof(double)*bitsize*num_resps_use);
XTCX=malloc(sizeof(double)*bitsize);

if(adjpreds==1)
{
Z2=malloc(sizeof(double)*num_samples_use*num_fixed*num_resps_use);
Z3=malloc(sizeof(double)*num_samples_use*num_fixed*num_resps_use);
ZTdata=malloc(sizeof(double)*num_fixed);
}
}

if(spatest==1){usedpreds=malloc(sizeof(int)*bitsize);}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//will divide thetas and SDs by prop non-missing (for spa, must divide test statistic by root)
for(m=0;m<num_resps_use;m++)
{Pscales[m]=pow((double)respcounts[m]/num_samples_use,-1);}

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y (maybe permuted)
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){Y[i+m*num_samples_use]=resp[order[i]+m*num_samples_use];}
}

if(pad==1)	//pad missing
{impute_matrix_missing(Y, num_samples_use, num_samples_use, num_resps_use, missingvalue);}

//fill (start of) Z (maybe permuted)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

for(m=0;m<num_resps_use;m++)	//solve null models
{
//get thetas (and maybe nullprobs, nullweights and W x adjusted response)
if(scoretest==0)
{
reg_covar_log(Y+m*num_samples_use, Z, num_samples_use, num_covars, num_tops, offsets, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, NULL, -9999, NULL, NULL, -9999, 0.001, 100);
}
else
{
reg_covar_log(Y+m*num_samples_use, Z, num_samples_use, num_covars, num_tops, offsets, thetas+m*num_fixed, thetasds+m*num_fixed, thetapvas+m*num_fixed, nullprobs+m*num_samples_use, 2, NULL, NULL, -9999, 0.001, 100);
for(i=0;i<num_samples_use;i++){nullweights[i+m*num_samples_use]=nullprobs[i+m*num_samples_use]*(1-nullprobs[i+m*num_samples_use]);}
for(i=0;i<num_samples_use;i++){Yadj[i+m*num_samples_use]=(Y[i+m*num_samples_use]-nullprobs[i+m*num_samples_use]);}
}

if(respcounts[m]<num_samples_use)	//adjust estimates for padding (can not have padded values for scoretest=0)
{
for(j=0;j<num_fixed;j++)
{
thetas[j+m*num_fixed]=thetas[j+m*num_fixed]*Pscales[m];
thetasds[j+m*num_fixed]=thetasds[j+m*num_fixed]*Pscales[m];
thetapvas[j+m*num_fixed]=erfc(fabs(thetas[j+m*num_fixed]/thetasds[j+m*num_fixed])*M_SQRT1_2);
}
}

//save coefficients
if(mpheno!=-1){sprintf(filename,"%s.coeff", outfile);}
else{sprintf(filename,"%s.pheno%d.coeff", outfile, m+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component\tLog_OR\tSE\tP\n");
fprintf(output, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0+m*num_fixed], thetasds[0+m*num_fixed], thetapvas[0+m*num_fixed]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[j+m*num_fixed], thetasds[j+m*num_fixed], thetapvas[j+m*num_fixed]);}
fclose(output);
}	//end of m loop

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

if(scoretest==1&&adjpreds==1)	//fill Z3 with WZ, then get Z2 = Z inv(ZTZ3) - do separately for each phenotype
{
for(m=0;m<num_resps_use;m++)
{
copy_matrix(num_samples_use, num_fixed, Z, Z3+m*num_samples_use*num_fixed, 1, nullweights+m*num_samples_use);
get_Z2(Z2+m*num_samples_use*num_fixed, Z, Z3+m*num_samples_use*num_fixed, num_samples_use, num_fixed);
}
}

if(spatest==0)	//set spastatus to zero
{
for(j=0;j<bitsize;j++){spastatus[j]=0;}
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
if(scoretest==0){fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tLog_OR\tSE\tMAF\tCallRate\tMachR2\tSPA_Status\n");}
else{fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tScore_Stat\tScore_P\tApprox_Log_OR\tApprox_SE\tMAF\tCallRate\tMachR2\tSPA_Status\n");}
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

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
fclose(output);
printf("Performing logistic regression for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing logistic regression for Chunk %d of %d\n", bit+1, bittotal);
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

if(scoretest==1)
{
//get XTY
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_resps_use, &num_samples_use, &alpha, data, &num_samples_use, Yadj, &num_samples_use, &beta, YTdata, &bitlength);
}

//ready to test
for(m=0;m<num_resps_use;m++)
{
if(scoretest==1)
{
if(adjpreds==2)	//regress out covariates (its ok if have already done this for previous phenotypes)
{reg_covar_weighted(data, Z, num_samples_use, bitlength, num_fixed, nullweights+m*num_samples_use);}

//get XTCX
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2)*nullweights[i+m*num_samples_use];}
}
}

for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//will be testing - remember to scale theta and SD by Pscales (not necessary for scoretest=0)
{
if(tindex[j]!=-9999)	//a top predictor, so have already tested (and scaled)
{
stats[0+mark]=thetas[num_covars+tindex[bitstart+j]+m*num_fixed];
stats[1+mark]=thetasds[num_covars+tindex[bitstart+j]+m*num_fixed];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=thetapvas[num_covars+tindex[bitstart+j]+m*num_fixed];
}
else	//not a top, so must test
{
if(scoretest==0)	//add test predictor to end of Z and test (reg_single_log returns statistics corresponding to final covariate)
{
for(i=0;i<num_samples_use;i++){Z[i+num_fixed*num_samples_use]=data[(size_t)j*num_samples_use+i];}
reg_covar_log_lite(Y, Z, num_samples_use, num_fixed+1, offsets, stats+mark, thetas, 0.001, 100);
}
else
{
//the score statistic is YTdata, with variance XTCX
//SAIGE paper explains how score/variance approx equals LogOR, and 1/var approx equals Var of LogOR
//this is because XT(Y-mu)/XTWX is 1-step NR estimate of LogOR
stats[0+mark]=YTdata[j+m*bitlength]/XTCX[j]*Pscales[m];
stats[1+mark]=pow(XTCX[j],-.5)*Pscales[m];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}
}	//end of testing
}	//end of j loop

if(scoretest==1&&adjpreds==1)	//revisit the most significant predictors
{
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//might consider
{
if(tindex[bitstart+j]==-9999&&stats[3+mark]<threshold)	//recalculate XTCX, then test statistics
{
//save original value
value4=XTCX[j];

//need to subtract Z2 x Z3TX from data (its ok if have already done this for previous phenotypes)
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_fixed, &alpha, Z3+m*num_samples_use*num_fixed, &num_samples_use, data+(size_t)j*num_samples_use, &one, &beta, ZTdata, &one);

alpha=-1.0;beta=1.0;
dgemv_("N", &num_samples_use, &num_fixed, &alpha, Z2+m*num_samples_use*num_fixed, &num_samples_use, ZTdata, &one, &beta, data+(size_t)j*num_samples_use, &one);

XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2)*nullweights[i+m*num_samples_use];}

stats[0+mark]=YTdata[j+m*bitlength]/XTCX[j];
stats[1+mark]=pow(XTCX[j],-.5);
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

//restore original value
XTCX[j]=value4;
}}	//end of revisiting
}	//end of j loop
}	//end of scoretest=1 and adjpreds=1

if(spatest==1)	//compute spa test statistic
{
//work out which predictors to test
count=0;
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//might test
{
if(tindex[bitstart+j]==-9999&&stats[3+mark]<spathresh)	//will compute SPA
{usedpreds[count]=j;count++;}
else	//will report non-SPA results
{spastatus[j]=0;}
}
else	//trivial, so will not test, but helps to set status to missing
{spastatus[j]=-9999;}
}

#pragma omp parallel for private(j2, j, mark) schedule(static)
for(j2=0;j2<count;j2++)
{
j=usedpreds[j2];
mark=4*j;

if(spaside==1){spastatus[j]=spa_logistic_one(YTdata[j+m*bitlength], data+(size_t)j*num_samples_use, num_samples_use, nullprobs, nullweights, stats+mark, 1);}
else{spastatus[j]=spa_logistic_two(YTdata[j+m*bitlength], XTCX[j], data+(size_t)j*num_samples_use, num_samples_use, nullprobs, nullweights, stats+mark, 1);}
}
}

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
mark=4*j;

if(mults[bitstart+j]!=-9999)	//tested predictor
{
//print assoc
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4f\t%.4e\t%.4e\t%.4e\t", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
if(spastatus[j]==0){fprintf(output2, "NOT_USED\n");}
if(spastatus[j]==1){fprintf(output2, "SUCCESS\n");}
if(spastatus[j]==2){fprintf(output2, "APPROX\n");}
if(spastatus[j]==-1){fprintf(output2, "FAILED\n");}

//print summaries, pvalues and scores
fprintf(output3, "%s\t%s\t%s\t%.4f\t%d\t%.4f\n", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j], stats[2+mark], (int)(rates[bitstart+j]*num_samples_use), centres[bitstart+j]/2);
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
if(spastatus[j]==0){fprintf(output2, "NOT_USED\n");}
}
}
}	//end of j loop

fclose(output2);
fclose(output3);
fclose(output4);
fclose(output5);
}	//end of m loop
}	//end of bit loop
printf("\n");

fclose(output);

count=0;for(j=0;j<num_preds_use;j++){count+=(mults[j]==-9999);}
if(count==num_preds_use)
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

free(data);
free(order);free(Y);free(Z);
free(Pscales);
free(thetas);free(thetasds);free(thetapvas);
free(tindex);free(stats);free(spastatus);
if(scoretest==1)
{
free(nullprobs);free(nullweights);free(Yadj);free(YTdata);free(XTCX);
if(adjpreds==1){free(Z2);free(Z3);free(ZTdata);}
}
if(spatest==1){free(usedpreds);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


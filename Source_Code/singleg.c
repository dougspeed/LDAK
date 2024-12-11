/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor linear regression with duos (no top-preds, enviro, spa or sample weights) - note that keepsamps might not be monotonic / unique
//the full model regresses the phenotype on the predictor and one set of parental genotypes

///////////////////////////

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//allocate variables

data_warn2(bitsize,4*num_samples_use);
data=malloc(sizeof(double)*2*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);

stats=malloc(sizeof(double)*4*bitsize);
stats2=malloc(sizeof(double)*8*bitsize);
stats3=malloc(sizeof(double)*4*bitsize);

YTdata=malloc(sizeof(double)*2*bitsize);
XTCX=malloc(sizeof(double)*4*bitsize);
XTCX2=malloc(sizeof(double)*4*bitsize);

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y (maybe permuted)
for(i=0;i<num_samples_use;i++){Y[i]=resp[order[i]];}

//fill Z (maybe permuted)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
}

//solve null model, get thetas and adjust response
reg_covar_lin(Y, Z, num_samples_use, num_covars, 0, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);

//save coefficients
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"ws"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component\tEffect\tSE\tP\n");
fprintf(output, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

////////

//deal with progress and on-the-fly files

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.basic",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tA1_Mean\tMAF\tCallRate\tMachR2\n");

if(duos==1)	//one parent
{
sprintf(filename3,"%s.duos",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Chromosome\tPredictor\tBasepair\tA1\tA2\tOffspring_Effect\tOffspring_SE\tOffspring_Z\tOffspring_P\tParent_Effect\tParent_SE\tParent_Z\tParent_P\n");
}
if(duos==2)	//father
{
sprintf(filename3,"%s.duos.fathers",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Chromosome\tPredictor\tBasepair\tA1\tA2\tOffspring_Effect\tOffspring_SE\tOffspring_Z\tOffspring_P\tFather_Effect\tFather_SE\tFather_Z\tFather_P\n");
}
if(duos==3)	//mother
{
sprintf(filename3,"%s.duos.mothers",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Chromosome\tPredictor\tBasepair\tA1\tA2\tOffspring_Effect\tOffspring_SE\tOffspring_Z\tOffspring_P\tMother_Effect\tMother_SE\tMother_Z\tMother_P\n");
}

if(duos==1)	//one parent
{
sprintf(filename4,"%s.swap",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tA1_Mean\tMAF\tCallRate\tMachR2\n");
}
if(duos==2)	//father
{
sprintf(filename4,"%s.swap.fathers",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tA1_Mean\tMAF\tCallRate\tMachR2\n");
}
if(duos==3)	//mother
{
sprintf(filename4,"%s.swap.mothers",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tA1_Mean\tMAF\tCallRate\tMachR2\n");
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
printf("Performing linear regression for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing linear regression for Chunk %d of %d\n", bit+1, bittotal);
}

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}

//read data (for all individuals), compute statistics, centre and set missing to zero
current=read_data_fly(datafile, dtype, data, NULL, 2*num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, 2*num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1);

//regress covariates out of both sets of data (necessary even if only one fixed effect)
reg_covar_matrix(data, Z, num_samples_use, 2*bitlength, num_fixed);

//get YTdata for both sets of data
token=2*bitlength;
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &token, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

//get XTXC and inverse for all predictors (set mults to -9999 if cholesky fails)
wcount=0;
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
//get XTCX
alpha=1.0;beta=0.0;
dgemm_("T", "N", &two, &two, &num_samples_use, &alpha, data+(size_t)2*j*num_samples_use, &num_samples_use, data+(size_t)2*j*num_samples_use, &num_samples_use, &beta, XTCX+4*j, &two);

//copy then invert
for(k=0;k<4;k++){XTCX2[k+4*j]=XTCX[k+4*j];}
(void)cholesky_invert(XTCX2+4*j, 2, -1, NULL, &info, 0);

if(info!=0)
{
if(wcount<5){printf("Warning, Predictor %s will be ignored due to colinearity with parental predictors\n", preds[bitstart+j]);}
mults[bitstart+j]=-9999;
wcount++;
}
}
}
if(wcount>5){printf("In total, %d predictors will be ignored due to colinearity with parental predictors\n", wcount);}
if(wcount>0){printf("\n");}

if(passqc==1)	//perform additional qc
{
#pragma omp parallel for private(j, maf) schedule(static)
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

////////

//ready to test
#pragma omp parallel for private(j, i, mark, mean, mean2, var, var2) schedule(static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//cholesky worked
{
//test basic model (just offspring)
mark=4*j;

//effect size estimate is inv(XTX) XTY
mean=YTdata[2*j]/XTCX[4*j];

//fill residuals
for(i=0;i<num_samples_use;i++){data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)2*j*num_samples_use+i];}

//compute variance
var=club_sandwich_one(XTCX[4*j], data+(size_t)2*j*num_samples_use, data2+(size_t)j*num_samples_use, num_fams, famindex, num_samples_use, num_fixed);

stats[0+mark]=mean;
stats[1+mark]=pow(var,.5);
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);

//test full model (offspring and parent)
mark=8*j;

//effect size estimates are inv(XTX) XTY
mean=XTCX2[0+4*j]*YTdata[2*j]+XTCX2[1+4*j]*YTdata[2*j+1];
mean2=XTCX2[2+4*j]*YTdata[2*j]+XTCX2[3+4*j]*YTdata[2*j+1];

//fill residuals
for(i=0;i<num_samples_use;i++)
{data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)2*j*num_samples_use+i]-mean2*data[(size_t)(2*j+1)*num_samples_use+i];}

//compute variances
club_sandwich_two(XTCX2+4*j, data+(size_t)2*j*num_samples_use, data2+(size_t)j*num_samples_use, &var, &var2, num_fams, famindex, num_samples_use, num_fixed, 0);

stats2[0+mark]=mean;
stats2[1+mark]=pow(var,.5);
stats2[2+mark]=stats2[0+mark]/stats2[1+mark];
stats2[3+mark]=erfc(fabs(stats2[2+mark])*M_SQRT1_2);

stats2[4+mark]=mean2;
stats2[5+mark]=pow(var2,.5);
stats2[6+mark]=stats2[4+mark]/stats2[5+mark];
stats2[7+mark]=erfc(fabs(stats2[6+mark])*M_SQRT1_2);

//test swap model (just parent)
mark=4*j;

//effect size estimate is inv(XTX) XTY
mean=YTdata[2*j+1]/XTCX[3+4*j];

//fill residuals
for(i=0;i<num_samples_use;i++){data2[(size_t)j*num_samples_use+i]=Yadj[i]-mean*data[(size_t)(2*j+1)*num_samples_use+i];}

//compute variance
var=club_sandwich_one(XTCX[3+4*j], data+(size_t)(2*j+1)*num_samples_use, data2+(size_t)j*num_samples_use, num_fams, famindex, num_samples_use, num_fixed);

stats3[0+mark]=mean;
stats3[1+mark]=pow(var,.5);
stats3[2+mark]=stats3[0+mark]/stats3[1+mark];
stats3[3+mark]=erfc(fabs(stats3[2+mark])*M_SQRT1_2);
}
}	//end of j loop

//save results
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)	//tested predictor
{
//basic model
mark=4*j;
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4f\t%.4e\t%.4e\t%.4e\t", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
fprintf(output2, "%.6f\t", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\n", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\n", rates[bitstart+j], infos[bitstart+j]);}

//full model
mark=8*j;
fprintf(output3, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output3, "%.4f\t%.4e\t%.4e\t%.4e\t", stats2[0+mark], stats2[1+mark], stats2[2+mark], stats2[3+mark]);
fprintf(output3, "%.4f\t%.4e\t%.4e\t%.4e\n", stats2[4+mark], stats2[5+mark], stats2[6+mark], stats2[7+mark]);

//swap model
mark=4*j;
fprintf(output4, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output4, "%.4f\t%.4e\t%.4e\t%.4e\t", stats3[2+mark], stats3[3+mark], stats3[0+mark], stats3[1+mark]);
fprintf(output4, "%.6f\t", centres[bitstart+j]);
if(nonsnp==0){fprintf(output4, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output4, "NA\t");}
if(genprobs<2){fprintf(output4, "%.4f\tNA\n", rates[bitstart+j]);}
else{fprintf(output4, "%.4f\t%.2f\n", rates[bitstart+j], infos[bitstart+j]);}
}
else	//did not test, but will include in results if not doing qc
{
if(passqc==0)
{
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "NA\tNA\tNA\tNA\t");
fprintf(output2, "%.6f\t", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\n", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\n", rates[bitstart+j], infos[bitstart+j]);}

fprintf(output3, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output3, "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n");

fprintf(output4, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output4, "NA\tNA\tNA\tNA\t");
fprintf(output4, "%.6f\t", centres[bitstart+j]);
if(nonsnp==0){fprintf(output4, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output4, "NA\t");}
if(genprobs<2){fprintf(output4, "%.4f\tNA\n", rates[bitstart+j]);}
else{fprintf(output4, "%.4f\t%.2f\n", rates[bitstart+j], infos[bitstart+j]);}
}
}
}	//end of j loop
}	//end of bit loop
if(wcount==0){printf("\n");}

fclose(output);
fclose(output2);
fclose(output3);
fclose(output4);

count=0;for(j=0;j<num_preds_use;j++){count+=(mults[j]==-9999);}
if(count==num_preds_use)
{
if(passqc==0)
{printf("Warning, all %d predictors were excluded due to colinearity with parental predictors or because they were trivial\n\n", count);}
else
{printf("Warning, all %d predictors were excluded due to colinearity with parental predictors, or because they were trivial or failed quality control\n\n", count);}
}
else
{
if(count>0)
{
if(passqc==0)
{printf("Warning, %d predictors were excluded due to colinearity with parental predictors or because they were trivial\n\n", count);}
else
{printf("Warning, %d predictors were excluded due to colinearity with parental predictors, or because they were trivial or failed quality control\n\n", count);}
}
}

//write details
sprintf(filename5,"%s.details", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
if(duos==1){fprintf(output5,"Analysis\tParent-Offspring-Duos\nNum_Samples\t%d\n", num_samples_use);}
if(duos==2){fprintf(output5,"Analysis\tFather-Offspring-Duos\nNum_Samples\t%d\n", num_samples_use);}
if(duos==3){fprintf(output5,"Analysis\tMother-Offspring-Duos\nNum_Samples\t%d\n", num_samples_use);}
fprintf(output5,"Num_Predictors\t%d\nNum_Tested\t%d\n", data_length, data_length-count);
fclose(output5);

printf("Main results saved in %s, %s and %s, with a summary of the analysis in %s\n\n", filename2, filename3, filename4, filename5);

////////

free(data);free(data2);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(stats);free(stats2);free(stats3);
free(YTdata);free(XTCX);free(XTCX2);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


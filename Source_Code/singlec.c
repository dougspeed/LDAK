/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor linear regression with one enviro - might permute (no top-preds, spa, or sample weights)

///////////////////////////

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//allocate variables

data_warn2(bitsize,2*num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
data2=malloc(sizeof(double)*num_samples_use*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);

stats=malloc(sizeof(double)*11*bitsize);

YTdata=malloc(sizeof(double)*bitsize);
YTdata2=malloc(sizeof(double)*bitsize);
XTCX=malloc(sizeof(double)*bitsize);
XTCX2=malloc(sizeof(double)*bitsize);
XTCX3=malloc(sizeof(double)*bitsize);

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

//solve null model, get thetas, adjust response and get variance
reg_covar_lin(Y, Z, num_samples_use, num_covars+1, 0, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/num_samples_use;varphen=sumsq/num_samples_use-pow(mean,2);

//save coefficients
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"ws"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component\tEffect\tSE\tP\n");
fprintf(output, "Intercept\t%.4e\t%.4e\t%.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d\t%.4e\t%.4e\t%.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fprintf(output, "Enviromental\t%.4e\t%.4e\t%.4e\n", thetas[num_covars], thetasds[num_covars], thetapvas[num_covars]);
fclose(output);

////////

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.interactions",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tMain_Effect\tMain_SE\tMain_P\tPred_Effect\tPred_SE\tPred_P\tInt_Effect\tInt_SE\tInt_P\tLRT_Stat\tLRT_P\tA1_Mean\tMAF\tCallRate\tMachR2\n");

sprintf(filename3,"%s.main.pvalues",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tP\n");

sprintf(filename4,"%s.int.pvalues",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor\tP\n");

sprintf(filename5,"%s.full.pvalues",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor\tP\n");

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
printf("Performing interaction tests for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing interaction tests for Chunk %d of %d\n", bit+1, bittotal);
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
fclose(output5);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//read data, compute statistics, centre and set missing to zero
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1);

//fill data2 (data multipled by environmental), then centre and check whether trivial
wcount=0;
#pragma omp parallel for private(j, i, sum, sumsq, mean, var) schedule(static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++)
{
data2[(size_t)j*num_samples_use+i]=data[(size_t)j*num_samples_use+i]*Z[i+num_covars*num_samples_use];
sum+=data2[(size_t)j*num_samples_use+i];
sumsq+=pow(data2[(size_t)j*num_samples_use+i],2);
}
mean=sum/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);
if(var>0)
{
for(i=0;i<num_samples_use;i++){data2[(size_t)j*num_samples_use+i]-=mean;}
}
else	//becomes trivial
{
#pragma omp critical
{
if(wcount<5){printf("Warning, Predictor %s becomes trivial when muliplied by the environmental, so will be ignored\n", preds[bitstart+j]);}
wcount++;
}
mults[bitstart+j]=-9999;
}
}}
if(wcount>5){printf("In total, %d predictors become trivial when multiplied by the environmental, so will be ignored\n", wcount);}

//regress covariates out of data, get YTdata and YTdata2 and XTCXs
reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);
reg_covar_matrix(data2, Z, num_samples_use, bitlength, num_fixed);

alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata, &one);
dgemv_("T", &num_samples_use, &bitlength, &alpha, data2, &num_samples_use, Yadj, &one, &beta, YTdata2, &one);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=ddot_(&num_samples_use, data+(size_t)j*num_samples_use, &one, data+(size_t)j*num_samples_use, &one);
XTCX2[j]=ddot_(&num_samples_use, data2+(size_t)j*num_samples_use, &one, data2+(size_t)j*num_samples_use, &one);
XTCX3[j]=ddot_(&num_samples_use, data+(size_t)j*num_samples_use, &one, data2+(size_t)j*num_samples_use, &one);
}

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

//ready to test
wcount=0;
for(j=0;j<bitlength;j++)
{
mark=11*j;

if(mults[bitstart+j]!=-9999)	//will be testing
{
//compute determinant of DTD
value=XTCX[j]*XTCX2[j]-pow(XTCX3[j],2);

if(value==0)	//must have X=XE, so can't test
{
if(wcount<5){printf("Warning, Predictor %s is unchanged when muliplied by the environmental, so will be ignored\n", preds[bitstart+j]);}
mults[bitstart+j]=-9999;
}
else
{
invDTD[0]=XTCX2[j]/value;
invDTD[1]=-XTCX3[j]/value;
invDTD[2]=invDTD[1];
invDTD[3]=XTCX[j]/value;

//get effect sizes
value=invDTD[0]*YTdata[j]+invDTD[2]*YTdata2[j];
value2=invDTD[1]*YTdata[j]+invDTD[3]*YTdata2[j];

//get residual sumsq and variances of effects
sumsq=varphen*num_samples_use-YTdata[j]*value-YTdata2[j]*value2;
var=sumsq/(num_samples_use-num_fixed-2)*invDTD[0];
var2=sumsq/(num_samples_use-num_fixed-2)*invDTD[3];

//save estimates for each effect
stats[0+mark]=value;
stats[1+mark]=pow(var,.5);
stats[2+mark]=erfc(fabs(stats[0+mark]/stats[1+mark])*M_SQRT1_2);
stats[3+mark]=value2;
stats[4+mark]=pow(var2,.5);
stats[5+mark]=erfc(fabs(stats[3+mark]/stats[4+mark])*M_SQRT1_2);

//perform LRT for joint model - pvalue is (1-Igamma(stat*beta, alpha)), with alpha=dof/2, beta=1/2
stats[6+mark]=num_samples_use*log(varphen*num_samples_use/sumsq);
stats[7+mark]=(1-gamain(stats[6+mark]/2, 1.0, &info));
if(stats[7+mark]==0||info>2)	//typically means lrt too large, so use pchisq(2*beta*stat,1)
{stats[7+mark]=erfc(pow(stats[6+mark],.5)*M_SQRT1_2);}

//get estimates for just predictor
value=YTdata[j]/XTCX[j];
sumsq=varphen*num_samples_use-YTdata[j]*value;
var=sumsq/(num_samples_use-num_fixed-1)/XTCX[j];

stats[8+mark]=value;
stats[9+mark]=pow(var,.5);
stats[10+mark]=erfc(fabs(stats[8+mark]/stats[9+mark])*M_SQRT1_2);
}
}	//end of testing
}	//end of j loop
if(wcount>5){printf("In total, %d predictors are unchanged when multipled by the environmental, so will be ignored\n", wcount);}
if(wcount>0){printf("\n");}

//save results
for(j=0;j<bitlength;j++)
{
mark=11*j;

if(mults[bitstart+j]!=-9999)	//tested predictor
{
//print interactions
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4f\t%.4e\t", stats[8+mark], stats[9+mark], stats[10+mark], stats[0+mark], stats[1+mark], stats[2+mark], stats[3+mark], stats[4+mark], stats[5+mark], stats[6+mark], stats[7+mark]);
fprintf(output2, "%.6f\t", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\n", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\n", rates[bitstart+j], infos[bitstart+j]);}

//print pvalues
fprintf(output3, "%s %.4e\n", preds[j], stats[10+mark]);
fprintf(output4, "%s %.4e\n", preds[j], stats[5+mark]);
fprintf(output5, "%s %.4e\n", preds[j], stats[7+mark]);
}
else	//did not test, but will include in interactions if not doing qc
{
if(passqc==0)
{
fprintf(output2, "%d\t%s\t%.0f\t%s\t%s\t", chr[bitstart+j], preds[bitstart+j], bp[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
fprintf(output2, "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t");
fprintf(output2, "%.6f\t", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\n", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\n", rates[bitstart+j], infos[bitstart+j]);}
}
}
}	//end of j loop

}	//end of bit loop
printf("\n");

fclose(output);
fclose(output2);
fclose(output3);
fclose(output4);
fclose(output5);

count=0;for(j=0;j<num_preds_use;j++){count+=(mults[j]==-9999);}
if(count==num_preds_use)
{
if(passqc==0)
{printf("Warning, all %d predictors were excluded when multiplied by the environmental or because they were trivial\n\n", count);}
else
{printf("Warning, all %d predictors were excluded when multiplied by the environmental, or because they were trivial or failed quality control\n\n", count);}
}
else
{
if(count>0)
{
if(passqc==0)
{printf("Warning, %d predictors were excluded when multiplied by the environmental or because they were trivial\n\n", count);}
else
{printf("Warning, %d predictors were excluded when multiplied by the environmental, or because they were trivial or failed quality control\n\n", count);}
}
}

printf("Main results saved in %s, with p-values in %s and %s\n\n", filename2, filename3, filename4);

free(data);free(data2);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(stats);
free(YTdata);free(YTdata2);free(XTCX);free(XTCX2);free(XTCX3);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Old-fashioned mixed-model linear regression - might permute (no top-preds, enviro, spa or sample weights)
//if exact==0, then have nk=1 and will regress UTY on UTZ using weights that are the inverse root of v1D+v0
//if permuting, we permute predictors when reading (whereas normally we permute phenotype and covariates)

///////////////////////////

//passqc indicates whether doing any qc
passqc=(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999||mininfo!=-9999);

//allocate variables

data_warn2(bitsize,2*num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
datasqs=malloc(sizeof(double)*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*(num_fixed+1));

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

UTY=malloc(sizeof(double)*num_samples_use);
UTZ=malloc(sizeof(double)*num_samples_use*(num_fixed+1));
UTdata=malloc(sizeof(double)*num_samples_use*bitsize);
vstarts=malloc(sizeof(double)*4);

stats=malloc(sizeof(double)*4*bitsize);

if(exact==0)
{
sweights=malloc(sizeof(double)*num_samples_use);
Yadj=malloc(sizeof(double)*num_samples_use);
YTdata=malloc(sizeof(double)*bitsize);
XTCX=malloc(sizeof(double)*bitsize);
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order (used only when reading predictors)
for(i=0;i<num_samples_use;i++){order[i]=keepsamps[i];}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y (not permuted)
for(i=0;i<num_samples_use;i++){Y[i]=resp[i];}

//fill start of Z (not permuted)
for(j=0;j<num_fixed;j++)
{
for(i=0;i<num_samples_use;i++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//get UTY and start of UTZ
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_samples_use, &alpha, U, &num_samples_use, Y, &one, &beta, UTY, &one);
dgemm_("T", "N", &num_samples_use, &num_fixed, &num_samples_use, &alpha, U, &num_samples_use, Z, &num_samples_use, &beta, UTZ, &num_samples_use);

//solve null model and get starting variances (and more things if exact=0)
printf("Solving Null Model\n");
linear_reml(num_samples_use, num_fixed, Y, Z, U, E, UTY, UTZ, kintraces, NULL, vstarts, thetas, thetasds, thetapvas, constrain, tol, maxiter, 1);
printf("\n");

if(dougvar2!=-9999)
{
printf("vstarts are %f %f and %f %f\n", vstarts[0], vstarts[1], vstarts[2], vstarts[3]);
sum=vstarts[2]+vstarts[3];
vstarts[3]=dougvar2*sum;
vstarts[2]=sum-vstarts[3];
printf("vstarts are %f %f and %f %f\n", vstarts[0], vstarts[1], vstarts[2], vstarts[3]);
}

if(exact==0)	//compute weights, then scale UTZ and UTY, adjust response and get variance
{
for(i=0;i<num_samples_use;i++)
{
if(vstarts[2]+vstarts[3]*E[i]<=0){printf("Error, variances %f %f for ind %d is %f\n", vstarts[2], vstarts[3], i+1, E[i]);exit(1);}
sweights[i]=pow(vstarts[2]+vstarts[3]*E[i],-.5);
}

for(i=0;i<num_samples_use;i++)
{
UTY[i]*=sweights[i];
for(j=0;j<num_fixed;j++){UTZ[i+j*num_samples_use]*=sweights[i];}
}

//the following is slightly redundant, because it recomputes thetas (values will not change)
reg_covar_lin(UTY, UTZ, num_samples_use, num_covars, 0, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/num_samples_use;varphen=sumsq/num_samples_use-pow(mean,2);
}

//save coefficients
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"w"))==NULL)
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

sprintf(filename2,"%s.assoc",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome\tPredictor\tBasepair\tA1\tA2\tWald_Stat\tWald_P\tEffect\tSE\tMAF\tCallRate\tMachR2\n");

sprintf(filename3,"%s.summaries",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tA1\tA2\tZ\tn\tA1Freq\n");

sprintf(filename4,"%s.pvalues",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor\tP\n");

sprintf(filename5,"%s.score",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor\tA1\tA2\tCentre\tALL\tP<0.1\tP<0.01\tP<0.001\tP<0.0001\tP<0.00001\tP<5e-8\n");

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
wcount=0;
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
fclose(output5);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//read data, compute statistics, centre and set missing to zero
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, order, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1);

//get UTdata
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_samples_use, &bitlength, &num_samples_use, &alpha, U, &num_samples_use, data, &num_samples_use, &beta, UTdata, &num_samples_use);

if(exact==0)	//scale UTdata by sweights, regress out covariates, get YTdata and XTCXs
{
copy_matrix(num_samples_use, bitlength, UTdata, UTdata, 1, sweights);
reg_covar_matrix(UTdata, UTZ, num_samples_use, bitlength, num_fixed);

alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, UTdata, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

for(j=0;j<bitlength;j++)
{
XTCX[j]=ddot_(&num_samples_use, UTdata+(size_t)j*num_samples_use, &one, UTdata+(size_t)j*num_samples_use, &one);
}
}

if(exact==0)	//do not test trivial predictors (have likely already been excluded)
{
#pragma omp parallel for private(j) schedule(static)
for(j=0;j<bitlength;j++)
{
if(XTCX[j]==0){mults[bitstart+j]=-9999;}
}
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
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//will be testing
{
if(exact==0)	//simple case
{
//compute YTCX/XTCX and YCT(YC-XCvalue)/XTCX(n-nf-1)
value=YTdata[j]/XTCX[j];
value2=(varphen*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);

stats[0+mark]=value;
stats[1+mark]=pow(value2,.5);
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
else	//complicated case
{
//add test predictor to end of Z, add UTdata for test predictor to end of UTZ, then test
for(i=0;i<num_samples_use;i++){Z[i+num_fixed*num_samples_use]=data[(size_t)j*num_samples_use+i];}
for(i=0;i<num_samples_use;i++){UTZ[i+num_fixed*num_samples_use]=UTdata[(size_t)j*num_samples_use+i];}
linear_reml(num_samples_use, num_fixed+1, Y, Z, U, E, UTY, UTZ, kintraces, stats+mark, vstarts, NULL, NULL, NULL, constrain, tol, maxiter, exact);
}
}	//end of testing
}	//end of j loop

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
fprintf(output2, "NOT_USED\n");

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
fprintf(output2, "NA NA NA NA ");
if(nonsnp==0){fprintf(output2, "%.6f\t", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA\t");}
if(genprobs<2){fprintf(output2, "%.4f\tNA\t", rates[bitstart+j]);}
else{fprintf(output2, "%.4f\t%.2f\t", rates[bitstart+j], infos[bitstart+j]);}
fprintf(output2, "NOT_USED\n");
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

printf("Main results saved in %s, with a summary version in %s, p-values in %s and score file in %s\n\n", filename2, filename3, filename4, filename5);

free(data);free(datasqs);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);
free(UTY);free(UTZ);free(UTdata);free(vstarts);
free(stats);
if(exact==0){free(sweights);free(Yadj);free(YTdata);free(XTCX);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


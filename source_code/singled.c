/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor logistic regression - might have top-preds and offsets (no sample weights, bootstrapping or enviro)
//can not have spatest=1 with scoretest=0

///////////////////////////

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*(num_fixed+1));

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4*bitsize);
spastatus=malloc(sizeof(int)*bitsize);

if(scoretest==1)
{
nullprobs=malloc(sizeof(double)*num_samples_use);
nullweights=malloc(sizeof(double)*num_samples_use);
Yadj=malloc(sizeof(double)*num_samples_use);
YTdata=malloc(sizeof(double)*bitsize);
XTCX=malloc(sizeof(double)*bitsize);
}

if(spatest==1){usedpreds=malloc(sizeof(int)*bitsize);}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y and start of Z
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use];}
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

//solve null model and get thetas (and maybe nullweights and W x adjusted response)
if(scoretest==0)
{reg_covar_log(Y, Z, num_samples_use, num_covars, num_tops, offsets, thetas, thetasds, thetapvas, NULL, -9999, NULL, NULL, -9999, 0.001, 100);}
else
{
reg_covar_log(Y, Z, num_samples_use, num_covars, num_tops, offsets, thetas, thetasds, thetapvas, nullprobs, 2, NULL, NULL, -9999, 0.001, 100);
for(i=0;i<num_samples_use;i++){nullweights[i]=nullprobs[i]*(1-nullprobs[i]);}
for(i=0;i<num_samples_use;i++){Yadj[i]=(Y[i]-nullprobs[i]);}
}

//save
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Log_OR SD P\n");
fprintf(output, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

if(spatest==0)	//set spastatus to zero
{
for(j=0;j<bitsize;j++){spastatus[j]=0;}
}

////////

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.assoc",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
if(scoretest==0){fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Log_OR SD A1_Mean MAF SPA_Status\n");}
else{fprintf(output2, "Chromosome Predictor Basepair A1 A2 Score_Stat Score_P Approx_Log_OR Approx_SD A1_Mean MAF SPA_Status\n");}

sprintf(filename3,"%s.summaries",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Z n\n");

sprintf(filename4,"%s.pvalues",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor P\n");

sprintf(filename5,"%s.score",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor A1 A2 Centre ALL P<0.1 P<0.01 P<0.001 P<0.0001 P<0.00001 P<5e-8\n");

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

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output5);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//read data for chunk and centre (and replace NAs with zero)
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1, preds+bitstart);

if(scoretest==1)	//maybe weighted regress covariates out of data, get YTdata and XTWXs
{
if(num_fixed>1){reg_covar_weighted(data, Z, num_samples_use, bitlength, num_fixed, nullweights);}

alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2)*nullweights[i];}
}
}

//ready to test
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//will be testing
{
if(tindex[j]!=-9999)	//a top predictor, so have already tested
{
stats[0+mark]=thetas[num_covars+tindex[bitstart+j]];
stats[1+mark]=thetasds[num_covars+tindex[bitstart+j]];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=thetapvas[num_covars+tindex[bitstart+j]];
}
else	//not a top, so must test
{
if(scoretest==0)	//add test predictor to end of Z and test
{
for(i=0;i<num_samples_use;i++){Z[i+num_fixed*num_samples_use]=data[(size_t)j*num_samples_use+i];}
reg_single_log(Y, Z, num_samples_use, num_fixed+1, offsets, stats+mark, thetas, 0.001, 100);
}
else	//have already computed score and variance (saved in YTdata and XTCX)
{
//the score statistic is YTdata, with variance XTCX
//SAIGE paper explains how score/variance approx equals LogOR, and 1/var approx equals Var of LogOR
//this is because XT(Y-mu)/XTWX is 1-step NR estimate of LogOR
stats[0+mark]=YTdata[j]/XTCX[j];
stats[1+mark]=pow(XTCX[j],-.5);

stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}	//end of not a top
}	//end of testing
}	//end of j loop

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

if(spaside==1){spastatus[j]=spa_logistic_one(YTdata[j], data+j*num_samples_use, num_samples_use, nullprobs, nullweights, stats+mark, 1.0);}
else{spastatus[j]=spa_logistic_two(YTdata[j], XTCX[j], data+j*num_samples_use, num_samples_use, nullprobs, nullweights, stats+mark, 1.0);}
}
}

//save results
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//include in all results
{
//print assoc
fprintf(output2, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output2, "%.4f %.4e %.4e %.4e ", stats[2+mark], stats[3+mark], stats[0+mark], stats[1+mark]);
fprintf(output2, "%.6f ", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f ", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA ");}
if(spastatus[j]==0){fprintf(output2, "NOT_USED\n");flag=1;}
if(spastatus[j]==1){fprintf(output2, "SUCCESS\n");flag=1;}
if(spastatus[j]==2){fprintf(output2, "APPROX\n");flag=1;}
if(spastatus[j]==-1){fprintf(output2, "FAILED\n");flag=1;}

//print summaries, pvalues and scores
fprintf(output3, "%s %c %c %.4f %d\n", preds[bitstart+j], al1[bitstart+j], al2[bitstart+j], stats[2+mark], num_samples_use);
fprintf(output4, "%s %.4e\n", preds[bitstart+j], stats[3+mark]);
fprintf(output5, "%s %c %c %.6f %.6f", preds[bitstart+j], al1[bitstart+j], al2[bitstart+j],  centres[bitstart+j], stats[0+mark]);
for(k=0;k<6;k++)
{
if(stats[3+mark]<cuts[k]){fprintf(output5, " %.6f", stats[0+mark]);}
else{fprintf(output5, " 0");}
}
fprintf(output5, "\n");
}
else	//trivial predictor - so include only in assoc
{
fprintf(output2, "%d %s %.0f %c %c ", chr[bitstart+j], preds[bitstart+j], cmbp[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
fprintf(output2, "NA NA NA NA ");
fprintf(output2, "%.6f ", centres[bitstart+j]);
if(nonsnp==0){fprintf(output2, "%.6f NOT_USED\n", centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]));}
else{fprintf(output2, "NA NOT_USED\n");}
}
}	//end of j loop
}	//end of bit loop
printf("\n");

fclose(output);
fclose(output2);
fclose(output3);
fclose(output4);
fclose(output5);

printf("Main results saved in %s, with a summary version in %s, p-values in %s and score file in %s\n\n", filename2, filename3, filename4, filename5);

free(data);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);
free(tindex);free(stats);free(spastatus);
if(scoretest==1){free(nullprobs);free(nullweights);free(Yadj);free(YTdata);free(XTCX);}
if(spatest==1){free(usedpreds);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


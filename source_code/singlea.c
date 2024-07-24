/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Single-predictor linear regression - might have top-preds, sample weights or bootstrapping (but no enviro)
//if exact==0, then have nk=1 and will regress UTY on UTZ using weights that are the inverse root of v1D+v0

///////////////////////////

//allocate variables

if(num_kins==0){data_warn2(bitsize,num_samples_use);}
else{data_warn2(bitsize,2*num_samples_use);}

data=malloc(sizeof(double)*num_samples_use*bitsize);
datasqs=malloc(sizeof(double)*bitsize);

sweights=malloc(sizeof(double)*num_samples_use);
order=malloc(sizeof(int)*num_samples_use);
keepsamps2=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*(num_fixed+1));

thetas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4*bitsize);
spastatus=malloc(sizeof(int)*bitsize);

if(num_kins==0||exact==0)
{
YTdata=malloc(sizeof(double)*bitsize);
XTCX=malloc(sizeof(double)*bitsize);
}

if(num_kins==1)
{
kweights=malloc(sizeof(double)*num_samples_use);
UTY=malloc(sizeof(double)*num_samples_use);
UTZ=malloc(sizeof(double)*num_samples_use*(num_fixed+1));
UTdata=malloc(sizeof(double)*num_samples_use*bitsize);
vstarts=malloc(sizeof(double)*4);
}

if(spatest==1)
{
knots=malloc(sizeof(double)*num_knots);
bins=malloc(sizeof(double)*num_bins);
CGF0=malloc(sizeof(double*)*num_knots);
CGF1=malloc(sizeof(double*)*num_knots);
CGF2=malloc(sizeof(double*)*num_knots);
CGF3=malloc(sizeof(double*)*num_knots);

for(j=0;j<num_knots;j++)
{
CGF0[j]=malloc(sizeof(double)*num_bins);
CGF1[j]=malloc(sizeof(double)*num_bins);
CGF2[j]=malloc(sizeof(double)*num_bins);
CGF3[j]=malloc(sizeof(double)*num_bins);
}

usedpreds=malloc(sizeof(int)*bitsize);
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

if(strcmp(sampwfile,"blank")!=0)	//get sample weights
{
read_sampwfile(sampwfile, sweights, num_samples_use, ids3, ids1, ids2);

//scale roots so the mean inverse is one
sum=0;for(i=0;i<num_samples_use;i++){sum+=pow(sweights[i],-1);}
mean=sum/num_samples_use;
for(i=0;i<num_samples_use;i++){sweights[i]*=mean;}

//will only use square root below
for(i=0;i<num_samples_use;i++){sweights[i]=pow(sweights[i],0.5);}
}
else
{
for(i=0;i<num_samples_use;i++){sweights[i]=1;}
}

//fill order and keepsamps2 assuming not permuting or bootstrapping
for(i=0;i<num_samples_use;i++){order[i]=i;}
for(i=0;i<num_samples_use;i++){keepsamps2[i]=keepsamps[i];}

if(permute==1)	//only change order
{
permute_int(order,num_samples_use);
}

if(booty==1)	//change order and keepsamps2
{
for(i=0;i<num_samples_use;i++){order[i]=rand()%num_samples_use;keepsamps2[i]=keepsamps[order[i]];}
}

//fill Y and start of Z
for(i=0;i<num_samples_use;i++)
{
Y[i]=resp[order[i]]*sweights[order[i]];
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[order[i]+j*num_samples_use]*sweights[order[i]];}
}

if(num_kins==1)	//get UTY and start of UTZ
{
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &num_samples_use, &alpha, U, &num_samples_use, Y, &one, &beta, UTY, &one);
dgemm_("T", "N", &num_samples_use, &num_fixed, &num_samples_use, &alpha, U, &num_samples_use, Z, &num_samples_use, &beta, UTZ, &num_samples_use);
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

if(num_kins==0)	//solve null model, get thetas, adjust response and get variance
{
reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
varphen=sumsq/num_samples_use-pow(sum/num_samples_use,2);
}
else	//solve null model and get starting variances (and more things if exact=0)
{
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
kweights[i]=pow(vstarts[2]+vstarts[3]*E[i],-.5);
}

for(i=0;i<num_samples_use;i++)
{
UTY[i]*=kweights[i];
for(j=0;j<num_fixed;j++){UTZ[i+j*num_samples_use]*=kweights[i];}
}

//the following is slightly redundant, because it recomputes thetas (values will not change)
reg_covar_lin(UTY, UTZ, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/num_samples_use;
varphen=sumsq/num_samples_use-pow(mean,2);
}
}

//save
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Component Effect SD P\n");
fprintf(output, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

if(spatest==0)	//set spastatus to zero
{
for(j=0;j<bitsize;j++){spastatus[j]=0;}
}
else	//prepare for spa
{
//set bins - evenly spaced from -2 to 2
for(k=0;k<num_bins;k++){bins[k]=-2+(double)k/(num_bins-1)*4;}

if(spamax==-9999)	//set spamax based on variance of Yadj (will already have variance)
{spamax=2000*pow(num_samples_use*varphen,-.5);}

printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n\n", num_knots, -spamax, spamax, num_bins);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
}

////////

//deal with progress and on-the-fly files

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.assoc",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Effect SD A1_Mean MAF SPA_Status\n");

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
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
fclose(output5);
if((output5=fopen(filename5,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename5);exit(1);}

//read data for chunk and centre (and replace NAs with zero)
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps2, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1, preds+bitstart);

if(num_kins==0)	//maybe scale data and regress out covariates, get YTdata and XTCXs
{
if(strcmp(sampwfile,"blank")!=0)
{copy_matrix(num_samples_use, bitlength, data, data, 1, sweights);}

if(num_fixed>1||strcmp(sampwfile,"blank")!=0)
{reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}

alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}
}
}
else	//get UTdata (and more things if exact=0)
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_samples_use, &bitlength, &num_samples_use, &alpha, U, &num_samples_use, data, &num_samples_use, &beta, UTdata, &num_samples_use);

if(exact==0)	//scale UTdata by kweights, regress out covariates, get YTdata and XTCXs
{
copy_matrix(num_samples_use, bitlength, UTdata, UTdata, 1, kweights);

reg_covar_matrix(UTdata, UTZ, num_samples_use, bitlength, num_fixed);

alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, UTdata, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(UTdata[(size_t)j*num_samples_use+i],2);}
}
}
}

//ready to test
for(j=0;j<bitlength;j++)
{
mark=4*j;

if(mults[bitstart+j]!=-9999)	//will be testing
{
if(tindex[bitstart+j]!=-9999)	//a top predictor, so have already tested
{
stats[0+mark]=thetas[num_covars+tindex[bitstart+j]];
stats[1+mark]=thetasds[num_covars+tindex[bitstart+j]];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=thetapvas[num_covars+tindex[bitstart+j]];
}
else	//not a top, so must test
{
if(num_kins==0)	//basic linear regression (maybe using sandwich estimates)
{
if(sandwich==0)	//use standard estimates
{
//compute YTCX/XTCX and YCT(YC-XCvalue)/XTCX(n-nf-1)
value=YTdata[j]/XTCX[j];
value2=(varphen*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);

stats[0+mark]=value;
stats[1+mark]=pow(value2,.5);
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
else	//use sandwich estimates
{
//effect size remains YTCX/XTCX
value=YTdata[j]/XTCX[j];

//variance is sum(Xi^2 ri^2)/XTCX^2, where ri are residuals
value3=0;
for(i=0;i<num_samples_use;i++)
{value3+=pow(data[(size_t)j*num_samples_use+i]*(Yadj[i]-value*data[(size_t)j*num_samples_use+i]),2);}
value2=value3*pow(XTCX[j],-2);

stats[0+mark]=value;
stats[1+mark]=pow(value2,.5);
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
}
else	//mixed model
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

while(1)
{
#pragma omp parallel for private(j2, j, mark) schedule(static)
for(j2=0;j2<count;j2++)
{
j=usedpreds[j2];
mark=4*j;

spastatus[j]=spa_test(YTdata[j], data+(size_t)j*num_samples_use, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, stats+mark, 1.0);
}

//see if necessary to increase range for any predictors (and recompute SPA)
count=0;
for(j=0;j<bitlength;j++)
{
if(spastatus[j]==-2){usedpreds[count]=j;count++;}
}
if(count==0){break;}

spamax*=5;
printf("Increasing SPA range to %f\n", spamax);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
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
else	//predictor not tested - so include only in assoc
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

free(data);free(datasqs);
free(sweights);free(order);free(keepsamps2);free(Y);free(Z);
free(thetas);free(thetasds);free(thetapvas);free(Yadj);
free(tindex);free(stats);free(spastatus);
if(num_kins==0||exact==0){free(YTdata);free(XTCX);}
if(num_kins==1){free(kweights);free(UTY);free(UTZ);free(UTdata);free(vstarts);}
if(spatest==1)
{
for(j=0;j<num_knots;j++){free(CGF0[j]);free(CGF1[j]);free(CGF2[j]);free(CGF3[j]);}
free(knots);free(bins);free(CGF0);free(CGF1);free(CGF2);free(CGF3);free(usedpreds);
}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//LOCO regression - might have top-preds (no enviro, sample weights or bootstrapping) - have set scal1, scal2
//performing linear regression (mode==131), quasi-logistic regression (mode==132&fastgwa==0) or logistic regression (mode==131&fastgwa!=0)

///////////////////////////

if(dougvar2!=-9999)
{
printf("Original scaling is %f divide by %f\n", scal1, 1.0/scal1);
scal1=dougvar2;
scal2=pow(scal1,.5);
printf("Revised scaling is %f divide by %f\n", scal1, 1.0/scal1);
}

//deal with prs file

wantids=malloc(sizeof(char*)*num_samples_use);
indexer=malloc(sizeof(int)*num_samples_use);

if(fastgwa==0)
{
chrindex=malloc(sizeof(int)*num_chr2);
prs=malloc(sizeof(double)*num_samples_use*num_chr2);
}
else
{
pedprs=malloc(sizeof(double)*num_samples_use);
}

//check ids
sprintf(filename,"%s.loco.prs", locofile);
read_ids(filename, NULL, NULL, wantids, num_samples_use, NULL, 1, 0);
count=find_strings(wantids, num_samples_use, ids3, num_samples_use, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count==0){printf("Error, none of the %d samples in %s are in %s\n\n", num_samples_use, filename, idsfile);exit(1);}
if(count<num_samples_use){printf("Error, only %d of the %d samples in %s are in %s\n\n", count, num_samples_use, blupfile, idsfile);exit(1);}

if(fastgwa==3)	//no need to read PRS (it is blank)
{
for(i=0;i<num_samples_use;i++){pedprs[i]=0;}
}
else	//read values
{
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fastgwa==0)
{
printf("Reading LOCO PRS from %s\n\n", filename);

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of %s\n\n", filename);exit(1);}

//read chr numbers
for(j=0;j<num_chr2;j++)
{
if(fscanf(input, "Chr%d ", chrindex+j)!=1)
{printf("Error reading Element %d of %s\n\n", j+3, filename);exit(1);}
}

//check all chr covered
usedpreds=malloc(sizeof(int)*(chrindex[num_chr2-1]+1));
for(j=0;j<chrindex[num_chr2-1]+1;j++){usedpreds[j]=0;}
for(j=0;j<num_chr2;j++){usedpreds[chrindex[j]]=1;}

for(j=0;j<data_length;j++)
{
if(chr[j]>chrindex[num_chr2-1])
{printf("Error, will be analyzing predictors on Chromosome %d, but this chromosome was not included when making the prediction model\n\n", chr[j]);exit(1);}
if(usedpreds[chr[j]]==0)
{printf("Error, will be analyzing predictors on Chromosome %d, but this chromosome was not included when making the prediction model\n\n", chr[j]);exit(1);}
}
free(usedpreds);

for(i=0;i<num_samples_use;i++)
{
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of Row %d of %s\n\n", i+2, filename);exit(1);}
for(j=0;j<num_chr2;j++)
{
if(fscanf(input, "%lf ", prs+indexer[i]+j*num_samples_use)!=1)
{printf("Error reading PRS %d from Row %d of %s\n\n", j+1, i+1, filename);exit(1);}
}
}
}
else
{
printf("Reading Pedigree PRS from %s\n\n", filename);

if(fscanf(input, "%s %s Pedigree", readstring, readstring2)!=2)
{printf("Error reading first row %s\n\n", filename);exit(1);}
	
for(i=0;i<num_samples_use;i++)
{
if(fscanf(input, "%s %s %lf ", readstring, readstring2, pedprs+indexer[i])!=3)
{printf("Error reading Row %d of %s\n\n", i+2, filename);exit(1);}
}
}
fclose(input);
}

for(i=0;i<num_samples_use;i++){free(wantids[i]);}free(wantids);free(indexer);

if(dougvar2==99)	//blank prs and set scaling to one
{
printf("blanking prs\n");
for(i=0;i<num_samples_use;i++)
{
if(fastgwa==0)
{
for(j=0;j<num_chr2;j++){prs[i+j*num_samples_use]=0;}
}
else{pedprs[i]=0;}
}
scal1=1;scal2=1;
printf("Revised scaling is one - this does not work for binary\n");
}

////////

//allocate variables

data_warn2(bitsize,num_samples_use);

data=malloc(sizeof(double)*num_samples_use*bitsize);
datasqs=malloc(sizeof(double)*bitsize);

order=malloc(sizeof(int)*num_samples_use);
Y=malloc(sizeof(double)*num_samples_use);
Z=malloc(sizeof(double)*num_samples_use*num_fixed);

thetas=malloc(sizeof(double)*num_fixed);
thetapvas=malloc(sizeof(double)*num_fixed);
thetasds=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*num_samples_use);
Yadj2=malloc(sizeof(double)*num_samples_use);

if(mode==132)
{
nullprobs=malloc(sizeof(double)*num_samples_use*2);
nullweights=malloc(sizeof(double)*num_samples_use*2);
}

tindex=malloc(sizeof(int)*data_length);
stats=malloc(sizeof(double)*4*bitsize);
spastatus=malloc(sizeof(int)*bitsize);

YTdata=malloc(sizeof(double)*bitsize);
XTCX=malloc(sizeof(double)*bitsize);

if(spatest==1)
{
usedpreds=malloc(sizeof(int)*bitsize);

if(mode==131||fastgwa==0)
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
}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with order
for(i=0;i<num_samples_use;i++){order[i]=i;}
if(permute==1){permute_int(order,num_samples_use);}

//fill Y and Z
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

if(mode==131)	//solve null model, get thetas and adjust response
{
reg_covar_lin(Y, Z, num_samples_use, num_covars, num_tops, thetas, thetasds, thetapvas, Yadj, 0, NULL, NULL);
}
else	//solve null model, get null weights, and compute W x adjusted response
{
reg_covar_log(Y, Z, num_samples_use, num_covars, num_tops, NULL, thetas, thetasds, thetapvas, nullprobs, 2, NULL, NULL, -9999, 0.001, 100);
for(i=0;i<num_samples_use;i++){nullweights[i]=nullprobs[i]*(1-nullprobs[i]);}
for(i=0;i<num_samples_use;i++){Yadj[i]=Y[i]-nullprobs[i];}
}

//save coefficients
sprintf(filename,"%s.coeff", outfile);
if((output=fopen(filename,"ws"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
if(mode==131){fprintf(output, "Component Effect SD P\n");}
else{fprintf(output, "Component LogOR SD P\n");}
fprintf(output, "Intercept %.6f %.6f %.4e\n", thetas[0], thetasds[0], thetapvas[0]);
for(j=1;j<num_covars;j++){fprintf(output, "Covariate_%d %.6f %.6f %.4e\n",j, thetas[j], thetasds[j], thetapvas[j]);}
fclose(output);

if(spatest==0)	//set spastatus to zero
{
for(j=0;j<bitsize;j++){spastatus[j]=0;}
}

if(spatest==1&&(mode==131||fastgwa==0))	//prepare for spa
{
//set bins - evenly spaced from -2 to 2
for(k=0;k<num_bins;k++){bins[k]=-2+(double)k/(num_bins-1)*4;}

if(spamax==-9999)	//set spamax based on variance of Yadj
{
sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
var=sumsq/num_samples_use-pow(sum/num_samples_use,2);
spamax=2000*pow(num_samples_use*var,-.5);
}
}

//store original Yadj (saves needing to recompute for each chromosome)
for(i=0;i<num_samples_use;i++){Yadj2[i]=Yadj[i];}

if(fastgwa==0)	//adjust response for first chromosome PRS, get variance and deal with SPA
{
mark2=0;
while(chr[0]!=chrindex[mark2]){mark2++;}

if(mode==131)	//adjusted phenotype is Yadj2 - PRS, get regular variance
{
for(i=0;i<num_samples_use;i++){Yadj[i]=Yadj2[i]-prs[order[i]+mark2*num_samples_use];}

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/num_samples_use;
varphen=sumsq/num_samples_use-pow(mean,2);

if(spatest==1)	//get cumulants
{
printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -spamax, spamax, num_bins);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
}
}
else	//adjusted phenotype is Yadj2 - W PRS, get variance of 1/root(W) Yadj
{
for(i=0;i<num_samples_use;i++){Yadj[i]=Yadj2[i]-prs[order[i]+mark2*num_samples_use]*nullweights[i];}

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i]*pow(nullweights[i],-.5);sumsq+=pow(Yadj[i],2)/nullweights[i];}
mean=sum/num_samples_use;
varphen=sumsq/num_samples_use-pow(mean,2);

if(spatest==1)	//get cumulants
{
printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -spamax, spamax, num_bins);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
}
}
}
else	//similar but with pedigree PRS and will use logistic regression
{
if(mode==131)
{
for(i=0;i<num_samples_use;i++){Yadj[i]=Yadj2[i]-pedprs[order[i]];}

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/num_samples_use;
varphen=sumsq/num_samples_use-pow(mean,2);

if(spatest==1)
{
printf("Computing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -spamax, spamax, num_bins);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
}
}
else
{
for(i=0;i<num_samples_use;i++){Yadj[i]=Yadj2[i]-pedprs[order[i]]*nullweights[i];}

if(spatest==1)	//convenient to get some extra probs and weights
{
for(i=0;i<num_samples_use;i++)
{
nullprobs[i+num_samples_use]=nullprobs[i]+pedprs[order[i]]*nullweights[i];
if(nullprobs[i+num_samples_use]<=0){printf("warn %f\n",nullprobs[i+num_samples_use]); nullprobs[i+num_samples_use]=1e-10;}
if(nullprobs[i+num_samples_use]>=1){printf("warn %f\n",nullprobs[i+num_samples_use]);nullprobs[i+num_samples_use]=1-1e-10;}
nullweights[i+num_samples_use]=nullprobs[i+num_samples_use]*(1-nullprobs[i+num_samples_use]);
}
}
}
}	//end of fastgwa!=0

////////

//deal with progress and on-the-fly files

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

sprintf(filename2,"%s.assoc",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
if(mode==131){fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Effect SD A1_Mean MAF SPA_Status\n");}
else{fprintf(output2, "Chromosome Predictor Basepair A1 A2 Wald_Stat Wald_P Approx_Log_OR Approx_SD A1_Mean MAF SPA_Status\n");}

sprintf(filename3,"%s.summaries",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Z n\n");

sprintf(filename4,"%s.pvalues",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4, "Predictor P\n");

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

if(fastgwa==0)	//see if changing chromosome
{
if(chr[bitstart]!=chrindex[mark2])	//must update phenotype, maybe variance and SPA
{
while(chr[bitstart]!=chrindex[mark2]){mark2++;}

if(mode==131)
{
for(i=0;i<num_samples_use;i++){Yadj[i]=Yadj2[i]-prs[order[i]+mark2*num_samples_use];}

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i];sumsq+=pow(Yadj[i],2);}
mean=sum/num_samples_use;
varphen=sumsq/num_samples_use-pow(mean,2);

if(spatest==1)
{
printf("Recomputing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -spamax, spamax, num_bins);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
}
}
else
{
for(i=0;i<num_samples_use;i++){Yadj[i]=Yadj2[i]-prs[order[i]+mark2*num_samples_use]*nullweights[i];}

sum=0;sumsq=0;for(i=0;i<num_samples_use;i++){sum+=Yadj[i]*pow(nullweights[i],-.5);sumsq+=pow(Yadj[i],2)/nullweights[i];}
mean=sum/num_samples_use;
varphen=sumsq/num_samples_use-pow(mean,2);

if(spatest==1)
{
if(fastgwa==0)
{
printf("Recomputing empirical CGF (%d knots between %f and %f, and %d bins between -2 and 2)\n", num_knots, -spamax, spamax, num_bins);
empirical_cumulants(Yadj, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, spamax);
}
else
{
for(i=0;i<num_samples_use;i++)
{
nullprobs[i+num_samples_use]=nullprobs[i]+prs[order[i]+mark2*num_samples_use]*nullweights[i];
if(nullprobs[i+num_samples_use]<=0){nullprobs[i+num_samples_use]=1e-10;}
if(nullprobs[i+num_samples_use]>=1){nullprobs[i+num_samples_use]=1-1e-10;}
nullweights[i+num_samples_use]=nullprobs[i+num_samples_use]*(1-nullprobs[i+num_samples_use]);
}
}
}
}
}}

if(bit%10==0)
{
if(mode==131)
{
printf("Performing linear regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing linear regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
}
else
{
printf("Performing quasi-logistic regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
fclose(output);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing quasi-logistic regression for Chunk %d of %d (Chromosome %d)\n", bit+1, bittotal, chr[bitstart]);
}
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

//read data for chunk and centre (and replace NAs with zero)
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 1, preds+bitstart);

if(num_fixed>1)	//(weighted) regress covariates out of data
{
if(mode==131){reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);}
else{reg_covar_weighted(data, Z, num_samples_use, bitlength, num_fixed, nullweights);}
}

//compute YTdata
alpha=1.0;beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, Yadj, &one, &beta, YTdata, &one);

if(mode==131)	//compute XTX
{
#pragma omp parallel for private(j, i) schedule(static)
for(j=0;j<bitlength;j++)
{
XTCX[j]=0;for(i=0;i<num_samples_use;i++){XTCX[j]+=pow(data[(size_t)j*num_samples_use+i],2);}
}
}
else	//compute XTWX
{
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
if(tindex[bitstart+j]!=-9999)	//a top predictor, so have already tested
{
stats[0+mark]=thetas[num_covars+tindex[bitstart+j]];
stats[1+mark]=thetasds[num_covars+tindex[bitstart+j]];
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
else	//not a top, so must test - remember to scale test statistic by scal2 (root lambda)
{
if(mode==131||fastgwa==0)
{
value=YTdata[j]/XTCX[j];
value2=(varphen*num_samples_use/XTCX[j]-pow(value,2))/(num_samples_use-num_fixed-1);
}
else
{
value=YTdata[j]/XTCX[j];
value2=pow(XTCX[j],-1);
}

stats[0+mark]=value*scal1;
stats[1+mark]=pow(value2,.5)*scal2;
stats[2+mark]=stats[0+mark]/stats[1+mark];
stats[3+mark]=erfc(fabs(stats[2+mark])*M_SQRT1_2);
}
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
else	//trivial, so will not test, but helps to set spastatus to missing
{spastatus[j]=-9999;}
}

while(1)
{
#pragma omp parallel for private(j2, j, mark) schedule(static)
for(j2=0;j2<count;j2++)
{
j=usedpreds[j2];
mark=4*j;

if(mode==131||fastgwa==0){spastatus[j]=spa_test(YTdata[j], data+j*num_samples_use, num_samples_use, num_knots, knots, num_bins, bins, CGF0, CGF1, CGF2, CGF3, stats+mark, scal2);}
else
{
if(spaside==1){spastatus[j]=spa_logistic_one(YTdata[j], data+j*num_samples_use, num_samples_use, nullprobs+num_samples_use, nullweights+num_samples_use, stats+mark, scal2);}
else{spastatus[j]=spa_logistic_two(YTdata[j], XTCX[j], data+j*num_samples_use, num_samples_use, nullprobs+num_samples_use, nullweights+num_samples_use, stats+mark, scal2);}
}
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

//print summaries and pvalues
fprintf(output3, "%s %c %c %.4f %.1f\n", preds[bitstart+j], al1[bitstart+j], al2[bitstart+j], stats[2+mark], neff);
fprintf(output4, "%s %.4e\n", preds[bitstart+j], stats[3+mark]);
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

bit++;
}	//end of while loop
printf("\n");

fclose(output);
fclose(output2);
fclose(output3);
fclose(output4);

printf("Main results saved in %s, with a summary version in %s and p-values in %s\n\n", filename2, filename3, filename4);

if(fastgwa==0){free(chrindex);free(prs);}
else{free(pedprs);}
free(data);free(datasqs);
free(order);free(Y);free(Z);
free(thetas);free(thetasds);;free(thetapvas);free(Yadj);free(Yadj2);
free(tindex);free(stats);free(spastatus);
free(YTdata);free(XTCX);
if(spatest==1)
{
free(usedpreds);
if(mode==131||fastgwa==0)
{
for(j=0;j<num_knots;j++){free(CGF0[j]);free(CGF1[j]);free(CGF2[j]);free(CGF3[j]);}
free(knots);free(bins);free(CGF0);free(CGF1);free(CGF2);free(CGF3);
}
}
if(mode==132){free(nullprobs);free(nullweights);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


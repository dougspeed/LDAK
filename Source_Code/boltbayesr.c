/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Perform ridge, bolt and bayesr
//dichot=0 - linear, dichot=1 - quasi-logistic (linear with weights nullweights)

///////////////////////////

//read in (all) the data - will be using either bed format or short SPEED format (dtype 1 or 4)
//will standardize so variance one
if(dtype==1)
{(void)read_bed_full(datafile, data_char, num_samples_use, keepsamps, data_length, keeppreds_use, num_samples, num_preds, maxthreads);}
else
{(void)read_speed_full(datafile, speedstarts, speedscales, data_char, NULL, num_samples_use, keepsamps, data_length, keeppreds_use, num_samples, num_preds, nonsnp, maxthreads);}

if(dtype==1)	//get look up for each predictor (have already computed means and variances)
{
printf("Constructing lookup table for each predictor\n\n");
#pragma omp parallel for private(j,k,readint) schedule (static)
for(j=0;j<data_length;j++)
{
for(k=0;k<256;k++)
{
readint=(k & 3);
if(readint!=3){bytelookup[j][k*4]=(readint-centres[j])*mults[j];}
else{bytelookup[j][k*4]=0;}
readint=((k >> 2) & 3);
if(readint!=3){bytelookup[j][k*4+1]=(readint-centres[j])*mults[j];}
else{bytelookup[j][k*4+1]=0;}
readint=((k >> 4) & 3);
if(readint!=3){bytelookup[j][k*4+2]=(readint-centres[j])*mults[j];}
else{bytelookup[j][k*4+2]=0;}
readint=((k >> 6) & 3);
if(readint!=3){bytelookup[j][k*4+3]=(readint-centres[j])*mults[j];}
else{bytelookup[j][k*4+3]=0;}
}
}
}

/*
//compute predictor-predictor covariances
printf("Computing predictor-predictor covariances\n\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Computing predictor-predictor covariances\n\n");
fclose(output);

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

//load standardized data
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

if(num_fixed>1&&adjpreds==1)	//adjust for covariates
{reg_covar_thetas(savethetas+bitstart*num_fixed, data, Z, num_samples_use, bitlength, num_fixed, NULL, 1);}

if(dichot==0||mpheno==-1)	//using non-weighted regression, so get XTX
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);
}
else	//using weighted regression, so get XTWX (must have one phenotype)
{
//put weighted version of data into data2
copy_matrix(num_samples_use, bitlength, data, data2, 1, nullweights);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data2, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);
}

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+j-bitstart];}
}	//end of bit loop
*/

////////

if(skipcv==0)	//make and test training models and maybe MCMC REML and/or top
{
total=num_small*num_resps_use;

//set effect sizes to zero, and fill residuals (will not use probs)
for(p=0;p<total;p++)
{
p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;
g=(p2-num_try-num_hers2-m2)/num_hers2;

if(p2<num_try)	//using real phenotypes with phenotypes of test samples set to zero
{
for(j=0;j<data_length;j++){effs[(size_t)p*data_length+j]=0;}
for(i=0;i<num_train;i++){residuals[(size_t)p*num_samples_use+keeptrain[i]]=Yadj[keeptrain[i]+m*num_samples_use];}
for(i=0;i<num_test;i++){residuals[(size_t)p*num_samples_use+keeptest[i]]=0;}
}
else	//MC REML terms
{
if(p2<num_try+num_hers2)	//using real phenotypes (but Yadj2 instead of Yadj)
{
for(j=0;j<data_length;j++){effs[(size_t)p*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)p*num_samples_use+i]=Yadj2[i+m*num_samples_use];}
}
else	//using random vectors
{
for(j=0;j<data_length;j++){effs[(size_t)p*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)p*num_samples_use+i]=gaussian[(size_t)g*num_samples_use+i];}
}
}
}

//screen and file print
if(mpheno!=-1){printf("Constructing %d PRS using training samples\n", num_try);}
else{printf("Constructing %d PRS for each phenotype using training samples\n", num_try);}
if(her==-9999&&revher==1){printf("Will also make %d MCMC REML models (using all samples)\n", num_hers2*(1+nmcmc)*num_resps_use);}
printf("\nIter\tNum_Con\tMax_Diff\tTarget\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output,"Constructing %d PRS using training samples\n", num_try);}
else{fprintf(output,"Constructing %d PRS for each phenotype using training samples\n", num_try);}
if(her==-9999&&revher==1){fprintf(output, "Will also make %d REML models (using all samples)\n", num_hers2*(1+nmcmc)*num_resps_use);}
fprintf(output,"Iter\tNum_Con\tMax_Diff\tTarget\n");
fclose(output);

//compute starting approx likelihoods (will use penalty zero)
#pragma omp parallel for private(p,p2,m,m2) schedule(static)
for(p=0;p<total;p++)
{
p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;

if(p2<num_try)	//using estimated heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], 0, dichot, nullweights+m*num_samples_use, -9999);}
else	//using fixed heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-tryhers[m2+m*num_hers2], 0, dichot, nullweights+m*num_samples_use, -9999);}
}

//iterate effect sizes
count=0;
while(1)
{
if(count==maxiter){printf("\nWarning, Variational Bayes failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)\n", maxiter, tol);break;}

//set pens to zero
for(p=0;p<total;p++){pens[p]=0;}

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

//load standardized data
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//will not adjust for covariates (can instead adjust residuals when necessary)

if(dichot==0||num_fixed==1)	//get t(X) residuals - can do all phenotypes at once
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitlength);
}
else	//get t(X) W residuals - first multiply columns of residuals by corresponding weights, then premultiply by t(X)
{
for(m=0;m<num_resps_use;m++)
{copy_matrix(num_samples_use, num_small, residuals+(size_t)m*num_small*num_samples_use, residuals2+(size_t)m*num_small*num_samples_use, 1, nullweights+m*num_samples_use);}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals2, &num_samples_use, &beta, YTdata, &bitlength);
}

#pragma omp parallel for private(p,p2,m,m2,j,j2,sum,value,value2,value3,value4,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;

for(j=bitstart;j<bitend;j++)
{
if(exps[j+m*data_length]>0)
{
//value4 is datasq
if(dichot==0){value4=datasqs[j];}
else{value4=datasqs[j+m*data_length];}

//get XjTresiduals
sum=YTdata[j-bitstart+p*bitlength]+effs[(size_t)p*data_length+j]*value4;

if(p2<num_try)	//model depends on mode, using estimated heritability
{
//value and value2 are how much to scale gaussian and laplace parameters, value3 is noise term
value=exps[j+m*data_length]*hers[m];
value2=pow(exps[j+m*data_length]*hers[m],-.5);
value3=1-hers[m];

if(mode==151)	//ridge
{postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, NULL);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[p2]*value, lambdas2[p2]*value, -9999, -9999, value4, value3, tryps[p2], tryp2s[p2], -9999, -9999, pens+p, 4, NULL);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[p2]*value, lambdas2[p2]*value, lambdas3[p2]*value, lambdas4[p2]*value, value4, value3, tryps[p2], tryp2s[p2], tryp3s[p2], tryp4s[p2], pens+p,  5+(pointmass==0), NULL);}
if(mode==154)	//elastic
{postmean=get_postmean(sum, lambdas[p2]*value2, lambdas2[p2]*value2, lambdas3[p2]*value, -9999, value4, value3, tryps[p2], tryp2s[p2], tryp3s[p2], -9999, pens+p, 7, NULL);}
}
else	//using ridge model and fixed heritability
{
//value and value3 are per-predictor and noise term
value=exps[j+m*data_length]*tryhers[m2+m*num_hers2];
value3=1-tryhers[m2+m*num_hers2];

postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, NULL);
}

if(pens[p]!=pens[p]){printf("Error p %d bit %d j %d; please tell doug\n", p+1, bit+1, j+1);exit(1);}
if(isinf(pens[p])){printf("Brror p %d bit %d j %d; please tell doug\n", p+1, bit+1, j+1);exit(1);}
}
else{postmean=0;}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitlength]=postmean-effs[(size_t)p*data_length+j];
effs[(size_t)p*data_length+j]=postmean;
if(dichot==0)
{
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)j*bitsize+j2-bitstart];}
}
else
{
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)(j+m*data_length)*bitsize+j2-bitstart];}
}
}	//end of j loop
}	//end of p loop

//update residuals - can do all phenotypes at once
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total, &bitlength, &alpha, data, &num_samples_use, changes, &bitlength, &beta, residuals, &num_samples_use);

if(num_fixed>1)	//regress out covariates
{
if(dichot==0)	//can do all phenotypes at once
{reg_covar_matrix(residuals, Z, num_samples_use, total, num_fixed);}
else	//do phenotypes separately
{
for(m=0;m<num_resps_use;m++)
{reg_covar_weighted(residuals+(size_t)m*num_small*num_samples_use, Z, num_samples_use, num_small, num_fixed, nullweights+m*num_samples_use);}
}
}
}	//end of bit loop

//save then update approx likelihoods
#pragma omp parallel for private(p,p2,m,m2) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];

p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;

if(p2<num_try)	//using estimated heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights+m*num_samples_use, -9999);}
else	//using fixed heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-tryhers[m2+m*num_hers2], pens[p], dichot, nullweights+m*num_samples_use, -9999);}
}

if(her==-9999&&revher==1)	//compute revised heritabilities
{
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_hers2;m2++)
{
p=num_try+m2+m*num_small;

if(dichot==0)	//numerator is t(R_Y) P_Y / t(R_Y) R_Y for real phenotypes
{
sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=residuals[(size_t)p*num_samples_use+i]*(Yadj[i+m*num_samples_use]-residuals[(size_t)p*num_samples_use+i]);
sumsq2+=pow(residuals[(size_t)p*num_samples_use+i],2);
}
value=sumsq/sumsq2;
}
else	//numerator is t(WR_Y) P_Y / t(WR_Y) R_Y for real phenotypes
{
sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=residuals[(size_t)p*num_samples_use+i]*nullweights[i+m*num_samples_use]*(Yadj[i+m*num_samples_use]-residuals[(size_t)p*num_samples_use+i]);
sumsq2+=pow(residuals[(size_t)p*num_samples_use+i],2)*nullweights[i+m*num_samples_use];
}
value=sumsq/sumsq2;
}

sum=0;
for(g=0;g<nmcmc;g++)
{
p=num_try+num_hers2+m2+g*num_hers2+m*num_small;

if(dichot==0)	//denominator is mean of t(r) P_r / t(r) R_r for random vectors
{
sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=gaussian[i+g*num_samples_use]*(gaussian[i+g*num_samples_use]-residuals[(size_t)p*num_samples_use+i]);
sumsq2+=gaussian[i+g*num_samples_use]*residuals[(size_t)p*num_samples_use+i];
}
sum+=sumsq/sumsq2;
}
else	//denominator is mean of t(r) P_r / t(r) R_r for random vectors
{
sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=gaussian[i+g*num_samples_use]*(gaussian[i+g*num_samples_use]-residuals[(size_t)p*num_samples_use+i]);
sumsq2+=gaussian[i+g*num_samples_use]*residuals[(size_t)p*num_samples_use+i];
}
sum+=sumsq/sumsq2;
}
}
value2=sum/nmcmc;

polates[m2]=value/value2;
}

//can now get the updated value
value=inter_hers(num_hers2, polates, tryhers+m*num_hers2, 1);

if(value!=-9999)	//worked
{
//if(mpheno!=-1){printf("The current MC REML heritability estimate is %.4f\n", value);}
hersold[m]=value;
}
else	//failed
{
if(mpheno!=-1){printf("Warning, MC REML has failed this iteration\n");}
else{printf("Warning, MC REML has failed this iteration for Phenotype %d\n", m+1);}
hersold[m]=hers[m];
}
}	//end of m loop
}	//end of revising heritability

if(count==0)	//can't stop after just one iteration, so just print an update
{
printf("1\t0\tNA\t%.4f\n", tol*num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "1\t0\tNA\t%.4f\n", tol*num_samples_use);
fclose(output);
}
else	//get number converged, largest difference, print update, then see if breaking
{
cflag=0;
diff=fabs(likes[0]-likesold[0]);
for(p=0;p<total;p++)
{
cflag+=(fabs(likes[p]-likesold[p])<tol*num_samples_use);
if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

printf("%d\t%d\t%.4f\t%.4f\n", count+1, cflag, diff, tol*num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol*num_samples_use);
fclose(output);

if(cflag==total){break;}
}

count++;
}	//end of while loop
printf("\n");

////////

//correct residuals for cv models
for(m=0;m<num_resps_use;m++)
{
for(p=0;p<num_try;p++)
{
for(i=0;i<num_test;i++){residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]]+=Yadj[keeptest[i]+m*num_samples_use];}
}
}

//measure accuracy, recording which has lowest MSE
if(mpheno!=-1){printf("Measuring accuracy of each model\n");}
else{printf("Finding the most accurate model for each phenotype\n");}

for(m=0;m<num_resps_use;m++)
{
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999)
{
if(mpheno!=-1){sprintf(filename3,"%s.accuracy", outfile);}
else{sprintf(filename3,"%s.pheno%d.accuracy", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

if(mode==151){fprintf(output3,"Model\tHeritability\tMean_Squared_Error\tCorrelation\n");}
if(mode==152){fprintf(output3,"Model\tHeritability\tp\tf2\tMean_Squared_Error\tCorrelation\n");}
if(mode==153){fprintf(output3,"Model\tHeritability\tp1\tp2\tp3\tp4\tMean_Squared_Error\tCorrelation\n");}
if(mode==154){fprintf(output3,"Model\tHeritability\tp\tf2\tMean_Squared_Error\tCorrelation\n");}
}

for(p=0;p<num_try;p++)
{
//compute relative mse (across non-missing individuals) and record lowest
sumsq=0;sumsq2=0;
for(i=0;i<num_test;i++)
{
if(resp[keeptest[i]+m*num_samples_use]!=missingvalue)
{
sumsq+=pow(residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]],2);
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2);
}
}
value=sumsq/sumsq2;

if(p==0){Mmses[m]=value;Mbests[m]=0;Mneffs[m]=respcounts[m]/value;}
if(value<Mmses[m]){Mmses[m]=value;Mbests[m]=p;Mneffs[m]=respcounts[m]/value;}

//compute correlation
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;indcount=0;
for(i=0;i<num_test;i++)
{
if(resp[keeptest[i]+m*num_samples_use]!=missingvalue)
{
sum+=Yadj[keeptest[i]+m*num_samples_use]-residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]];
sum2+=Yadj[keeptest[i]+m*num_samples_use];
sumsq+=pow(Yadj[keeptest[i]+m*num_samples_use]-residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]],2);
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2);
sumsq3+=(Yadj[keeptest[i]+m*num_samples_use]-residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]])*Yadj[keeptest[i]+m*num_samples_use];
indcount++;
}
}
mean=sum/indcount;mean2=sum2/indcount;
value3=(sumsq3-indcount*mean*mean2)*pow(sumsq-indcount*mean*mean,-.5)*pow(sumsq2-indcount*mean2*mean2,-.5);

if(mpheno!=-1)	//screen print all models - must have m=0
{
if(mode==151){printf("Model %d: heritability %.4f, mean squared error %.4f\n", p+1, hers[0], value);}
if(mode==152){printf("Model %d: heritability %.4f, p %.4f, f2 %.4f, mean squared error %.4f\n", p+1, hers[0], tryps[p], tryf2s[p], value);}
if(mode==153){printf("Model %d: heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f, mean squared error %.4f\n", p+1, hers[0], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], value);}
if(mode==154){printf("Model %d: heritability %.4f, p %.4f, f2 %.4f, mean squared error %.4f\n", p+1, hers[0], 1-tryp3s[p], tryf2s[p], value);}
}

if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999)	//save all models
{
if(mode==151){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], value, value3);}
if(mode==152){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], tryps[p], tryf2s[p], value, value3);}
if(mode==153){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], value, value3);}
if(mode==154){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], 1-tryp3s[p], tryf2s[p], value, value3);}
}
}	//end of p loop
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999){fclose(output3);}

if(mpheno==-1)	//screen print best model
{
if(mode==151){printf("Phenotype %d: best model %d, heritability %.4f, mean squared error %.4f\n", m+1, Mbests[m]+1, hers[m], Mmses[m]);}
if(mode==152){printf("Phenotype %d: best model %d, heritability %.4f, p %.4f, f2 %.4f, mean squared error %.4f\n", m+1, Mbests[m]+1, hers[m], tryps[Mbests[m]], tryf2s[Mbests[m]], Mmses[m]);}
if(mode==153){printf("Phenotype %d: best model %d, heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f, mean squared error %.4f\n", m+1, Mbests[m]+1, hers[m], tryps[Mbests[m]], tryp2s[Mbests[m]], tryp3s[Mbests[m]], tryp4s[Mbests[m]], Mmses[m]);}
if(mode==154){printf("Phenotype %d: best model %d, heritability %.4f, p %.4f, f2 %.4f, mean squared error %.4f\n", m+1, Mbests[m]+1, hers[m], 1-tryp3s[Mbests[m]], tryf2s[Mbests[m]], Mmses[m]);}
}
}	//end of m loop
printf("\n");

//assume we are using PRS
for(m=0;m<num_resps_use;m++){Mincs[m]=1;}

if(fprs==0)	//see if excluding the PRS
{
for(m=0;m<num_resps_use;m++)
{
if(Mmses[m]>0.995)
{
printf("Warning, the MSE for Phenotype %d is only %.4f%%, so will exclude the polygenic contribution\n\n", m+1, Mmses[m]);
Mincs[m]=0;Mneffs[m]=respcounts[m];
}
}
}

////////

if(mpheno==-1)	//compute multi-phenotype calibrations
{
for(m=0;m<num_resps_use;m++)	//estimating effects for phenotype m
{
//load best PRS into rows of PT and target phenotype into first column of Q
indcount=0;
for(i=0;i<num_test;i++)
{
if(resp[keeptest[i]+m*num_samples_use]!=missingvalue)
{
for(m2=0;m2<num_resps_use;m2++)
{
PT[m2+indcount*num_resps_use]=Yadj[keeptest[i]+m2*num_samples_use]-residuals[(size_t)(Mbests[m2]+m2*num_small)*num_samples_use+keeptest[i]];
}
Q[indcount]=Yadj[keeptest[i]+m*num_samples_use];
indcount++;
}
}

//centre rows of PT
for(m2=0;m2<num_resps_use;m2++)
{
sum=0;for(i=0;i<indcount;i++){sum+=PT[m2+i*num_resps_use];}
mean=sum/indcount;
for(i=0;i<indcount;i++){PT[m2+i*num_resps_use]-=mean;}
}

//centre first column of Q
sum=0;for(i=0;i<indcount;i++){sum+=Q[i];}
mean=sum/indcount;
for(i=0;i<indcount;i++){Q[i]-=mean;}

//compute PTP and PTQ (saved in mth column of Mcals)
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &indcount, &alpha, PT, &num_resps_use, PT, &num_resps_use, &beta, PTP, &num_resps_use);
dgemv_("N", &num_resps_use, &indcount, &alpha, PT, &num_resps_use, Q, &one, &beta, Mcals+m*num_resps_use, &one);

if(num_resps_use>1)	//add small amount of shrinkage
{
for(m2=0;m2<num_resps_use;m2++){PTP[m2+m2*num_resps_use]*=1.01;}
}

//get weights
(void)eigen_invert(PTP, num_resps_use, PTP2, 1, Mcals+m*num_resps_use, 1);

//make combined PRS and save in second column of Q
alpha=1.0;beta=0.0;
dgemv_("T", &num_resps_use, &indcount, &alpha, PT, &num_resps_use, Mcals+m*num_resps_use, &one, &beta, Q+num_test, &one);

//compute relative mse
sumsq=0;sumsq2=0;
for(i=0;i<indcount;i++)
{
sumsq+=pow(Q[i]-Q[i+num_test],2);
sumsq2+=pow(Q[i],2);
}
value=sumsq/sumsq2;
Mcombs[m]=respcounts[m]/value;
}
}	//end of mpheno==-1

if(her==-9999&&revher==1)	//update hers and print out / save
{
for(m=0;m<num_resps_use;m++)
{
if(hersold[m]<0.01){hersold[m]=0.01;}
if(hersold[m]>maxher){hersold[m]=maxher;}
hers[m]=hersold[m];

if(mpheno!=-1){printf("The revised estimate of heritability is %.4f\n", hers[m]);}
else{printf("Phenotype %d: revised heritability %.4f\n", m+1, hers[m]);}
}
printf("\n");

if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999)	//save
{
sprintf(filename3,"%s.hers", outfile);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++){fprintf(output3,"MCREML %d %.4f %.4f NA\n", m+1, powers[Mtops[m]], hers[m]);}
fclose(output3);
}
}
}	//end of making training models

////////

if(skipcv==0)	//make final models (will have best) and maybe loco and calibration models - note that have already set cors
{
total=num_full*num_resps_use;

//save best models from cv stage
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<data_length;j++){effs2[(size_t)m*data_length+j]=effs[(size_t)(Mbests[m]+m*num_small)*data_length+j];}
for(i=0;i<num_samples_use;i++){residuals2[(size_t)m*num_samples_use+i]=residuals[(size_t)(Mbests[m]+m*num_small)*num_samples_use+i];}
}

for(m=0;m<num_resps_use;m++)	//set effects and residuals (and set probs to zero)
{
//first model for each phenotype is final model from cv stage
for(j=0;j<data_length;j++){effs[(size_t)m*num_full*data_length+j]=effs2[(size_t)m*data_length+j];probs[(size_t)m*num_full*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)m*num_full*num_samples_use+i]=residuals2[(size_t)m*num_samples_use+i];}

if(loco==1)
{
//also put the final model into the next num_chr positions
for(p=0;p<num_chr;p++)
{
p2=1+p+m*num_full;
for(j=0;j<data_length;j++){effs[(size_t)p2*data_length+j]=effs2[(size_t)m*data_length+j];probs[(size_t)p2*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)p2*num_samples_use+i]=residuals2[(size_t)m*num_samples_use+i];}
}

if(usecal==1)	//add in calibration models
{
for(p=0;p<ncal;p++)
{
p2=1+num_chr+p+m*num_full;
for(j=0;j<data_length;j++){effs[(size_t)p2*data_length+j]=0;probs[(size_t)p2*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)p2*num_samples_use+i]=cdata[i+(p+m*ncal)*num_samples_use];}
}
}
}
}

////////

//screen and file print
if(mpheno!=-1)	//must have m=0
{
if(mode==151){printf("Constructing final PRS (heritability %.4f) using all samples\n", hers[0]);}
if(mode==152){printf("Constructing final PRS (heritability %.4f, p %.4f, f2 %.4f) using all samples\n", hers[0], tryps[Mbests[0]], tryf2s[Mbests[0]]);}
if(mode==153){printf("Constructing final PRS (heritability %.4f, p1 %.4f, p2 %.4f, p3 %.4f, p4 %.4f) using all samples\n", hers[0], tryps[Mbests[0]], tryp2s[Mbests[0]], tryp3s[Mbests[0]], tryp4s[Mbests[0]]);}
if(mode==154){printf("Constructing final PRS (heritability %.4f, p %.4f, f2 %.4f) using all samples\n", hers[0], 1-tryp3s[Mbests[0]], tryf2s[Mbests[0]]);}
if(loco==1)
{
printf("Will also make %d LOCO PRS", num_chr);
if(usecal==1){printf(", and %d calibration models", ncal);}
printf("\n");
}
}
else
{
printf("Constructing final PRS for each phenotype\n");
if(loco==1)
{
printf("Will also make %d LOCO PRS", num_chr);
if(usecal==1){printf(", and %d calibration models", ncal);}
printf("\n");
}
}
printf("\nIter\tNum_Con\tMax_Diff\tTarget\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output,"Constructing final PRS using all samples\n");}
else{fprintf(output,"Constructing final PRS for %d phenotypes using all samples\n", num_resps_use);}
if(loco==1){fprintf(output, "Will also make %d LOCO PRS\n", num_chr);}
fprintf(output,"Iter\tNum_Con\tMax_Diff\tTarget\n");
fclose(output);

//iterate effect sizes
count=0;
while(1)
{
if(count==maxiter){printf("\nWarning, Variational Bayes failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)\n", maxiter, tol);break;}

//set pens to zero
for(p=0;p<total;p++){pens[p]=0;}

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

//load standardized data
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//will not adjust for covariates (can instead adjust residuals when necessary)

if(dichot==0||num_fixed==1)	//get t(X) residuals - can do all phenotypes at once
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitlength);
}
else	//get t(X) W residuals - first multiply columns of residuals by corresponding weights, then premultiply by t(X)
{
for(m=0;m<num_resps_use;m++)
{copy_matrix(num_samples_use, num_full, residuals+(size_t)m*num_full*num_samples_use, residuals2+(size_t)m*num_full*num_samples_use, 1, nullweights+m*num_samples_use);}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals2, &num_samples_use, &beta, YTdata, &bitlength);
}

#pragma omp parallel for private(p,p2,m,j,j2,sum,value,value2,value3,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
p2=p%num_full;
m=(p-p2)/num_full;

for(j=bitstart;j<bitend;j++)
{
if(exps[j+m*data_length]>0&&chr[j]!=chrindex[p2])
{
//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j+m*data_length]*hers[m]/chrprops[p];
value2=pow(exps[j+m*data_length]*hers[m]/chrprops[p],-.5);

//value3 is noise term
if(dichot==0||dichot==1){value3=1-hers[m];}
else{value3=1;}

//value4 is datasq
if(dichot==0){value4=datasqs[j];}
else{value4=datasqs[j+m*data_length];}

//get XjTresiduals (or XjTWresiduals)
sum=YTdata[j-bitstart+p*bitlength]+effs[(size_t)p*data_length+j]*value4;

if(mode==151)	//ridge
{postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, probs+(size_t)p*data_length+j);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[Mbests[m]]*value, lambdas2[Mbests[m]]*value, -9999, -9999, value4, value3, tryps[Mbests[m]], tryp2s[Mbests[m]], -9999, -9999, pens+p, 4, probs+(size_t)p*data_length+j);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[Mbests[m]]*value, lambdas2[Mbests[m]]*value, lambdas3[Mbests[m]]*value, lambdas4[Mbests[m]]*value, value4, value3, tryps[Mbests[m]], tryp2s[Mbests[m]], tryp3s[Mbests[m]], tryp4s[Mbests[m]], pens+p,  5+(pointmass==0), probs+(size_t)p*data_length+j);}
if(mode==154)	//elastic
{postmean=get_postmean(sum, lambdas[Mbests[m]]*value2, lambdas2[Mbests[m]]*value2, lambdas3[Mbests[m]]*value, -9999, value4, value3, tryps[Mbests[m]], tryp2s[Mbests[m]], tryp3s[Mbests[m]], -9999, pens+p, 7, probs+(size_t)p*data_length+j);}

if(pens[p]!=pens[p]){printf("Error p %d bit %d j %d; please tell doug\n", p+1, bit+1, j+1);exit(1);}
if(isinf(pens[p])){printf("Brror p %d bit %d j %d\n; please tell doug", p+1, bit+1, j+1);exit(1);}
}
else{postmean=0;}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitlength]=postmean-effs[(size_t)p*data_length+j];
effs[(size_t)p*data_length+j]=postmean;
if(dichot==0)
{
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)j*bitsize+j2-bitstart];}
}
else
{
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)(j+m*data_length)*bitsize+j2-bitstart];}
}
}	//end of j loop
}	//end of p loop

//update residuals - can do all phenotypes at once
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total, &bitlength, &alpha, data, &num_samples_use, changes, &bitlength, &beta, residuals, &num_samples_use);

if(num_fixed>1)	//regress out covariates
{
if(dichot==0)	//can do all phenotypes at once
{reg_covar_matrix(residuals, Z, num_samples_use, total, num_fixed);}
else	//do phenotypes separately
{
for(m=0;m<num_resps_use;m++)
{reg_covar_weighted(residuals+(size_t)m*num_full*num_samples_use, Z, num_samples_use, num_full, num_fixed, nullweights+m*num_samples_use);}
}
}
}	//end of bit loop

//save then update approx likelihoods
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];

p2=p%num_full;
m=(p-p2)/num_full;

likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights+m*num_samples_use, -9999);
}

if(count==0)	//can't stop after just one iteration, so just print an update
{
printf("1\t0\tNA\t%.4f\n", tol*num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "1\t0\tNA\t%.4f\n", tol*num_samples_use);
fclose(output);
}
else	//get number converged, largest difference, print update, then see if breaking
{
cflag=0;
diff=fabs(likes[0]-likesold[0]);
for(p=0;p<total;p++)
{
cflag+=(fabs(likes[p]-likesold[p])<tol*num_samples_use);
if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

printf("%d\t%d\t%.4f\t%.4f\n", count+1, cflag, diff, tol*num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol*num_samples_use);
fclose(output);

if(cflag==total){break;}
}

count++;
}	//end of while loop
printf("\n");

////////

if(usecal==1)	//get "ridge lambda"
{
for(m=0;m<num_resps_use;m++)
{
if(Mincs[m]==0){cgammas[m]=1;csds[m]=0;}
else
{
//ridge scaling is GTG/(GT invV G) x var(invV Y) / t(Y) invV Y x n = t(G) G / t(G) R_G x nvar(R_Y) / t(Y) R_Y
//for dichot=1, it is GTWG/(GT invV G) x var(1/sqrt(W) invV Y) / t(Y) invV Y x n = t(G) WG / t(G) WR_G x nvar(sqrt(W) R_Y) / t(Y) WR_Y
//for dichot=2, it is simply GTWG/(GT invV G) = t(G) WG / t(G) WR_G

sum=0;sum2=0;count=0;
for(j=0;j<ncal;j++)
{
if(cmults[j]!=-9999)
{
//p and p2 index residuals corresponding to R_Y and R_G
p=1;while(chr[cindex[j]]>chrindex[p]){p++;}
p2=1+num_chr+j;

if(dichot==0)	//sumsq is t(G) G; sumsq2 is t(G) R_G; sumsq3 is t(R_Y) R_Y; sumsq4 is t(Y) R_Y - (R_Y has mean zero)
{
sumsq=0;sumsq2=0;sumsq3=0;sumsq4=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=pow(cdata[i+(j+m*ncal)*num_samples_use],2);
sumsq2+=cdata[i+(j+m*ncal)*num_samples_use]*residuals[(size_t)(p2+m*num_full)*num_samples_use+i];
sumsq3+=pow(residuals[(size_t)(p+m*num_full)*num_samples_use+i],2);
sumsq4+=Yadj[i+m*num_samples_use]*residuals[(size_t)(p+m*num_full)*num_samples_use+i];
}
value=sumsq/sumsq2*sumsq3/sumsq4;
}
if(dichot==1)	//sumsq is t(G) WG; sumsq2 is t(G) WR_G; sum3 + sumsq3 give var(sqrt(W) R_Y); sumsq4 is t(Y) WR_Y
{
sum3=0;sumsq=0;sumsq2=0;sumsq3=0;sumsq4=0;
for(i=0;i<num_samples_use;i++)
{
sum3+=residuals[(size_t)(p+m*num_full)*num_samples_use+i]*pow(nullweights[i+m*num_samples_use],.5);
sumsq+=pow(cdata[i+(j+m*ncal)*num_samples_use],2)*nullweights[i+m*num_samples_use];
sumsq2+=cdata[i+(j+m*ncal)*num_samples_use]*residuals[(size_t)(p2+m*num_full)*num_samples_use+i]*nullweights[i+m*num_samples_use];
sumsq3+=pow(residuals[(size_t)(p+m*num_full)*num_samples_use+i],2)*nullweights[i+m*num_samples_use];
sumsq4+=Yadj[i+m*num_samples_use]*residuals[(size_t)(p+m*num_full)*num_samples_use+i]*nullweights[i+m*num_samples_use];
}
value=sumsq/sumsq2*(sumsq3-sum3/num_samples_use*sum3)/sumsq4;
}
if(dichot==2)	//sumsq is t(G) WG; sumsq2 is t(G) WR_G;
{
sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=pow(cdata[i+(j+m*ncal)*num_samples_use],2)*nullweights[i+m*num_samples_use];
sumsq2+=cdata[i+(j+m*ncal)*num_samples_use]*residuals[(size_t)(p2+m*num_full)*num_samples_use+i]*nullweights[i+m*num_samples_use];
}
value=sumsq/sumsq2;
}

sum+=value;
sum2+=pow(value,2);
count++;
}
}
mean=sum/count;
var=sum2/count-pow(mean,2);

cgammas[m]=mean;
csds[m]=pow(var/ncal,.5);

if(mpheno!=-1){printf("Test statistic scaling factor is %.4f (SE %.4f)\n", cgammas[m], csds[m]);}
else{printf("Phenotype %d: test statistic scaling factor is %.4f (SE %.4f)\n", m+1, cgammas[m], csds[m]);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output, "Test statistic scaling factor is %.4f (SE %.4f)\n", cgammas[m], csds[m]);}
else{fprintf(output, "Phenotype %d: test statistic scaling factor is %.4f (SE %.4f)\n", m+1, cgammas[m], csds[m]);}
fclose(output);
}
}	//end of m loop
printf("\n");
}

////////

if(loco==1)	//save root file
{
sprintf(filename3,"%s.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Analysis KVIK\n");
if(dichot==0){fprintf(output3,"Regression Linear\n");}
else{fprintf(output3,"Regression Logistic\n");}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Phenotypes %s\n", respfile);
if(strcmp(covarfile,"blank")!=0){fprintf(output3,"Covariates %s\n", covarfile);}
else{fprintf(output3,"Covariates none\n");}
if(strcmp(topfile,"blank")!=0){fprintf(output3,"Top_Predictors %s\n", topfile);}
else{fprintf(output3,"Top_Predictors none\n");}
if(strcmp(factorfile,"blank")!=0){fprintf(output3,"Factors %s\n", factorfile);}
else{fprintf(output3,"Factors none\n");}
if(mpheno!=-1){fprintf(output3,"Phenotype_Number %d\n", mpheno);}
else{fprintf(output3,"Num_Phenotypes %d\n", num_resps_use);}
fprintf(output3,"Num_Samples_Used %d\n", num_samples_use);
fprintf(output3,"Num_Predictors_Used %d\n", data_length);
fprintf(output3,"Num_Chromosomes %d\n", num_chr);
fclose(output3);
}

for(m=0;m<num_resps_use;m++)	//save - remember that we scaled Yadj to have variance one
{
//save effects for final models
if(mpheno!=-1){sprintf(filename3,"%s.effects", outfile);}
else{sprintf(filename3,"%s.pheno%d.effects", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor A1 A2 Centre Effect\n");

for(j=0;j<num_tops;j++)
{fprintf(output3, "%s %c %c %.6f %.4e\n", tpreds[j], tal1[j], tal2[j], tcentres[j], thetas[num_covars+j+m*num_fixed]);}
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){fprintf(output3, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[(size_t)m*num_full*data_length+j]*mults[j]/Mscales[m]);}
}
fclose(output3);

if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999)	//save some more things
{
//save probabilities for final models
if(mpheno!=-1){sprintf(filename3,"%s.probs", outfile);}
else{sprintf(filename3,"%s.pheno%d.probs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Probability\n");

for(j=0;j<num_tops;j++){fprintf(output3, "%s 1\n", tpreds[j]);}
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){fprintf(output3, "%s %.4f\n", preds[j], probs[(size_t)m*num_full*data_length+j]);}
}
fclose(output3);

//save final PRS (for non-missing individuals)
if(mpheno!=-1){sprintf(filename3,"%s.prs", outfile);}
else{sprintf(filename3,"%s.pheno%d.prs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"FID IID Adjusted_Phenotype PRS\n");

for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]!=missingvalue)
{fprintf(output3, "%s %s %.4f %.4f\n", ids1[i], ids2[i], Yadj[i+m*num_samples_use]/Mscales[m], (Yadj[i+m*num_samples_use]-residuals[(size_t)m*num_full*num_samples_use+i])/Mscales[m]);}
}
fclose(output3);
}

if(loco==1)	//save LOCO PRS and details
{
if(mpheno!=-1){sprintf(filename3,"%s.loco.prs", outfile);}
else{sprintf(filename3,"%s.pheno%d.loco.prs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

fprintf(output3,"FID IID");
for(p=1;p<1+num_chr;p++){fprintf(output3," Chr%d", chrindex[p]);}
fprintf(output3,"\n");

for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]!=missingvalue)
{
fprintf(output3, "%s %s", ids1[i], ids2[i]);
for(p=1;p<1+num_chr;p++)
{
if(Mincs[m]==1){fprintf(output3," %.4f", (Yadj[i+m*num_samples_use]-residuals[(size_t)(p+m*num_full)*num_samples_use+i])/Mscales[m]);}
else{fprintf(output3, " 0");}
}
fprintf(output3,"\n");
}
}
fclose(output3);

if(mpheno!=-1){sprintf(filename3,"%s.loco.details", outfile);}
else{sprintf(filename3,"%s.pheno%d.loco.details", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Actual_Sample_Size %d\n", respcounts[m]);
fprintf(output3,"Approx_Effective_Sample_Size %.1f\n", Mneffs[m]);
if(usecal==1)
{
fprintf(output3,"Scaling_Estimate %.4f\n", cgammas[m]);
fprintf(output3,"Scaling_SE %.4f\n", csds[m]);
}
else
{
fprintf(output3,"Scaling_Estimate 1\n");
fprintf(output3,"Scaling_SE 0\n");
}
fprintf(output3,"Power %.4f\n", powers[Mtops[m]]);
fprintf(output3,"Heritability %.4f\n", hers[m]);
fclose(output3);
}

if(mpheno==-1&&kvikstep==-9999&&gctastep==-9999&&faststep==-9999)	//also save multi-phen versions
{
//save calibrations
sprintf(filename3,"%s.pheno%d.multi.weights", outfile, m+1);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Phenotype Effect\n");
for(m2=0;m2<num_resps_use;m2++){fprintf(output3,"%d %.4f\n", m2+1, Mcals[m2+m*num_resps_use]*Mscales[m2]/Mscales[m]);}
fclose(output3);

if(loco==1)	//save LOCO PRS and root
{
sprintf(filename3,"%s.pheno%d.multi.loco.prs", outfile, m+1);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"FID IID");
for(p=1;p<1+num_chr;p++){fprintf(output3," Chr%d", chrindex[p]);}
fprintf(output3,"\n");

for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]!=missingvalue)
{
fprintf(output3, "%s %s", ids1[i], ids2[i]);

for(p=1;p<1+num_chr;p++)
{
sum=0;
for(m2=0;m2<num_resps_use;m2++){sum+=Mcals[m2+m*num_resps_use]*(Yadj[i+m2*num_samples_use]-residuals[(size_t)(p+m2*num_full)*num_samples_use+i]);}
fprintf(output3," %.4f", sum/Mscales[m]);
}
fprintf(output3,"\n");
}
}
fclose(output3);

if(mpheno!=-1){sprintf(filename3,"%s.multi.loco.details", outfile);}
else{sprintf(filename3,"%s.pheno%d.multi.loco.details", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Phenotypes %s\n", respfile);
fprintf(output3,"Analysis KVIK_Multi\n");
if(dichot==0){fprintf(output3,"Regression Linear\n");}
else{fprintf(output3,"Regression Logistic\n");}
if(mpheno!=-1){fprintf(output3,"Phenotype_Number %d\n", mpheno+1);}
else{fprintf(output3,"Num_Phenotypes ALL\n");}
if(strcmp(covarfile,"blank")!=0){fprintf(output3,"Covariates %s\n", covarfile);}
else{fprintf(output3,"Covariates none\n");}
if(strcmp(topfile,"blank")!=0){fprintf(output3,"Top_Predictors %s\n", topfile);}
else{fprintf(output3,"Top_Predictors none\n");}
if(strcmp(factorfile,"blank")!=0){fprintf(output3,"Factors %s\n", factorfile);}
else{fprintf(output3,"Factors none\n");}
fprintf(output3,"Num_Samples_Used %d\n", respcounts[m]);
fprintf(output3,"Num_Predictors_Used %d\n", data_length);
fprintf(output3,"Num_Chromosomes %d\n", num_chr);
fprintf(output3,"Approx_Effective_Sample_Size %.1f\n", Mcombs[m]);
if(usecal==1)
{
fprintf(output3,"Scaling_Estimate %.4f\n", cgammas[m]);
fprintf(output3,"Scaling_SE %.4f\n", csds[m]);
}
else
{
fprintf(output3,"Scaling_Estimate 1\n");
fprintf(output3,"Scaling_SE 0\n");
}
fprintf(output3,"Power %.4f\n", powers[Mtops[m]]);
fprintf(output3,"Heritability %.4f\n", hers[m]);
fclose(output3);
}
}	//end of saving multi-phen versions
}	//end of m loop

if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999)
{
if(mpheno!=-1){printf("Best-fitting model saved in %s.effects, with posterior probabilities in %s.probs and in-sample PRS in %s.prs\n\n", outfile, outfile, outfile);}
else{printf("Best-fitting models saved in %s.phenoX.effects, with posterior probabilities in %s.phenoX.probs, and in-sample PRS in %s.phenoX.prs, where X is the phenotype number\n\n", outfile, outfile, outfile);}
}

if(loco==1)
{
if(mpheno!=-1){printf("LOCO results saved in %s.loco.prs and %s.loco.details\n\n", outfile, outfile);}
else{printf("LOCO results saved in %s.phenoX.loco.prs and %s.phenoX.loco.prs, where X is the phenotype number\n\n", outfile, outfile);}
}
}	//end of skipcv=0

////////

if(skipcv==1)	//make all models using all samples - will not be performing loco
{
total=num_try*num_resps_use;

//set effect sizes and probs to zero, and fill residuals
for(p=0;p<total;p++)
{
p2=p%num_try;
m=(p-p2)/num_try;

for(j=0;j<data_length;j++){effs[(size_t)p*data_length+j]=0;probs[(size_t)p*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)p*num_samples_use+i]=Yadj[i+m*num_samples_use];}
}

//screen and file print
if(mpheno!=-1){printf("Constructing %d PRS using all samples\n\n", num_try);}
else{printf("Constructing %d PRS for each phenotype using all samples\n\n", num_try);}
printf("Iter\tNum_Con\tMax_Diff\tTarget\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output,"Constructing %d PRS using all samples\n", total);}
else{fprintf(output,"Constructing %d PRS for each phenotype using all samples\n", num_try);}
fprintf(output,"Iter\tNum_Con\tMax_Diff\tTarget\n");
fclose(output);

//compute starting approx likelihoods (will use penalty zero)
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
p2=p%num_try;
m=(p-p2)/num_try;

likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], 0, dichot, nullweights+m*num_samples_use, -9999);
}

//iterate effect sizes
count=0;
while(1)
{
if(count==maxiter){printf("\nWarning, Variational Bayes failed to converge; consider using \"--max-iter\" and/or \"--tolerance\" to increase the iteration limit and tolerance (currently %d and %.4e)\n", maxiter, tol);break;}

//set pens to zero
for(p=0;p<total;p++){pens[p]=0;}

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

//load standardized data
load_data_chunk(data_char, dtype, data, num_samples_use, bitstart, bitend, NULL, bytelookup, speedstarts, speedscales, centres);

//will not adjust for covariates (can instead adjust residuals when necessary)

if(dichot==0||num_fixed==1)	//get t(X) residuals - can do all phenotypes at once
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitlength);
}
else	//get t(X) W residuals - first multiply columns of residuals by corresponding weights, then premultiply by t(X)
{
for(m=0;m<num_resps_use;m++)
{copy_matrix(num_samples_use, num_try, residuals+(size_t)m*num_try*num_samples_use, residuals2+(size_t)m*num_try*num_samples_use, 1, nullweights+m*num_samples_use);}

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals2, &num_samples_use, &beta, YTdata, &bitlength);
}

#pragma omp parallel for private(p,p2,m,j,j2,sum,value,value2,value3,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
p2=p%num_try;
m=(p-p2)/num_try;

for(j=bitstart;j<bitend;j++)
{
if(exps[j+m*data_length]>0)
{
//value4 is datasq
if(dichot==0){value4=datasqs[j];}
else{value4=datasqs[j+m*data_length];}

//get XjTresiduals (or XjTWresiduals)
sum=YTdata[j-bitstart+p*bitlength]+effs[(size_t)p*data_length+j]*value4;

//value and value2 are how much to scale gaussian and laplace parameters, value3 is noise term
value=exps[j+m*data_length]*hers[m];
value2=pow(exps[j+m*data_length]*hers[m],-.5);
value3=1-hers[m];

if(mode==151)	//ridge
{postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, probs+(size_t)p*data_length+j);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[p2]*value, lambdas2[p2]*value, -9999, -9999, value4, value3, tryps[p2], tryp2s[p2], -9999, -9999, pens+p, 4, probs+(size_t)p*data_length+j);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[p2]*value, lambdas2[p2]*value, lambdas3[p2]*value, lambdas4[p2]*value, value4, value3, tryps[p2], tryp2s[p2], tryp3s[p2], tryp4s[p2], pens+p,  5+(pointmass==0), probs+(size_t)p*data_length+j);}
if(mode==154)	//elastic
{postmean=get_postmean(sum, lambdas[p2]*value2, lambdas2[p2]*value2, lambdas3[p2]*value, -9999, value4, value3, tryps[p2], tryp2s[p2], tryp3s[p2], -9999, pens+p, 7, probs+(size_t)p*data_length+j);}
}
else{postmean=0;}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitlength]=postmean-effs[(size_t)p*data_length+j];
effs[(size_t)p*data_length+j]=postmean;
if(dichot==0)
{
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)j*bitsize+j2-bitstart];}
}
else
{
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)(j+m*data_length)*bitsize+j2-bitstart];}
}
}	//end of j loop
}	//end of p loop

//update residuals - can do all phenotypes at once
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total, &bitlength, &alpha, data, &num_samples_use, changes, &bitlength, &beta, residuals, &num_samples_use);

if(num_fixed>1)	//regress out covariates
{
if(dichot==0)	//can do all phenotypes at once
{reg_covar_matrix(residuals, Z, num_samples_use, total, num_fixed);}
else	//do phenotypes separately
{
for(m=0;m<num_resps_use;m++)
{reg_covar_weighted(residuals+(size_t)m*num_try*num_samples_use, Z, num_samples_use, num_try, num_fixed, nullweights+m*num_samples_use);}
}
}
}	//end of bit loop

//save then update approx likelihoods
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];

p2=p%num_try;
m=(p-p2)/num_try;

likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights, -9999);
}

if(count==0)	//can't stop after just one iteration, so just print an update
{
printf("1\t0\tNA\t%.4f\n", tol*num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "1\t0\tNA\t%.4f\n", tol*num_samples_use);
fclose(output);
}
else	//get number converged, largest difference, print update, then see if breaking
{
cflag=0;
diff=fabs(likes[0]-likesold[0]);
for(p=0;p<total;p++)
{
cflag+=(fabs(likes[p]-likesold[p])<tol*num_samples_use);
if(fabs(likes[p]-likesold[p])>diff){diff=fabs(likes[p]-likesold[p]);}
}

printf("%d\t%d\t%.4f\t%.4f\n", count+1, cflag, diff, tol*num_samples_use);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"%d\t%d\t%.4f\t%.4f\n", count, cflag, diff, tol*num_samples_use);
fclose(output);

if(cflag==total){break;}
}

count++;
}	//end of while loop
printf("\n");

////////

for(m=0;m<num_resps_use;m++)	//save - remember that we may have scaled Yadj to have variance one
{
//save effects and probabilities
if(mpheno!=-1){sprintf(filename3,"%s.effects", outfile);}
else{sprintf(filename3,"%s.pheno%d.effects", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor A1 A2 Centre");
for(p=0;p<num_try;p++){fprintf(output3," Effect%d", p+1);}
fprintf(output3,"\n");

for(j=0;j<num_tops;j++)
{
fprintf(output3, "%s %c %c %.6f", tpreds[j], tal1[j], tal2[j], tcentres[j]);
for(p=0;p<num_try;p++){fprintf(output3," %.4e", thetas[num_covars+j+m*num_fixed]);}
fprintf(output3,"\n");
}
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output3, "%s %c %c %.6f ", preds[j], al1[j], al2[j], centres[j]);
for(p=0;p<num_try;p++){fprintf(output3," %.4e", effs[(size_t)(p+m*num_try)*data_length+j]*mults[j]/Mscales[m]);}
fprintf(output3,"\n");
}
}
fclose(output3);

if(mpheno!=-1){sprintf(filename3,"%s.probs", outfile);}
else{sprintf(filename3,"%s.pheno%d.probs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor");
for(p=0;p<num_try;p++){fprintf(output3," Probability%d", p+1);}
fprintf(output3,"\n");

for(j=0;j<num_tops;j++)
{
fprintf(output3, "%s", tpreds[j]);
for(p=0;p<num_try;p++){fprintf(output3," 1");}
fprintf(output3,"\n");
}
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output3, "%s", preds[j]);
for(p=0;p<num_try;p++){fprintf(output3," %.4f", probs[(size_t)(p+m*num_try)*data_length+j]);}
fprintf(output3,"\n");
}
}
fclose(output3);

//print out final PRS (for non-missing individuals)
if(mpheno!=-1){sprintf(filename3,"%s.prs", outfile);}
else{sprintf(filename3,"%s.pheno%d.prs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"FID IID Adjusted_Phenotype");
for(p=0;p<num_try;p++){fprintf(output3," PRS%d", p+1);}
fprintf(output3,"\n");

for(i=0;i<num_samples_use;i++)
{
if(resp[i+m*num_samples_use]!=missingvalue)
{
fprintf(output3, "%s %s %.4f ", ids1[i], ids2[i], Yadj[i+m*num_samples_use]/Mscales[m]);
for(p=0;p<num_try;p++){fprintf(output3,"  %.4f", (Yadj[i+m*num_samples_use]-residuals[(size_t)(p+m*num_try)*num_samples_use+i])/Mscales[m]);}
fprintf(output3,"\n");
}
}
fclose(output3);
}	//end of m loop

if(mpheno!=-1){printf("Models saved in %s.effects, with posterior probabilities in %s.probs, and in-sample PRS in %s.prs\n\n", outfile, outfile, outfile);}
else{printf("Models saved in %s.phenoX.effects, with posterior probabilities in %s.phenoX.probs, and in-sample PRS in %s.phenoX.prs, where X is the phenotype number\n\n", outfile, outfile, outfile);}
}	//end of skipcv=1

///////////////////////////


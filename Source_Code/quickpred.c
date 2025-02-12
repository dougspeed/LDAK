/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//A quick version of ridge, bolt, bayesr and elastic (data must be in binary format)
//dichot=0 - linear, dichot=1 - quasi-logistic (linear with weights nullweights)

///////////////////////////

if(skipcv==0)	//make and test training models and maybe MCMC REML and/or comparison models
{
total=num_small*num_resps_use;

/*
//save likelihood difference per bit
sprintf(filename4,"%s.likes.train", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Scan Bit Iterations Chr");
for(p=0;p<total;p++){fprintf(output4, " Model%d", p+1);}
fprintf(output4,"\n");
fclose(output4);
*/

//set effect sizes to zero, and fill residuals (will not use probs)
for(p=0;p<total;p++)
{
p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;
g=(p2-num_try-num_hers2-m2)/num_hers2;
k=p2-num_try-revher*num_hers2*(1+nmcmc);

if(p2<num_try)	//cv models - using real phenotypes with phenotypes of test samples set to zero
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
if(revher==1){printf("Will also make %d MCMC REML models (using all samples)\n", revher*num_hers2*(1+nmcmc)*num_resps_use);}
printf("\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output,"Constructing %d PRS using training samples\n", num_try);}
else{fprintf(output,"Constructing %d PRS for each phenotype using training samples\n", num_try);}
if(revher==1){fprintf(output, "Will also make %d REML models (using all samples)\n", revher*num_hers2*(1+nmcmc)*num_resps_use);}
fclose(output);

//set bit penalties to zero (not sure if correct, but we need a value)
for(p=0;p<total;p++)
{
for(bit=0;bit<bittotal;bit++){bitpens[bit+p*bittotal]=0;}
}

//bitrun records how many models running for each phenotypes
for(m=0;m<num_resps_use;m++)
{
for(bit=0;bit<bittotal;bit++){bitrun[bit+m*bittotal]=num_try;}
}

//bitdo records how many phenotypes running for each bit
for(bit=0;bit<bittotal;bit++){bitdo[bit]=num_resps_use;}

//total2 records how many bits to visit
total2=bittotal;

//set MSEs to one
for(m=0;m<num_resps_use;m++){Mmses[m]=1;}

//save current estimates of heritability
for(m=0;m<num_resps_use;m++){hersold[m]=hers[m];}

count4=0;	//number of scans
while(total2>0)
{
if(count4==nscan){printf("Warning, Variational Bayes did not converge after %d scans (this is not normally a problem)\n\n", nscan);break;}

//ready for bit loop
count2=0;	//number of bits tested
count3=0;	//total number of iterations
for(bit=0;bit<bittotal;bit++)
{
if(bitdo[bit]>0)	//will consider bit
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

bitdet2[bit]=chr[bitstart];

if(count2%200==0)
{
printf("Scan %d: estimating training effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Scan %d: estimating training effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);
fclose(output);
}

//compute current (partial) approx likelihood and save minus this in bitdiffs
#pragma omp parallel for private(p,p2,m,m2) schedule(static)
for(p=0;p<total;p++)
{
p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
if(p2<num_try)	//using estimated heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], bitpens[bit+p*bittotal], dichot, nullweights+m*num_samples_use, -9999);}
else	//using fixed heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-tryhers[m2+m*num_hers2], bitpens[bit+p*bittotal], dichot, nullweights+m*num_samples_use, -9999);}
}

bitdiffs[bit+p*bittotal]=-likes[p];
}

//read data, standardize and set missing to zero - already have statistics
if(dtype==1)	//fast way
{(void)read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, NULL, NULL, NULL, 0, maxthreads);}
else	//slow way
{
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs,bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, NULL, NULL, NULL, num_samples_use, bitlength, missingvalue, -9999, -9999, -9999, NULL, 2);
}

//will not adjust for covariates (can instead adjust residuals when necessary)

//iterate effect sizes
count=0;
while(count<maxiter)
{
count++;

//reset pens to zero for each iteration
for(p=0;p<total;p++){pens[p]=0;}

for(m=0;m<num_resps_use;m++)
{
if(bitrun[bit+m*bittotal]>0)
{
if(dichot==0||num_fixed==1)	//get t(X) residuals (if dichot=1 and no covariates, weights will be one)
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_small, &num_samples_use, &alpha, data, &num_samples_use, residuals+(size_t)m*num_small*num_samples_use, &num_samples_use, &beta, YTdata+m*num_small*bitlength, &bitlength);
}
else	//get t(X) W residuals - first multiply columns of residuals by corresponding weights, then premultiply by t(X)
{
copy_matrix(num_samples_use, num_small, residuals+(size_t)m*num_small*num_samples_use, residuals2+(size_t)m*num_small*num_samples_use, 1, nullweights+m*num_samples_use);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_small, &num_samples_use, &alpha, data, &num_samples_use, residuals2+(size_t)m*num_small*num_samples_use, &num_samples_use, &beta, YTdata+m*num_small*bitlength, &bitlength);
}
}
}

#pragma omp parallel for private(p,p2,m,m2,j,j2,sum,value,value2,value3,value4,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j+m*data_length]>0)
{
//value4 is datasq
if(dichot==0){value4=datasqs[j];}
else{value4=datasqs[j+m*data_length];}

//get XjTresiduals
sum=YTdata[j-bitstart+p*bitlength]+effs[(size_t)p*data_length+j]*value4;

if(p2<num_try)	//training models - model depends on mode, using estimated heritability
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
else	//MC REML models - using ridge model and fixed heritability
{
//value and value3 are per-predictor and noise term
value=exps[j+m*data_length]*tryhers[m2+m*num_hers2];
value3=1-tryhers[m2+m*num_hers2];

postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, NULL);
}

if(pens[p]!=pens[p]){printf("Error p %d bit %d j %d, m %d; please tell doug\n", p+1, bit+1, j+1, m+1);exit(1);}
if(isinf(pens[p])){printf("Error p %d bit %d j %d; please tell doug\n", p+1, bit+1, j+1);exit(1);}
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
}
}	//end of p loop

for(m=0;m<num_resps_use;m++)
{
if(bitrun[bit+m*bittotal]>0)	//update residuals
{
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_small, &bitlength, &alpha, data, &num_samples_use, changes+m*num_small*bitlength, &bitlength, &beta, residuals+(size_t)m*num_small*num_samples_use, &num_samples_use);

if(num_fixed>1)	//regress out covariates
{
if(dichot==0){reg_covar_matrix(residuals+(size_t)m*num_small*num_samples_use, Z, num_samples_use, num_small, num_fixed);}
else{reg_covar_weighted(residuals+(size_t)m*num_small*num_samples_use, Z, num_samples_use, num_small, num_fixed, nullweights+m*num_samples_use);}
}
}
}

//save and update current (partial) approx likelihood
#pragma omp parallel for private(p,p2,m,m2) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];

p2=p%num_small;
m=(p-p2)/num_small;
m2=(p2-num_try)%num_hers2;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
if(p2<num_try)	//using estimated heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights+m*num_samples_use, -9999);}
else	//using fixed heritability
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-tryhers[m2+m*num_hers2], bitpens[p], dichot, nullweights+m*num_samples_use, -9999);}
}
}

//see which phenotypes have converged
for(m=0;m<num_resps_use;m++)
{
bitrun[bit+m*bittotal]=0;
for(p=0;p<num_small;p++){bitrun[bit+m*bittotal]+=(fabs(likes[p+m*num_small]-likesold[p+m*num_small])>tol*num_samples_use);}
}

//see if breaking
cflag=0;for(p=0;p<total;p++){cflag+=(fabs(likes[p]-likesold[p])>tol*num_samples_use);}
if(cflag==0){break;}
}	//end of inner loop
count3+=count;
bitdet1[bit]=count;

//add current likelihood onto bitdiffs (so that it now stores difference)
for(p=0;p<total;p++){bitdiffs[bit+p*bittotal]+=likes[p];}

//update bitrun
for(m=0;m<num_resps_use;m++)
{
bitrun[bit+m*bittotal]=0;
for(p=0;p<num_small;p++){bitrun[bit+m*bittotal]+=(fabs(bitdiffs[bit+(p+m*num_small)*bittotal])>tol*num_samples_use);}
}

//update bitdo
bitdo[bit]=0;for(m=0;m<num_resps_use;m++){bitdo[bit]+=(bitrun[bit+m*bittotal]>0);}

//store penalties in bitpens
for(p=0;p<total;p++){bitpens[bit+p*bittotal]=pens[p];}

count2++;
}	//end of testing bit
else
{
bitdet1[bit]=0;
for(p=0;p<total;p++){bitdiffs[bit+p*bittotal]=0;}
}
}	//end of bit loop

printf("Average number of iterations per chunk: %.2f\n\n", (double)count3/total2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Average number of iterations per chunk: %.2f\n", (double)count3/total2);
fclose(output);

//update total2
total2=0;for(bit=0;bit<bittotal;bit++){total2+=(bitdo[bit]>0);}

/*
//append likelihood differences per bit
sprintf(filename4,"%s.likes.train", outfile);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
printf("print to %s for %d \n", filename4, bittotal);
for(bit=0;bit<bittotal;bit++)
{
fprintf(output4,"%d %d %d %d", count4+1, bit+1, bitdet1[bit], bitdet2[bit]);
for(p=0;p<total;p++){fprintf(output4," %f", bitdiffs[bit+p*bittotal]);}
fprintf(output4,"\n");
}
fclose(output4);
*/

////////

//eflag indicates whether smallest mse has converged (will not adjust for padding, but should affect all models equally)
eflag=1;

for(m=0;m<num_resps_use;m++)
{
for(p=0;p<num_try;p++)	//when computing mse, remember we started residuals for test samples at zero
{
sumsq=0;sumsq2=0;
for(i=0;i<num_test;i++)
{
if(Y[keeptest[i]+m*num_samples_use]!=missingvalue||multi==1)
{
if(dichot==0)	//compute regular MSE
{
sumsq+=pow(Yadj[keeptest[i]+m*num_samples_use]+residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]],2);
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2);
}
else	//compute weighted MSE
{
sumsq+=pow(Yadj[keeptest[i]+m*num_samples_use]+residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]],2)*nullweights[i+m*num_samples_use];
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2)*nullweights[i+m*num_samples_use];
}
}
}
value=sumsq/sumsq2;

if(p==0){value2=value;}
if(value<value2){value2=value;}
}
//if(mpheno!=-1){printf("The current lowest MSE is %.4f\n", value2);}

//set eflag=0 if not converged, then save
if(fabs(value2-Mmses[m])>0.005){eflag=0;}
Mmses[m]=value2;
}

//hflag indicates whether heritabilities have converged
hflag=1;

if(revher==1)	//compute revised heritabilities and set hflag=0 if not converged
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
/*
if(dichot==2)	//numerator is t(WR_Y) P_Y for real phenotypes
{
sumsq=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=residuals[(size_t)p*num_samples_use+i]*nullweights[i+m*num_samples_use]*(Yadj[i+m*num_samples_use]-residuals[(size_t)p*num_samples_use+i]);
}
value=sumsq;
}
*/

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

/*
if(dichot==2)	//denominator is mean of t(r) P_r for random vectors
{
sumsq=0;
for(i=0;i<num_samples_use;i++)
{
sumsq+=gaussian[i+g*num_samples_use]*(gaussian[i+g*num_samples_use]-residuals[(size_t)p*num_samples_use+i]);
}
sum+=sumsq;
}
*/
}
value2=sum/nmcmc;

polates[m2]=value/value2;
}

//can now get the updated value
value=inter_hers(num_hers2, polates, tryhers+m*num_hers2, 1);

if(value!=-9999)	//worked - set hflag=0 if not converged, then save
{
//if(mpheno!=-1){printf("The current MC REML heritability estimate is %.4f\n", value);}

if(fabs(value-hersold[m])>0.005){hflag=0;}
hersold[m]=value;
}
else	//failed - set hflag=1 and reset estimate
{
if(mpheno!=-1){printf("Warning, MC REML has failed this iteration\n");}
else{printf("Warning, MC REML has failed this iteration for Phenotype %d\n", m+1);}

hflag=1;
hersold[m]=hers[m];
}
}	//end of m loop
}	//end of revising heritability

if(eflag==1&&hflag==1){break;}

count4++;
}	//end of outer loop

////////

//correct residuals of cv models
for(m=0;m<num_resps_use;m++)
{
for(p=0;p<num_try;p++)
{
for(i=0;i<num_test;i++){residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]]+=Yadj[keeptest[i]+m*num_samples_use];}
}
}

if(multi==0)	//correct residuals of cv models for padded values (will reverse this after testing accuracy)
{
for(m=0;m<num_resps_use;m++)
{
if(respcounts[m]<num_samples_use)	//need to go from R=Y-Xbeta to R=Y-Xbeta/p - means adding on (1/p-1)(R-Y)
{
value=(double)num_samples_use/respcounts[m]-1;
for(p=0;p<num_try;p++)
{
for(i=0;i<num_samples_use;i++){residuals[(size_t)(p+m*num_small)*num_samples_use+i]+=value*(residuals[(size_t)(p+m*num_small)*num_samples_use+i]-Yadj[i+m*num_samples_use]);}
}
}
}
}

//measure accuracy, recording which has lowest MSE
if(mpheno!=-1){printf("Measuring accuracy of each model\n");}
else{printf("Finding the most accurate model for each phenotype\n");}

for(m=0;m<num_resps_use;m++)
{
if(verbose==1)
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
if(Y[keeptest[i]+m*num_samples_use]!=missingvalue||multi==1)
{
sumsq+=pow(residuals[(size_t)(p+m*num_small)*num_samples_use+keeptest[i]],2);
sumsq2+=pow(Yadj[keeptest[i]+m*num_samples_use],2);
}
}
value=sumsq/sumsq2;

if(p==0){Mmses[m]=value;Mmses2[m]=value;Mbests[m]=0;Mneffs[m]=respcounts[m]/value;}
if(value<Mmses[m]){Mmses[m]=value;Mbests[m]=p;Mneffs[m]=respcounts[m]/value;}

//compute correlation
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;indcount=0;
for(i=0;i<num_test;i++)
{
if(Y[keeptest[i]+m*num_samples_use]!=missingvalue||multi==1)
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

if(verbose==1)	//save all models
{
if(mode==151){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], value, value3);}
if(mode==152){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], tryps[p], tryf2s[p], value, value3);}
if(mode==153){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], value, value3);}
if(mode==154){fprintf(output3,"%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, hers[m], 1-tryp3s[p], tryf2s[p], value, value3);}
}
}	//end of p loop
if(verbose==1){fclose(output3);}

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
printf("Warning, the MSE for Phenotype %d is only %.4f, so will exclude the polygenic contribution\n\n", m+1, Mmses[m]);
Mincs[m]=0;Mneffs[m]=respcounts[m];
}
}
}

//at bottom is code for regressing multivariate traits (did not help)

if(multi==0)	//uncorrect residuals of cv models for padded values
{
for(m=0;m<num_resps_use;m++)
{
if(respcounts[m]<num_samples_use)	//need to go from R=Y-Xbeta/p to R=Y-Xbeta - means adding on (p-1)(R-Y)
{
value=(double)respcounts[m]/num_samples_use-1;
for(p=0;p<num_try;p++)
{
for(i=0;i<num_samples_use;i++){residuals[(size_t)(p+m*num_small)*num_samples_use+i]+=value*(residuals[(size_t)(p+m*num_small)*num_samples_use+i]-Yadj[i+m*num_samples_use]);}
}
}
}
}

if(revher==1)	//update hers and print out / save
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

if(verbose==1)	//save
{
sprintf(filename3,"%s.hers", outfile);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++){fprintf(output3,"MCREML %d %.4f %.4f NA\n", m+1, powers[Mtops[m]], hers[m]);}
fclose(output3);
}
}

time(&midtime);
printf("Time check: have so far spent %.2f hours\n\n", (double)(midtime-starttime)/60/60);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Time check: have so far spent %.2f hours\n\n", (double)(midtime-starttime)/60/60);
fclose(output);
}	//end of making training models

////////

if(skipcv==0)	//make final models (will have best) and maybe loco and calibration models - note that have already set cors
{
total=num_full*num_resps_use;

/*
//save likelihood difference per bit
sprintf(filename4,"%s.likes.final", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Scan Bit Iterations Chr");
for(p=0;p<total;p++){fprintf(output4, " Model%d", p+1);}
fprintf(output4,"\n");
fclose(output4);
*/

//save best models from cv stage
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<data_length;j++){effs2[(size_t)m*data_length+j]=effs[(size_t)(Mbests[m]+m*num_small)*data_length+j];}
for(i=0;i<num_samples_use;i++){residuals2[(size_t)m*num_samples_use+i]=residuals[(size_t)(Mbests[m]+m*num_small)*num_samples_use+i];}
}

//save ridge models from cv stage
for(m=0;m<num_resps_use;m++)
{
for(j=0;j<data_length;j++){effs3[(size_t)m*data_length+j]=effs[(size_t)(m*num_small)*data_length+j];}
for(i=0;i<num_samples_use;i++){residuals3[(size_t)m*num_samples_use+i]=residuals[(size_t)(m*num_small)*num_samples_use+i];}
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

if(usecomp==1)	//add in ridge models
{
for(p=0;p<num_chr;p++)
{
p2=1+num_chr+ncal+p+m*num_full;
for(j=0;j<data_length;j++){effs[(size_t)p2*data_length+j]=effs3[(size_t)m*data_length+j];probs[(size_t)p2*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)p2*num_samples_use+i]=residuals3[(size_t)m*num_samples_use+i];}
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
printf("Will also make %d LOCO PRS", (1+usecomp)*num_chr);
if(usecal==1){printf(", and %d calibration models", ncal);}
printf("\n");
}
printf("\n");
}
else
{
printf("Constructing final PRS for each phenotype\n");
if(loco==1)
{
printf("Will also make %d LOCO PRS", (1+usecomp)*num_chr);
if(usecal==1){printf(", and %d calibration models", ncal);}
printf("\n");
}
printf("\n");
}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output,"Constructing final PRS using all samples\n");}
else{fprintf(output,"Constructing final PRS for %d phenotypes using all samples\n", num_resps_use);}
if(loco==1){fprintf(output, "Will also make %d LOCO PRS", (1+usecomp)*num_chr);}
if(usecal==1){fprintf(output, ", and %d calibration models", ncal);}
fprintf(output, "\n");
fclose(output);

//set bit penalties to zero (will be incorrect, but we need a value)
for(p=0;p<total;p++)
{
for(bit=0;bit<bittotal;bit++){bitpens[bit+p*bittotal]=0;}
}

//bitrun records how many models running for each phenotypes
for(m=0;m<num_resps_use;m++)
{
for(bit=0;bit<bittotal;bit++){bitrun[bit+m*bittotal]=num_full;}
}

//bitdo records how many phenotypes running for each bit
for(bit=0;bit<bittotal;bit++){bitdo[bit]=num_resps_use;}

//total2 records how many bits to visit
total2=bittotal;

count4=0;	//number of scans
while(total2>0)
{
if(count4==nscan){printf("Warning, Variational Bayes did not converge after %d scans (this is not normally a problem)\n\n", nscan);break;}

//ready for bit loop
count2=0;	//number of bits tested
count3=0;	//total number of iterations
for(bit=0;bit<bittotal;bit++)
{
if(bitdo[bit]>0)	//will consider bit
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

bitdet2[bit]=chr[bitstart];

if(count2%200==0)
{
printf("Scan %d: estimating final effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Scan %d: estimating final effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);
fclose(output);
}

//compute current (partial) approx likelihood and save minus this in bitdiffs
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
p2=p%num_full;
m=(p-p2)/num_full;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], bitpens[bit+p*bittotal], dichot, nullweights+m*num_samples_use, -9999);
}

bitdiffs[bit+p*bittotal]=-likes[p];
}

//read data, standardize and set missing to zero - already have statistics
if(dtype==1)	//fast way
{
(void)read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, NULL, NULL, NULL, 0, maxthreads);}
else	//slow way
{
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, NULL, NULL, NULL, num_samples_use, bitlength, missingvalue, -9999, -9999, -9999, NULL, 2);
}

//will not adjust for covariates (can instead adjust residuals when necessary)

//iterate effect sizes
count=0;
while(count<maxiter)
{
count++;

//reset pens to zero for each iteration
for(p=0;p<total;p++){pens[p]=0;}

for(m=0;m<num_resps_use;m++)
{
if(bitrun[bit+m*bittotal]>0)
{
if(dichot==0||num_fixed==1)	//get t(X) residuals (if dichot=1 and no covariates, weights will be one)
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_full, &num_samples_use, &alpha, data, &num_samples_use, residuals+(size_t)m*num_full*num_samples_use, &num_samples_use, &beta, YTdata+m*num_full*bitlength, &bitlength);
}
else	//get t(X) W residuals - first multiply columns of residuals by corresponding weights, then premultiply by t(X)
{
copy_matrix(num_samples_use, num_full, residuals+(size_t)m*num_full*num_samples_use, residuals2+(size_t)m*num_full*num_samples_use, 1, nullweights+m*num_samples_use);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_full, &num_samples_use, &alpha, data, &num_samples_use, residuals2+(size_t)m*num_full*num_samples_use, &num_samples_use, &beta, YTdata+m*num_full*bitlength, &bitlength);
}
}
}

#pragma omp parallel for private(p,p2,m,j,j2,sum,value,value2,value3,value4,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
//bit active records if any loco predictors
bitactive[p]=0;

p2=p%num_full;
m=(p-p2)/num_full;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
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

if(p2<1+num_chr)	//model depends on model
{
if(mode==151)	//ridge
{postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, probs+(size_t)p*data_length+j);}
if(mode==152)	//bolt
{postmean=get_postmean(sum, lambdas[Mbests[m]]*value, lambdas2[Mbests[m]]*value, -9999, -9999, value4, value3, tryps[Mbests[m]], tryp2s[Mbests[m]], -9999, -9999, pens+p, 4, probs+(size_t)p*data_length+j);}
if(mode==153)	//bayesr
{postmean=get_postmean(sum, lambdas[Mbests[m]]*value, lambdas2[Mbests[m]]*value, lambdas3[Mbests[m]]*value, lambdas4[Mbests[m]]*value, value4, value3, tryps[Mbests[m]], tryp2s[Mbests[m]], tryp3s[Mbests[m]], tryp4s[Mbests[m]], pens+p,  5+(pointmass==0), probs+(size_t)p*data_length+j);}
if(mode==154)	//elastic
{postmean=get_postmean(sum, lambdas[Mbests[m]]*value2, lambdas2[Mbests[m]]*value2, lambdas3[Mbests[m]]*value, -9999, value4, value3, tryps[Mbests[m]], tryp2s[Mbests[m]], tryp3s[Mbests[m]], -9999, pens+p, 7, probs+(size_t)p*data_length+j);}
}
else	//using ridge model
{
postmean=get_postmean(sum, value, -9999, -9999, -9999, value4, value3, -9999, -9999, -9999, -9999, pens+p, 3, NULL);
}

if(pens[p]!=pens[p]){printf("Error p %d bit %d j %d; please tell doug\n", p+1, bit+1, j+1);exit(1);}
if(isinf(pens[p])){printf("Error p %d bit %d j %d\n; please tell doug", p+1, bit+1, j+1);exit(1);}

bitactive[p]=1;
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
}
}	//end of p loop

for(m=0;m<num_resps_use;m++)
{
if(bitrun[bit+m*bittotal]>0)	//update residuals
{
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_full, &bitlength, &alpha, data, &num_samples_use, changes+m*num_full*bitlength, &bitlength, &beta, residuals+(size_t)m*num_full*num_samples_use, &num_samples_use);

if(num_fixed>1)	//regress out covariates
{
if(dichot==0){reg_covar_matrix(residuals+(size_t)m*num_full*num_samples_use, Z, num_samples_use, num_full, num_fixed);}
else{reg_covar_weighted(residuals+(size_t)m*num_full*num_samples_use, Z, num_samples_use, num_full, num_fixed, nullweights+m*num_samples_use);}
}
}
}

//save and update current (partial) approx likelihood
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];

p2=p%num_full;
m=(p-p2)/num_full;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
if(bitactive[p]==1)	//only need to update if some active predictors
{likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights+m*num_samples_use, -9999);}
}
}

//see which phenotypes have converged
for(m=0;m<num_resps_use;m++)
{
bitrun[bit+m*bittotal]=0;
for(p=0;p<num_full;p++){bitrun[bit+m*bittotal]+=(fabs(likes[p+m*num_full]-likesold[p+m*num_full])>tol*num_samples_use);}
}

//see if breaking
cflag=0;for(p=0;p<total;p++){cflag+=(fabs(likes[p]-likesold[p])>tol*num_samples_use);}
if(cflag==0){break;}
}	//end of inner loop
count3+=count;
bitdet1[bit]=count;

//add current likelihood onto bitdiffs (so that it now stores difference)
for(p=0;p<total;p++){bitdiffs[bit+p*bittotal]+=likes[p];}

//update bitrun
for(m=0;m<num_resps_use;m++)
{
bitrun[bit+m*bittotal]=0;
for(p=0;p<num_full;p++){bitrun[bit+m*bittotal]+=(fabs(bitdiffs[bit+(p+m*num_full)*bittotal])>tol*num_samples_use);}
}

//update bitdo
bitdo[bit]=0;for(m=0;m<num_resps_use;m++){bitdo[bit]+=(bitrun[bit+m*bittotal]>0);}

//store penalties in bitpens
for(p=0;p<total;p++){bitpens[bit+p*bittotal]=pens[p];}

count2++;
}	//end of testing bit
else
{
bitdet1[bit]=0;
for(p=0;p<total;p++){bitdiffs[bit+p*bittotal]=0;}
}
}	//end of bit loop

printf("Average number of iterations per chunk: %.2f\n\n", (double)count3/total2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Average number of iterations per chunk: %.2f\n", (double)count3/total2);
fclose(output);

//update total2
total2=0;for(bit=0;bit<bittotal;bit++){total2+=(bitdo[bit]>0);}

/*
//append likelihood differences per bit
sprintf(filename4,"%s.likes.final", outfile);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
for(bit=0;bit<bittotal;bit++)
{
fprintf(output4,"%d %d %d %d", count4+1, bit+1, bitdet1[bit], bitdet2[bit]);
for(p=0;p<total;p++){fprintf(output4," %f", bitdiffs[bit+p*bittotal]);}
fprintf(output4,"\n");
}
fclose(output4);
*/

count4++;
}	//end of outer loop

////////

//see if excluding any prs
for(m=0;m<num_resps_use;m++)
{
if(Mincs[m]==0)	//blank the corresponding effect sizes and residuals
{
for(p=0;p<1+num_chr;p++)
{
for(j=0;j<data_length;j++){effs[(size_t)(p+m*num_full)*data_length+j]=0;}
for(i=0;i<num_samples_use;i++){residuals[(size_t)(p+m*num_full)*num_samples_use+i]=Yadj[i+m*num_samples_use];}
}
}
}

if(usecal==0)	//set lambda to one, and work out caleff
{
for(m=0;m<num_resps_use;m++)
{
cgammas[m]=1;
csds[m]=0;

if(dichot==0)	//caleff is root(MSE / var(R_Y)) - R_Y has mean zero
{
sumsq=0;
for(i=0;i<num_samples_use;i++){sumsq+=pow(residuals[(size_t)m*num_full*num_samples_use+i],2);}
var=sumsq/num_samples_use;
ceffs[m]=pow(Mmses[m]/var,.5);
}
else	//caleff is root (var(root(W) Y) MSE / var(root(R_Y)))
{
sum=0;sum2=0;sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
sum+=Yadj[i+m*num_samples_use]*pow(nullweights[i+m*num_samples_use],.5);
sum2+=residuals[(size_t)m*num_full*num_samples_use+i]*pow(nullweights[i+m*num_samples_use],.5);
sumsq+=pow(Yadj[i+m*num_samples_use],2)*nullweights[i+m*num_samples_use];
sumsq2+=pow(residuals[(size_t)m*num_full*num_samples_use+i],2)*nullweights[i+m*num_samples_use];
}
mean=sum/num_samples_use;
mean2=sum2/num_samples_use;
var=sumsq/num_samples_use-pow(mean,2);
var2=sumsq2/num_samples_use-pow(mean2,2);
ceffs[m]=pow(var*Mmses[m]/var2,.5);
}
}
}
else	//compute using grammar-gamma formula (or one if not using PRS)
{
for(m=0;m<num_resps_use;m++)
{
if(Mincs[m]==1)
{
//ridge scaling is GTG/(GT invV G) x var(invV Y) / t(Y) invV Y x n = t(G) G / t(G) R_G x nvar(R_Y) / t(Y) R_Y
//for dichot=1, it is GTWG/(GT invV G) x var(1/sqrt(W) invV Y) / t(Y) invV Y x n = t(G) WG / t(G) WR_G x nvar(sqrt(W) R_Y) / t(Y) WR_Y
//for dichot=2, it is simply GTWG/(GT invV G) = t(G) WG / t(G) WR_G
//effect calibration always first ratio

sum=0;sum2=0;value3=0;count=0;
for(j=0;j<ncal;j++)
{
if(cmults[j]!=-9999)
{
//p and p2 index residuals corresponding to R_Y and R_G
if(mode==151){p=1;}
else{p=1+num_chr+ncal;}
while(chr[cindex[j]]>chrindex[p]){p++;}
p2=1+num_chr+j;

if(dichot==0)	//sumsq is t(G) G; sumsq2 is t(G) R_G; sumsq3 is t(R_Y) R_Y; sumsq4 is t(Y) R_Y - (R_Y has mean zero)
{
sumsq=0;sumsq2=0;sumsq3=0;sumsq4=0;
for(i=0;i<num_samples_use;i++)
{
if(Y[i+m*num_samples_use]!=missingvalue||multi==1)
{
sumsq+=pow(cdata[i+(j+m*ncal)*num_samples_use],2);
sumsq2+=cdata[i+(j+m*ncal)*num_samples_use]*residuals[(size_t)(p2+m*num_full)*num_samples_use+i];
sumsq3+=pow(residuals[(size_t)(p+m*num_full)*num_samples_use+i],2);
sumsq4+=Yadj[i+m*num_samples_use]*residuals[(size_t)(p+m*num_full)*num_samples_use+i];
}}
value=sumsq/sumsq2*sumsq3/sumsq4;
value2=sumsq/sumsq2;
}
if(dichot==1)	//sumsq is t(G) WG; sumsq2 is t(G) WR_G; sum3 + sumsq3 give var(sqrt(W) R_Y); sumsq4 is t(Y) WR_Y
{
sum3=0;sumsq=0;sumsq2=0;sumsq3=0;sumsq4=0;
for(i=0;i<num_samples_use;i++)
{
if(Y[i+m*num_samples_use]!=missingvalue||multi==1)
{
sum3+=residuals[(size_t)(p+m*num_full)*num_samples_use+i]*pow(nullweights[i+m*num_samples_use],.5);
sumsq+=pow(cdata[i+(j+m*ncal)*num_samples_use],2)*nullweights[i+m*num_samples_use];
sumsq2+=cdata[i+(j+m*ncal)*num_samples_use]*residuals[(size_t)(p2+m*num_full)*num_samples_use+i]*nullweights[i+m*num_samples_use];
sumsq3+=pow(residuals[(size_t)(p+m*num_full)*num_samples_use+i],2)*nullweights[i+m*num_samples_use];
sumsq4+=Yadj[i+m*num_samples_use]*residuals[(size_t)(p+m*num_full)*num_samples_use+i]*nullweights[i+m*num_samples_use];
}}
value=sumsq/sumsq2*(sumsq3-sum3/num_samples_use*sum3)/sumsq4;
value2=sumsq/sumsq2;
}
if(dichot==2)	//sumsq is t(G) WG; sumsq2 is t(G) WR_G;
{
sumsq=0;sumsq2=0;
for(i=0;i<num_samples_use;i++)
{
if(Y[i+m*num_samples_use]!=missingvalue||multi==1)
{
sumsq+=pow(cdata[i+(j+m*ncal)*num_samples_use],2)*nullweights[i+m*num_samples_use];
sumsq2+=cdata[i+(j+m*ncal)*num_samples_use]*residuals[(size_t)(p2+m*num_full)*num_samples_use+i]*nullweights[i+m*num_samples_use];
}}
value=sumsq/sumsq2;
value2=sumsq/sumsq2;
}

//if(dichot==2){printf("lambda for %d value is %f\n", j+1, value);}
//else{printf("lambda for %d value is %f, var is %f sums %f %f %f\n", j+1, value, sumsq4/num_samples_use/(1-hers[m]), sumsq, sumsq2, sumsq4);}

sum+=value;
sum2+=pow(value,2);
value3+=value2;
count++;
}
}
mean=sum/count;
var=sum2/count-pow(mean,2); 

cgammas[m]=mean;
csds[m]=pow(var/ncal,.5);
ceffs[m]=value3/count;

if(mpheno!=-1){printf("The ridge regression test statistic scaling factor is %.4f (SE %.4f)\n", cgammas[m], csds[m]);}
else{printf("Phenotype %d: the ridge regression test statistic scaling factor is %.4f (SE %.4f)\n", m+1, cgammas[m], csds[m]);}

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output, "The ridge regression test statistic scaling factor is %.4f (SE %.4f)\n", cgammas[m], csds[m]);}
else{fprintf(output, "Phenotype %d: the ridge regression test statistic scaling factor is %.4f (SE %.4f)\n", m+1, cgammas[m], csds[m]);}
fclose(output);
}
else{cgammas[m]=1;csds[m]=0;ceffs[m]=1;}
}	//end of m loop
printf("\n");
}	//end of usecal=1

if(usecomp==1)	//revise scaling factors - for simplicity, will always assume dichot=0
{
//no need to revise if all prs ridge or excluded (can not be using fracfile, so first model always ridge)
flag=0;
for(m=0;m<num_resps_use;m++){flag+=(Mbests[m]!=0&&Mincs[m]==1);}

if(flag>0)
{
printf("Comparing test statistics for %d predictors\n", ncomp);

//pick ncomp predictors at random
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){eindex[count]=j;count++;}
}

if(ncomp>count)
{
ncomp=count;
printf("Warning, the number of comparison predictors has been reduced to %d (the total number of non-trivial predictors)\n", ncomp);
}

permute_int(eindex,count);
qsort(eindex,ncomp,sizeof(int), compare_int);
for(j=0;j<ncomp;j++){eindex2[j]=keeppreds_use[eindex[j]];}

if(num_blocks==-9999){num_blocks=100;}
if(num_blocks>ncomp){num_blocks=ncomp;}

stats=malloc(sizeof(double)*ncomp*2*num_resps_use);
stats2=malloc(sizeof(double)*ncomp*2*num_resps_use);
stats3=malloc(sizeof(double)*num_blocks);

bittotal2=(ncomp-1)/bitsize+1;

for(bit=0;bit<bittotal2;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>ncomp){bitend=ncomp;}
bitlength=bitend-bitstart;

if(bit%200==0)
{
printf("Computing test statistics for Chunk %d of %d\n", bit+1, bittotal2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Computing test statistics for Chunk %d of %d\n", bit+1, bittotal2);
fclose(output);
}

//read data and set missing to mean
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, eindex2, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, ecentres+bitstart, emults+bitstart, esqdevs+bitstart, erates+bitstart, einfos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 3);

//adjust predictors for covariates and standardize
reg_covar_matrix(data, Z, num_samples_use, bitlength, num_fixed);
stand_matrix_nomiss(data, num_samples_use, num_samples_use, bitlength);

for(m=0;m<num_resps_use;m++)	//get two sets of test statistics
{
value3=pow(ceffs[m],2)/cgammas[m];
for(j=bitstart;j<bitend;j++)
{
//p and p2 index residuals corresponding to mixture and ridge PRS
p=1;while(chr[eindex[j]]>chrindex[p]){p++;}
p2=1+num_chr+ncal;while(chr[eindex[j]]>chrindex[p2]){p2++;}

//get uncalibrated test statistic for mixture residuals (sumsq is XTR, sumsq2 is RTR)
sumsq=ddot_(&num_samples_use, data+(size_t)(j-bitstart)*num_samples_use, &one, residuals+(size_t)(p+m*num_full)*num_samples_use, &one);
sumsq2=ddot_(&num_samples_use, residuals+(size_t)(p+m*num_full)*num_samples_use, &one, residuals+(size_t)(p+m*num_full)*num_samples_use, &one);

value=sumsq/num_samples_use;
value2=(sumsq2/num_samples_use-pow(value,2))/(num_samples_use-num_fixed-1);
stats[j+m*2*ncomp]=pow(value,2)/value2;
stats2[j+m*2*ncomp]=value2;

//get something proportional to test statistic for ridge residuals (sumsq is XTR, sumsq2 is RTR)
sumsq=ddot_(&num_samples_use, data+(size_t)(j-bitstart)*num_samples_use, &one, residuals+(size_t)(p2+m*num_full)*num_samples_use, &one);
sumsq2=ddot_(&num_samples_use, residuals+(size_t)(p2+m*num_full)*num_samples_use, &one, residuals+(size_t)(p2+m*num_full)*num_samples_use, &one);

value=sumsq/num_samples_use;
value2=(sumsq2/num_samples_use-pow(value,2))/(num_samples_use-num_fixed-1);
stats[j+ncomp+m*2*ncomp]=cgammas[m]*pow(value,2)/value2;
stats2[j+ncomp+m*2*ncomp]=value3*value2;
}
}
}	//end of bit loop
printf("\n");

//save values
/*
for(m=0;m<num_resps_use;m++)
{
if(mpheno!=-1){sprintf(filename3,"%s.stats", outfile);}
else{sprintf(filename3,"%s.pheno%d.stats", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Mixture Ridge\n");

for(j=0;j<ncomp;j++){fprintf(output3,"%s %.6f %.6f\n", preds[eindex[j]], stats[j+m*2*ncomp], stats[j+ncomp+m*2*ncomp]);}
fclose(output3);
}
*/

for(m=0;m<num_resps_use;m++)
{
if(Mbests[m]!=0&&Mincs[m]==1)	//revise scaling factor so that sums are equal
{
//value is how much to scale current gamma - starting estimate is one
value=1;
count2=0;
while(count2<100)
{
//consider predictors where both calibrated stats under the threshold
sum=0;sum2=0;count=0;
for(j=0;j<ncomp;j++)
{
if(value*cgammas[m]*stats[j+m*2*ncomp]<=cthresh&stats[j+ncomp+m*2*ncomp]<=cthresh){sum+=stats[j+m*2*ncomp];sum2+=stats[j+ncomp+m*2*ncomp];count++;}
}
if(count==0){printf("Error, no stats under the value (phen %d ncomp %d thresh %f); please tell Doug\n\n", m+1, ncomp, cthresh);exit(1);}

value2=value;
value=sum2/sum/cgammas[m];

if(fabs(value-value2)<.001){break;}
count2++;
}

//get a jackknife variance (using final threshold)
for(p=0;p<num_blocks;p++)
{
start=(double)(p)/num_blocks*count;
end=(double)(p+1)/num_blocks*count;

sum=0;sum2=0;count2=0;
for(j=0;j<ncomp;j++)
{
if(value*cgammas[m]*stats[j+m*2*ncomp]<=cthresh&stats[j+ncomp+m*2*ncomp]<=cthresh)
{
if(count2<start||count2>=end){sum+=stats[j+m*2*ncomp];sum2+=stats[j+ncomp+m*2*ncomp];}
count2++;
}
}
stats3[p]=sum2/sum/cgammas[m];
}

sum=0;sumsq=0;
for(p=0;p<num_blocks;p++){sum+=stats3[p];sumsq+=pow(stats3[p],2);}
value2=(num_blocks-1)*(sumsq/num_blocks-pow(sum/num_blocks,2));

//revised scaling is cgammas[m] x value
value3=cgammas[m];
value4=pow(csds[m],2);
cgammas[m]=value*value3;
csds[m]=pow(value4*value2+value4*pow(value,2)+value2*pow(value3,2),.5);

//revise effect calibration so ceffs^2 = cgammas var(ridge) / var(mixture) times improvement in mse
sum=0;sum2=0;
for(j=0;j<ncomp;j++){sum+=stats2[j+m*2*ncomp];sum2+=stats2[j+ncomp+m*2*ncomp];}
ceffs[m]=pow(cgammas[m]*sum2/sum*Mmses[m]/Mmses2[m],.5);


/*
if(dougvar==1)	//alt approach - based on expected increase relative to ridge - similar to quickdraws
{
printf("Current estimate is %f %f - try to change target %f\n", cgammas[m], csds[m], Mmses2[m]/Mmses[m]);
sum=0;sum2=0;
for(j=0;j<ncomp;j++){sum+=cgammas[m]*stats[j+m*2*ncomp];sum2+=stats[j+ncomp+m*2*ncomp]-1;}
value=(Mmses2[m]/Mmses[m]*sum2+ncomp)/sum;

printf("New estimates is %f ratio %f\n", cgammas[m]*value, value);
cgammas[m]*=value;

sum=0;sum2=0;
for(j=0;j<ncomp;j++){sum+=cgammas[m]*stats[j+m*2*ncomp]-1;sum2+=stats[j+ncomp+m*2*ncomp]-1;}
printf("check %f vs %f\n", sum/sum2, Mmses2[m]/Mmses[m]);
}
*/

if(mpheno!=-1){printf("The revised test statistic scaling factor is %.4f (SE %.4f)\n", cgammas[m], csds[m]);}
else{printf("Phenotype %d: the revised test statistic scaling factor is %.4f (SE %.4f)\n", m+1, cgammas[m], csds[m]);}

if(csds[m]>0.02)	//very imprecise
{
if(mode==152){printf("Warning, this is very imprecise, it may be better to repeat using \"--ridge\" instead of \"--bolt\"\n");}
if(mode==153){printf("Warning, this is very imprecise, it may be better to repeat using \"--ridge\" instead of \"--bayesr\"\n");}
if(mode==154){printf("Warning, this is very imprecise, it may be better to repeat using \"--ridge\" instead of \"--elastic\"\n");}
}
}}	//end of updating scaling and m loop
printf("\n");

free(stats);
free(stats2);
free(stats3);
}
}	//end of usecomp=1

////////

if(multi==0)	//correct effects and residuals of gwide and loco models for padded values
{
for(m=0;m<num_resps_use;m++)
{
if(respcounts[m]<num_samples_use)
{
//divide effects by p
value=(double)num_samples_use/respcounts[m];
for(p=0;p<num_full;p++)
{
for(j=0;j<data_length;j++){effs[(size_t)(p+m*num_full)*data_length+j]*=value;}
}

//residuals need to go from R=Y-Xbeta to R=Y-Xbeta/p - means adding on (1/p-1)(R-Y)
value=(double)num_samples_use/respcounts[m]-1;
for(i=0;i<num_samples_use;i++){residuals[(size_t)m*num_full*num_samples_use+i]+=value*(residuals[(size_t)m*num_full*num_samples_use+i]-Yadj[i+m*num_samples_use]);}
for(p=0;p<num_full;p++)
{
for(i=0;i<num_samples_use;i++){residuals[(size_t)(p+m*num_full)*num_samples_use+i]+=value*(residuals[(size_t)(p+m*num_full)*num_samples_use+i]-Yadj[i+m*num_samples_use]);}
}
}
}
}

if(multi==1)	//untransform genome-wide prs and effects
{
for(m=0;m<num_resps_use;m++)
{
for(i=0;i<num_samples_use;i++){residuals2[(size_t)m*num_samples_use+i]=residuals[(size_t)m*num_full*num_samples_use+i];}
}

token=num_samples_use*num_full;
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, residuals2, &num_samples_use, Umat2, &num_resps_use, &beta, residuals, &token);

for(m=0;m<num_resps_use;m++)
{
for(j=0;j<data_length;j++){effs2[(size_t)m*data_length+j]=effs[(size_t)m*num_full*data_length+j];}
}

token=data_length*num_full;
alpha=1.0;beta=0.0;
dgemm_("N", "T", &data_length, &num_resps_use, &num_resps_use, &alpha, effs2, &data_length, Umat2, &num_resps_use, &beta, effs, &token);
}

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

if(verbose==1)	//save some more things
{
//save probabilities for final models (do not make sense for multi=1)
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
if(Y[i+m*num_samples_use]!=missingvalue)
{fprintf(output3, "%s %s %.4f %.4f\n", ids1[i], ids2[i], Yadj[i+m*num_samples_use]/Mscales[m], (Yadj[i+m*num_samples_use]-residuals[(size_t)m*num_full*num_samples_use+i])/Mscales[m]);}
}
fclose(output3);
}

if(loco==1)	//save LOCO PRS and details (saving depends on value of multi)
{
if(multi==0){value=Mscales[m];}
else{value=1;}

if(mpheno!=-1){sprintf(filename3,"%s.loco.prs", outfile);}
else{sprintf(filename3,"%s.pheno%d.loco.prs", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

fprintf(output3,"FID IID");
for(p=1;p<1+num_chr;p++){fprintf(output3," Chr%d", chrindex[p]);}
fprintf(output3,"\n");

if(multi==0)	//ignore missing phenotypes and scale
{
for(i=0;i<num_samples_use;i++)
{
if(Y[i+m*num_samples_use]!=missingvalue)
{
fprintf(output3, "%s %s", ids1[i], ids2[i]);
for(p=1;p<1+num_chr;p++)
{fprintf(output3," %.4f", (Yadj[i+m*num_samples_use]-residuals[(size_t)(p+m*num_full)*num_samples_use+i])/Mscales[m]);}
fprintf(output3,"\n");
}
}
}
else	//do not ignore missing nor scale
{
for(i=0;i<num_samples_use;i++)
{
fprintf(output3, "%s %s", ids1[i], ids2[i]);
for(p=1;p<1+num_chr;p++)
{fprintf(output3," %.4f", (Yadj[i+m*num_samples_use]-residuals[(size_t)(p+m*num_full)*num_samples_use+i]));}
fprintf(output3,"\n");
}
}
fclose(output3);

if(mpheno!=-1){sprintf(filename3,"%s.loco.details", outfile);}
else{sprintf(filename3,"%s.pheno%d.loco.details", outfile, m+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(multi==0)
{
fprintf(output3,"Actual_Sample_Size %d\n", respcounts[m]);
fprintf(output3,"Approx_Effective_Sample_Size %.1f\n", Mneffs[m]);
}
else
{
fprintf(output3,"Actual_Sample_Size %d\n", num_samples_use);
fprintf(output3,"Approx_Effective_Sample_Size %.1f\n", Mneffs[m]*num_samples_use/respcounts[m]);
}
if(Mincs[m]==1)
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
if(Mincs[m]==1){fprintf(output3,"Effect_Calibration %.4f\n", ceffs[m]);}
else{fprintf(output3,"Effect_Calibration 1\n");}
fclose(output3);
}
}	//end of m loop

if(verbose==0)
{
if(mpheno!=-1){printf("Best-fitting model saved in %s.effects\n\n", outfile);}
else{printf("Best-fitting models saved in %s.phenoX.effects, where X is the phenotype number\n\n", outfile);}
}
else
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

if(skipcv==1)	//make all models using all samples - will not be performing loco or using multi=1
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

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mpheno!=-1){fprintf(output,"Constructing %d PRS using all samples\n", total);}
else{fprintf(output,"Constructing %d PRS for each phenotype using all samples\n", num_try);}
fclose(output);

//set bit penalties to zero
for(p=0;p<total;p++)
{
for(bit=0;bit<bittotal;bit++){bitpens[bit+p*bittotal]=0;}
}

//bitdo says how many models not converged for each bit
for(bit=0;bit<bittotal;bit++){bitdo[bit]=total;}

//total2 records how many bits to visit
total2=bittotal;

count4=0;	//number of scans
while(total2>0)
{
if(count4==nscan)
{
printf("Warning, Variational Bayes did not converge after %d scans (this is not normally a problem)\n\n", nscan);
break;
}

//ready for bit loop
count2=0;	//number of bits tested
count3=0;	//total number of iterations
for(bit=0;bit<bittotal;bit++)
{
if(bitdo[bit]>0)	//will consider bit
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

bitdet2[bit]=chr[bitstart];

if(count2%200==0)
{
printf("Scan %d: estimating effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Scan %d: estimating effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);
fclose(output);
}

//compute current (partial) approx likelihood and save minus this in bitdiffs
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
p2=p%num_try;
m=(p-p2)/num_try;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], bitpens[bit+p*bittotal], dichot, nullweights+m*num_samples_use, -9999);
}

bitdiffs[bit+p*bittotal]=-likes[p];
}

//read data, standardize and set missing to zero - already have statistics
if(dtype==1)	//fast way
{(void)read_bed_wrapper(datafile, data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, keepsamps, bitlength, keeppreds_use+bitstart, num_samples, num_preds, missingvalue, NULL, NULL, NULL, 0, maxthreads);}
else	//slow way
{
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, NULL, NULL, NULL, num_samples_use, bitlength, missingvalue, -9999, -9999, -9999, NULL, 2);
}

//will not adjust for covariates (can instead adjust residuals when necessary)

//iterate effect sizes
count=0;
while(count<maxiter)
{
count++;

//reset pens to zero for each iteration
for(p=0;p<total;p++){pens[p]=0;}

for(m=0;m<num_resps_use;m++)
{
if(bitrun[bit+m*bittotal]>0)
{
if(dichot==0||num_fixed==1)	//get t(X) residuals (if dichot=1 and no covariates, weights will be one)
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_full, &num_samples_use, &alpha, data, &num_samples_use, residuals+(size_t)m*num_full*num_samples_use, &num_samples_use, &beta, YTdata+m*num_full*bitlength, &bitlength);
}
else	//get t(X) W residuals - first multiply columns of residuals by corresponding weights, then premultiply by t(X)
{
copy_matrix(num_samples_use, num_full, residuals+(size_t)m*num_full*num_samples_use, residuals2+(size_t)m*num_full*num_samples_use, 1, nullweights+m*num_samples_use);

alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &num_full, &num_samples_use, &alpha, data, &num_samples_use, residuals2+(size_t)m*num_full*num_samples_use, &num_samples_use, &beta, YTdata+m*num_full*bitlength, &bitlength);
}
}
}

#pragma omp parallel for private(p,p2,m,j,j2,sum,value,value2,value3,value4,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
p2=p%num_try;
m=(p-p2)/num_try;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
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

if(pens[p]!=pens[p]){printf("Error p %d bit %d j %d, m %d; please tell doug\n", p+1, bit+1, j+1, m+1);exit(1);}
if(isinf(pens[p])){printf("Error p %d bit %d j %d; please tell doug\n", p+1, bit+1, j+1);exit(1);}

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
}
}	//end of p loop

for(m=0;m<num_resps_use;m++)
{
if(bitrun[bit+m*bittotal]>0)	//update residuals
{
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &num_try, &bitlength, &alpha, data, &num_samples_use, changes+m*num_try*bitlength, &bitlength, &beta, residuals+(size_t)m*num_try*num_samples_use, &num_samples_use);

if(num_fixed>1)	//regress out covariates
{
if(dichot==0){reg_covar_matrix(residuals+(size_t)m*num_try*num_samples_use, Z, num_samples_use, num_try, num_fixed);}
else{reg_covar_weighted(residuals+(size_t)m*num_try*num_samples_use, Z, num_samples_use, num_try, num_fixed, nullweights+m*num_samples_use);}
}
}
}

//save and update current (partial) approx likelihood
#pragma omp parallel for private(p,p2,m) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];

p2=p%num_try;
m=(p-p2)/num_try;

if(bitrun[bit+m*bittotal]>0)	//still using this phenotype
{
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights, -9999);
}
}

//see which phenotypes have converged
for(m=0;m<num_resps_use;m++)
{
bitrun[bit+m*bittotal]=0;
for(p=0;p<num_try;p++){bitrun[bit+m*bittotal]+=(fabs(likes[p+m*num_try]-likesold[p+m*num_try])>tol*num_samples_use);}
}

//see if breaking
cflag=0;for(p=0;p<total;p++){cflag+=(fabs(likes[p]-likesold[p])>tol*num_samples_use);}
if(cflag==0){break;}
}	//end of inner loop
count3+=count;
bitdet1[bit]=count;

//add current likelihood onto bitdiffs (so that it now stores difference)
for(p=0;p<total;p++){bitdiffs[bit+p*bittotal]+=likes[p];}

//update bitrun
for(m=0;m<num_resps_use;m++)
{
bitrun[bit+m*bittotal]=0;
for(p=0;p<num_try;p++){bitrun[bit+m*bittotal]+=(fabs(bitdiffs[bit+(p+m*num_try)*bittotal])>tol*num_samples_use);}
}

//update bitdo
bitdo[bit]=0;for(m=0;m<num_resps_use;m++){bitdo[bit]+=(bitrun[bit+m*bittotal]>0);}

//store penalties in bitpens
for(p=0;p<total;p++){bitpens[bit+p*bittotal]=pens[p];}

count2++;
}	//end of testing bit
}	//end of bit loop

printf("Average number of iterations per chunk: %.2f\n\n", (double)count3/total2);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Average number of iterations per chunk: %.2f\n", (double)count3/total2);
fclose(output);

//update total2
total2=0;for(bit=0;bit<bittotal;bit++){total2+=(bitdo[bit]>0);}

/*
//compute current (full) approx likelihood - necessary if comparing with fast=0
#pragma omp parallel for private(p,p2,m,bit) schedule(static)
for(p=0;p<total;p++)
{
p2=p%num_try;
m=(p-p2)/num_try;

pens[p]=0;for(bit=0;bit<bittotal;bit++){pens[p]+=bitpens[bit+p*bittotal];}
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[m], pens[p], dichot, nullweights, -9999);
}
*/

count4++;
}	//end of outer loop

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
if(Y[i+m*num_samples_use]!=missingvalue)
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

//old code for multiphens
/*
if(multi==1)	//compute multi-phenotype calibrations
{
for(m=0;m<num_resps_use;m++)	//estimating effects for phenotype m
{
//load best PRS into rows of PT and target phenotype into first column of Q
indcount=0;
for(i=0;i<num_test;i++)
{
if(Y[keeptest[i]+m*num_samples_use]!=missingvalue)
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

//compute PTP and PTQ (saved latter in mth column of Mcals)
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_resps_use, &num_resps_use, &indcount, &alpha, PT, &num_resps_use, PT, &num_resps_use, &beta, PTP, &num_resps_use);
dgemv_("N", &num_resps_use, &indcount, &alpha, PT, &num_resps_use, Q, &one, &beta, Mcals+m*num_resps_use, &one);

if(num_resps_use>1)	//add small amount of shrinkage
{
for(m2=0;m2<num_resps_use;m2++){PTP[m2+m2*num_resps_use]*=1.01;}
}

//get weights
(void)eigen_invert(PTP, num_resps_use, PTP2, 1, Mcals+m*num_resps_use, 1);
}

printf("Here is the regression matrix:\n");
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){printf("%.2f ", Mcals[m+m2*num_resps_use]);}
printf("\n");
}
printf("\n");

//save transformation, remembering scaling

sprintf(filename3,"%s.regressed.matrix", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(m=0;m<num_resps_use;m++)
{
for(m2=0;m2<num_resps_use;m2++){fprintf(output3, "%.6f ", Mcals[m+m2*num_resps_use]*Mscales[m2]);}
fprintf(output3, "\n");
}
fclose(output3);

//compute and save regressed (adjusted) phenotypes

alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_samples_use, &num_resps_use, &num_resps_use, &alpha, Yadj, &num_samples_use, Mcals, &num_resps_use, &beta, Ymat, &num_samples_use);

sprintf(filename3,"%s.regressed.pheno", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output3,"%s %s", ids1[i], ids2[i]);
for(m=0;m<num_resps_use;m++){fprintf(output3, " %.6f", Ymat[i+m*num_samples_use]);}
fprintf(output3,"\n");
}
fclose(output3);

printf("Regressed phenotypes saved in %s\n\n", filename3);
}	//end of multi=1

if(mpheno==-1&&verbose==1)	//also save multi-phen version
{
sprintf(filename3,"%s.multi.root", outfile);
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
fprintf(output3,"Num_Phenotypes %d\n", num_resps_use);
fprintf(output3,"Num_Samples_Used %d\n", num_samples_use);
fprintf(output3,"Num_Predictors_Used %d\n", data_length);
fprintf(output3,"Num_Chromosomes %d\n", num_chr);
fclose(output3);
}


if(mpheno==-1&&verbose==1)	//also save multi-phen versions
{
//save calibrations
sprintf(filename3,"%s.multi.pheno%d.weights", outfile, m+1);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Phenotype Effect\n");
for(m2=0;m2<num_resps_use;m2++){fprintf(output3,"%d %.4f\n", m2+1, Mcals[m2+m*num_resps_use]*Mscales[m2]/Mscales[m]);}
fclose(output3);

if(loco==1)	//save LOCO PRS and details
{
sprintf(filename3,"%s.multi.pheno%d.loco.prs", outfile, m+1);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"FID IID");
for(p=1;p<1+num_chr;p++){fprintf(output3," Chr%d", chrindex[p]);}
fprintf(output3,"\n");

for(i=0;i<num_samples_use;i++)
{
if(Y[i+m*num_samples_use]!=missingvalue)
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

sprintf(filename3,"%s.multi.pheno%d.loco.details", outfile, m+1);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Actual_Sample_Size %d\n", respcounts[m]);
fprintf(output3,"Approx_Effective_Sample_Size %.1f\n", Mcombs[m]);
if(Mincs[m]==1)
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
*/


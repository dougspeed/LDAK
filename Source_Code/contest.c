/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Obtain inverses by conjugate gradient descent, then variational bayes

///////////////////////////

//only have (use) one phenotype
total=1+ncal;

//squeeze down chrindex, shares and chrprops
for(j=0;j<ncal;j++)
{
chrindex[1+j]=chrindex[1+num_chr+nridge+j];
shares[1+j]=shares[1+num_chr+nridge+j];
chrprops[1+j]=chrprops[1+num_chr+nridge+j];
}

//test time to just read data two times

sprintf(filename2,"%s.start",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output2,"now\n");
fclose(output2);

count=0;
while(count<2)
{
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%200==0){printf("Reading data for Chunk %d of %d\n", bit+1, bittotal);}

//read data for chunk and standardize (and replace NAs with zero)
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, hwestand, NULL, 1);
}
count++;
}
printf("\n");

sprintf(filename2,"%s.mid",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output2,"now\n");
fclose(output2);

////////

cY=malloc(sizeof(double)*num_samples_use*total);
cX=malloc(sizeof(double)*num_samples_use*total);
cR=malloc(sizeof(double)*num_samples_use*total);
cP=malloc(sizeof(double)*num_samples_use*total);
cHP=malloc(sizeof(double)*bitsize*total);
cKP=malloc(sizeof(double)*num_samples_use*total);
cVP=malloc(sizeof(double)*num_samples_use*total);
cP=malloc(sizeof(double)*num_samples_use*total);
ctemps=malloc(sizeof(double)*total);
calphas=malloc(sizeof(double)*total);
cbetas=malloc(sizeof(double)*total);

results1a=malloc(sizeof(double)*100*total);
results1b=malloc(sizeof(double)*100*total);
results2a=malloc(sizeof(double)*100*total);
results2b=malloc(sizeof(double)*100*total);

//set cY to outcomes
for(i=0;i<num_samples_use;i++){cY[i]=Yadj[i];}
for(j=0;j<ncal;j++)
{
for(i=0;i<num_samples_use;i++){cY[(size_t)(1+j)*num_samples_use+i]=data2[(size_t)j*num_samples_use+i];}
}

sprintf(filename2,"%s.outcomes",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output2,"%s %s ", ids1[i], ids2[i]);
for(p=0;p<total;p++){fprintf(output2,"%f ", cY[(size_t)p*num_samples_use+i]);}
fprintf(output2, "\n");
}
fclose(output2);

////////

//compute inverses by cgd

//set cX to zero
for(p=0;p<total;p++)
{
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]=0;}
}

//set cR to cY
for(p=0;p<total;p++)
{
for(i=0;i<num_samples_use;i++){cR[(size_t)p*num_samples_use+i]=cY[(size_t)p*num_samples_use+i];}
}

//set ctemps = RTR
for(p=0;p<total;p++)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(cR[(size_t)p*num_samples_use+i],2);}
ctemps[p]=sumsq;
printf("start %d is %f\n", p+1, ctemps[p]/num_samples_use);
}

//set cP=R
for(p=0;p<total;p++)
{
for(i=0;i<num_samples_use;i++){cP[(size_t)p*num_samples_use+i]=cR[(size_t)p*num_samples_use+i];}
}

count=0;
while(count<50)
{
//compute K cP (multiplied by chrprops)

//set cKP to zero
for(p=0;p<total;p++)
{
#pragma omp parallel for private(i) schedule(static)
for(i=0;i<num_samples_use;i++){cKP[(size_t)p*num_samples_use+i]=0;}
}

for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%200==0){printf("Analyzing data for Chunk %d of %d\n", bit+1, bittotal);}

//read data for chunk and standardize (and replace NAs with zero)
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, hwestand, NULL, 1);

//premultiply cP by t(data)
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, cP, &num_samples_use, &beta, cHP, &bitlength);

//scale values by exps*her
for(p=0;p<total;p++)
{
#pragma omp parallel for private(j) schedule(static)
for(j=0;j<bitlength;j++)
{
if(chrindex[p]!=chr[bitstart+j]){cHP[j+p*bitlength]*=exps[bitstart+j];}
else{cHP[j+p*bitlength]=0;}
}
}

//premultiply cHP by data, and add on to cKP
alpha=1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total, &bitlength, &alpha, data, &num_samples_use, cHP, &bitlength, &beta, cKP, &num_samples_use);
}

//set cVP = h2 cKP + (1-h2) P (must divide K by chrprops)
for(p=0;p<total;p++)
{
#pragma omp parallel for private(i) schedule(static)
for(i=0;i<num_samples_use;i++){cVP[(size_t)p*num_samples_use+i]=shares[p]/chrprops[p]*cKP[(size_t)p*num_samples_use+i]+(1-shares[p])*cP[(size_t)p*num_samples_use+i];}
}

//calphas = ctemps / P cVP
for(p=0;p<total;p++)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=cP[(size_t)p*num_samples_use+i]*cVP[(size_t)p*num_samples_use+i];}
calphas[p]=ctemps[p]/sumsq;
}

//cX = cX + calphas cP
for(p=0;p<total;p++)
{
#pragma omp parallel for private(i) schedule(static)
for(i=0;i<num_samples_use;i++){cX[(size_t)p*num_samples_use+i]+=calphas[p]*cP[(size_t)p*num_samples_use+i];}
}

//cR = cR - calphas cVP
for(p=0;p<total;p++)
{
#pragma omp parallel for private(i) schedule(static)
for(i=0;i<num_samples_use;i++){cR[(size_t)p*num_samples_use+i]-=calphas[p]*cVP[(size_t)p*num_samples_use+i];}
}

//cbetas = RTR / ctemps
for(p=0;p<total;p++)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(cR[(size_t)p*num_samples_use+i],2);}
cbetas[p]=sumsq/ctemps[p];
}

//set ctemps = RTR
for(p=0;p<total;p++)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(cR[(size_t)p*num_samples_use+i],2);}
ctemps[p]=sumsq;
if(p<6){printf("now %d is %f\n", p+1, ctemps[p]/num_samples_use);}
}

//cP = cR + cbetas P
for(p=0;p<total;p++)
{
#pragma omp parallel for private(i) schedule(static)
for(i=0;i<num_samples_use;i++){cP[(size_t)p*num_samples_use+i]=cR[(size_t)p*num_samples_use+i]+cbetas[p]*cP[(size_t)p*num_samples_use+i];}
}

//save current estimates of cY cX/n and var(cX)
for(p=0;p<total;p++)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=cX[(size_t)p*num_samples_use+i]*cY[(size_t)p*num_samples_use+i];}
results1a[p+count*total]=sumsq/num_samples_use;
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(cX[(size_t)p*num_samples_use+i],2);}
results1b[p+count*total]=sumsq/num_samples_use;
}

diff=ctemps[0];
for(p=1;p<total;p++)
{
if(ctemps[p]>diff){diff=ctemps[p];}
}
printf("Iteration %d - difference is %f\n", count+1, diff/num_samples_use);
if(diff<tol){break;}

count++;
}
printf("\n");

printf("took count %d\n", count);

sprintf(filename2,"%s.res1a",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(j=0;j<count;j++)
{
for(p=0;p<total;p++){fprintf(output2,"%f ", results1a[p+j*total]);}
fprintf(output2,"\n");
}
fclose(output2);

sprintf(filename2,"%s.res1b",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(j=0;j<count;j++)
{
for(p=0;p<total;p++){fprintf(output2,"%f ", results1b[p+j*total]);}
fprintf(output2,"\n");
}
fclose(output2);

sprintf(filename2,"%s.res1c",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output2,"%s %s ", ids1[i], ids2[i]);
for(p=0;p<total;p++){fprintf(output2,"%f ", cX[(size_t)p*num_samples_use+i]);}
fprintf(output2, "\n");
}
fclose(output2);

////////

//now use vb

//set effs to zero
for(p=0;p<total;p++)
{
for(j=0;j<data_length;j++){effs[(size_t)p*data_length+j]=0;}
}

//set residuals to outcomes
for(i=0;i<num_samples_use;i++){residuals[i]=Yadj[i];}
for(j=0;j<ncal;j++)
{
for(i=0;i<num_samples_use;i++){residuals[(size_t)(1+j)*num_samples_use+i]=data2[(size_t)j*num_samples_use+i];}
}

//set bit penalties to zero (will be incorrect, but we need a value)
for(p=0;p<total;p++)
{
for(bit=0;bit<bittotal;bit++){bitpens[bit+p*bittotal]=0;}
}

//bitdo indicates whether to test each bit
for(bit=0;bit<bittotal;bit++){bitdo[bit]=total;}

//total2 counts how many bits to test
total2=bittotal;

count4=0;
while(total2>0)
{
if(count4==nscan){printf("Warning, Variational Bayes did not converge after %d scans (this is not normally a problem)\n\n", nscan);break;}

//ready for bit loop
count2=0;
count3=0;
for(bit=0;bit<bittotal;bit++)
{
if(bitdo[bit]>0)	//will test
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(count2%200==0)
{
printf("Scan %d: estimating final effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Scan %d: estimating final effect sizes for Chunk %d of %d\n", count4+1, count2+1, total2);
fclose(output);
}

//read data for chunk and standardize (and replace NAs with zero)
(void)read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, -9999, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, hwestand, NULL, 1);

if(count4==0)	//compute predictor-predictor covariances, from which can extract XTX
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors+(size_t)bitstart*bitsize, &bitsize);

for(j=bitstart;j<bitend;j++){datasqs[j]=cors[(size_t)j*bitsize+j-bitstart];}
}

//compute current (partial) approx likelihood and save minus this in bitdiffs
#pragma omp parallel for private(p) schedule(static)
for(p=0;p<total;p++)
{
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-shares[p], bitpens[bit+p*bittotal], 0, NULL, -9999);
bitdiffs[bit+p*bittotal]=-likes[p];
}

//iterate effect sizes
count=0;
while(count<maxiter)
{
count++;

//reset pens to zero for each iteration
for(p=0;p<total;p++){pens[p]=0;}

//get t(X) residuals
alpha=1.0;beta=0.0;
dgemm_("T", "N", &bitlength, &total, &num_samples_use, &alpha, data, &num_samples_use, residuals, &num_samples_use, &beta, YTdata, &bitlength);

#pragma omp parallel for private(p,m,value,value2,j,j2,sum,postmean) schedule(dynamic, 1)
for(p=0;p<total;p++)
{
for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0&&chr[j]!=chrindex[p])
{
//get XjTresiduals
sum=YTdata[j-bitstart+p*bitlength]+effs[(size_t)p*data_length+j]*datasqs[j];

postmean=get_postmean(sum, exps[j]*shares[p]/chrprops[p], -9999, -9999, -9999, datasqs[j], 1-shares[p], -9999, -9999, -9999, -9999, pens+p, 3, probs+(size_t)p*data_length+j);
}
else{postmean=0;}

//get difference, then update effects and YTdata for remaining predictors in bit
changes[j-bitstart+p*bitlength]=postmean-effs[(size_t)p*data_length+j];
effs[(size_t)p*data_length+j]=postmean;
for(j2=j+1;j2<bitend;j2++){YTdata[j2-bitstart+p*bitlength]-=changes[j-bitstart+p*bitlength]*cors[(size_t)j*bitsize+j2-bitstart];}
}	//end of j loop
}	//end of p loop

//update residuals
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &num_samples_use, &total, &bitlength, &alpha, data, &num_samples_use, changes, &bitlength, &beta, residuals, &num_samples_use);

//save and update current (partial) approx likelihood
#pragma omp parallel for private(p) schedule(static)
for(p=0;p<total;p++)
{
likesold[p]=likes[p];
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-hers[0], pens[p], 0, NULL, -9999);
}

//see if breaking
cflag=0;
for(p=0;p<total;p++){cflag+=(fabs(likes[p]-likesold[p])<tol);}
if(cflag==total){break;}
}	//end of inner loop
count3+=count;

//add current likelihood onto bitdiffs (so that it now stores difference)
for(p=0;p<total;p++){bitdiffs[bit+p*bittotal]+=likes[p];}

//store penalties in bitpens
for(p=0;p<total;p++){bitpens[bit+p*bittotal]=pens[p];}

count2++;
}	//end of testing bit
}	//end of bit loop

printf("Average number of iterations per chunk: %.2f\n\n", (double)count3/total2);

if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Average number of iterations per chunk: %.2f\n", (double)count3/total2);
fclose(output);

//compute current (full) approx likelihood
#pragma omp parallel for private(p,p2,m,bit) schedule(static)
for(p=0;p<total;p++)
{
pens[p]=0;for(bit=0;bit<bittotal;bit++){pens[p]+=bitpens[bit+p*bittotal];}
likes[p]=comp_like(num_samples_use, residuals+(size_t)p*num_samples_use, 1-shares[p], pens[p], 0, NULL, -9999);
}

//update bitdo and total2
total2=0;
for(bit=0;bit<bittotal;bit++)
{
bitdo[bit]=0;
for(p=0;p<total;p++){bitdo[bit]+=(fabs(bitdiffs[bit+p*bittotal])>tol);}
total2+=(bitdo[bit]>0);
}

//save current estimates of cY cX/n and var(cX)
for(p=0;p<total;p++)
{
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=residuals[(size_t)p*num_samples_use+i]/(1-shares[p])*cY[(size_t)p*num_samples_use+i];}
results2a[p+count4*total]=sumsq/num_samples_use;
sumsq=0;for(i=0;i<num_samples_use;i++){sumsq+=pow(residuals[(size_t)p*num_samples_use+i]/(1-shares[p]),2);}
results2b[p+count4*total]=sumsq/num_samples_use;
}

count4++;
}	//end of outer loop

printf("took count4 %d\n", count4);

sprintf(filename2,"%s.res2a",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(j=0;j<count4;j++)
{
for(p=0;p<total;p++){fprintf(output2,"%f ", results2a[p+j*total]);}
fprintf(output2,"\n");
}
fclose(output2);

sprintf(filename2,"%s.res2b",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(j=0;j<count4;j++)
{
for(p=0;p<total;p++){fprintf(output2,"%f ", results2b[p+j*total]);}
fprintf(output2,"\n");
}
fclose(output2);

sprintf(filename2,"%s.res2c",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
for(i=0;i<num_samples_use;i++)
{
fprintf(output2,"%s %s ", ids1[i], ids2[i]);
for(p=0;p<total;p++){fprintf(output2,"%f ", residuals[(size_t)p*num_samples_use+i]/(1-shares[p]));}
fprintf(output2, "\n");
}
fclose(output2);

////////

free(cY);free(cX);free(cR);free(cP);free(cHP);free(cKP);free(cVP);free(ctemps);free(calphas);free(cbetas);
free(results1a);free(results1b);free(results2a);free(results2b);

///////////////////////////


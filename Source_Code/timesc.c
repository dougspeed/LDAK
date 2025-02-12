/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Do multiplications for calculating scores - think it might be possible to estimate SE of correlation (predsds)

///////////////////////////

//sort out step
step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//allocate variables

data_warn2(bitsize,(1+savecounts)*num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);
if(savecounts==1){data2=malloc(sizeof(double)*num_samples_use*bitsize);}

blupmids=malloc(sizeof(double)*bitsize*num_scores);
blupprojs=malloc(sizeof(double)*num_samples_use*num_scores);

guesses=malloc(sizeof(double*)*(num_scores+1));	//last one is covs
for(k=0;k<num_scores+1;k++){guesses[k]=malloc(sizeof(double)*num_samples_use);}

if(loco==1)
{
chrindex=malloc(sizeof(int)*num_chr);
guesses2=malloc(sizeof(double**)*num_chr);
for(p=0;p<num_chr;p++)
{
guesses2[p]=malloc(sizeof(double*)*num_scores);
for(k=0;k<num_scores;k++){guesses2[p][k]=malloc(sizeof(double)*num_samples_use);}
}
}

if(prsvar==1){indsds=malloc(sizeof(double)*num_samples_use);}

if(savecounts==1)	//saving counts
{
nums=malloc(sizeof(int*)*num_scores);
for(k=0;k<num_scores;k++){nums[k]=malloc(sizeof(int)*num_samples_use);}
}

if(strcmp(respfile,"blank")!=0||strcmp(sumsfile,"blank")!=0)	//computing accuracy
{
predmeans=malloc(sizeof(double)*num_scores);
predvars=malloc(sizeof(double)*num_scores);
predcors=malloc(sizeof(double)*num_scores);
predsds=malloc(sizeof(double)*num_scores);
}

//set guesses to zero
for(k=0;k<num_scores+1;k++)
{
for(i=0;i<num_samples_use;i++){guesses[k][i]=0;}
}

if(loco==1)	//fill chrindex and set guesses2 to zero
{
chrindex[0]=chr[0];
count=1;
for(j=1;j<data_length;j++)	
{
if(chr[j]>chr[j-1]){chrindex[count]=chr[j];count++;}
}
if(count!=num_chr){printf("Doug error, 66HB3\n\n");exit(1);}

for(p=0;p<num_chr;p++)
{
for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){guesses2[p][k][i]=0;}
}
}
}

if(savecounts==1)	//set nums to zero
{
for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){nums[k][i]=0;}
}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

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
p=0;
while(bitend<data_length)
{
bitstart=bitend;
bitend=bitstart+bitsize;
if(bitend>data_length){bitend=data_length;}
while(chr[bitend-1]>chr[bitstart]){bitend--;}
bitlength=bitend-bitstart;

if(bit%step==0)
{
printf("Calculating scores for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating scores for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read data and get statistics (but do not standardize)
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, power, 0, hwestand, NULL, 0);

//update centres if available (maybe should also update mults, but will only matter when power!=0, which almost never happens)
for(j=0;j<bitlength;j++)
{
if(blupcentres[0][bitstart+j]!=-9999){centres[bitstart+j]=blupcentres[0][bitstart+j];}
}

if(savecounts==1)	//must copy data into data2
{copy_matrix(num_samples_use, bitlength, data, data2, 0, NULL);}

//now standardize data
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -9999, -9999, -9999, NULL, 2);

//load effect sizes into blupmids
for(k=0;k<num_scores;k++)
{
for(j=0;j<bitlength;j++){blupmids[j+k*bitlength]=blupfactors[k][bitstart+j];}
}

//get contribution to projections
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_samples_use, &num_scores, &bitlength, &alpha, data, &num_samples_use, blupmids, &bitlength, &beta, blupprojs, &num_samples_use);

//add these to guesses
for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){guesses[k][i]+=blupprojs[i+k*num_samples_use];}
}

if(loco==1)	//subtract from corresponding guesses2
{
while(chr[bitstart]!=chrindex[p]){p++;}

for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){guesses2[p][k][i]-=blupprojs[i+k*num_samples_use];}
}
}

if(savecounts==1)	//deal with nums
{
//first increase nums assuming all relevant values are present for each score
for(k=0;k<num_scores;k++)
{
count=0;for(j=0;j<bitlength;j++){count+=(blupfactors[k][bitstart+j]!=0&&mults[bitstart+j]!=-9999);}
for(i=0;i<num_samples_use;i++){nums[k][i]+=count;}
}

//now loop through values subtracting from nums for each missing
#pragma omp parallel for private(j,i,k) schedule (static)
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999&&rates[bitstart+j]<1)	//there are some missing values
{
for(i=0;i<num_samples_use;i++)
{
if(data2[(size_t)j*num_samples_use+i]==missingvalue)
{
for(k=0;k<num_scores;k++){nums[k][i]-=(blupfactors[k][bitstart+j]!=0);}
}}
}}
}
bit++;
}	//end of while loop
printf("\n");

if(strcmp(cofile,"blank")!=0)	//get contribution of covariates
{
alpha=1.0;beta=0.0;
dgemv_("N", &num_samples_use, &num_covars, &alpha, covar, &num_samples_use, thetas, &one, &beta, guesses[num_scores], &one);
}
//else will remain zero

if(prsvar==0)	//save
{
write_scores(guesses, indsds, nums, num_samples_use, num_scores, ids1, ids2, outfile, resp, missingvalue, savecounts);

if(loco==1)
{
sprintf(filename2,"%s.loco",outfile);

if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"FID IID");
for(k=0;k<num_scores;k++)
{
fprintf(output2," Profile");
for(p=0;p<num_chr;p++){fprintf(output2," Chr%d", chrindex[p]);}
}
fprintf(output2,"\n");

for(i=0;i<num_samples_use;i++)
{
fprintf(output2, "%s %s", ids1[i], ids2[i]);
for(k=0;k<num_scores;k++)
{
fprintf(output2," %.4f", guesses[k][i]);
for(p=0;p<num_chr;p++){fprintf(output2," %.4f", guesses[k][i]+guesses2[p][k][i]);}
}
fprintf(output2,"\n");
}
fclose(output2);

printf("LOCO profiles saved in %s\n\n", filename2);
}
}
else	//get sds then save
{
for(i=0;i<num_samples_use;i++)
{
sum=0;sumsq=0;
for(p=0;p<num_blocks;p++)
{sum+=guesses[1+p][i];sumsq+=pow(guesses[1+p][i],2);}
mean=sum/num_blocks;
var=(1.0-jackprop)/jackprop*(sumsq/num_blocks-pow(mean,2));
indsds[i]=pow(var,.5);
}

write_scores(guesses, indsds, nums, num_samples_use, num_scores, ids1, ids2, outfile, resp, missingvalue, 2);
}

////////

if(strcmp(respfile,"blank")!=0||strcmp(sumsfile,"blank")!=0)	//get variances and correlations
{
if(strcmp(respfile,"blank")!=0)	//use only individuals with phenotypes
{
//get mean, variance and indcount for phenotype (after subtracting covariates) - have checked not trivial
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(resp[i]!=missingvalue)
{sum+=resp[i]-guesses[num_scores][i];sumsq+=pow(resp[i]-guesses[num_scores][i],2);indcount++;}
}
mean=sum/indcount;
var=sumsq/indcount-pow(mean,2);

//now get vars of score and cors / sds
for(k=0;k<num_scores;k++)
{
sum2=0;sumsq2=0;sumsq3=0;
for(i=0;i<num_samples_use;i++)
{
if(resp[i]!=missingvalue)
{sum2+=guesses[k][i];sumsq2+=pow(guesses[k][i],2);sumsq3+=guesses[k][i]*(resp[i]-guesses[num_scores][i]);
}
}

predmeans[k]=sum2/indcount;
predvars[k]=sumsq2/indcount-pow(predmeans[k],2);
if(predvars[k]>0)	//score not trivial
{
predcors[k]=(sumsq3/indcount-mean*predmeans[k])*pow(var,-.5)*pow(predvars[k],-.5);
value=(1-pow(predcors[k],2))/indcount;
predsds[k]=pow(value,.5);
}
else{predcors[k]=-9999;predsds[k]=-9999;}
}
}
else	//so have summaries
{
//get neff - average number of individuals
sum=0;for(j=0;j<data_length;j++){sum+=nss[j];}
neff=sum/data_length;

for(k=0;k<num_scores;k++)
{
//get vars of score
sum=0;sumsq=0;
for(i=0;i<num_samples_use;i++){sum+=guesses[k][i];sumsq+=pow(guesses[k][i],2);}
predmeans[k]=sum/num_samples_use;
predvars[k]=sumsq/num_samples_use-pow(predmeans[k],2);

//get cov(score,Y)/SE(Y) = weighted sum of cov(Xj,Y)/SE(Y) - cov(Xj,Y)/SE(Y) = rhos[j]*SE(Xj)
value2=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){value2+=blupfactors[k][j]*rhos[j]*pow(sqdevs[j],.5)*mults[j];}
}
if(predvars[k]>0)
{
predcors[k]=value2*pow(predvars[k],-.5);
value=(1-pow(predcors[k],2))/neff;
predsds[k]=pow(value,.5);
}
else{predcors[k]=-9999;predsds[k]=-9999;}
}
}

//get best, min and max
best=-1;
for(k=0;k<num_scores;k++)
{
if(predcors[k]!=-9999)
{
if(best==-1){best=k;min=predcors[k];max=predcors[k];}
if(predcors[k]<min){min=predcors[k];}
if(predcors[k]>max){best=k;max=predcors[k];}
}
}

if(best==-1)	//not possible to compute a correlation for any models
{printf("Error, it was not possible to compute a correlation for any of the scores (suggesting they all have zero effect sizes)\n\n");exit(1);}

if(prsvar==0)	//save all
{
sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Profile\tCorrelation\n");
for(k=0;k<num_scores;k++)
{
if(predcors[k]!=-9999){fprintf(output2,"Score_%d\t%.6f\n", k+1, predcors[k]);}
else{fprintf(output2,"Score_%d\tNA\n", k+1);}
}
fclose(output2);

if(num_scores==1)
{printf("Correlation between score and phenotype is %.4f, saved in %s\n\n", predcors[0], filename2);}
else
{printf("Correlations between %d scores and phenotype range from %.4f to %.4f, saved in %s\n\n", num_scores, min, max, filename2);}
}
else	//save only first
{
sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Profile\tCorrelation\n");
if(predcors[0]!=-9999){fprintf(output2,"Score_%d\t%.6f\n", k+1, predcors[0]);}
else{fprintf(output2,"Score_%d\tNA\n", k+1);}
fclose(output2);

printf("Correlation between score and phenotype is %.4f, saved in %s\n\n", predcors[0], filename2);
}
}	//end of getting correlations for scores

////////

if(strcmp(finalfile,"blank")!=0)	//get new model
{
count=countrows(finalfile)-1;
printf("Extracting effects for %d predictors from %s\n\n", count, finalfile);

if((input=fopen(finalfile,"r"))==NULL)
{printf("Error opening %s\n\n",finalfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

sprintf(filename5,"%s.effects.best",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor A1 A2 Centre Model%d\n", best+1);

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %s %s %s ", readstring, readstring2, readstring3, readstring4)!=4)
{printf("Error reading first four values of Row %d of %s\n\n", j+1, finalfile);exit(1);}

for(k=0;k<num_scores;k++)
{
if(fscanf(input, "%lf ", &readdouble)!=1)
{printf("Error reading Effect %d from Row %d of %s\n\n", k+1, j+1, finalfile);exit(1);}
if(k==best){value=readdouble;}
}

fprintf(output5, "%s %s %s %s %.4e\n", readstring, readstring2, readstring3, readstring4, value);
}
fclose(input);
fclose(output5);

printf("The best-fitting model is saved in %s\n\n", filename5);
}

////////

//free allocations from setdl.c
for(k=0;k<num_scores;k++){free(blupcentres[k]);free(blupfactors[k]);}free(blupcentres);free(blupfactors);

//frees from above
free(data);
if(savecounts==1){free(data2);}
free(blupmids);free(blupprojs);
for(k=0;k<num_scores+1;k++){free(guesses[k]);}free(guesses);
if(loco==1)
{
free(chrindex);
for(p=0;p<num_chr;p++)
{
for(k=0;k<num_scores;k++){free(guesses2[p][k]);}
free(guesses2[p]);
}
free(guesses2);
}
if(prsvar==1){free(indsds);}
if(savecounts==1)
{
for(k=0;k<num_scores;k++){free(nums[k]);}free(nums);
}
if(strcmp(respfile,"blank")!=0||strcmp(sumsfile,"blank")!=0){free(predmeans);free(predvars);free(predsds);free(predcors);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Do multiplications for finding the best prediction model - will not have covariates, nor use savecount or power

///////////////////////////

//allocate variables

data_warn2(bitsize,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitsize);

blupmids=malloc(sizeof(double)*bitsize*num_scores);
blupprojs=malloc(sizeof(double)*num_samples_use*num_scores);

guesses=malloc(sizeof(double*)*num_scores);	//last one is covs
for(k=0;k<num_scores;k++){guesses[k]=malloc(sizeof(double)*num_samples_use);}

predcors=malloc(sizeof(double)*num_scores);

//set guesses
for(k=0;k<num_scores;k++)
{
for(i=0;i<num_samples_use;i++){guesses[k][i]=0;}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

printf("Calculating scores for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating scores for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

//read data
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);

//centre - can use provided mean by copying blupcentre into centres
for(j=bitstart;j<bitend;j++){centres[j]=blupcentres[0][j];}
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, num_samples_use, bitlength, missingvalue, 0, 1, -9999, NULL, 1, preds+bitstart);

//load scaled effect sizes into blupmids
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
}	//end of bit loop
printf("\n");

////////

//get correlations with phenotype
for(r=0;r<num_scores;r++)
{
sum=0;sum2=0;sumsq=0;sumsq2=0;sumsq3=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(resp[i]!=missingvalue)
{
sum+=guesses[r][i];sum2+=resp[i];
sumsq+=pow(guesses[r][i],2);sumsq2+=pow(resp[i],2);
sumsq3+=guesses[r][i]*resp[i];
indcount++;
}
}

if(sumsq>sum*sum/indcount)	//score might be trivial (but have checked phenotype is not, so indcount>1)
{
mean=sum/indcount;mean2=sum2/indcount;
predcors[r]=(sumsq3-indcount*mean*mean2)/pow(sumsq-indcount*mean*mean,.5)/pow(sumsq2-indcount*mean2*mean2,.5);
}
else{predcors[r]=-9999;}
}	//end of r loop

sprintf(filename2,"%s.cors",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
for(r=0;r<num_scores;r++)
{
if(predcors[r]!=-9999){fprintf(output2,"Score_%d %.6f\n", r+1, predcors[r]);}
else{fprintf(output2,"Score_%d NA\n", r+1);}
}
fclose(output2);

min=predcors[0];max=predcors[0];
for(r=1;r<num_scores;r++)
{
if(predcors[r]<min){min=predcors[r];}
if(predcors[r]>max){max=predcors[r];}
}
printf("Correlations between %d scores and phenotype range from %.4f to %.4f, saved in %s\n\n", num_scores, min, max, filename2);

////////

//first get best (checking at least one non-trivial score)
best=-9999;value=-9999;
for(r=0;r<num_scores;r++)
{
if(predcors[r]>value){best=r;value=predcors[r];}
}
if(best==-9999){printf("Error, all %d scores are trivial\n\n", num_scores);exit(1);}

//create new effects

count=countrows(scorefile)-1;
printf("Extracting effects for %d predictors from %s\n\n", count, scorefile);

if((input=fopen(scorefile,"r"))==NULL)
{printf("Error re-opening %s\n\n",scorefile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

sprintf(filename3,"%s.best.effects",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor A1 A2 Centre Model%d\n", best+1);

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %s %s %s ", readstring, readstring2, readstring3, readstring4)!=4)
{printf("Error reading first four values of Row %d of %s\n\n", j+1, finalfile);exit(1);}

for(r=0;r<num_scores;r++)
{
if(fscanf(input, "%lf ", &readdouble)!=1)
{printf("Error reading Effect %d from Row %d of %s\n\n", r+1, j+1, finalfile);exit(1);}
if(r==best){value=readdouble;}
}

fprintf(output3, "%s %s %s %s %.4e\n", readstring, readstring2, readstring3, readstring4, value);
}
fclose(input);
fclose(output3);

printf("The best-fitting model is saved in %s\n\n", filename3);

////////

//free allocations from setdl.c
for(k=0;k<num_scores;k++){free(blupcentres[k]);free(blupfactors[k]);}free(blupcentres);free(blupfactors);

//frees from above
free(data);
free(blupmids);free(blupprojs);
for(k=0;k<num_scores;k++){free(guesses[k]);}free(guesses);
free(predcors);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


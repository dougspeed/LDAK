/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculating correlations, saving those above a threshold - ignore self correlations (and maybe near duplicates)

///////////////////////////

if(mincor==-9999)	//find threshold corresponding to P=0.01
{
mincor=1-exp(-6.634897/num_samples_use);
printf("Will record pairs of predictors with correlation squared at least %.4e (this corresponds to P<0.01)\n\n", mincor);
}

if(bitsize==-9999)	//get optimal bitsize (average number of neighbours)
{
scount=0;
k=1;
for(j=0;j<data_length;j++)
{
for(k=j+1;k<data_length;k++)
{
if(chr[k]!=chr[j]||cmbp[k]-cmbp[j]>1000*window_kb){break;}
}
scount+=k-j-1;
}
bitsize=(int)((scount-1)/data_length)+1;

if(bitsize<20){bitsize=20;}
if(bitsize>8000){bitsize=8000;}
printf("The bit-size will be set to %d (you can change this using \"--bit-size\")\n\n", bitsize);
}

step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//work out bitmax
bitmax=bitsize;
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
if(bitend2-bitstart>bitmax){bitmax=bitend2-bitstart;}
}

//work out max number of correlations for each predictor
maxnums=malloc(sizeof(int)*data_length);
for(j=0;j<data_length;j++){maxnums[j]=0;}

//first go to left
k=data_length-2;
for(j=data_length-1;j>=0;j--)
{
while(1)
{
if(k==-1){break;}
if(chr[k]!=chr[j]||cmbp[j]-cmbp[k]>1000*window_kb){break;}
k--;
}
maxnums[j]+=j-k-1;
}

//now right
k=1;
for(j=0;j<data_length;j++)
{
while(1)
{
if(k==data_length){break;}
if(chr[k]!=chr[j]||cmbp[k]-cmbp[j]>1000*window_kb){break;}
k++;
}
maxnums[j]+=k-j-1;
}

//work out max allocations within any bit (as multiple of bitmax)
bitmax=bitsize;
smax=0;
bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
if(bitend2-bitstart>bitmax){bitmax=bitend2-bitstart;}

scount=0;for(j=bitstart;j<bitend2;j++){scount+=maxnums[j];}
if(scount/bitmax>smax){smax=scount/bitmax;}
}

////////

//allocate variables
data_warn3(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);

Z=malloc(sizeof(double)*num_samples_use*num_fixed);

anal_warn(bitmax, bitmax+smax);
cors=malloc(sizeof(double)*bitmax*bitmax);
bigs=malloc(sizeof(int*)*data_length);
rjks=malloc(sizeof(float*)*data_length);

actnums=malloc(sizeof(int)*data_length);
rjksums=malloc(sizeof(double)*data_length);
rjkaves=malloc(sizeof(double)*data_length);
rjktemp=malloc(sizeof(double)*data_length);

randnorms=malloc(sizeof(double)*num_samples_use);
datarands=malloc(sizeof(double)*data_length);

//fill covariates
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//set actnums to zero and rjksums to one
for(j=0;j<data_length;j++){actnums[j]=0;rjksums[j]=1;}

//set randnorms
for(i=0;i<num_samples_use;i++){randnorms[i]=rnorm_safe();}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//open bin file
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//write placeholders for counts, means, scalings and variances (latter three arrays were initialized in ldak.c)
fwrite(actnums, sizeof(int), data_length, output2);
fwrite(centres, sizeof(double), data_length, output2);
fwrite(mults, sizeof(double), data_length, output2);
fwrite(sqdevs, sizeof(double), data_length, output2);

//open high ld file
sprintf(filename3,"%s.cors.ldpairs",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Kept_Predictor Kept_A1 Kept_A2 Removed_Predictor Removed_A1 Removed_A2\n");

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
bitlength=bitend2-bitstart;

if(bit%step==0)
{
printf("Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating correlations for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}

shuffle=0;
for(j=0;j<end-bitstart;j++)	//using values already in data, so shuffle back
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(bitstart-start+j)*num_samples_use+i];}
shuffle++;
}

//allocate for predictors not in previous chunk
for(j=bitstart+shuffle;j<bitend2;j++)
{
if(maxnums[j]>0)
{
bigs[j]=malloc(sizeof(int)*maxnums[j]);
rjks[j]=malloc(sizeof(float)*maxnums[j]);
}
}

//read data and standardize  (maybe regressing out covariates)
current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, bitstart+shuffle, bitend2, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data+(size_t)shuffle*num_samples_use, centres+bitstart+shuffle, mults+bitstart+shuffle, sqdevs+bitstart+shuffle, num_samples_use, bitlength-shuffle, missingvalue, -1, 0, 0, NULL, 1, preds+bitstart+shuffle);

if(num_covars>1)	//regress out covariates - must then standardize again
{
reg_covar_matrix(data+(size_t)shuffle*num_samples_use, Z, num_samples_use, bitlength-shuffle, num_covars);
stand_matrix_nomiss(data+(size_t)shuffle*num_samples_use, num_samples_use, num_samples_use, bitlength-shuffle);
}

//get correlation
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors, &bitlength);

//now loop through predictors in chunk (have already set actnums to zero and rjksums to one)
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)	//non trivial
{
for(k=j+1;k<bitend2;k++)
{
if(chr[k]!=chr[j]){break;}
if(cmbp[k]-cmbp[j]>1000*window_kb){break;}

if(mults[k]!=-9999)	//non-trivial and within range
{
value=cors[(size_t)(k-bitstart)*bitlength+(j-bitstart)];
value2=pow(value,2);

if(value2>=maxcor)	//near duplicate - will ignore second predictor
{
mults[k]=-9999;
fprintf(output3, "%s\t%c\t%c\t%s\t%c\t%c\n", preds[j], al1[j], al2[j], preds[k], al1[k], al2[k]);

if(wcount<5){printf("Warning, Predictor %s is a near-duplicate (has correlation %.4f with %s) and will be ignored\n", preds[k], value2, preds[j]);}
wcount++;
}
else
{
if(value2>=mincor)	//above threshold
{
bigs[j][actnums[j]]=k;
rjks[j][actnums[j]]=value;
actnums[j]++;
rjksums[j]+=value2;

bigs[k][actnums[k]]=j;
rjks[k][actnums[k]]=value;
actnums[k]++;
rjksums[k]+=value2;
}
}
}}	//end of using j and k loops

if(actnums[j]>0)	//print correlations for this predictor
{
fwrite(bigs[j], sizeof(int), actnums[j], output2);
fwrite(rjks[j], sizeof(float), actnums[j], output2);
}
}}	//end of using j and j loop

//free predictors that will not be used again
for(j=bitstart;j<bitend;j++)
{
if(maxnums[j]>0){free(bigs[j]);free(rjks[j]);}
}

//get t(data) * randnorms * root(1/num_samples_use)
token=bitlength-shuffle;
alpha=pow(1.0/num_samples_use,.5);beta=0.0;
dgemv_("T", &num_samples_use, &token, &alpha, data+(size_t)shuffle*num_samples_use, &num_samples_use, randnorms, &one, &beta, datarands+bitstart+shuffle, &one);

start=bitstart;if(bitend2>end){end=bitend2;}
}	//end of bit loop
printf("\n");
if(wcount>5){printf("In total, %d predictors are near duplicates\n", wcount);}
if(wcount>0){printf("\n");}

fclose(output2);
fclose(output3);

//compute average ldscore for each window
bitstart=0;
while(bitstart<data_length)
{
for(bitend=bitstart+1;bitend<data_length;bitend++)
{
if(cmbp[bitend]-cmbp[bitstart]>1000*window_kb||chr[bitend]!=chr[bitstart]){break;}
}
bitlength=bitend-bitstart;

sum=0;
for(j=bitstart;j<bitend;j++){sum+=rjksums[j];}
mean=sum/bitlength;
for(j=bitstart;j<bitend;j++){rjkaves[j]=mean;}

bitstart=bitend;
}

//work out median
for(j=0;j<data_length;j++){rjktemp[j]=rjkaves[j];}
qsort(rjktemp, data_length, sizeof(double), compare_double);

//re-open cors file
if((output2=fopen(filename2,"rb+"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
fseeko(output2, 0, SEEK_SET);	//probably unnecessary, but not sure

//update counts, means, scalings and variances
fwrite(actnums, sizeof(int), data_length, output2);
fwrite(centres, sizeof(double), data_length, output2);
fwrite(mults, sizeof(double), data_length, output2);
fwrite(sqdevs, sizeof(double), data_length, output2);

fseeko(output2, 0, SEEK_END);
fclose(output2);

//save root
count=0;
scount=0;
for(j=0;j<data_length;j++){count+=(mults[j]!=-9999);scount+=actnums[j];}

if(count==data_length){printf("For each of the %d predictors, there are on average %.2f other predictors with correlation squared at least %.6f\n\n", data_length, (double)scount/data_length, mincor);}
else{printf("For each of the %d remaining predictors, there are on average %.2f other predictors with correlation squared at least %.6f\n\n", count, (double)scount/count, mincor);}

sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Num_Samples %d\nNum_Predictors %d\n", num_samples, num_preds);
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, data_length);
fprintf(output3,"Num_Pairs %jd\n", scount);
fprintf(output3,"Threshold %.4e\n", mincor);
if(window_cm!=-9999){fprintf(output3,"Window_cM %.4f\n", window_cm);}
else{fprintf(output3,"Window_kb %.4f\n", window_kb);}
fclose(output3);

//save bim
sprintf(filename4,"%s.cors.bim", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<data_length;j++)
{
fprintf(output4, "%d\t%s\t", chr[j], preds[j]);
if(cm[j]==0){fprintf(output4, "0\t");}
else{fprintf(output4, "%.6f\t", cm[j]);}
fprintf(output4, "%ld\t%c\t%c\n", (long int)bp[j], al1[j], al2[j]);
}
fclose(output4);

//save ldscores
sprintf(filename5,"%s.cors.ldscores",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor Individual Region\n");
for(j=0;j<data_length;j++){fprintf(output5, "%s %.4f %.4f\n", preds[j], rjksums[j], rjkaves[j]);}
fclose(output5);

//save highld predictors
sprintf(filename6,"%s.cors.highld",outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
count=0;
for(j=0;j<data_length;j++)
{
if(rjkaves[j]>5*rjktemp[data_length/2]){fprintf(output6, "%s\n", preds[j]);count++;}
}
if(count==0){fprintf(output6,"None\n");}
fclose(output6);

//save noise
sprintf(filename7,"%s.cors.noise",outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(j=0;j<data_length;j++){fprintf(output7, "%.10f\n", datarands[j]);}
fclose(output7);

printf("The correlations are saved in files with prefix %s\n\n", outfile);

free(maxnums);
free(data);
free(Z);
free(cors);free(bigs);free(rjks);
free(actnums);free(rjksums);free(rjkaves);free(rjktemp);
free(randnorms);free(datarands);
if(binary==0){gzclose(datainputgz);}


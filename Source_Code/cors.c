/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculating correlations, saving those above a threshold - ignore self correlations (and maybe near duplicates)

///////////////////////////

//divide predictors into blocks
blockstarts=malloc(sizeof(int)*data_length);
blockends=malloc(sizeof(int)*data_length);
indexer=malloc(sizeof(int)*data_length);

if(strcmp(blockfile,"blank")==0)	//divide based on window_kb (or window_cm)
{
blockstarts[0]=0;
bittotal=1;
for(j=1;j<data_length;j++)
{
if(chr[j]>chr[blockstarts[bittotal-1]]||cmbp[j]-cmbp[blockstarts[bittotal-1]]>=1000*window_kb)
{
blockends[bittotal-1]=j;
bittotal++;
blockstarts[bittotal-1]=j;
}
}
blockends[bittotal-1]=data_length;
}
else	//use blockfile
{
count=countrows(blockfile)-1;

//read in breakpoints, then add a final one
readints=malloc(sizeof(int)*(count+1));
readdoubles=malloc(sizeof(double)*(count+1));
read_integers(blockfile, readints, count, NULL, 1, 1, 0);
read_values(blockfile, readdoubles, count, NULL, 2, 1, 0);
readints[count]=chr[data_length-1]+1;
readdoubles[count]=1;

if(chr[0]>readints[count-1]||(chr[0]==readints[count-1]&&bp[0]>readdoubles[count-1]))
{printf("Warning, the location of the first predictor (Chr %d, BP %.2f) is after the last breakpoint in %s (Chr %d, BP %.2f)\n\n", chr[0], bp[0], blockfile, readints[count-1], readdoubles[count-1]);}

if(chr[data_length-1]<readints[0]||(chr[data_length-1]==readints[0]&&bp[data_length-1]<=readdoubles[0]))
{printf("Warning, the location of the last predictor (Chr %d, BP %.2f) is before the first breakpoint in %s (Chr %d, BP %.2f)\n\n", chr[data_length-1], bp[data_length-1], blockfile, readints[0], readdoubles[0]);}

//find first breakpoint after first predictor
for(count2=0;count2<count;count2++)
{
if(readints[count2]>chr[0]||(readints[count2]==chr[0]&&readdoubles[count2]>bp[0])){break;}
}

blockstarts[0]=0;
bittotal=1;
for(j=1;j<data_length;j++)
{
if(chr[j]>readints[count2]||(chr[j]==readints[count2]&&bp[j]>readdoubles[count2]))
{
blockends[bittotal-1]=j;
bittotal++;
blockstarts[bittotal-1]=j;
count2++;
}
else
{
if(chr[j]>chr[blockstarts[bittotal-1]])
{
blockends[bittotal-1]=j;
bittotal++;
blockstarts[bittotal-1]=j;
}
}
}
blockends[bittotal-1]=data_length;

free(readints);free(readdoubles);
}

//get bitmax
bitmax=blockends[0]-blockstarts[0];
for(bit=1;bit<bittotal;bit++)
{
if(blockends[bit]-blockstarts[bit]>bitmax){bitmax=blockends[bit]-blockstarts[bit];}
}

printf("The predictors have been divided into %d windows, with the largest containing %d predictors\n\n", bittotal, bitmax);

////////

//allocate variables

data_warn(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);

Z=malloc(sizeof(double)*num_samples_use*num_fixed);

anal_warn(bitmax, bitmax);
cors=malloc(sizeof(double)*bitmax*bitmax);
writefloats=malloc(sizeof(float)*data_length);

rjksums=malloc(sizeof(double)*data_length);
rjkaves=malloc(sizeof(double)*data_length);
rjktemp=malloc(sizeof(double)*data_length);

randnorms=malloc(sizeof(double)*num_samples_use);
datarands=malloc(sizeof(double)*data_length);

randnorms2=malloc(sizeof(double)*num_samples_use*num_blocks);
datarands2=malloc(sizeof(double)*data_length*num_blocks);

//fill covariates
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

//set rjksums to 1
for(j=0;j<data_length;j++){rjksums[j]=1;}

//set randnorms
for(i=0;i<num_samples_use;i++){randnorms[i]=rnorm_safe();}

//and randnorms2
for(p=0;p<num_blocks;p++)
{
for(i=0;i<num_samples_use;i++){randnorms2[i+p*num_samples_use]=rnorm_safe();}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress file
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//open temporary bin file
sprintf(filename2,"%s.cors.temp", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//ready for bit loop
scount=0;
ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=bitend-bitstart;

if(bit%10==0)
{
printf("Calculating correlations for Window %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating correlations for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read data and standardize
current=read_data_fly(datafile, dtype, data, NULL, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, -1, 0, 0, NULL, 1);

if(minmaf>0)	//check for low maf predictors
{
for(j=bitstart;j<bitend;j++)
{
maf=centres[j]/2+(centres[j]>1)*(1-centres[j]);
if(maf<minmaf)
{
mults[j]=-9999;
if(ecount<5){printf("Warning, Predictor %s has MAF %.6f and will be ignored\n", preds[j], maf);}
ecount++;
}
}
}

if(num_covars>1)	//regress out covariates - must then standardize again
{
reg_covar_matrix(data, Z, num_samples_use, bitlength, num_covars);
stand_matrix_nomiss(data, num_samples_use, num_samples_use, bitlength);
}

//get correlations
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors, &bitlength);

if(maxcor<1)	//remove near duplicates
{
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)	//not trivial or near duplicate 
{
for(k=j+1;k<bitend;k++)
{
if(mults[k]!=-9999)	//not trivial or near duplicate
{
value=pow(cors[(size_t)(j-bitstart)*bitlength+(k-bitstart)],2);

if(value>maxcor)	//near duplicate - will ignore second predictor
{
mults[k]=-9999;
if(wcount<5){printf("Warning, Predictor %s is a near-duplicate (has correlation squared %.4f with %s) and will be ignored\n", preds[k], value, preds[j]);}
wcount++;
}
}}	//end of using k and k loop
}}	//end of using j and j loop
}

if(strip>0){eigen_strip_mults(cors, bitlength, mults+bitstart, strip);}

//how many predictors remain?
count=0;for(j=bitstart;j<bitend;j++){count+=(mults[j]!=-9999);}
scount+=count*(count-1)/2;

if(count>0)	//save off-diagonal correlations for predictors not excluded (and add contributions to rjksums, maybe truncate)
{
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)
{
count2=0;
for(k=j+1;k<bitend;k++)
{
if(mults[k]!=-9999)
{
value=cors[(size_t)(j-bitstart)*bitlength+k-bitstart];
writefloats[count2]=value;
rjksums[j]+=pow(value,2);
rjksums[k]+=pow(value,2);
count2++;
}
}

if(mincor>0)	//truncate
{
for(k=0;k<count2;k++)
{
if(pow(writefloats[k],2)<mincor){writefloats[k]=0;}
}
}

fwrite(writefloats, sizeof(float), count2, output2);
}
}
}

//get t(data) * randnorms * root(1/num_samples_use)
alpha=pow(1.0/num_samples_use,.5);beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, randnorms, &one, &beta, datarands+bitstart, &one);

//and same for jacknifes (except a matrix)
alpha=pow(1.0/num_samples_use,.5);beta=0.0;
dgemm_("T", "N", &bitlength, &num_blocks, &num_samples_use, &alpha, data, &num_samples_use, randnorms2, &num_samples_use, &beta, datarands2+bitstart, &data_length);
}	//end of bit loop
printf("\n");
if(ecount>5){printf("In total, %d predictors have MAF below %.6f\n", ecount, minmaf);}
if(wcount>5){printf("In total, %d predictors are near duplicates\n", wcount);}
if(ecount>0||wcount>0){printf("\n");}

fclose(output2);

//work out number of predictors that remain and their indexes
total=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){indexer[total]=j;total++;}
}
if(total==0){printf("Error, after filtering, no predictors remain\n\n");exit(1);}

//open bin file and save counts, means, scalings, variances and rjksums, for predictors that remain
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

for(j=0;j<total;j++){fwrite(centres+indexer[j], sizeof(double), 1, output2);}
for(j=0;j<total;j++){fwrite(mults+indexer[j], sizeof(double), 1, output2);}
for(j=0;j<total;j++){fwrite(sqdevs+indexer[j], sizeof(double), 1, output2);}
for(j=0;j<total;j++){fwrite(rjksums+indexer[j], sizeof(double), 1, output2);}

//now add temp file to end of bin file
readfloats=malloc(sizeof(float)*1000000);
printf("Moving contents of %s.cors.temp to end of %s.cors.bin\n\n", outfile, outfile);

sprintf(filename,"%s.cors.temp", outfile);
if((input=fopen(filename,"r"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_SET);

smax=scount;
while(smax>1000000)
{
if(fread(readfloats, sizeof(float), 1000000, input)!=1000000)
{printf("Error reading correlations from %s\n\n", filename);exit(1);}
fwrite(readfloats, sizeof(float), 1000000, output2);
smax-=1000000;
}

if(smax>0)
{
if(fread(readfloats, sizeof(float), smax, input)!=smax)
{printf("Error reading correlations from %s\n\n", filename);exit(1);}
fwrite(readfloats, sizeof(float), smax, output2);
}

free(readfloats);
fclose(input);
fclose(output2);

//delete temp file
sprintf(cmd, "rm %s.cors.temp", outfile);
system(cmd);

////////

//compute average ldscore for each block
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];

sum=0;
count=0;
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999){sum+=rjksums[j];count++;}
}
if(count>0){mean=sum/count;}
else{mean=0;}
for(j=bitstart;j<bitend;j++){rjkaves[j]=mean;}
}

//work out median average ldscore (across predictors that remain)
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){rjktemp[count]=rjkaves[j];count++;}
}
qsort(rjktemp, count, sizeof(double), compare_double);
med=rjktemp[count/2];

//adjust blockstarts and blockends to take into account missing predictors
count=0;
for(bit=0;bit<bittotal;bit++)
{
for(j=blockstarts[bit];j<blockends[bit];j++){count+=(mults[j]==-9999);}
blockends[bit]-=count;
}
for(bit=1;bit<bittotal;bit++){blockstarts[bit]=blockends[bit-1];}

//only retain blocks that contain some predictors
count=0;
for(bit=0;bit<bittotal;bit++)
{
if(blockends[bit]>blockstarts[bit])
{
if(count!=bit){blockstarts[count]=blockstarts[bit];blockends[count]=blockends[bit];}
count++;
}
}
if(count<bittotal){printf("Due to predictor filtering, the number of windows is reduced to %d\n\n", count);}
bittotal=count;

//save root
sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Num_Samples %d\nNum_Predictors %d\n", num_samples, num_preds);
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, total);
fprintf(output3,"Num_Windows %d\n", bittotal);
fprintf(output3,"Num_Pairs %jd\n", scount);
fprintf(output3,"Num_Jackknifes %d\n", num_blocks);
fclose(output3);

//save bim
sprintf(filename4,"%s.cors.bim", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output4, "%d\t%s\t", chr[j], preds[j]);
if(cm[j]==0){fprintf(output4, "0\t");}
else{fprintf(output4, "%.6f\t", cm[j]);}
fprintf(output4, "%ld\t%c\t%c\n", (long int)bp[j], al1[j], al2[j]);
}
}
fclose(output4);

//save windows
sprintf(filename5,"%s.cors.windows", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5,"Window Start_Predictor End_Predictor Start_Location End_Location\n");
for(bit=0;bit<bittotal;bit++){fprintf(output5, "%d %d %d %d:%.0f %d:%.0f\n", bit+1, blockstarts[bit]+1, blockends[bit], chr[indexer[blockstarts[bit]]], bp[indexer[blockstarts[bit]]], chr[indexer[blockends[bit]-1]], bp[indexer[blockends[bit]-1]]);}
fclose(output5);

//save highld predictors
sprintf(filename6,"%s.cors.highld",outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999&&rjkaves[j]>2*med){fprintf(output6, "%s\n", preds[j]);count++;}
}
if(count==0){fprintf(output6,"None\n");}
fclose(output6);

//save noise
sprintf(filename7,"%s.cors.noise",outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){fprintf(output7, "%.6f\n", datarands[j]);}
}
fclose(output7);

//and jackknifes
sprintf(filename8,"%s.cors.jackknifes",outfile);
if((output8=fopen(filename8,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename8);exit(1);}
for(p=0;p<num_blocks;p++)
{
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){writefloats[count]=(float)datarands2[(size_t)p*data_length+j];count++;}
}
fwrite(writefloats, sizeof(float), count, output8);
}
fclose(output8);

printf("The correlations are saved in files with prefix %s\n\n", outfile);

free(blockstarts);free(blockends);free(indexer);
free(data);
free(Z);
free(cors);free(writefloats);
free(rjksums);free(rjkaves);free(rjktemp);
free(randnorms);free(datarands);free(randnorms2);free(datarands2);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


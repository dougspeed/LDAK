/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Calculating correlations

///////////////////////////

if(strcmp(betafile,"blank")!=0)    //will be making fake summaries
{
//allocate variables
gaussian=malloc(sizeof(double)*num_samples_use*num_phenos);
effs=malloc(sizeof(double)*data_length*num_phenos);
effs2=malloc(sizeof(double)*data_length*num_phenos);

//noise should be standard gaussian divided by root (num_samples_use x num_inds)
value=pow(num_samples_use,-.5)*pow(num_inds,-.5);
for(m=0;m<num_phenos;m++)
{
for(i=0;i<num_samples_use;i++){gaussian[i+m*num_samples_use]=rnorm_safe()*value;}
}

//read effects
read_betas(betafile, effs, num_phenos, data_length, preds, al1, al2, bimfile);
for(m=0;m<num_phenos;m++)
{
printf("First few effect sizes for Phenotype %d are;", m+1);
count=0;
for(j=0;j<data_length;j++)
{
if(effs[j+m*data_length]!=0)
{
if(count==0){printf(" %s %.4e", preds[j], effs[j+m*data_length]);}
else{printf(" | %s %.4e", preds[j], effs[j+m*data_length]);}
count++;
}
if(count==3){break;}
}
printf("\n");
if(m==5){break;}
}
printf("\n");
}

//will use five powers
num_pows=5;
powers=malloc(sizeof(double)*5);
powers[0]=-1.0;powers[1]=-0.75;powers[2]=-0.5;powers[3]=-0.25;powers[4]=0;

//get the heritability model (must have parttype=0, so addpart will always be set to one)
model_warn(data_length*3/2, num_parts+1);
pindexes=malloc(sizeof(int *)*(num_parts+1));
pweights=malloc(sizeof(double *)*(num_parts+1));
addpart=get_her_model(num_parts, partpref, pindexes, pweights, data_length, keeppreds_use, num_preds, allpreds, predorder, parttype, backpart, allone, -9999);

//get labels (will always have labfile if num_parts>0)
catlabels=malloc(sizeof(char *)*(num_parts+1));
for(q=0;q<num_parts+1;q++){catlabels[q]=malloc(sizeof(char)*500);}
strcpy(catlabels[num_parts],"Base");

if(num_parts>0) //read names from labfile
{
if((input=fopen(labfile,"r"))==NULL)
{printf("Error opening %s\n\n",labfile);exit(1);}
for(q=0;q<num_parts;q++)
{
catlabels[q]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", catlabels[q])!=1)
{printf("Error reading Row %d of %s\n\n", q+1, labfile);exit(1);}

if(strcmp(catlabels[q],"Base")==0||strcmp(catlabels[q],"Background")==0)
{printf("Error reading Row %d of %s; labels can not be either \"Base\" or \"Background\"\n\n", q+1, labfile);exit(1);}

for(q2=0;q2<q;q2++)	//check different to previous
{
if(strcmp(catlabels[q],catlabels[q2])==0)
{printf("Error reading %s; the label %s appears twice\n\n", labfile, catlabels[q]);exit(1);}
}
}
fclose(input);
}

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
if(chr[j]>readints[count2]||(chr[j]==readints[count2]&&bp[j]>readdoubles[count2]))	//finished with current block
{
blockends[bittotal-1]=j;
bittotal++;
blockstarts[bittotal-1]=j;
count2++;

//skip any empty blocks (must be at least one non empty at the end)
while(chr[j]>readints[count2]||(chr[j]==readints[count2]&&bp[j]>readdoubles[count2])){count2++;}
}
else
{
if(chr[j]>chr[blockstarts[bittotal-1]])	//this can only happen if there is a chromosome that does not contain any blocks
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

anal_warn(bitmax, bitmax);

Z=malloc(sizeof(double)*num_samples_use*num_fixed);
cors=malloc(sizeof(double)*bitmax*bitmax);
writeshorts=malloc(sizeof(unsigned short)*data_length);
writefloats=malloc(sizeof(float)*data_length);
exps=malloc(sizeof(double)*bitmax*num_pows);

rjksums=malloc(sizeof(double)*data_length*num_pows);
rjkaves=malloc(sizeof(double)*data_length);
rjktemp=malloc(sizeof(double)*data_length);

randnorms=malloc(sizeof(double)*num_samples_use);
datarands=malloc(sizeof(double)*data_length);

highlds=malloc(sizeof(int)*data_length);

//fill covariates
for(i=0;i<num_samples_use;i++)
{
for(j=0;j<num_fixed;j++){Z[i+j*num_samples_use]=covar[i+j*num_samples_use];}
}

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

if(strcmp(betafile,"blank")==0) //open temporary bin file
{
sprintf(filename2,"%s.cors.temp", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
}

//ready for bit loop
scount=0;
xcount=0;
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

//check for trivial, low maf or callrate predictors
for(j=bitstart;j<bitend;j++)
{
if(mults[j]==-9999)
{
if(xcount<5){printf("Warning, Predictor %s is trivial and will be ignored\n", preds[j]);}
xcount++;
}
else
{
maf=centres[j]/2+(centres[j]>1)*(1-centres[j]);
if(maf<minmaf||rates[j]<minobs)
{
mults[j]=-9999;
if(ecount<5){printf("Warning, Predictor %s has MAF %.6f and callrate %.4f and will be ignored\n", preds[j], maf, rates[j]);}
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

if(dougvar==99) //save as kinship matrix then exit
{
sprintf(filename3,"%s.grm.raw", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

for(k=0;k<bitlength;k++)
{
for(j=0;j<bitlength;j++){fprintf(output3,"%.6f ", cors[j+k*bitlength]);}
fprintf(output3,"\n");
}
fclose(output3);

sprintf(filename3,"%s.grm.id", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

for(j=0;j<bitlength;j++){fprintf(output3,"%s %s\n", preds[bitstart+j], preds[bitstart+j]);}
fclose(output3);

printf("have saved kinships for %d snps\n\n", bitlength);
exit(1);
}

if(strcmp(betafile,"blank")==0) //normal case - saving correlations
{
if(maxcor<1)	//remove near duplicates
{
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)	//not trivial, low maf or near duplicate 
{
for(k=j+1;k<bitend;k++)
{
if(mults[k]!=-9999)	//not trivial, low maf or near duplicate (yet)
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

if(strip>0) //strip and scale so diagonals are one
{
eigen_strip_mults(cors, bitlength, mults+bitstart, strip);
#pragma omp parallel for private(j,k) schedule(static)
for(j=0;j<bitlength;j++)
{
for(k=0;k<bitlength;k++)
{
if(k!=j){cors[(size_t)j*bitlength+k]*=pow(cors[(size_t)j*bitlength+j]*cors[(size_t)k*bitlength+k],-.5);}
}
}
for(j=0;j<bitlength;j++){cors[(size_t)j*bitlength+j]=1;}
}

//how many predictors remain?
count=0;for(j=bitstart;j<bitend;j++){count+=(mults[j]!=-9999);}

if(count>0)	//some predictors remain
{
//compute exps for the different powers (do not scale, to ensure exp is 1 for first power)
for(k=0;k<num_pows;k++)
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999)
{
if(hwestand==1){exps[j+k*bitlength]=pow(centres[bitstart+j]*(1-centres[bitstart+j]/2),1+powers[k]);}
else{exps[j+k*bitlength]=pow(sqdevs[bitstart+j],1+powers[k]);}
}
else{exps[j+k*bitlength]=0;}
}
}

for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)
{
//extract values
count2=0;
for(k=j+1;k<bitend;k++)
{
if(mults[k]!=-9999)
{
writefloats[count2]=cors[(size_t)(j-bitstart)*bitlength+k-bitstart];
count2++;
}
}

if(mincor>0)	//truncate before printing
{
for(k=0;k<count2;k++)
{
if(pow(writefloats[k],2)<mincor){writefloats[k]=0;}
}
}

//convert to a short then save
for(k=0;k<count2;k++){writeshorts[k]=(int)(20000*writefloats[k]+20000.5);}
fwrite(writeshorts, sizeof(unsigned short), count2, output2);
}
}

//replace correlations with either their square (corrected for sample size) or zero
value=(double)(num_samples_use-1)/(num_samples_use-2);
value2=pow(num_samples_use-2,-1);
#pragma omp parallel for private(j,k) schedule(static)
for(j=0;j<bitlength;j++)
{
for(k=0;k<bitlength;k++)
{
if(mults[bitstart+j]!=-9999&&mults[bitstart+k]!=-9999){cors[(size_t)j*bitlength+k]=pow(cors[(size_t)j*bitlength+k],2)*value-value2;}
else{cors[(size_t)j*bitlength+k]=0;}
}
}

//compute rjksums
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_pows, &bitlength, &alpha, cors, &bitlength, exps, &bitlength, &beta, rjksums+bitstart, &data_length);

//ensure rjksums are not below exps
for(k=0;k<num_pows;k++)
{
for(j=0;j<bitlength;j++)
{
if(mults[bitstart+j]!=-9999&&rjksums[bitstart+j+k*data_length]<exps[j+k*bitlength]){rjksums[bitstart+j+k*data_length]=exps[j+k*bitlength];}
}
}

//get t(data) * randnorms * root(1/num_samples_use)
alpha=pow(1.0/num_samples_use,.5);beta=0.0;
dgemv_("T", &num_samples_use, &bitlength, &alpha, data, &num_samples_use, randnorms, &one, &beta, datarands+bitstart, &one);

//update number of pairs
scount+=count*(count-1)/2;
}	//end of count>0
}
else    //making fake summaries
{
//scale effect sizes (or set to zero for trivial predictors)
for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)
{
for(m=0;m<num_phenos;m++){effs[j+m*data_length]=effs[j+m*data_length]/mults[j];}
}
else
{
for(m=0;m<num_phenos;m++){effs[j+m*data_length]=0;}
}
}

//compute cor beta + t(data) gaussian
alpha=1.0;beta=0.0;
dgemm_("N", "N", &bitlength, &num_phenos, &bitlength, &alpha, cors, &bitlength, effs+bitstart, &data_length, &beta, effs2+bitstart, &data_length);
alpha=1.0;beta=1.0;
dgemm_("T", "N", &bitlength, &num_phenos, &num_samples_use, &alpha, data, &num_samples_use, gaussian, &num_samples_use, &beta, effs2+bitstart, &data_length);
}
}	//end of bit loop

if(xcount>5){printf("In total, %d predictors are trivial\n", xcount);}
if(ecount>5){printf("In total, %d predictors have MAF below %.6f or callrate below %.4f\n", ecount, minmaf, minobs);}
if(wcount>5){printf("In total, %d predictors are near duplicates\n", wcount);}
printf("\n");

if(strcmp(betafile,"blank")==0){fclose(output2);}

////////

if(strcmp(betafile,"blank")==0) //finish saving correlations
{
//work out number of predictors that remain and their indexes
total=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999){indexer[total]=j;total++;}
}
if(total==0){printf("Error, after filtering, no predictors remain\n\n");exit(1);}

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

//determine highld predictors
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999&&rjkaves[j]>2*med){highlds[j]=1;}
else{highlds[j]=0;}
}

//adjust blockstarts and blockends to account for missing predictors
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

////////

//open bin file and save counts, means, scalings, variances, rjksums, datarands and highlds, for predictors that remain
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

for(j=0;j<total;j++){fwrite(centres+indexer[j], sizeof(double), 1, output2);}
for(j=0;j<total;j++){fwrite(mults+indexer[j], sizeof(double), 1, output2);}
for(j=0;j<total;j++){fwrite(sqdevs+indexer[j], sizeof(double), 1, output2);}
for(k=0;k<num_pows;k++)
{
for(j=0;j<total;j++){fwrite(rjksums+indexer[j]+k*data_length, sizeof(double), 1, output2);}
}
for(j=0;j<total;j++){fwrite(datarands+indexer[j], sizeof(double), 1, output2);}
for(j=0;j<total;j++){fwrite(highlds+indexer[j], sizeof(int), 1, output2);}

//now add temp file to end of bin file
printf("Moving contents of %s.cors.temp to end of %s.cors.bin\n\n", outfile, outfile);

sprintf(filename,"%s.cors.temp", outfile);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_SET);

readshorts=malloc(sizeof(unsigned short)*1000000);

scount2=scount;
while(scount2>1000000)
{
if(fread(readshorts, sizeof(unsigned short), 1000000, input)!=1000000)
{printf("Error reading correlations from %s\n\n", filename);exit(1);}
fwrite(readshorts, sizeof(unsigned short), 1000000, output2);
scount2-=1000000;
}

if(scount2>0)
{
if(fread(readshorts, sizeof(unsigned short), scount2, input)!=scount2)
{printf("Error reading correlations from %s\n\n", filename);exit(1);}
fwrite(readshorts, sizeof(unsigned short), scount2, output2);
}

free(readshorts);

fclose(input);
fclose(output2);

//delete temp file
sprintf(cmd, "rm %s.cors.temp", outfile);
system(cmd);

////////

//save root
sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Datafile %s\n", datafile);
fprintf(output3,"Num_Samples %d\nNum_Predictors %d\n", num_samples, num_preds);
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, total);
fprintf(output3,"Num_Windows %d\n", bittotal);
fprintf(output3,"Num_Pairs %jd\n", scount);
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

if(num_parts>0) //save annotations
{
sprintf(filename6,"%s.cors.annotations", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6,"Predictor");
for(q=0;q<num_parts+addpart;q++){fprintf(output6," %s", catlabels[q]);}
fprintf(output6,"\n");
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output6, "%s", preds[j]);
for(q=0;q<num_parts+addpart;q++)
{
if(pindexes[q][j]!=0)
{
if(pweights[q][j]==1){fprintf(output6, " 1");}
else{fprintf(output6, " %.4f", pweights[q][j]);}
}
else{fprintf(output6, " 0");}
}
fprintf(output6,"\n");
}
}
fclose(output6);
}

printf("The correlations are saved in files with stem %s\n\n", outfile);
}
else    //save fake summaries
{
//variance of estimate is 1/num_inds
value=pow(num_inds,.5);

for(m=0;m<num_phenos;m++)
{
sprintf(filename3,"%s.pheno%d.summaries", outfile, m+1);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tA1\tA2\tZ\tn\tA1Freq\tP\n");

for(j=0;j<data_length;j++){fprintf(output3, "%s\t%c\t%c\t%.4f\t%d\t%.4f\t%.4e\n", preds[j], al1[j], al2[j], effs2[j+m*data_length]*value, num_inds, centres[j]/2, erfc(fabs(effs2[j+m*data_length]*value)*M_SQRT1_2));}
fclose(output3);
}

printf("The fake summary statistics are saved in %s.phenoX.summaries\n\n", outfile);
}

////////

if(strcmp(betafile,"blank")!=0){free(gaussian);free(effs);free(effs2);}

free(powers);
for(q=0;q<num_parts+1;q++){free(catlabels[q]);}free(catlabels);
for(q=0;q<num_parts+1;q++){free(pindexes[q]);free(pweights[q]);}free(pindexes);free(pweights);
free(blockstarts);free(blockends);free(indexer);
free(data);
free(Z);free(cors);free(exps);free(writeshorts);free(writefloats);
free(rjksums);free(rjkaves);free(rjktemp);
free(randnorms);free(datarands);
free(highlds);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


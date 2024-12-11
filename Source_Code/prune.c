/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Thinning - flag=0: 1-pass pruning, flag=1: pass1, flag=2: pass2

///////////////////////////

if(bitsize==-9999)	//get optimal bitsize
{
if(window_kb!=-9999)	//set to average number of predictors within window
{
scount=0;
k=1;
for(j=0;j<data_length;j++)
{
while(1)	//move along until out of range, tallying each data_length predictors encountered
{
if(k==data_length){break;}
if(chr[k]!=chr[j]||cmbp[k]-cmbp[j]>1000*window_kb){break;}
k++;
}
scount+=k-j-1;
}
bitsize=(int)((scount-1)/data_length)+1;
}
else	//set to window_length
{bitsize=window_length;}

if(bitsize<20){bitsize=20;}
if(bitsize>8000){bitsize=8000;}
if(flag!=2){printf("The bit-size will be set to %d (you can change this using \"--bit-size\")\n\n", bitsize);}
else{printf("The bit-size has now been set to %d\n\n", bitsize);}
}

step=(50000/bitsize);if(step<20){step=20;}
if(step>20){step=10*(step/10);}if(step>50){step=20*(step/20);}
if(step>100){step=50*(step/50);}if(step>300){step=100*(step/100);}
if(step>1000){step=500*(step/500);}

//work out bitmax
if(window_kb!=-9999)
{
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
}
else{bitmax=bitsize+window_length;}

//can now allocate variables

data_warn3(bitmax,num_samples_use);
data=malloc(sizeof(double)*num_samples_use*bitmax);

anal_warn(bitmax, bitmax);
cors=malloc(sizeof(double)*bitmax*bitmax);
if(mode==101){replace=malloc(sizeof(int)*num_preds);}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress file
if(mode==101){sprintf(filename,"%sthin.progress",folder);}
else{sprintf(filename,"%s.progress",outfile);}
if(flag==0||flag==1)	//only pass or first pass
{
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);
}

//ready for bit loop
bittotal=(data_length-1)/bitsize+1;
count=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
if(window_kb!=-9999)
{
for(bitend2=bitend;bitend2<data_length;bitend2++)
{
if(cmbp[bitend2]-cmbp[bitend-1]>1000*window_kb||chr[bitend2]!=chr[bitend-1]){break;}
}
}
else	//using window_length
{
bitend2=bitend+window_length;
if(bitend2>data_length){bitend2=data_length;}
}
bitlength=bitend2-bitstart;

if(flag==0)
{
if(bit%step==0)
{
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(bit==0)
{
printf("Thinning for Chunk %d of %d; stay tuned for updates ;)\n", bit+1, bittotal);
fprintf(output, "Thinning for Chunk %d of %d; stay tuned for updates ;)\n", bit+1, bittotal);
}
else
{
printf("Thinning for Chunk %d of %d; kept %d out of %d predictors (%.2f%%)\n", bit+1, bittotal, count, bitstart, 100.0*count/bitstart);
fprintf(output, "Thinning for Chunk %d of %d; kept %d out of %d predictors (%.2f%%)\n", bit+1, bittotal, count, bitstart, 100.0*count/bitstart);
}
fclose(output);
}
}
else
{
if(bit%step==0)
{
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(bit==0)
{
printf("Pass %d: Thinning for Chunk %d of %d; stay tuned for updates ;)\n", flag, bit+1, bittotal);
fprintf(output, "Pass %d: Thinning for Chunk %d of %d; stay tuned for updates ;)\n", flag, bit+1, bittotal);}
else
{
printf("Pass %d: Thinning for Chunk %d of %d; kept %d out of %d predictors (%.2f%%)\n", flag, bit+1, bittotal, count, bitstart, 100.0*count/bitstart);
fprintf(output, "Pass %d: Thinning for Chunk %d of %d; kept %d out of %d predictors (%.2f%%)\n", flag, bit+1, bittotal, count, bitstart, 100.0*count/bitstart);
}
fclose(output);
}
}

shuffle=0;
for(j=0;j<end-bitstart;j++)	//using values already in data, so shuffle back
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(bitstart-start+j)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, bitstart+shuffle, bitend2, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data+(size_t)shuffle*num_samples_use, centres+bitstart+shuffle, mults+bitstart+shuffle, sqdevs+bitstart+shuffle, rates+bitstart+shuffle, infos+bitstart+shuffle, num_samples_use, bitlength-shuffle, missingvalue, -1, 0, 0, NULL, 1);

for(j=bitstart+shuffle;j<bitend2;j++)	//trivial should get usedpreds-3
{
if(mults[j]==-9999){usedpreds[keeppreds_use[j]]=-3;}
}

if(wprune>0)	//get correlation
{
alpha=1.0/num_samples_use;beta=0.0;
dgemm_("T", "N", &bitlength, &bitlength, &num_samples_use, &alpha, data, &num_samples_use, data, &num_samples_use, &beta, cors, &bitlength);
}
else	//can set all correlations to one
{
for(j=0;j<bitlength;j++)
{
for(j2=0;j<bitlength;j++){cors[(size_t)j+j2]=1;}
}
}

for(j=bitstart;j<bitend;j++)
{
if(mults[j]!=-9999)
{
for(k=j+1;k<bitend2;k++)
{
if(chr[k]!=chr[j]){break;}
if(window_kb!=-9999&&cmbp[k]-cmbp[j]>1000*window_kb){break;}
if(window_length!=-9999&&k-j>window_length){break;}

if(mults[k]!=-9999)	//check correlation
{
if(pow(cors[(size_t)(j-bitstart)*bitlength+(k-bitstart)],2)>=wprune)	//will lose i, keep i2
{
i=j;i2=k;
if(j+k%2==1){i=k;i2=j;}
if(pvalues[j]<pvalues[k]){i=k;i2=j;}
if(pvalues[k]<pvalues[j]){i=j;i2=k;}
mults[i]=-9999;
usedpreds[keeppreds_use[i]]=keeppreds_use[i2];
}
if(mults[j]==-9999){break;}
}}
}
if(mults[j]!=-9999){count++;}
}
start=bitstart;if(bitend2>end){end=bitend2;}
}	//end of bit loop
printf("\n");

////////

//squeeze down (not usedpreds)
count=0;
for(j=0;j<data_length;j++)
{
if(mults[j]!=-9999)
{
if(count!=j)
{
chr[count]=chr[j];cmbp[count]=cmbp[j];al1[count]=al1[j];al2[count]=al2[j];
free(preds[count]);copy_string(preds,count,preds[j]);
free(along1[count]);copy_string(along1,count,along1[j]);free(along2[count]);copy_string(along2,count,along2[j]);
keeppreds_use[count]=keeppreds_use[j];
centres[count]=centres[j];mults[count]=mults[j];sqdevs[count]=sqdevs[j];rates[count]=rates[j];infos[count]=infos[j];
weights[count]=weights[j];pvalues[count]=pvalues[j];
}
count++;
}}
for(j=count;j<data_length;j++){free(preds[j]);free(along1[j]);free(along2[j]);}
data_length=count;

////////

if(flag==0||flag==2)	//final loop, so save - usedpreds -3/-2/-1/>=0 - trivial/filtered/in/pruned
{
if(mode==101)	//work out the replacements
{
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==-3){replace[j]=-1;}	//trivial
if(usedpreds[j]==-2){replace[j]=-2;}	//filtered out
if(usedpreds[j]==-1){replace[j]=j;}	//remained
if(usedpreds[j]>=0)	//has a replacement
{
found=usedpreds[j];
while(usedpreds[found]>=0){found=usedpreds[found];}
replace[j]=found;
}
}
}

if(mode==101){sprintf(filename2,"%sthin.in",folder);}
else{sprintf(filename2,"%s.in",outfile);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(mode==101){sprintf(filename3,"%sthin.out",folder);}
else{sprintf(filename3,"%s.out",outfile);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

if(mode!=110)
{
if(mode==101){sprintf(filename4,"%sthin.trivial",folder);}
else{sprintf(filename4,"%s.trivial",outfile);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
}

count=0;count2=0;count3=0;
for(j=0;j<num_preds;j++)
{
if(usedpreds[j]==-1){fprintf(output2,"%s\n", allpreds[j]);count++;}
if(usedpreds[j]>=0){fprintf(output3,"%s\n", allpreds[j]);count2++;}
if(usedpreds[j]==-3){fprintf(output4,"%s\n", allpreds[j]);count3++;}
}
fclose(output2);
fclose(output3);
if(mode!=110){fclose(output4);}

if(mode==101)	//print dup indexes
{
sprintf(filename5,"%sthin.dups",folder);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Replacement_Index\n"); 
for(j=0;j<num_preds;j++)
{
if(replace[j]<0){fprintf(output5, "-1\n");}	//trivial or thinned out
else{fprintf(output5, "%d\n", replace[j]+1);}	//write replacement (could be itself)
}
fclose(output5);
}

if(mode!=110){printf("Thinning complete: %d predictors kept (saved in %s), %d lost (%s) and %d trivial (%s)\n", count, filename2, count2, filename3, count3, filename4);}
else{printf("Thinning complete: %d predictors kept (saved in %s) and %d lost (%s)\n", count, filename2, count2, filename3);}

if(mode==101)
{printf("If the subsequent cutting of weights fails, or if you wish to re-cut (with different parameters) you can avoid having to thin again by adding \"--no-thin DONE\"\n");}
printf("\n");
}

free(data);
free(cors);
if(mode==101){free(replace);}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


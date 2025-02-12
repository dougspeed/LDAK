/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Code for quickly merging bim files

///////////////////////////

//work out from which file each predictor comes
indexer=malloc(sizeof(int)*num_preds_use);
indexer2=malloc(sizeof(int)*num_preds_use);
for(k=0;k<num_files;k++)
{
for(j=0;j<XNuse[k];j++){indexer[Xkp2[k][j]]=k;indexer2[Xkp2[k][j]]=Xkp[k][j];}
}

if(bitsize>num_preds_use){bitsize=num_preds_use;}

//some allocations

data_warn2(bitsize,num_samples_use/32);

rowlength=(num_samples_use-1)/4+1;
rowchars=malloc(sizeof(unsigned char)*rowlength*bitsize);

centres=malloc(sizeof(double)*num_preds_use);
mults=malloc(sizeof(double)*num_preds_use);
sqdevs=malloc(sizeof(double)*num_preds_use);
rates=malloc(sizeof(double)*num_preds_use);
infos=malloc(sizeof(double)*num_preds_use);

nums=malloc(sizeof(int*)*6);
for(k=0;k<6;k++)
{
nums[k]=malloc(sizeof(int)*num_preds_use);
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}

Xinput=malloc(sizeof(FILE*)*num_files);
Xcurrent=malloc(sizeof(int)*num_files);

//set nums to zero
for(k=0;k<6;k++)
{
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}

//open all bim files, check their length and first three digits, and move to first predictor
for(k=0;k<num_files;k++)
{
if((Xinput[k]=fopen(datastems[k],"rb"))==NULL)
{printf("Error opening %s\n\n",datastems[k]);exit(1);}

fseeko(Xinput[k], 0, SEEK_END);
if(ftello(Xinput[k])!=(off_t)sizeof(unsigned char)*rowlength*XNall[k]+sizeof(unsigned char)*3)
{printf("Error reading %s; should have size %jd (%d ind x %d predictors), but instead has size %jd\n\n", datastems[k], (off_t)sizeof(unsigned char)*rowlength*XNall[k]+sizeof(unsigned char)*3, num_samples_use, XNall[k], ftello(Xinput[k]));exit(1);}

fseeko(Xinput[k], 0, SEEK_SET);
if(fread(startchars, sizeof(unsigned char), 3, Xinput[k])!=3)
{printf("Error reading first three values of %s\n\n", datastems[k]);exit(1);}
if(startchars[0]!=108||startchars[1]!=27)
{printf("Error reading %s; does not appear to be in binary PLINK format\n\n", datastems[k]);exit(1);}
if(startchars[2]!=1)
{printf("Error reading %s; can only read in SNP-major mode\n\n", datastems[k]);exit(1);}

Xcurrent[k]=0;
}

//open output files
sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.stats",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor\tA1\tA2\tA1_Mean\tMAF\tCall_Rate\tInfo\n");

sprintf(filename3,"%s.bed", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
onechar=108;fwrite(&onechar, sizeof(unsigned char), 1, output3);
onechar=27;fwrite(&onechar, sizeof(unsigned char), 1, output3);
onechar=1;fwrite(&onechar, sizeof(unsigned char), 1, output3);

////////

//ready for bit loop
bittotal=(num_preds_use-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_preds_use){bitend=num_preds_use;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
printf("Making data for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Making data for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
}

for(j=0;j<bitlength;j++)
{
k=indexer[bitstart+j];j2=indexer2[bitstart+j];
if(j2!=Xcurrent[k])
{
if(fseeko(Xinput[k], (off_t)sizeof(unsigned char)*rowlength*j2+sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d\n\n", datastems[k], j2+1);exit(1);}
}

if(fread(rowchars+(size_t)j*rowlength, sizeof(unsigned char), rowlength, Xinput[k])!=rowlength)
{printf("Error reading values for Predictor %d from %s\n\n", j2+1, datastems[k]);exit(1);}
Xcurrent[k]=j2+1;
}	//end of j loop

//calculate statistics
#pragma omp parallel for private(j, i, c0, c1, c2, rowchars2, indcount, mean, var, value) schedule(static)
for(j=0;j<bitlength;j++)
{
//tally numbers of zeros, ones and twos
i=0;
c0=0;c1=0;c2=0;
rowchars2=rowchars+(size_t)j*rowlength;

//process blocks of four samples
while(i+4<=num_samples_use)
{
c0+=bedzeros[rowchars2[0]];
c1+=bedones[rowchars2[0]];
c2+=bedtwos[rowchars2[0]];
i+=4;
rowchars2++;
}

//now any remaining samples (will be at most three, so no need to update rowchars2)
for(i2=i;i2<num_samples_use;i2++)
{
switch ((rowchars2[0]>>(2*(i2%4))) & 3)
{
case 3: c0++;break;
case 2: c1++;break;
case 0: c2++;
}
}

//now compute stats (assuming power=0)
indcount=c0+c1+c2;
if(indcount>0){mean=(double)(c1+2*c2)/indcount;var=(double)(c1+4*c2)/indcount-pow(mean,2);}
else{mean=0;var=0;}
centres[bitstart+j]=mean;
if(var>0){mults[bitstart+j]=1;}
else{mults[bitstart+j]=-9999;}
sqdevs[bitstart+j]=var*indcount/num_samples_use;
rates[bitstart+j]=(double)indcount/num_samples_use;
value=centres[bitstart+j]*(1-centres[bitstart+j]/2);
if(value>0){infos[bitstart+j]=sqdevs[bitstart+j]/value;}
else{infos[bitstart+j]=0;}
}

if(passqc==1)	//perform QC
{
for(j=0;j<bitlength;j++)
{
maf=centres[bitstart+j]/2+(centres[bitstart+j]>1)*(1-centres[bitstart+j]);
if(minmaf!=-9999&&maf<minmaf){nums[1][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(maxmaf!=-9999&&maf>maxmaf){nums[2][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(minvar!=-9999&&sqdevs[bitstart+j]<minvar){nums[3][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(minobs!=-9999&&rates[bitstart+j]<minobs){nums[4][bitstart+j]=1;nums[0][bitstart+j]=1;}
if(mininfo!=-9999&&infos[bitstart+j]<mininfo){nums[5][bitstart+j]=1;nums[0][bitstart+j]=1;}
}
}

//save predictor stats (allowing for QC)
for(j=0;j<bitlength;j++)
{
if(nums[0][bitstart+j]==0)
{
fprintf(output2, "%s\t%s\t%s\t", preds[bitstart+j], along1[bitstart+j], along2[bitstart+j]);
if(rates[bitstart+j]>0)
{
mean=centres[bitstart+j];
maf=mean/2+(mean>1)*(1-mean);
fprintf(output2,"%.6f\t",mean);
if(nonsnp==0){fprintf(output2, "%.6f\t", maf);}
else{fprintf(output2, "NA\t");}
}
else{fprintf(output2, "NA\tNA\t");}
if(rates[bitstart+j]==0||rates[bitstart+j]==1){fprintf(output2, "%d\t", (int)rates[bitstart+j]);}
else{fprintf(output2, "%.6f\t", rates[bitstart+j]);}
if(genprobs<2){fprintf(output2, "NA\n");}
else{fprintf(output2, "\t%.4f\n", infos[bitstart+j]);}
}
}

//save data (allowing for QC)
for(j=0;j<bitlength;j++)
{
if(nums[0][bitstart+j]==0){fwrite(rowchars+(size_t)j*rowlength, sizeof(unsigned char), rowlength, output3);}
}
}	//end of bit loop
printf("\n");

for(k=0;k<num_files;k++){fclose(Xinput[k]);}
fclose(output2);
fclose(output3);

//see how many passed qc
count=0;for(j=0;j<num_preds_use;j++){count+=(nums[0][j]==0);}
if(count==0)
{
sprintf(cmd,"rm %s", filename3);
system(cmd);
printf("Error, no predictors passed quality control, so no data were saved\n\n");

if((output=fopen(filename,"w"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fprintf(output, "Error, all %d chunks were read, but no predictors passed quality control (so no files were written)\n\n", bittotal);
fclose(output);
exit(1);
}

//save bim
sprintf(filename4,"%s.bim", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<num_preds_use;j++)
{
if(nums[0][j]==0)
{
fprintf(output4, "%d %s ", chr[j], preds[j]);
if(cm[j]==0){fprintf(output4, "0 ");}
else{fprintf(output4, "%.6f ", cm[j]);}
fprintf(output4, "%ld %s %s\n", (long int)bp[j], along1[j], along2[j]);
}
}
fclose(output4);

//and fam
sprintf(filename5,"%s.fam", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output5, "%s %s %s %s %s %s\n", ids1[i], ids2[i], pid[i], mid[i], schar[i], pchar[i]);}
fclose(output5);

if(passqc==1)	//save qc info
{
sprintf(filename6,"%s.qc", outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
fprintf(output6, "Predictor\tFail_Any\tFail_Min_MAF\tFail_Max_MAF\tFail_Min_Var\tFail_Min_Obs\tFail_Min_Info\n");
for(j=0;j<num_preds_use;j++)
{
fprintf(output6, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", preds[j], nums[0][j], nums[1][j], nums[2][j], nums[3][j], nums[4][j], nums[5][j]);}
fclose(output6);
}

free(indexer);free(indexer2);
free(rowchars);
free(centres);free(mults);free(sqdevs);free(rates);free(infos);
for(k=0;k<6;k++){free(nums[k]);}free(nums);
free(Xinput);free(Xcurrent);

printf("Data for %d samples and %d predictors saved in %s, %s and %s, with statistics in %s", num_samples_use, count, filename3, filename4, filename5, filename2);
if(passqc==1){printf("; %s indicates which predictors failed quality control",filename6);}
printf("\n\n");

///////////////////////////


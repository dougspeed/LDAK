/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Get predictor and sample statistics (maybe because making data)
//for consistency with merge.c, use num_preds_use and keeppreds instead of data_length and keeppreds_use

///////////////////////////

//allocate variables

if(genprobs<2){data_warn2(bitsize,num_samples_use);}
else{data_warn2(bitsize,2*num_samples_use);}

data=malloc(sizeof(double)*num_samples_use*bitsize);

if(genprobs>=2)
{
ps=malloc(sizeof(float*)*2);
ps[0]=malloc(sizeof(float)*num_samples_use*bitsize);
ps[1]=malloc(sizeof(float)*num_samples_use*bitsize);
}
else{ps=NULL;}

nums=malloc(sizeof(int*)*6);
for(k=0;k<6;k++)
{
nums[k]=malloc(sizeof(int)*num_preds_use);
}

presents=malloc(sizeof(int)*num_samples_use);
hets=malloc(sizeof(int)*num_samples_use);

pid=malloc(sizeof(char*)*num_samples_use);
mid=malloc(sizeof(char*)*num_samples_use);
schar=malloc(sizeof(char*)*num_samples_use);
pchar=malloc(sizeof(char*)*num_samples_use);

//set nums, presents and hets to zero
for(k=0;k<6;k++)
{
for(j=0;j<num_preds_use;j++){nums[k][j]=0;}
}
for(i=0;i<num_samples_use;i++){presents[i]=0;hets[i]=0;}

//see if possible to read parental ids
flag=0;
if(famhead==0)
{
count=countcols(famfile);
if(count<6)
{printf("Warning, %s has only %d columns; famfiles should generally have (at least) 6 columns\n\n", famfile, count);}
else
{
read_strings(famfile, pid, num_samples_use, keepsamps, 3, 0);
read_strings(famfile, mid, num_samples_use, keepsamps, 4, 0);
read_strings(famfile, schar, num_samples_use, keepsamps, 5, 0);
read_strings(famfile, pchar, num_samples_use, keepsamps, 6, 0);
flag=1;
}
}

if(flag==0)
{
for(i=0;i<num_samples_use;i++)
{
pid[i]=malloc(sizeof(char)*2);strcpy(pid[i],"0");
mid[i]=malloc(sizeof(char)*2);strcpy(mid[i],"0");
schar[i]=malloc(sizeof(char)*2);strcpy(schar[i],"0");
pchar[i]=malloc(sizeof(char)*3);strcpy(pchar[i],"NA");
}
}

//prepare for reading data

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.stats",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor\tA1\tA2\tA1_Mean\tMAF\tCall_Rate\tInfo\n");

if(mode==181)	//bed
{
sprintf(filename3,"%s.bed", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
onechar=108;fwrite(&onechar, sizeof(unsigned char), 1, output3);
onechar=27;fwrite(&onechar, sizeof(unsigned char), 1, output3);
onechar=1;fwrite(&onechar, sizeof(unsigned char), 1, output3);
}
if(mode==182)	//sp
{
sprintf(filename3,"%s.sp", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}
if(mode==183)	//sped
{
sprintf(filename3,"%s.sped", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}
if(mode==184)	//speed
{
sprintf(filename3,"%s.speed", outfile);
if((output3=fopen(filename3,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}
if(mode==185)	//gen
{
sprintf(filename3,"%s.gen", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
}

////////

bittotal=(num_preds_use-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>num_preds_use){bitend=num_preds_use;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
if(mode==171){printf("Calculating statistics for Chunk %d of %d\n", bit+1, bittotal);}
else{printf("Making data for Chunk %d of %d\n", bit+1, bittotal);}
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
if(mode==171){fprintf(output, "Calculating statistics for Chunk %d of %d\n", bit+1, bittotal);}
else{fprintf(output, "Making data for Chunk %d of %d\n", bit+1, bittotal);}
fclose(output);
}

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read data
current=read_data_fly(datafile, dtype, data, ps, num_samples_use, keepsamps, bitstart, bitend, keeppreds, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, threshold, minprob, nonsnp, maxthreads);

if(mode!=171&&encoding!=1)	//change coding - and maybe also alleles
{change_coding(data, al1+bitstart, al2+bitstart, num_samples_use, bitlength, encoding, missingvalue);}

//compute statistics
stand_data(data, centres+bitstart, mults+bitstart, sqdevs+bitstart, rates+bitstart, infos+bitstart, num_samples_use, bitlength, missingvalue, 0, 0, -9999, NULL, 0);

if(minprob==0)	//info scores are not valid
{
for(j=0;j<bitlength;j++){infos[bitstart+j]=1;}
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
}	//end of j loop

//tally individual stats (allowing for QC)
for(j=0;j<bitlength;j++)
{
if(nums[0][bitstart+j]==0)
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{presents[i]++;hets[i]+=(data[(size_t)j*num_samples_use+i]==1);}
}
}
}

if(mode!=171)	//save data (allowing for QC)
{
for(j=0;j<bitlength;j++)
{
if(nums[0][bitstart+j]==0)
{
#include "savedata.c"
}
}
}
}	//end of bit loop
printf("\n");

fclose(output2);
if(mode!=171){fclose(output3);}

if(mode!=171)	//write bimfile, samplefile/famfile and maybe qc
{
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

if(mode==185)	//sample
{
sprintf(filename5,"%s.sample", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "ID_1 ID_2 missing sex\n");
fprintf(output5, "0 0 0 D\n");
for(i=0;i<num_samples_use;i++)
{fprintf(output5, "%s %s 0 %s\n", ids1[i], ids2[i], schar[i]);}
fclose(output5);
}
else
{
sprintf(filename5,"%s.fam", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output5, "%s %s %s %s %s %s\n", ids1[i], ids2[i], pid[i], mid[i], schar[i], pchar[i]);}
fclose(output5);
}

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
}

//save individual statistics - will only print heterozygosity if have 0, 1, 2 values
sprintf(filename7,"%s.missing",outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
fprintf(output7, "FID\tIID\tMissing_Rate\tHeterozygosity_Rate\n");
for(i=0;i<num_samples_use;i++)
{
fprintf(output7, "%s\t%s\t%.6f\t", ids1[i], ids2[i], 1-(double)presents[i]/count);
if(nonsnp==0&&(dtype==1||threshold!=-9999||minprob!=-9999)){fprintf(output7, "%.6f\n", (double)hets[i]/presents[i]);}
else{fprintf(output7, "NA\n");}
}
fclose(output7);

if(mode==171){printf("Statistics saved in %s and %s\n\n", filename2, filename7);}
else
{
printf("Data for %d samples and %d predictors saved in %s, %s and %s, with statistics in %s and %s", num_samples_use, count, filename3, filename4, filename5, filename2, filename7);
if(passqc==1){printf("; %s indicates which predictors failed quality control",filename6);}
printf("\n\n");
}

free(data);
if(genprobs>=2){free(ps[0]);free(ps[1]);free(ps);}
for(k=0;k<6;k++){free(nums[k]);}free(nums);
free(presents);free(hets);
for(i=0;i<num_samples_use;i++){free(pid[i]);free(mid[i]);free(schar[i]);free(pchar[i]);}
free(pid);free(mid);free(schar);free(pchar);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


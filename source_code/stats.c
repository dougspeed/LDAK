/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Get predictor and sample statistics

///////////////////////////

//allocate variables

data_warn2(bitsize,2*num_samples_use);

data=malloc(sizeof(double)*num_samples_use*bitsize);
ps=malloc(sizeof(float*)*2);
ps[0]=malloc(sizeof(float)*num_samples_use*bitsize);
ps[1]=malloc(sizeof(float)*num_samples_use*bitsize);

presents=malloc(sizeof(double)*num_samples_use);
hets=malloc(sizeof(double)*num_samples_use);
for(i=0;i<num_samples_use;i++){presents[i]=0;hets[i]=0;}

if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

sprintf(filename2,"%s.stats",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor\tA1\tA2\tA1_Mean\tMAF\tCall_Rate\tInfo");
for(s=0;s<num_subs;s++){fprintf(output2,"\tA1_Mean_S%d\tMAF_S%d\tCall_Rate_S%d\tInfo_S%d", s+1, s+1, s+1, s+1);}
fprintf(output2, "\n");

bittotal=(data_length-1)/bitsize+1;
for(bit=0;bit<bittotal;bit++)
{
bitstart=bit*bitsize;
bitend=(bit+1)*bitsize;
if(bitend>data_length){bitend=data_length;}
bitlength=bitend-bitstart;

if(bit%10==0)
{
printf("Calculating statistics for Chunk %d of %d\n", bit+1, bittotal);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Calculating statistics for Chunk %d of %d\n", bit+1, bittotal);
fclose(output);
}

fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

current=read_data_fly(datafile, dtype, data, ps, num_samples_use, keepsamps, bitstart, bitend, keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);

for(j=0;j<bitlength;j++)
{
//first all samples
sum=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)j*num_samples_use+i]!=missingvalue)
{
sum+=data[(size_t)j*num_samples_use+i];
indcount++;
presents[i]++;
hets[i]+=fabs(data[(size_t)j*num_samples_use+i]-1);
}
}
if(indcount>0){mean=sum/indcount;}
if(genprobs>1){value=compute_info(data+(size_t)j*num_samples_use, ps[0]+j*num_samples_use, ps[1]+j*num_samples_use, NULL, num_samples_use, missingvalue);}

fprintf(output2, "%s\t%c\t%c", preds[bitstart+j], al1[bitstart+j], al2[bitstart+j]);
if(indcount>0)
{
fprintf(output2,"\t%.6f",mean);
if(nonsnp==0){fprintf(output2, "\t%.6f", mean/2+(mean>1)*(1-mean));}
else{fprintf(output2, "\t-1");}
}
else{fprintf(output2, "\t-1\t-1");}
if(indcount==0||indcount==num_samples_use){fprintf(output2, "\t%d", (indcount==num_samples_use));}
else{fprintf(output2, "\t%.6f", (double)indcount/num_samples_use);}
if(genprobs>1){fprintf(output2, "\t%.4f", value);}
else{fprintf(output2, "\t-1");}

for(s=0;s<num_subs;s++)
{
sum=0;indcount=0;
for(i=0;i<subindex[s][0];i++)
{
i2=subindex[s][1+i];
if(data[(size_t)j*num_samples_use+i2]!=missingvalue){sum+=data[(size_t)j*num_samples_use+i2];indcount++;}
}
if(indcount>0){mean=sum/indcount;}
if(genprobs>1){value=compute_info(data+(size_t)j*num_samples_use, ps[0]+j*num_samples_use, ps[1]+j*num_samples_use, subindex[s]+1, num_samples_use, missingvalue);}

if(indcount>0)
{
fprintf(output2,"\t%.6f",mean);
if(nonsnp==0){fprintf(output2, "\t%.6f", mean/2+(mean>1)*(1-mean));}
else{fprintf(output2, "\t-1");}
}
else{fprintf(output2, "\t-1\t-1");}
if(indcount==0||indcount==subindex[s][0]){fprintf(output2, "\t%d", (indcount==subindex[s][0]));}
else{fprintf(output2, "\t%.6f", (double)indcount/subindex[s][0]);}
if(genprobs>1){fprintf(output2, "\t%.4f", value);}
else{fprintf(output2, "\t-1");}
}	//end of s loop

fprintf(output2, "\n");
}	//end of j loop
}	//end of bit loop
fclose(output2);

sprintf(filename3,"%s.missing",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "FID\tIID\tMissing_Rate\tHeterozygosity_Rate\n");
for(i=0;i<num_samples_use;i++)
{
fprintf(output3, "%s\t%s\t%.6f\t", ids1[i], ids2[i], 1-(double)presents[i]/data_length);
if(nonsnp==0){fprintf(output3, "%.6f\n", 1-hets[i]/presents[i]);}
else{fprintf(output3, "-1\n");}
}
fclose(output3);

printf("\nStatistics saved in %s and %s\n\n", filename2, filename3);

free(data);free(ps[0]);free(ps[1]);free(ps);
if(binary==0){gzclose(datainputgz);}

///////////////////////////


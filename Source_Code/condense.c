/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Condensing predictors

///////////////////////////

//get an upper limit on number of genes/chunks - copied from ldak.c
if(strcmp(genefile,"blank")!=0)	//using genefile - easy
{num_genes=countrows(genefile)-check_head(genefile,"Gene","Name",1);}

if(chunks!=-9999)	//using weights
{
sum=weights[0];count=1;
for(j=1;j<data_length;j++)
{
sum+=weights[j];
if(chr[j]!=chr[j-1]){count++;}
}
num_genes=sum/chunks*(1+overlap)+2*count;
}

if(chunksbp==1)	//each basepair observed becomes a gene - easy
{num_genes=data_length;}

if(chunksbp>1)
{
sum=0;count=1;
for(j=1;j<data_length;j++)
{
if(chr[j]!=chr[j-1]){sum+=bp[j-1];count++;}
}
sum+=bp[data_length-1];

if(sum/chunksbp*(1+overlap)+2*count>pow(2,31))
{printf("Error, there are too many chunks; to continue, you must increase chunksbp (currently %d)\n\n", chunksbp);exit(1);}
num_genes=sum/chunksbp*(1+overlap)+2*count;
}

gchr=malloc(sizeof(int)*num_genes);
gbp1=malloc(sizeof(double)*num_genes);
gbp2=malloc(sizeof(double)*num_genes);
gnames=malloc(sizeof(char*)*num_genes);
gstrand=malloc(sizeof(int)*num_genes);
gstarts=malloc(sizeof(int)*num_genes);
gends=malloc(sizeof(int)*num_genes);
for(g=0;g<num_genes;g++)
{
if(chunks!=-9999||chunksbp!=-9999){gnames[g]=malloc(sizeof(char)*100);}
gstarts[g]=-9999;gends[g]=-9999;
}

cut_genes(gchr, gbp1, gbp2, gnames, gstrand, gstarts, gends, data_length, chr, preds, cmbp, weights, genefile, chunks, chunksbp, up_buffer, down_buffer, minweight, overlap, pvafile, pvalues, -9999, binary, 1, NULL, datafile, bimfile, extract, 0);

//count number of genes and set genemax to size of longest gene
count=0;genemax=0;
for(g=0;g<num_genes;g++)
{
if(gstarts[g]!=-9999)
{
count++;
if(gends[g]-gstarts[g]>genemax){genemax=gends[g]-gstarts[g];}
}
}
if(count==0){printf("Doug Error, there are no genes\n\n");exit(1);}
printf("The longest gene/chunk contains %d predictors\n\n", genemax);

//allocate variables
data_warn(num_samples_use, genemax);
data=malloc(sizeof(double)*num_samples_use*genemax);
data2=malloc(sizeof(double)*num_samples_use);

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

if(mode==186)	//bed
{
sprintf(filename2,"%s.bed", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
onechar=108;fwrite(&onechar, sizeof(unsigned char), 1, output2);
onechar=27;fwrite(&onechar, sizeof(unsigned char), 1, output2);
onechar=1;fwrite(&onechar, sizeof(unsigned char), 1, output2);
}
if(mode==187)	//sp
{
sprintf(filename2,"%s.sp", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
}
if(mode==188)	//sped
{
sprintf(filename2,"%s.sped", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
}
if(mode==189)	//speed
{
sprintf(filename2,"%s.speed", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
}

sprintf(filename3,"%s.bim",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

sprintf(filename4,"%s.details",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Gene_Name Length Weight Minimum Maximum Mean Gene_Chr Gene_Start Gene_End Start_Predictor_BP End_Predictor_BP Min_Pvalue Min_Pvalue_Bonferroni\n");

////////

found=0;xcount=0;ecount=0;
for(g=0;g<num_genes;g++)
{
total=gends[g]-gstarts[g];

if(g%500==0)
{
printf("Condensing data for Gene/Chunk %d of %d\n", g+1, num_genes);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Condensing data for Gene/Chunk %d of %d\n", g+1, num_genes);
fclose(output);
}

if(gstarts[g]!=-9999)
{
shuffle=0;	//this is how many predictors we already have data for (can be more than total)
for(j=gstarts[g];j<end;j++)
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(j-start)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, gstarts[g]+shuffle, gends[g], keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bgen_indexes, missingvalue, -9999, -9999, nonsnp, maxthreads);

//get mults, centre, do qg, set missing to mean and maybe flip
wcount=0;
for(j=gstarts[g]+shuffle;j<gends[g];j++)
{
sum=0;sumsq=0;indcount=0;
for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)(j-gstarts[g])*num_samples_use+i];
if(value!=missingvalue){sum+=value;sumsq+=pow(value,2);indcount++;}
}

if(indcount>0){mean=sum/indcount;var=sumsq/indcount-pow(mean,2);}
else{mean=0;var=0;}
maf=mean/2+(mean>1)*(1-mean);	//will not use when indcount=0 (or if nonsnp=1)

centres[j]=mean;
if(var>0)
{
if(hwestand==1){mults[j]=pow(mean*(1-mean/2),power/2);}
else{mults[j]=pow(var*indcount/num_samples_use,power/2);}
}
else	//trivial
{
mults[j]=-9999;
if(wcount<5){printf("Warning, Predictor %s is trivial (takes at most one non-missing value) and will be ignored\n", preds[j]);}
wcount++;
}
sqdevs[j]=var*indcount/num_samples_use;

//predictors with indcount=0 will already have mults -9999
if(minmaf!=-9999&&maf<minmaf){mults[j]=-9999;}
if(maxmaf!=-9999&&maf>maxmaf){mults[j]=-9999;}
if(minvar!=-9999&&var<minvar){mults[j]=-9999;}
if(minobs!=-9999&&indcount<minobs*num_samples_use){mults[j]=-9999;}

if(mults[j]!=-9999&&indcount<num_samples_use)	//set missing values to mean
{
for(i=0;i<num_samples_use;i++)
{
if(data[(size_t)(j-gstarts[g])*num_samples_use+i]==missingvalue)
{data[(size_t)(j-gstarts[g])*num_samples_use+i]=centres[j];}
}}

if(mults[j]!=-9999&&centres[j]>1&&useminor==1)	//flip data
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)(j-gstarts[g])*num_samples_use+i]=2-data[(size_t)(j-gstarts[g])*num_samples_use+i];}
}
}	//end of j loop

if(wcount>5){printf("In total, %d predictors are trivial\n", wcount);}
if(wcount>0){printf("\n");}

count2=0;for(j=gstarts[g];j<gends[g];j++){count2+=(mults[j]!=-9999);}
if(count2==0)
{
if(xcount<10)
{
printf("Warning, all %d predictors in Gene/Chunk %d are trivial", gends[g]-gstarts[g], g+1);
if(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999){printf(" or failed QC");}
printf("\n");
}
xcount++;
}
else	//will use this gene/chunk
{
for(i=0;i<num_samples_use;i++){data2[i]=0;}

weightsum=0;minpvalue=1;count3=0;
for(j=gstarts[g];j<gends[g];j++)
{
if(mults[j]!=-9999&&weights[j]>0)
{
value=mults[j]*pow(weights[j],.5);
for(i=0;i<num_samples_use;i++){data2[i]+=data[(size_t)(j-gstarts[g])*num_samples_use+i]*value;}
weightsum+=weights[j];
}

if(pvalues[j]!=2){count3++;}
if(pvalues[j]<minpvalue){minpvalue=pvalues[j];}
}
minpvalue2=1-pow(1-minpvalue,count3);
if(minpvalue<1e-10){minpvalue2=minpvalue*count3;}

if(mode==186)	//convert to 0, 1, 2s, warning if values above 2
{
flag=0;
for(i=0;i<num_samples_use;i++)
{
data2[i]=round(data2[i]);
if(data2[i]>2){flag=1;data2[i]=2;}
}
if(flag==1)
{
if(ecount<10){printf("Warning, Gene/Chunk %d has values above 2\n\n", g+1);}
ecount++;
}
}

//get min and max, and sum
minfloat=data2[0];maxfloat=data2[0];sum=0;
for(i=0;i<num_samples_use;i++)
{
if(data2[i]<minfloat){minfloat=data2[i];}
if(data2[i]>maxfloat){maxfloat=data2[i];}
sum+=data2[i];
}

if(mode==189)	//print headers
{
if(speedlong==0){writefloat=1;}
else{writefloat=2;}
fwrite(&writefloat, sizeof(float), 1, output2);
fwrite(&minfloat, sizeof(float), 1, output2);
fwrite(&maxfloat, sizeof(float), 1, output2);
writefloat=0;
for(k=0;k<13;k++){fwrite(&writefloat, sizeof(float), 1, output2);}
}

for(i=0;i<num_samples_use;i++)
{
value=data2[i];	//can not be missing	

if(mode==186)
{
if(i%4==0){onechar=0;}
if(value==0){onechar+=(3<<(2*(i%4)));}
if(value==1){onechar+=(2<<(2*(i%4)));}
//if(value==2){onechar+=(0<<(2*(i%4)));}
//if(value==missingvalue){onechar+=(1<<(2*(i%4)));}
if(i%4==3||i==num_samples_use-1){fwrite(&onechar, sizeof(unsigned char), 1, output2);}
}
if(mode==187)
{
if(fabs(value-round(value))<0.00005){fprintf(output2,"%d ", (int)round(value));}
else{fprintf(output2,"%.4f ", value);}
}
if(mode==188)
{
writefloat=(float)value;
fwrite(&writefloat, sizeof(float), 1, output2);
}
if(mode==189)
{
if(speedlong==0)	//use 256
{
onechar=round((value-minfloat)/(maxfloat-minfloat)*254);
fwrite(&onechar, sizeof(unsigned char), 1, output2);
}
else	//use 65536
{
oneshort=round((value-minfloat)/(maxfloat-minfloat)*65534);
fwrite(&oneshort, sizeof(unsigned short), 1, output3);
}
}
}	//end of i loop
if(mode==187){fprintf(output2,"\n");}

if(chunks!=-9999||chunksbp!=-9999)
{
sprintf(gnames[g],"Chunk_%d",found+1);
fprintf(output3, "%d %s %.0f %.0f A B\n", gchr[g], gnames[g], gbp2[g], gbp1[g]);
fprintf(output4, "%s %d %.4f %.6f %.6f %.6f %d %.0f %.0f %.0f %.0f ", gnames[g], gends[g]-gstarts[g], weightsum, minfloat, maxfloat, sum/num_samples_use, gchr[g], gbp1[g], gbp2[g], cmbp[gstarts[g]], cmbp[gends[g]-1]);
}
else
{
fprintf(output3, "%d %s_%d %.0f %.0f A B\n", gchr[g], gnames[g], found+1, gbp2[g], gbp1[g]);
fprintf(output4, "%s_%d %d %.4f %.6f %.6f %.6f %d %.0f %.0f %.0f %.0f ", gnames[g], found+1, gends[g]-gstarts[g], weightsum, minfloat, maxfloat, sum/num_samples_use, gchr[g], gbp1[g], gbp2[g], cmbp[gstarts[g]], cmbp[gends[g]-1]);
}

if(strcmp(pvafile,"blank")==0){fprintf(output4, "NA NA\n");}
else{fprintf(output4, "%.2e %.2e\n", minpvalue, minpvalue2);}

found++;
}	//end of count2>0

start=gstarts[g];if(gends[g]>end){end=gends[g];}
}}	//end of using g and g loop

fclose(output2);
fclose(output3);
fclose(output4);

if(found==0)
{
printf("\nError, all of the predictors are trivial");
if(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999){printf(" or failed QC");}
printf("\n\n");exit(1);
}
if(xcount>10)
{
printf("In total, %d genes/chunks are trivial", xcount);
if(minmaf!=-9999||maxmaf!=-9999||minvar!=-9999||minobs!=-9999){printf(" or failed QC");}
printf("\n");
}
if(ecount>10){printf("In total, %d genes/chunks had values above 2\n", ecount);}
printf("\n");

////////

pid=malloc(sizeof(char*)*num_samples_use);
mid=malloc(sizeof(char*)*num_samples_use);
schar=malloc(sizeof(char*)*num_samples_use);
pchar=malloc(sizeof(char*)*num_samples_use);

if(famhead==0)
{
read_strings(famfile, pid, num_samples_use, keepsamps, 3, 0);
read_strings(famfile, mid, num_samples_use, keepsamps, 4, 0);
read_strings(famfile, schar, num_samples_use, keepsamps, 5, 0);
read_strings(famfile, pchar, num_samples_use, keepsamps, 6, 0);
}
else
{
for(i=0;i<num_samples_use;i++)
{
pid[i]=malloc(sizeof(char)*2);strcpy(pid[i],"0");
mid[i]=malloc(sizeof(char)*2);strcpy(mid[i],"0");
schar[i]=malloc(sizeof(char)*2);strcpy(schar[i],"0");
pchar[i]=malloc(sizeof(char)*3);strcpy(pchar[i],"NA");
}
}

sprintf(filename5,"%s.fam", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<num_samples_use;i++)
{fprintf(output5, "%s %s %s %s %s %s\n", ids1[i], ids2[i], pid[i], mid[i], schar[i], pchar[i]);}
fclose(output5);

for(i=0;i<num_samples_use;i++){free(pid[i]);free(mid[i]);free(schar[i]);free(pchar[i]);}
free(pid);free(mid);free(schar);free(pchar);

////////

if(strcmp(genefile,"blank")!=0){printf("Condensed data for %d genes saved in %s, %s and %s, with details in %s\n\n", found, filename2, filename3, filename5, filename4);}
else{printf("Condensed data for %d chunks saved in %s, %s and %s, with details in %s\n\n", found, filename2, filename3, filename5, filename4);}

free(data);free(data2);
if(binary==0){gzclose(datainputgz);}
for(g=0;g<num_genes;g++){free(gnames[g]);}free(gnames);
free(gchr);free(gbp1);free(gbp2);free(gstrand);free(gstarts);free(gends);

///////////////////////////


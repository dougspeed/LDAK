/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Joining correlations

///////////////////////////

//read first root file, which gets datafile, num_samples, num_preds, num_samples_use, num_preds_use (readint), num_windows (readint2) and num_blocks

sprintf(filename, "%s.cors.root", corstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, %s should begin \"Datafile\" or \"Multiple\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readstring);exit(1);}

//flag indicates if using multiple datafiles
if(strcmp(readstring,"Datafile")==0){flag=0;strcpy(datafile,readstring2);}
else{flag=1;}

if(fscanf(input, "Num_Samples %d ", &num_samples)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Samples\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors %d ", &num_preds)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Samples_Used %d ", &num_samples_use)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors_Used %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Windows %d ", &readint2)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Windows\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 7 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Jackknifes %d ", &num_blocks)!=1)
{printf("Error reading Row 8 of %s (should begin \"Num_Jackknifes\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
fclose(input);

//now check other root files
for(k=1;k<num_cors;k++)
{
sprintf(filename, "%s.cors.root", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, %s should begin \"Datafile\" or \"Multiple\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readstring);exit(1);}

if(strcmp(readstring,"Multiple")==0){flag=1;}
if(flag==0)	//see if name matches previous
{
if(strcmp(datafile,readstring2)!=0)
{flag=1;}
}

if(fscanf(input, "Num_Samples %d ", &readint)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Samples\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(readint!=num_samples)
{printf("Error, the correlations with stems %s and %s correspond to data files with different numbers of samples (%d and %d)\n\n", corstems[0], corstems[k], num_samples, readint);exit(1);}

if(fscanf(input, "Num_Predictors %d ", &readint)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(flag==0&&readint!=num_preds)
{printf("Error, the correlations with stems %s and %s correspond to the same data file (%s) but have different numbers of predictors (%d and %d)\n\n", corstems[0], corstems[k], datafile, num_preds, readint);exit(1);}

if(fscanf(input, "Num_Samples_Used %d ", &readint)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(readint!=num_samples_use)
{printf("Error, the correlations with stems %s and %s were made using different numbers of samples (%d and %d)\n\n", corstems[0], corstems[k], num_samples_use, readint);exit(1);}

if(fscanf(input, "Num_Predictors_Used %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Windows %d ", &readint2)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Windows\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 7 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Jackknifes %d ", &readint3)!=1)
{printf("Error reading Row 8 of %s (should begin \"Num_Jackknifes\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(flag==0&&readint3!=num_blocks)
{printf("Error, the correlations with stems %s and %s were constructed using different numbers of jackknifes (%d and %d)\n\n", corstems[0], corstems[k], num_blocks, readint3);exit(1);}

fclose(input);
}

//check sizes of bin, bim, windows, noise and jackknife files - also get total numbers of preds, windows and pairs
num_preds_use=0;
bittotal=0;
smax=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use, num_windows and num_pairs for file from rows 5, 6 and 7 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);
read_integers(filename, &readint2, 1, NULL, 2, 5, 0);
read_long(filename, &scount, 1, NULL, 2, 6, 0);

//check size of bin file
sprintf(filename,"%s.cors.bin", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(double)*readint*4+sizeof(float)*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)sizeof(double)*readint*4+sizeof(float)*scount, ftello(input));exit(1);}
fclose(input);

//check size of bim file
sprintf(filename,"%s.cors.bim", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

//check size of windows file
sprintf(filename,"%s.cors.windows", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint2+1)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint2+1, count);exit(1);}

//check size of noise file
sprintf(filename,"%s.cors.noise", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

//check size of jackknife file
sprintf(filename,"%s.cors.jackknifes", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*num_blocks*readint)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)sizeof(float)*num_blocks*readint, ftello(input));exit(1);}
fclose(input);

num_preds_use+=readint;
bittotal+=readint2;
smax+=scount;
}

////////

//allocate variables

rs=malloc(sizeof(char)*10000000);
chr=malloc(sizeof(int)*num_preds_use); 
cm=malloc(sizeof(double)*num_preds_use);
bp=malloc(sizeof(double)*num_preds_use);
preds=malloc(sizeof(char*)*num_preds_use);
al1=malloc(sizeof(char)*num_preds_use);
al2=malloc(sizeof(char)*num_preds_use);

centres=malloc(sizeof(double)*num_preds_use);
mults=malloc(sizeof(double)*num_preds_use);
sqdevs=malloc(sizeof(double)*num_preds_use);
rjksums=malloc(sizeof(double)*num_preds_use);
rjkaves=malloc(sizeof(double)*num_preds_use);
rjktemp=malloc(sizeof(double)*num_preds_use);

blockstarts=malloc(sizeof(int)*bittotal);
blockends=malloc(sizeof(int)*bittotal);

datarands=malloc(sizeof(double)*num_preds_use);
datarands2=malloc(sizeof(double)*num_blocks*num_preds_use);

//read bim details
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//get predictor details from bim file
sprintf(filename,"%s.cors.bim", corstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
for(j=0;j<readint;j++)
{
if(fscanf(input, "%d %s %lf %lf %c %c ", chr+count, rs, cm+count, bp+count, al1+count, al2+count)!=6)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-cors\"\n\n", j+1, filename);exit(1);}
copy_string(preds,count,rs);
count++;
}
fclose(input);
}

//check the files are ordered
for(j=1;j<num_preds_use;j++)
{
if(chr[j]<chr[j-1])
{
printf("Error, chromosome for %s (%d) is lower than that for %s (%d); the correlations in %s must be in order (e.g., the predictors used to construct the correlations with stem %s must be before those used to construct the correlations with stem %s)\n\n", preds[j], chr[j], preds[j-1], chr[j-1], corslist, corstems[0], corstems[1]);exit(1);
}
if(chr[j]==chr[j-1])
{
if(cm[j]<cm[j-1])
{printf("Error, genetic distance for %s (%.2f) is lower than that for %s (%.2f); the correlations in %s must be in order (e.g., the predictors used to construct the correlations with stem %s must be before those used to construct the correlations with stem %s)\n\n", preds[j], cm[j], preds[j-1], cm[j-1], corslist, corstems[0], corstems[1]);exit(1);}

if(bp[j]<bp[j-1])
{printf("Error, basepair for %s (%.2f) is lower than that for %s (%.2f); the correlations in %s must be in order (e.g., the predictors used to construct the correlations with stem %s must be before those used to construct the correlations with stem %s)\n\n", preds[j], bp[j], preds[j-1], bp[j-1], corslist, corstems[0], corstems[1]);exit(1);}
}
}

//read windows (adjusting indexes)
count=0;
count2=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use and num_windows for file from rows 5 and 6 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);
read_integers(filename, &readint2, 1, NULL, 2, 5, 0);

//read values (then add count and maybe subtract one)
sprintf(filename, "%s.cors.windows", corstems[k]);
read_integers(filename, blockstarts+count2, readint2, NULL, 2, 1, 0);
read_integers(filename, blockends+count2, readint2, NULL, 3, 1, 0);
for(j=0;j<readint2;j++){blockstarts[count2+j]+=count-1;blockends[count2+j]+=count;}

count+=readint;
count2+=readint2;
}

//read noises
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//read values
sprintf(filename, "%s.cors.noise", corstems[k]);
read_values(filename, datarands+count, readint, NULL, 1, 0, 0);

count+=readint;
}

//read jackknifes
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//read values
sprintf(filename, "%s.cors.jackknifes", corstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readfloats=malloc(sizeof(float)*readint);

for(p=0;p<num_blocks;p++)
{
if(fread(readfloats, sizeof(float), readint, input)!=readint)
{printf("Error reading jackknifes for Block %d from %s\n\n", p+1, filename);exit(1);}

for(j=0;j<readint;j++){datarands2[(size_t)p*num_preds_use+count+j]=(double)readfloats[j];}
}

fclose(input);
free(readfloats);

count+=readint;
}

////////

//read means, scalings, variances and rjksums
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//get counts, means, scales, variances and rjksums from bin file
sprintf(filename,"%s.cors.bin", corstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fread(centres+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor means from %s\n\n", filename);exit(1);}
if(fread(mults+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor scalings from %s\n\n", filename);exit(1);}
if(fread(sqdevs+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor variances from %s\n\n", filename);exit(1);}
if(fread(rjksums+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor variances from %s\n\n", filename);exit(1);}

fclose(input);
count+=readint;
}

//compute average ldscore for each block - will be using all predictors, and no empty blocks
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];

sum=0;for(j=bitstart;j<bitend;j++){sum+=rjksums[j];}
mean=sum/(bitend-bitstart);
for(j=bitstart;j<bitend;j++){rjkaves[j]=mean;}
}

//work out median window ldscore
for(j=0;j<num_preds_use;j++){rjktemp[j]=rjkaves[j];}
qsort(rjktemp, num_preds_use, sizeof(double), compare_double);
med=rjktemp[num_preds_use/2];

//open output bin and save means, scalings, variances and rjksums
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fwrite(centres, sizeof(double), num_preds_use, output2);
fwrite(mults, sizeof(double), num_preds_use, output2);
fwrite(sqdevs, sizeof(double), num_preds_use, output2);
fwrite(rjksums, sizeof(double), num_preds_use, output2);

//now copy correlations to output bin
readfloats=malloc(sizeof(float)*1000000);

count=0;
for(k=0;k<num_cors;k++)
{
printf("Processing correlations from %s.cors.bin (File %d out of %d)\n", corstems[k], k+1, num_cors);

//get num_preds_use and num_pairs for file from rows 5 and 7 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);
read_long(filename, &scount, 1, NULL, 2, 6, 0);

//get correlations from bin file
sprintf(filename,"%s.cors.bin", corstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, sizeof(double)*4*readint, SEEK_SET);

while(scount>1000000)
{
if(fread(readfloats, sizeof(float), 1000000, input)!=1000000)
{printf("Error reading correlations from %s\n\n", filename);exit(1);}
fwrite(readfloats, sizeof(float), 1000000, output2);
scount-=1000000;
}

if(scount>0)
{
if(fread(readfloats, sizeof(float), scount, input)!=scount)
{printf("Error reading correlations from %s\n\n", filename);exit(1);}
fwrite(readfloats, sizeof(float), scount, output2);
}

fclose(input);
count+=readint;
}	//end of k loop
printf("\n");

free(readfloats);
fclose(output2);

//save windows
sprintf(filename2,"%s.cors.windows", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Window Start_Predictor End_Predictor Start_Location End_Location\n");
for(bit=0;bit<bittotal;bit++){fprintf(output2, "%d %d %d %d:%.0f %d:%.0f\n", bit+1, blockstarts[bit]+1, blockends[bit], chr[blockstarts[bit]], bp[blockstarts[bit]], chr[blockends[bit]-1], bp[blockends[bit]-1]);}
fclose(output2);

//save root
sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(flag==0){fprintf(output3,"Datafile %s\n", datafile);}
else{fprintf(output3,"Multiple Datafiles\n");}
fprintf(output3,"Num_Samples %d\n", num_samples);
if(flag==0){fprintf(output3,"Num_Predictors %d\n", num_preds);}
else{fprintf(output3,"Num_Predictors NA\n");}
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, num_preds_use);
fprintf(output3,"Num_Windows %d\n", bittotal);
fprintf(output3,"Num_Pairs %jd\n", smax);
fprintf(output3,"Num_Jackknifes %d\n", num_blocks);
fclose(output3);

//save bim
sprintf(filename4,"%s.cors.bim", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(j=0;j<num_preds_use;j++)
{
fprintf(output4, "%d\t%s\t", chr[j], preds[j]);
if(cm[j]==0){fprintf(output4, "0\t");}
else{fprintf(output4, "%.6f\t", cm[j]);}
fprintf(output4, "%ld\t%c\t%c\n", (long int)bp[j], al1[j], al2[j]);
}
fclose(output4);

//save highld predictors
sprintf(filename6,"%s.cors.highld",outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
count=0;
for(j=0;j<num_preds_use;j++)
{
if(rjkaves[j]>5*med){fprintf(output6, "%s\n", preds[j]);count++;}
}
if(count==0){fprintf(output6,"None\n");}
fclose(output6);

//save noise
sprintf(filename7,"%s.cors.noise", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(j=0;j<num_preds_use;j++){fprintf(output7, "%.6f\n", datarands[j]);}
fclose(output7);

//save jackknife
writefloats=malloc(sizeof(float)*num_preds_use);

sprintf(filename8,"%s.cors.jackknifes",outfile);
if((output8=fopen(filename8,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename8);exit(1);}
for(p=0;p<num_blocks;p++)
{
for(j=0;j<num_preds_use;j++){writefloats[j]=(float)datarands2[(size_t)p*num_preds_use+j];}
fwrite(writefloats, sizeof(float), num_preds_use, output8);
}
fclose(output8);

free(writefloats);

printf("The joined correlations are saved in files with prefix %s\n\n", outfile);

free(rs);free(chr);free(cm);free(bp);free(al1);free(al2);
for(j=0;j<num_preds_use;j++){free(preds[j]);}free(preds);
free(centres);free(mults);free(sqdevs);free(rjksums);free(rjkaves);free(rjktemp);
free(blockstarts);free(blockends);
free(datarands);free(datarands2);

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Joining correlations

///////////////////////////

//read first root file, which gets num_samples, num_samples_use, maybe datafile and num_preds, cutoff, window_kb/window_cm
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
if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Threshold %lf ", &mincor)!=1)
{printf("Error reading Row 7 of %s (should begin \"Threshold\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 8 of %s, suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(strcmp(readstring,"Window_kb")!=0&&strcmp(readstring,"Window_cM")!=0)
{printf("Error, Row 8 of %s should begin with \"Window_kb\" or \"Window_cM\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n", filename, readstring);exit(1);}
if(strcmp(readstring,"Window_kb")==0){window_kb=atof(readstring2);}
else{window_cm=atof(readstring2);}
fclose(input);

//check size of first bin file
sprintf(filename,"%s.cors.bin", corstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)(sizeof(int)+sizeof(double)*3)*readint+(sizeof(int)+sizeof(float))*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)(sizeof(int)+sizeof(double)*3)*readint+(sizeof(int)+sizeof(float))*scount, ftello(input));exit(1);}
fclose(input);

//check size of first bim file
sprintf(filename,"%s.cors.bim", corstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

////////

//now read remaining roots, bins and bims, checking consistent and sizes, and also get num_preds_use
num_preds_use=readint;
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
if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Threshold %lf ", &readdouble)!=1)
{printf("Error reading Row 7 of %s (should begin \"Threshold\"), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(readdouble!=mincor)
{printf("Error, the correlations with stems %s and %s were constructed using different thresholds (%f and %f)\n\n", corstems[0], corstems[k], mincor, readdouble);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 8 of %s, suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename);exit(1);}
if(strcmp(readstring,"Window_kb")!=0&&strcmp(readstring,"Window_cM")!=0)
{printf("Error, Row 8 of %s should begin with \"Window_kb\" or \"Window_cM\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\"\n\n", filename, readstring);exit(1);}
if(atof(readstring2)!=window_kb&&atof(readstring2)!=window_cm)
{printf("Error, the correlations with stems %s and %s were constructed using different window sizes\n\n", corstems[0], corstems[k]);exit(1);}

fclose(input);

sprintf(filename,"%s.cors.bin", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)(sizeof(int)+sizeof(double)*3)*readint+(sizeof(int)+sizeof(float))*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)(sizeof(int)+sizeof(double)*3)*readint+(sizeof(int)+sizeof(float))*scount, ftello(input));exit(1);}
fclose(input);

sprintf(filename,"%s.cors.bim", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

sprintf(filename,"%s.cors.ldscores", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint+1)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint+1, count);exit(1);}

sprintf(filename,"%s.cors.noise", corstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\"\n\n", filename);exit(1);}
count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\"\n\n",filename, readint, count);exit(1);}

num_preds_use+=readint;
}	//end of k loop

////////

//allocate variables
rs=malloc(sizeof(char)*10000000);
chr=malloc(sizeof(int)*num_preds_use); 
cm=malloc(sizeof(double)*num_preds_use);
bp=malloc(sizeof(double)*num_preds_use);
preds=malloc(sizeof(char*)*num_preds_use);
al1=malloc(sizeof(char)*num_preds_use);
al2=malloc(sizeof(char)*num_preds_use);

actnums=malloc(sizeof(int)*num_preds_use);
centres=malloc(sizeof(double)*num_preds_use);
mults=malloc(sizeof(double)*num_preds_use);
sqdevs=malloc(sizeof(double)*num_preds);

bigs=malloc(sizeof(int*)*num_preds_use);
rjks=malloc(sizeof(float*)*num_preds_use);
rjksums=malloc(sizeof(double)*num_preds_use);
rjkaves=malloc(sizeof(double)*num_preds_use);
rjktemp=malloc(sizeof(double)*num_preds_use);

datarands=malloc(sizeof(double)*num_preds_use);

//read bim details, checking in order
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//get predictor details from bim file
sprintf(filename,"%s.cors.bim", corstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
for(j=0;j<readint;j++)
{
if(fscanf(input, "%d %s %lf %lf %c %c ", chr+count, rs, cm+count, bp+count, al1+count, al2+count)!=6)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-cors\"\n\n", j+1, filename);exit(1);}
copy_string(preds,count,rs);
count++;
}
fclose(input);
}

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

//read ldscores
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//read ldscores
sprintf(filename, "%s.cors.ldscores", corstems[k]);
read_values(filename, rjksums+count, readdouble, NULL, 2, 1, 0);
read_values(filename, rjkaves+count, readdouble, NULL, 3, 1, 0);

count+=readint;
}

//read noises
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//read noise
sprintf(filename, "%s.cors.noise", corstems[k]);
read_values(filename, datarands+count, readint, NULL, 1, 0, 0);

count+=readint;
}

//work out median region ldscore
for(j=0;j<num_preds_use;j++){rjktemp[j]=rjkaves[j];}
qsort(rjktemp, num_preds_use, sizeof(double), compare_double);

////////

//read counts, means, scales and variances
count=0;
for(k=0;k<num_cors;k++)
{
//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//get counts, means, scales and variances from bin file
sprintf(filename,"%s.cors.bin", corstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}

if(fread(actnums+count, sizeof(int), readint, input)!=readint)
{printf("Error reading predictor counts from %s\n\n", filename);exit(1);}
if(fread(centres+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor means from %s\n\n", filename);exit(1);}
if(fread(mults+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor scalings from %s\n\n", filename);exit(1);}
if(fread(sqdevs+count, sizeof(double), readint, input)!=readint)
{printf("Error reading predictor variances from %s\n\n", filename);exit(1);}

fclose(input);
count+=readint;
}

//open output bin and save counts, means, scalings and variances
sprintf(filename2,"%s.cors.bin", outfile);
if((output2=fopen(filename2,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fwrite(actnums, sizeof(int), num_preds_use, output2);
fwrite(centres, sizeof(double), num_preds_use, output2);
fwrite(mults, sizeof(double), num_preds_use, output2);
fwrite(sqdevs, sizeof(double), num_preds_use, output2);

//now copy correlations to output bin, adjusting indexes
count=0;
for(k=0;k<num_cors;k++)
{
printf("Processing File %d out of %d\n", k+1, num_cors);

//get num_preds_use for file from row 5 of root file
sprintf(filename, "%s.cors.root", corstems[k]);
read_integers(filename, &readint, 1, NULL, 2, 4, 0);

//get correlations from bin file
sprintf(filename,"%s.cors.bin", corstems[k]);
if((input=fopen(filename,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fseeko(input, (sizeof(int)+sizeof(double)*3)*readint, SEEK_SET);

for(j=count;j<count+readint;j++)
{
if(actnums[j]>0)
{
bigs[j]=malloc(sizeof(int)*actnums[j]);
rjks[j]=malloc(sizeof(float)*actnums[j]);

if(fread(bigs[j], sizeof(int), actnums[j], input)!=actnums[j])
{printf("Error reading predictor indexes from %s\n\n", filename);exit(1);}
if(fread(rjks[j], sizeof(float), actnums[j], input)!=actnums[j])
{printf("Error reading predictor correlations from %s\n\n", filename);exit(1);}

//add on count to indexes
for(j2=0;j2<actnums[j];j2++){bigs[j][j2]+=count;}

fwrite(bigs[j], sizeof(int), actnums[j], output2);
fwrite(rjks[j], sizeof(float), actnums[j], output2);

free(bigs[j]);
free(rjks[j]);
}
}

fclose(input);
count+=readint;
}	//end of k loop
printf("\n");
fclose(output2);

//save root
scount=0;for(j=0;j<num_preds_use;j++){scount+=actnums[j];}

sprintf(filename3,"%s.cors.root", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
if(flag==0){fprintf(output3,"Datafile %s\n", datafile);}
else{fprintf(output3,"Multiple Datafiles\n");}
fprintf(output3,"Num_Samples %d\n", num_samples);
if(flag==0){fprintf(output3,"Num_Predictors %d\n", num_preds);}
else{fprintf(output3,"Num_Predictors NA\n");}
fprintf(output3,"Num_Samples_Used %d\nNum_Predictors_Used %d\n", num_samples_use, num_preds_use);
fprintf(output3,"Num_Pairs %jd\n", scount);
fprintf(output3,"Threshold %.4e\n", mincor);
if(window_cm!=-9999){fprintf(output3,"Window_cM %.4f\n", window_cm);}
else{fprintf(output3,"Window_kb %.4f\n", window_kb);}
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

//save ldscores
sprintf(filename5,"%s.cors.ldscores",outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
fprintf(output5, "Predictor Individual Region\n");
for(j=0;j<num_preds_use;j++){fprintf(output5, "%s %.4f %.4f\n", preds[j], rjksums[j], rjkaves[j]);}
fclose(output5);

//save highld predictors
sprintf(filename6,"%s.cors.highld",outfile);
if((output6=fopen(filename6,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename6);exit(1);}
count=0;
for(j=0;j<num_preds_use;j++)
{
if(rjkaves[j]>5*rjktemp[num_preds_use/2]){fprintf(output6, "%s\n", preds[j]);count++;}
}
if(count==0){fprintf(output6,"None\n");}
fclose(output6);

//save noise
sprintf(filename7,"%s.cors.noise", outfile);
if((output7=fopen(filename7,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename7);exit(1);}
for(j=0;j<num_preds_use;j++){fprintf(output7, "%.10f\n", datarands[j]);}
fclose(output7);

printf("The joined correlations are saved in files with prefix %s\n\n", outfile);

free(rs);free(chr);free(cm);free(bp);free(al1);free(al2);
for(j=0;j<num_preds_use;j++){free(preds[j]);}free(preds);
free(actnums);free(centres);free(mults);free(sqdevs);
free(bigs);free(rjks);free(rjksums);free(rjkaves);free(rjktemp);

///////////////////////////


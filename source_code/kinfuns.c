/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Kinship functions

///////////////////////////

int cut_partitions(int length, int *chr, char **preds, int part_length, int bychr, int num_parts, char *partpref, int checkpart, char *folder, char *datafile, char *bimfile, int extract)
{
int j, q, count, count2, count3;
int *pstarts, *pends, *usedpreds, *indexer;

char **wantpreds, **preds2;

char filename[500], filename2[500], filename3[500];
FILE *output, *output2, *output3;


if(part_length!=-9999)	//divide evenly
{
num_parts=(length-1)/part_length+1;
count=(length-1)/num_parts+1;

pstarts=malloc(sizeof(int)*num_parts);
pends=malloc(sizeof(int)*num_parts);
for(q=0;q<num_parts;q++)
{
pstarts[q]=q*count;
pends[q]=(q+1)*count;
if(pends[q]>length){pends[q]=length;}
}
}

if(bychr==1)	//divide by chromosomes
{
pstarts=malloc(sizeof(int)*length);
pends=malloc(sizeof(int)*length);

num_parts=0;count=0;
pstarts[0]=0;
for(j=1;j<length;j++)
{
if(chr[j]!=chr[j-1])	//new chromosome
{
pends[num_parts]=j;
printf("Chromosome %d contains %d predictors\n", chr[j-1], pends[num_parts]-pstarts[num_parts]);
if(pends[num_parts]-pstarts[num_parts]>count){count=pends[num_parts]-pstarts[num_parts];}
num_parts++;
pstarts[num_parts]=j;

if(chr[j]==chr[j-1]+2){printf("Warning, no predictors on Chromosome %d\n", chr[j-1]+1);}
if(chr[j]>chr[j-1]+2){printf("Warning, no predictors on Chromosomes %d-%d\n", chr[j-1]+1, chr[j]-1);}
}
}
pends[num_parts]=length;
if(pends[num_parts]-pstarts[num_parts]>count){count=pends[num_parts]-pstarts[num_parts];}
num_parts++;
}

if(strcmp(partpref,"blank")!=0&&checkpart==1)	//partitions provided - just need to check
{
//will speed up if we use a sorted version of preds
preds2=malloc(sizeof(char*)*length);
for(j=0;j<length;j++){preds2[j]=preds[j];}
qsort(preds2, length, sizeof(char *), compare_string);

usedpreds=malloc(sizeof(int)*length);
for(j=0;j<length;j++){usedpreds[j]=0;}

count3=0;
for(q=0;q<num_parts;q++)
{
sprintf(filename, "%s%d", partpref, q+1);
count=countrows(filename);
printf("Reading %d predictors from %s\n", count, filename);
wantpreds=malloc(sizeof(char*)*count);
read_strings(filename, wantpreds, count, NULL, 1, 0);
indexer=malloc(sizeof(int)*count);
count2=find_strings(preds2, length, wantpreds, count,  indexer, NULL, NULL, filename, NULL, NULL, 2);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
if(count2<count){printf("Warning, only %d of these are in the data\n", count2);}
for(j=0;j<count2;j++){usedpreds[indexer[j]]++;}
if(count2>count3){count3=count2;}
for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}	//end of q loop

if(num_parts>1)
{
count=0;
for(j=0;j<length;j++){count+=(usedpreds[j]>=1);}
if(count==length)
{printf("\nAll %d predictors are in at least one partition", length);}
else
{printf("\n%d of the %d predictors are in at least one partition",count, length);}

count2=0;
for(j=0;j<length;j++){count2+=(usedpreds[j]>1);}
if(count2>0)
{printf("\nWarning, %d predictors are in more than one partition; partitions are usually (but not necessarily) non-overlapping\n", count2);}
else{printf(" (no predictor is in more than one partition)\n");}
}
printf("\n");

free(preds2);free(usedpreds);
}	//end of partitions provided

////////

//open file which stores partition details and write header lines
sprintf(filename,"%spartition.details",folder);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Datafiles %s %s\n", datafile, bimfile);
if(extract==1){fprintf(output, "Using Filtered Predictors\n");}
else{fprintf(output, "Using All Predictors\n");}
fprintf(output,"Partition Start_Predictor End_Predictor\n");

for(q=0;q<num_parts;q++)
{
if(strcmp(partpref,"blank")!=0){fprintf(output, "%d -1 %s%d\n", q+1, partpref, q+1);}
else{fprintf(output, "%d %d %d\n", q+1, pstarts[q]+1, pends[q]);}
}
fclose(output);

sprintf(filename2,"%spartition.number",folder);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"%d\n", num_parts);
fclose(output2);

sprintf(filename3,"%spartition.list",folder);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
for(q=0;q<num_parts;q++){fprintf(output3,"%skinships.%d\n", folder, q+1);}
fclose(output3);

if(part_length!=-9999)
{printf("The %d predictors have been split into %d partitions with (approx) length %d; ", length, num_parts, part_length);}
if(bychr==1)
{printf("The %d predictors have been split by chromosome into %d partitions of average (maximum) length %d (%d); ", length, num_parts, length/num_parts, count);}
if(strcmp(partpref,"blank")!=0&&checkpart==1)
{printf("The %d predictors have been split into %d partitions, with average (maximum) length %d (%d); ", count, num_parts, count/num_parts, count3);}

printf("details saved in %spartition.details\n\n", folder);

if(strcmp(partpref,"blank")==0){free(pstarts);free(pends);}

return(0);
}	//end of cut_partitions

///////////////////////////

int manip_kins(char *outfile, int num_kins, char **kinstems, double *kinsums, char **ids1, char **ids2, char **ids3, int ns, int kingz, int kinraw, int type, double *scales, int maxthreads)
//type=0 - join-kins, type=1 - add, type=2 - subtract, type=3 - solve null
{
int j, j2, k, count, count2, count3, flag, samedata;
double value;

int length, *indexer, *usedpreds, *keeppreds_use;
char **preds, **wantpreds;

double power;

int *kindex, *kinrecord;
double *centres, *mults, *weights, *exps, *kcentres, *kmults, *kweights, *kexps;
char **kpreds, *al1, *al2, *kal1, *kal2;

float *kins_single;

int readint;
char readchar, *rs, readstring[500], readstring2[500], datafile[500];

char filename[500], filename2[500];
FILE *input;

rs=malloc(sizeof(char)*10000000);


//read adjust files recording datafile and power
flag=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.adjust", kinstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(k==0&&strcmp(readstring,"Datafile")==0){strcpy(datafile,readstring2);}
if(strcmp(readstring,"Multiple")==0){flag=1;}
if(flag==0)	//see if name matches previous
{
if(strcmp(datafile,readstring2)!=0)
{
if(type==0)	//should not mismatch when using join-kins
{printf("Error, the kinship matrices %s and %s were constructed from different data files (%s and %s); re-examine your commands, or override this error by manually editing the \".grm.adjust\" files\n\n", kinstems[0], kinstems[k], datafile, readstring2);exit(1);}
flag=1;
}
}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(k==0&&strcmp(readstring,"Power")==0){power=atof(readstring2);}
if(strcmp(readstring,"Multiple")==0){power=-9999;}
if(power!=-9999)	//check power matches previous
{
if(power!=atof(readstring2)){power=-9999;}
}

fclose(input);
}	//end of k loop

////////

//now sort out details - have already checked details exist
printf("Reading kinship matrix details\n");

//first read in all predictor names and indexes
count=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.details", kinstems[k]);
count+=countrows(filename)-1;
}

wantpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*count);
count3=0;
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.details", kinstems[k]);
count2=countrows(filename)-1;
printf("%s was constructed from %d predictors\n", filename, count2);

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}
for(j=0;j<count2;j++)
{
if(fscanf(input, "%s %d ", rs, &readint)!=2)
{printf("Error reading Row %d of %s\n\n", j+2, filename);exit(1);}
if(readint<=0)
{printf("Error reading Row %d of %s; index is not positive (%d), suggesting the file has been changed since creation\n\n", j+2, filename, readint);exit(1);}
copy_string(wantpreds, count3, rs);
indexer[count3]=readint-1;
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}
count3++;
}
fclose(input);
}
printf("\n");

//see if indexes valid and make a list of predictors
length=0;
for(j=0;j<count;j++)
{
if(indexer[j]+1>length){length=indexer[j]+1;}
}

preds=malloc(sizeof(char*)*length);
usedpreds=malloc(sizeof(int)*length);
for(j=0;j<length;j++){usedpreds[j]=0;}

samedata=1;
for(j=0;j<count;j++)
{
j2=indexer[j];
if(usedpreds[j2]==0)	//new index
{copy_string(preds, j2, wantpreds[j]);usedpreds[j2]=1;}
else	//has already appeared
{
if(strcmp(wantpreds[j],preds[j2])!=0)	//mismatch
{
if(type==0)	//should not mismatch when using join-kins
{printf("Error reading %s; index %d refers to %s, but it previously referred to %s, indicating that the kinship matrices were not all created from the same datafile\"\n\n", filename, j2+1, wantpreds[j], preds[j2]);exit(1);}
samedata=0;
break;
}
}
}

if(samedata==0)	//cant use indexes, so fill preds with sorted version of wantpreds
{
printf("Finding new predictor indexes; if there are very many predictors, this can take a while\n\n");
//reset preds
for(j=0;j<length;j++)
{
if(usedpreds[j]==1){free(preds[j]);}
}
free(preds);free(usedpreds);
preds=malloc(sizeof(char*)*count);
usedpreds=malloc(sizeof(int)*count);

qsort(wantpreds, count, sizeof(char *), compare_string);
copy_string(preds,0,wantpreds[0]);
length=1;
for(j=1;j<count;j++)
{
if(strcmp(wantpreds[j],preds[length-1])!=0)	//its new
{copy_string(preds,length,wantpreds[j]);length++;}
}
for(j=0;j<length;j++){usedpreds[j]=1;}
}

for(j=0;j<count;j++){free(wantpreds[j]);}
free(wantpreds);
free(indexer);

////////

//now read in details
centres=malloc(sizeof(double)*length);
mults=malloc(sizeof(double)*length);
weights=malloc(sizeof(double)*length);
al1=malloc(sizeof(char)*length);
al2=malloc(sizeof(char)*length);
exps=malloc(sizeof(double)*length);

//kinrecord indicates in which kinship matrix each predictor first found
kinrecord=malloc(sizeof(int)*length);
for(j=0;j<length;j++){kinrecord[j]=-1;}

//read in details for each kinship
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.details", kinstems[k]);
count=countrows(filename)-1;
kpreds=malloc(sizeof(char*)*count);
kindex=malloc(sizeof(int)*count);
kcentres=malloc(sizeof(double)*count);
kmults=malloc(sizeof(double)*count);
kweights=malloc(sizeof(double)*count);
kal1=malloc(sizeof(char)*count);
kal2=malloc(sizeof(char)*count);
kexps=malloc(sizeof(double)*count);
read_details(filename, kpreds, kindex, kcentres, kmults, kweights, kal1, kal2, kexps, count);

if(samedata==0)	//update kindex
{(void)find_strings(kpreds, count, preds, length,  NULL, kindex, NULL, filename, NULL, NULL, 1);}

//now can load up
value=1;if(scales!=NULL){value=scales[k];}
for(j=0;j<count;j++)
{
j2=kindex[j];
if(kinrecord[j2]==-1)	//it's new
{
if(type==2&&k!=0)
{printf("Error reading %s; can not subtract Predictor %s as it was not in %s.grm.details\n\n", filename, preds[j2], kinstems[0]);exit(1);}

centres[j2]=kcentres[j];mults[j2]=kmults[j];
weights[j2]=kweights[j]*value;
al1[j2]=kal1[j];al2[j2]=kal2[j];
exps[j2]=kexps[j]*value;
kinrecord[j2]=k;
}
else	//already used (so can't have k=0)
{
if(kal1[j]!=al1[j2]&&kal1[j]!=al2[j2])
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s.grm.details (%c and %c)\n\n", filename, kpreds[j], kal1[j], kal2[j], kinstems[kinrecord[j2]], al1[j2], al2[j2]);exit(1);}
if(kal2[j]!=al1[j2]&&kal2[j]!=al2[j2])
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s.grm.details (%c and %c)\n\n", filename, kpreds[j], kal1[j], kal2[j], kinstems[kinrecord[j2]], al1[j2], al2[j2]);exit(1);}

if(kal1[j]!=al1[j2])	//turn
{kcentres[j]=2-kcentres[j];}

if(fabs(kcentres[j]-centres[j2])>0.00001||fabs(kmults[j]-mults[j2])>0.00001)	//should match
{printf("Error reading %s; centres and scaling for %s different to those in %s.grm.details; if a predictor features in more than one kinship file, the same samples and predictor scalings must be used when constructing each\n\n", filename, kpreds[j], kinstems[kinrecord[j2]]);exit(1);}

if(type==0||type==1||type==3)
{weights[j2]+=kweights[j]*value;exps[j2]+=kexps[j]*value;}
else
{
weights[j2]-=kweights[j]*value;
exps[j2]-=kexps[j]*value;
if(weights[j2]<0||exps[j2]<0)
{printf("Error reading %s, predictor %s has been subtracted \"too much\"\n\n", filename2, preds[j2]);exit(1);}
}
}	//end of already used
}	//end of j loop

for(j=0;j<count;j++){free(kpreds[j]);}free(kpreds);
free(kindex);free(kcentres);free(kmults);free(kweights);free(kal1);free(kal2);free(kexps);
}	//end of k loop

count=length;
for(j=0;j<length;j++)
{
if(kinrecord[j]==-1||weights[j]==0){mults[j]=-9999;count--;}
}
if(count==0)
{
sprintf(filename, "%s.grm.details", kinstems[0]);
printf("Error, all %d predictors in the first kinship matrix have been subtracted\n\n", countrows(filename)-1);
exit(1);
}

//for printing, will require keeppreds_use (an index for each SNP)
keeppreds_use=malloc(sizeof(int)*length);
for(j=0;j<length;j++){keeppreds_use[j]=j;}

////////

//allocate kins and set to zero
kin_warn(1, ns, 0, 1);
kins_single=calloc((size_t)ns*ns,sizeof(float));

for(k=0;k<num_kins;k++)	//read through, adding or subtracting
{
value=1;if(scales!=NULL){value=scales[k];}
if(type==2&&k>0){value=-value;}
if(value!=0){read_kins(kinstems[k], NULL, kins_single, kinsums[k]*value, ns, ids3, 1, maxthreads);}
}

//now write
if(flag==0)
{write_kins(outfile, NULL, kins_single, ns, ids1, ids2, 1, preds, keeppreds_use, centres, mults, weights, al1, al2, exps, length, datafile, power, kingz, kinraw, 1+(type==3));}
else
{write_kins(outfile, NULL, kins_single, ns, ids1, ids2, 1, preds, keeppreds_use, centres, mults, weights, al1, al2, exps, length, NULL, power, kingz, kinraw, 1+(type==3));}

////////

for(j=0;j<length;j++){if(usedpreds[j]==1){free(preds[j]);}}
free(preds);free(usedpreds);
free(centres);free(mults);free(weights);free(al1);free(al2);free(exps);
free(kinrecord);free(keeppreds_use);
free(kins_single);

free(rs);
return(0);
}	//end of manip_kins

///////////////////////////


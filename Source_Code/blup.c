/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//BLUP, pca loads and scores (and find tags)

///////////////////////////

int read_effects(double **blupcentres, double **blupfactors, int *keeppreds_use, int num_kins, char **kinstems, int num_scores, char *efffile, int num_preds, char **allpreds, char *allal1, char *allal2, int *predorder, char *bimfile)
{
int j, j2, j3, k, r, count, count2, flag;
int *indexer, *indexer2;
double value, value2, **allcentres, **allfactors, *centres, *mults, *weights, *exps, *factors;
char **preds, *al1, *al2;

char readchar, readstring[500], *rs, *ra1, *ra2;

char filename[500];
FILE *input;

rs=malloc(sizeof(char)*10000000);
ra1=malloc(sizeof(char)*10000000);ra2=malloc(sizeof(char)*10000000);


//allocate allcentres and allfactors
value=(double)num_preds/1024/1024/1024*8*(num_kins+num_scores)*2;
if(value>1){printf("Warning, to read the centres and scales for all predictors requires %.1f Gb\n\n", value);}

allcentres=malloc(sizeof(double*)*(num_kins+num_scores));
allfactors=malloc(sizeof(double*)*(num_kins+num_scores));
for(k=0;k<num_kins+num_scores;k++)
{allcentres[k]=malloc(sizeof(double)*num_preds);allfactors[k]=malloc(sizeof(double)*num_preds);}

//set allcentres and allfactors to zero
for(k=0;k<num_kins+num_scores;k++)
{
for(j=0;j<num_preds;j++){allcentres[k][j]=0;allfactors[k][j]=0;}
}

//read in kinship details - all predictors must be in data
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.details", kinstems[k]);
count=countrows(filename)-1;
printf("Reading details for %d predictors from %s\n", count, filename);

preds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*count);
centres=malloc(sizeof(double)*count);
mults=malloc(sizeof(double)*count);
weights=malloc(sizeof(double)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);
exps=malloc(sizeof(double)*count);
read_details(filename, preds, indexer, centres, mults, weights, al1, al2, exps, count);

//see if indexer correct
flag=0;
for(j=0;j<count;j++)
{
j2=indexer[j];
if(j2<0||j2>=num_preds){flag=1;break;}
if(strcmp(allpreds[j2],preds[j])!=0){flag=1;break;}
}
if(flag==1)	//update indexer
{
count2=find_strings(preds, count, allpreds, num_preds, NULL, indexer, filename, NULL, NULL, predorder, 3);
if(count2!=count){printf("Error, only %d of these are in the data\n\n", count2);exit(1);}
}

//now load up - using kinships so alleles will be different and should be consistent
for(j=0;j<count;j++)
{
j2=indexer[j];
if((al1[j]!=allal1[j2]&&al1[j]!=allal2[j2])||(al2[j]!=allal1[j2]&&al2[j]!=allal2[j2]))
{printf("Error reading %s; alleles for %s (%c and %c) not consistent with those in %s (%c and %c), indicating that the kinship matrix was computed using a different datafile\n\n", filename, allpreds[j2], al1[j], al2[j], bimfile, allal1[j2], allal2[j2]);exit(1);}

//set assuming alleles align, then turn if necessary (reading from kinship details, so will always have centre)
allcentres[k][j2]=centres[j];allfactors[k][j2]=mults[j]*pow(weights[j],.5);
if(al1[j]!=allal1[j2])
{allcentres[k][j2]=2-allcentres[k][j2];allfactors[k][j2]=-allfactors[k][j2];}
}	//end of j loop

for(j=0;j<count;j++){free(preds[j]);}free(preds);
free(indexer);free(centres);free(mults);free(weights);free(al1);free(al2);free(exps);
}	//end of k loop

////////

if(num_scores>0)	//read in regions (must be type 0)
{
count=countrows(efffile)-1;

if(count==0)	//will not change any of allfactors (or allcentres)
{printf("Warning, there are no predictors in %s\n", efffile);}
else
{
printf("Reading details for %d predictors from %s\n", count, efffile);

preds=malloc(sizeof(char*)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);
centres=malloc(sizeof(double)*count);
factors=malloc(sizeof(double)*count*num_scores);
indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);

//open, check then skip header line
if((input=fopen(efffile,"r"))==NULL)
{printf("Error opening %s\n\n", efffile);exit(1);}
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading first element of %s\n\n", efffile);exit(1);}
if(strcmp(readstring,"Predictor")!=0&&strcmp(readstring,"SNP")!=0)
{printf("Error reading %s; first element should be \"Predictor\" or \"SNP\" (not %s)\n\n", efffile, readstring);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %s %s ", rs, ra1, ra2)!=3)
{printf("Error reading first three elements of Row %d of %s\n\n", j+2, efffile);exit(1);}
copy_string(preds,j,rs);
if(strlen(ra1)>1||strlen(ra2)>1)
{printf("Error reading Row %d of %s; alleles must be single characters (not %s and %s)\n\n", j+2, efffile, ra1, ra2);exit(1);}
al1[j]=ra1[0];al2[j]=ra2[0];
if(al1[j]==al2[j])
{printf("Error reading %s; both alleles for %s are the same (%c)\n\n", efffile, preds[j], al1[j]);exit(1);}

/*
if(fscanf(input, "%[0-9eE.+-] ", rc)==1){centres[j]=atof(rc);}
else
{
if(fscanf(input, "%s ", rc)!=1)
{printf("Error reading centre for Predictor %s in %s\n\n", preds[j], efffile);exit(1);}
if(strcmp(rc,"NA")==0)
{printf("Error reading %s; centre for Predictor %s is NA\n\n", efffile, preds[j]);exit(1);}
else
{printf("Error reading %s; centre for Predictor %s is unrecognisable (%s)\n\n", efffile, preds[j], rc);exit(1);}
}
*/

if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading centre for Predictor %s in %s\n\n", preds[j], efffile);exit(1);}
if(strcmp(readstring,"NA")==0)
{printf("Error reading %s; centre for Predictor %s is NA\n\n", efffile, preds[j]);exit(1);}
else
{
if(sscanf(readstring, "%lf%c", centres+j, &readchar)!=1)
{printf("Error reading %s; centre for Predictor %s does not appear to be numeric (%s)\n\n", efffile, preds[j], readstring);exit(1);}
}

for(r=0;r<num_scores;r++)
{
/*
if(fscanf(input, "%[0-9eE.+-] ", re)==1){factors[j+r*count]=atof(re);}
else
{
if(fscanf(input, "%s ", re)!=1)
{printf("Error reading Effect %d for Predictor %s in %s\n\n", r+1, preds[j], efffile);exit(1);}
if(strcmp(re,"NA")==0)
{printf("Error reading %s; Effect %d for Predictor %s is NA\n\n", efffile, r+1, preds[j]);exit(1);}
printf("Error reading %s; Effect %d for Predictor %s is unrecognisable (%s)\n\n", efffile, r+1, preds[j], re);exit(1);
}
*/

if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Effect %d for Predictor %s in %s\n\n", r+1, preds[j], efffile);exit(1);}
if(strcmp(readstring,"NA")==0)
{printf("Error reading %s; Effect %d for Predictor %s is NA\n\n", efffile, r+1, preds[j]);exit(1);}
else
{
if(sscanf(readstring, "%lf%c", factors+j+r*count, &readchar)!=1)
{printf("Error reading %s; Effect %d for Predictor %s does not appear to be numeric (%s)\n\n", efffile, r+1, preds[j], readstring);exit(1);}
}
}
}	//end of j loop
fclose(input);

//set indexer and indexer2
count2=find_strings(allpreds, num_preds, preds, count, indexer, indexer2, NULL, efffile, predorder, NULL, 3);
if(count2<count)
{
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}
printf("Error, only %d of these are in the data\n\n", count2);exit(1);
}

//now load up
for(j=0;j<count2;j++)
{
j2=indexer[j];
j3=indexer2[j];

if((al1[j3]!=allal1[j2]&&al1[j3]!=allal2[j2])||(al2[j3]!=allal1[j2]&&al2[j3]!=allal2[j2]))
{printf("Error, alleles for %s (%c and %c) not consistent with those in %s (%c and %c)\n", allpreds[j2], al1[j3], al2[j3], bimfile, allal1[j2], allal2[j2]);exit(1);}

for(r=0;r<num_scores;r++)
{
//set centres and effects assuming alleles align, then turn if necessary (note that we might be missing centre)
allcentres[num_kins+r][j2]=centres[j3];
allfactors[num_kins+r][j2]=factors[j3+r*count];
if(al1[j3]!=allal1[j2])
{
if(allcentres[num_kins+r][j2]!=-9999){allcentres[num_kins+r][j2]=2-allcentres[num_kins+r][j2];}
allfactors[num_kins+r][j2]=-allfactors[num_kins+r][j2];
}
}
}	//end of j loop

for(j=0;j<count;j++){free(preds[j]);}free(preds);
free(al1);free(al2);free(centres);free(factors);free(indexer);free(indexer2);
}	//end of count
}	//end of num_scores>0

printf("\n");

////////

//store, reducing to remove gaps - will always have at least one predictor
count=0;
for(j=0;j<num_preds;j++)
{
flag=0;
for(k=0;k<num_kins+num_scores;k++){flag+=(allfactors[k][j]!=0);}
if(flag>0){keeppreds_use[count]=j;count++;}
}

for(k=0;k<num_kins+num_scores;k++)
{
blupcentres[k]=malloc(sizeof(double)*count);
blupfactors[k]=malloc(sizeof(double)*count);
for(j=0;j<count;j++){blupcentres[k][j]=allcentres[k][keeppreds_use[j]];blupfactors[k][j]=allfactors[k][keeppreds_use[j]];}
free(allcentres[k]);
free(allfactors[k]);
}
free(allcentres);free(allfactors);

value2=(double)count/1024/1024/1024*(num_kins+num_scores)*2;
if(value-value2>1){printf("Good news, the memory required has been reduced to %.1f Gb\n\n", value2);}

free(rs);free(ra1);free(ra2);
return(count);
}	//end of read_effects

////////

int read_scores(double **blupcentres, double **blupfactors, int *keeppreds_use, int num_scores, char *efffile, int num_preds, char **allpreds, char *allal1, char *allal2, int *predorder, char *bimfile, int num_preds_use, int *keeppreds, int type)
//type=0 - scores, type=1 - find tags
{
int j, j2, j3, r, count, count2, count3, flag;
int *usedpreds, *indexer;
double value, *centres, **factors;
char **preds, *al1, *al2;

char readchar, readstring[500], *rs, *ra1, *ra2;

FILE *input;

rs=malloc(sizeof(char)*10000000);
ra1=malloc(sizeof(char)*10000000);ra2=malloc(sizeof(char)*10000000);


count=countrows(efffile)-1;

if(count==0)	//add one fake predictor 
{
printf("Warning, there are no predictors in %s\n", efffile);
for(r=0;r<num_scores;r++)
{
blupcentres[r]=malloc(sizeof(double));
blupfactors[r]=malloc(sizeof(double));
blupcentres[r][0]=0;
blupfactors[r][0]=0;
}
keeppreds_use[0]=0;
count=1;
}
else
{
printf("Reading details for %d predictors from %s\n", count, efffile);
value=(double)count/1024/1024/1024*8*num_scores*2;
if(value>1){printf("Warning, this requires %.1f Gb\n\n", value);}

preds=malloc(sizeof(char*)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);
centres=malloc(sizeof(double)*count);
factors=malloc(sizeof(double*)*num_scores);
indexer=malloc(sizeof(int)*count);
for(r=0;r<num_scores;r++){factors[r]=malloc(sizeof(double)*count);}

//open, check then skip header line
if((input=fopen(efffile,"r"))==NULL)
{printf("Error opening %s\n\n", efffile);exit(1);}
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading first element of %s\n\n", efffile);exit(1);}
if(strcmp(readstring,"Predictor")!=0&&strcmp(readstring,"SNP")!=0)
{printf("Error reading %s; first element should be \"Predictor\" or \"SNP\" (not %s)\n\n", efffile, readstring);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%s %s %s ", rs, ra1, ra2)!=3)
{printf("Error reading first three elements of Row %d of %s\n\n", j+2, efffile);exit(1);}
copy_string(preds,j,rs);
if(strlen(ra1)>1||strlen(ra2)>1)
{printf("Error reading Row %d of %s; alleles must be single characters (not %s and %s)\n\n", j+2, efffile, ra1, ra2);exit(1);}
al1[j]=ra1[0];al2[j]=ra2[0];
if(al1[j]==al2[j])
{printf("Error reading %s; both alleles for %s are the same (%c)\n\n", efffile, preds[j], al1[j]);exit(1);}

/*
if(fscanf(input, "%[0-9eE.+-] ", rc)==1){centres[j]=atof(rc);}
else
{
if(fscanf(input, "%s ", rc)!=1)
{printf("Error reading centre for Predictor %s in %s\n\n", preds[j], efffile);exit(1);}
if(strcmp(rc,"NA")==0)	//NAs allowed
else
{printf("Error reading %s; centre for Predictor %s is unrecognisable (%s)\n\n", efffile, preds[j], rc);exit(1);}
}
*/

if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading centre for Predictor %s in %s\n\n", preds[j], efffile);exit(1);}
if(strcmp(readstring,"NA")==0)
{centres[j]=-9999;}
else
{
if(sscanf(readstring, "%lf%c", centres+j, &readchar)!=1)
{printf("Error reading %s; centre for Predictor %s does not appear to be numeric (%s)\n\n", efffile, preds[j], readstring);exit(1);}
}

for(r=0;r<num_scores;r++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Effect %d for Predictor %s in %s\n\n", r+1, preds[j], efffile);exit(1);}
if(strcmp(readstring,"NA")==0)
{printf("Error reading %s; Effect %d for Predictor %s is NA\n\n", efffile, r+1, preds[j]);exit(1);}
else
{
if(sscanf(readstring, "%lf%c", factors[r]+j, &readchar)!=1)
{printf("Error reading %s; Effect %d for Predictor %s does not appear to be numeric (%s)\n\n", efffile, r+1, preds[j], readstring);exit(1);}
}
}
}	//end of j loop
fclose(input);

//set keeppreds_use and indexer
count2=find_strings(allpreds, num_preds, preds, count, keeppreds_use, indexer, NULL, efffile, predorder, NULL, 3);
if(count2==0){printf("Error, none of these are in the data\n\n");exit(1);}

if(type==0&&num_preds_use<num_preds)	//allow for filtering
{
usedpreds=malloc(sizeof(int)*num_preds);
for(j=0;j<num_preds;j++){usedpreds[j]=0;}
for(j=0;j<num_preds_use;j++){usedpreds[keeppreds[j]]=1;}

count3=0;
for(j=0;j<count2;j++)
{
if(usedpreds[keeppreds_use[j]]==1)
{
keeppreds_use[count3]=keeppreds_use[j];
indexer[count3]=indexer[j];
count3++;
}
}
if(count3==0){printf("Error, after filtering predictors, none remain\n\n");exit(1);}
free(usedpreds);
}
else{count3=count2;}

if(count3<count){printf("Warning, only %d of these are in the data\n", count3);}

//now load up - save a bit of memory by doing it one score at a time
for(r=0;r<num_scores;r++)
{
blupcentres[r]=malloc(sizeof(double)*count3);
blupfactors[r]=malloc(sizeof(double)*count3);

for(j=0;j<count3;j++)
{
j2=keeppreds_use[j];
j3=indexer[j];

if((al1[j3]!=allal1[j2]&&al1[j3]!=allal2[j2])||(al2[j3]!=allal1[j2]&&al2[j3]!=allal2[j2]))
{printf("Error, alleles for %s (%c and %c) not consistent with those in %s (%c and %c)\n", allpreds[j2], al1[j3], al2[j3], bimfile, allal1[j2], allal2[j2]);exit(1);}

//set centres and effects assuming alleles align, then turn if necessary (note that we might be missing centre)
blupcentres[r][j]=centres[j3];
blupfactors[r][j]=factors[r][j3];
if(al1[j3]!=allal1[j2])
{
if(blupcentres[r][j]!=-9999){blupcentres[r][j]=2-blupcentres[r][j];}
blupfactors[r][j]=-blupfactors[r][j];
}
}	//end of j loop

free(factors[r]);
}	//end of r loop

if(type==0)	//check centres
{
count2=0;for(j=0;j<count3;j++){count2+=(blupcentres[0][j]==-9999);}
if(count2>0)
{
if(count2==count3){printf("Warning, %s does not provide centres for any of the %d predictors, so these will be estimated from the data; these estimates will only be accurate if you have a reasonable sample size\n", efffile, count3);}
else{printf("Warning, %s is missing centres for %d of the %d predictors, so these will be estimated from the data; these estimates will only be accurate if you have a reasonable sample size\n\n", efffile, count2, count3);}
}
}

for(j=0;j<count;j++){free(preds[j]);}free(preds);
free(al1);free(al2);free(centres);free(factors);free(indexer);

//see if can reduce further
count=0;
for(j=0;j<count3;j++)
{
flag=0;
for(r=0;r<num_scores;r++){flag+=(blupfactors[r][j]!=0);}
if(flag>0)
{
for(r=0;r<num_scores;r++){blupcentres[r][count]=blupcentres[r][j];blupfactors[r][count]=blupfactors[r][j];}
keeppreds_use[count]=keeppreds_use[j];
count++;
}
}

if(count<count3){printf("After excluding predictors with zero effect sizes, %d remain\n", count);}

if(count==0)	//add one fake predictor back in
{
for(r=0;r<num_scores;r++){blupcentres[r][0]=0;blupfactors[r][0]=0;}
keeppreds_use[0]=0;
count=1;
}
}	//end of count
printf("\n");

free(rs);free(ra1);free(ra2);
return(count);
}	//end of read_scores

///////////////////////////

void read_rand(char *blupfile, double *bluprands, int num_kins, int ns, char **ids3)
{
int i, k, count, count2;

int *indexer;
char **wantids;

char readchar, readstring[500], *rs1, *rs2, *re;

FILE *input;

rs1=malloc(sizeof(char)*10000000);rs2=malloc(sizeof(char)*10000000);
re=malloc(sizeof(char)*10000000);


count=countrows(blupfile);
printf("Reading random effects for %d samples from %s\n", count, blupfile);
wantids=malloc(sizeof(char*)*count);
read_ids(blupfile, NULL, NULL, wantids, count, NULL, 0, 0);

//have already checked data available for all
indexer=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, NULL, indexer, blupfile, NULL, NULL, NULL, 3);
if(count2<count){printf("Doug error D422, %d\n\n", count2);exit(1);}

//set values to zero
for(k=0;k<num_kins;k++)
{
for(i=0;i<ns;i++){bluprands[i+k*ns]=0;}
}

if((input=fopen(blupfile,"r"))==NULL)
{printf("Error opening %s\n\n", blupfile);exit(1);}
for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", rs1, rs2)!=2)
{printf("Error reading IDs of Row %d of %s\n\n", i+2, blupfile);exit(1);}

for(k=0;k<num_kins;k++)	//save first effect, skip second
{
/*
if(fscanf(input, "%[0-9eE.+-] ", re)!=1)
{
if(fscanf(input, "%s ", re)!=1)
{printf("Error reading Column %d for Sample %s %s in %s\n\n", 3+2*k, rs1, rs2, blupfile);exit(1);}
if(strcmp(re,"NA")==0)
{printf("Error reading %s; Column %d for Sample %s %s is NA\n\n", blupfile, 3+2*k, rs1, rs2);exit(1);}
printf("Error reading %s; Column %d for Sample %s %s is unrecognisable (%s)\n\n", blupfile, 3+2*k, rs1, rs2, re);exit(1);
}
bluprands[indexer[i]+k*ns]=atof(re);
*/

if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Column %d for Sample %s %s in %s\n\n", 3+2*k, rs1, rs2, blupfile);exit(1);}
if(strcmp(readstring,"NA")==0)
{printf("Error reading %s; Column %d for Sample %s %s is NA\n\n", blupfile, 3+2*k, rs1, rs2);exit(1);}
else
{
if(sscanf(readstring, "%lf%c", bluprands+indexer[i]+k*ns, &readchar)!=1)
{printf("Error reading %s; Column %d for Sample %s %s does not appear to be numeric (%s)\n\n", blupfile, 3+2*k, rs1, rs2, readstring);exit(1);}
}

if(fscanf(input, "%s%c", readstring, &readchar)!=2)
{printf("Error reading Column %d for Sample %s %s in %s\n\n", 4+2*k, rs1, rs2, blupfile);exit(1);}
}

//skip to end of row
while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}	//end of i loop
fclose(input);
printf("\n");

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);

free(rs1);free(rs2);free(re);
}	//end of read_rand

////////

void read_pca(char *pcastem, double *bluprands, double *blupvects, int axes, int ns, char **ids3)
{
int i, k, count, count2;

int *indexer;
char **wantids;

char *rs1, *rs2;

char filename[500];
FILE *input;

rs1=malloc(sizeof(char)*10000000);rs2=malloc(sizeof(char)*10000000);


sprintf(filename,"%s.vect",pcastem);
count=countrows(filename);
printf("Reading eigen vectors for %d samples from %s\n", count, filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

//see which are available
indexer=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, NULL, indexer, filename, NULL, NULL, NULL, 3);
if(count2<count){printf("Error, only %d of these are in the data\n\n", count2);exit(1);}

//set values to zero
for(k=0;k<axes;k++)
{
for(i=0;i<ns;i++){bluprands[i+k*ns]=0;}
}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
for(i=0;i<count;i++)
{
if(fscanf(input, "%s %s ", rs1, rs2)!=2)
{printf("Error reading IDs of Row %d of %s\n\n", i+2, filename);exit(1);}
for(k=0;k<axes;k++)
{
if(fscanf(input, "%lf ", bluprands+indexer[i]+k*ns)!=1)
{printf("Error reading Column %d for Sample %s %s in %s, suggesting the file has been changed since creatino with \"--pca\"\n\n", 3+k, rs1, rs2, filename);exit(1);}
}
}	//end of i loop
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);free(indexer);

//read eigen values
sprintf(filename,"%s.values",pcastem);
printf("Reading eigen values from %s\n", filename);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
for(k=0;k<axes;k++)
{
if(fscanf(input, "%lf ", blupvects+k)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, filename);exit(1);}
}
fclose(input);
printf("\n");

free(rs1);free(rs2);
}	//end of read_pca

///////////////////////////

void write_blups(double **effects, int length, double **guesses, int ns, int num_kins, int num_scores, int num_tops, char **preds, char *al1, char *al2, char **ids1, char **ids2, char *outfile, double *resp, double missingvalue)
{
int i, j, k, r;

FILE *output, *output2, *output3, *output4;
char filename[500], filename2[500], filename3[500], filename4[500];


sprintf(filename,"%s.blup",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\tEffect\n");
for(j=0;j<length;j++)
{
if(effects[num_kins+num_scores+(num_tops>0)+2][j]==1)
{fprintf(output, "%s\t%c\t%c\t%.6f\t%.6e\n", preds[j], al1[j], al2[j], effects[num_kins+num_scores+(num_tops>0)+1][j], effects[num_kins+num_scores+(num_tops>0)][j]);}
}
fclose(output);

sprintf(filename2,"%s.pred",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"ID1\tID2\tPhenotype\tCovariates\tGenetics\tTotal\n");
for(i=0;i<ns;i++)
{
fprintf(output2, "%s\t%s\t", ids1[i], ids2[i]);
if(resp[i]!=missingvalue){fprintf(output2, "%.6f\t", resp[i]);}
else{fprintf(output2, "NA\t");}
fprintf(output2, "%.6f\t%.6f\t%.6f\n", guesses[num_kins+num_scores+(num_tops>0)+1][i], guesses[num_kins+num_scores+(num_tops>0)][i], guesses[num_kins+num_scores+(num_tops>0)+2][i]);
}
fclose(output2);

sprintf(filename3,"%s.blup.full",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Predictor\tA1\tA2\tCentre");
for(k=0;k<num_kins;k++){fprintf(output3,"\tKinship_%d", k+1);}
for(r=0;r<num_scores;r++){fprintf(output3,"\tRegion_%d", r+1);}
if(num_tops>0){fprintf(output3,"\tTops_Effect");}
fprintf(output3, "\n");
for(j=0;j<length;j++)
{
if(effects[num_kins+num_scores+(num_tops>0)+2][j]==1)
{
fprintf(output3, "%s\t%c\t%c\t%.6f", preds[j], al1[j], al2[j], effects[num_kins+num_scores+(num_tops>0)+1][j]);
for(k=0;k<num_kins+num_scores+(num_tops>0);k++){fprintf(output3, "\t%.6e", effects[k][j]);}
fprintf(output3, "\n");
}
}
fclose(output3);

sprintf(filename4,"%s.pred.full",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"ID1\tID2");
for(k=0;k<num_kins;k++){fprintf(output4,"\tKinship_%d", k+1);}
for(r=0;r<num_scores;r++){fprintf(output4,"\tRegion_%d", r+1);}
if(num_tops>0){fprintf(output4,"\tTop_Preds");}
fprintf(output4,"\n");
for(i=0;i<ns;i++)
{
fprintf(output4, "%s\t%s", ids1[i], ids2[i]);
for(k=0;k<num_kins+num_scores+(num_tops>0);k++){fprintf(output4, "\t%.6f", guesses[k][i]);}
fprintf(output4, "\n");
}
fclose(output4);

printf("Effect sizes saved in %s (and %s), with predictions in %s (and %s)\n\n", filename, filename3, filename2, filename4);
}	//end of write_blups

////////

void write_loads(double **effects, int length, double **guesses, int ns, int axes, double **blupcentres, char **preds, char *al1, char *al2, char **ids1, char **ids2, char *outfile)
{
int i, j, k;

FILE *output, *output2;
char filename[500], filename2[500];


sprintf(filename,"%s.load",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fprintf(output, "Predictor\tA1\tA2\tCentre\t");
for(k=0;k<axes;k++){fprintf(output, "Axes_%d\t", k+1);}
fprintf(output, "\n");
for(j=0;j<length;j++)
{
fprintf(output, "%s\t%c\t%c\t%.6e", preds[j], al1[j], al2[j], blupcentres[0][j]);
for(k=0;k<axes;k++){fprintf(output, "\t%e", effects[k][j]);}
fprintf(output,"\n");
}
fclose(output);

sprintf(filename2,"%s.proj",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"ID1\tID2");
for(k=0;k<axes;k++){fprintf(output2, "\tPCA_%d", k+1);}
fprintf(output2, "\n");
for(i=0;i<ns;i++)
{
fprintf(output2, "%s\t%s", ids1[i], ids2[i ]);
for(k=0;k<axes;k++){fprintf(output2, "\t%.6f", guesses[k][i]);}
fprintf(output2, "\n");
}
fclose(output2);

printf("Loadings saved in %s, with projections in %s\n\n", filename, filename2);
}	//end of write_loads

////////

void write_scores(double **guesses, double *indsds, int **nums, int ns, int num_scores, char **ids1, char **ids2, char *outfile, double *resp, double missingvalue, int type)
//type=0, save only scores, type=1, save also counts, type=2, save sds
{
int i, r;

FILE *output;
char filename[500];

sprintf(filename,"%s.profile",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

if(type==0)	//save just scores
{
fprintf(output,"ID1\tID2\tPhenotype\tCovariates");
for(r=0;r<num_scores;r++){fprintf(output,"\tProfile_%d\tBLANK", r+1);}
fprintf(output,"\n");

for(i=0;i<ns;i++)
{
fprintf(output, "%s\t%s\t", ids1[i], ids2[i]);
if(resp[i]!=missingvalue){fprintf(output, "%.6f\t", resp[i]);}
else{fprintf(output, "NA\t");}
fprintf(output, "%.6f", guesses[num_scores][i]);
for(r=0;r<num_scores;r++){fprintf(output, "\t%.6f\tNA", guesses[r][i]);}
fprintf(output, "\n");
}
}

if(type==1)	//save scores and counts
{
fprintf(output,"ID1\tID2\tPhenotype\tCovariates");
for(r=0;r<num_scores;r++){fprintf(output,"\tProfile_%d\tCount_%d", r+1, r+1);}
fprintf(output,"\n");

for(i=0;i<ns;i++)
{
fprintf(output, "%s\t%s\t", ids1[i], ids2[i]);
if(resp[i]!=missingvalue){fprintf(output, "%.6f\t", resp[i]);}
else{fprintf(output, "NA\t");}
fprintf(output, "%.6f", guesses[num_scores][i]);
for(r=0;r<num_scores;r++){fprintf(output, "\t%.6f\t%d", guesses[r][i], nums[r][i]);}
fprintf(output, "\n");
}
}

if(type==2)	//save scores and sds (only one score)
{
fprintf(output,"ID1\tID2\tPhenotype\tCovariates\tProfile\tSE\n");

for(i=0;i<ns;i++)
{
fprintf(output, "%s\t%s\t", ids1[i], ids2[i]);
if(resp[i]!=missingvalue){fprintf(output, "%.6f\t", resp[i]);}
else{fprintf(output, "NA\t");}
fprintf(output, "%.6f\t", guesses[num_scores][i]);
fprintf(output, "%.6f\t%.6f\n", guesses[0][i], indsds[i]);}
}

fclose(output);

if(num_scores==1){printf("Profile saved in %s\n\n", filename);}
else{printf("Profiles saved in %s\n\n", filename);}
}	//end of write_scores

///////////////////////////


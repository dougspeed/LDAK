/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for reading and writing kinships

///////////////////////////

int read_kins(char *kinstem, double *kins, float *kins_single, double scale, int ns, char **ids3, int type, int maxthreads)
//type=0 - read kins, type=1 - add kins, type=2 - read quiet, type=3 - add quiet, type=4 - read upper, type=5 - add upper
//in all cases, actually using kins*scale
{
int i, i2, count, count2, flag;
int thread, threadstart, threadend, *Mcurrent;
int *indexer, *indexer2;
float **Mdatatemp;
char **wantids;

char filename[500];
FILE **Minput;


Mcurrent=malloc(sizeof(int)*maxthreads);
Minput=malloc(sizeof(FILE *)*maxthreads);
Mdatatemp=malloc(sizeof(float*)*maxthreads);

if(type==0||type==1){printf("Reading kinship matrix with stem %s\n",kinstem);}

//first get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*ns);
indexer2=malloc(sizeof(int)*ns);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count2!=ns){printf("Error finding indexes read_kin, please tell Doug\n\n");exit(1);}

//see if indexer2 is in order (when considering only upper triangle, computation will be faster if it is)
flag=0;
for(i=1;i<ns;i++)
{
if(indexer2[i]<indexer2[i-1]){flag=1;break;}
}

//check size
sprintf(filename, "%s.grm.bin", kinstem);
if((Minput[0]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

fseeko(Minput[0], 0, SEEK_END);
if(ftello(Minput[0])!=(off_t)sizeof(float)*count*(count+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*count*(count+1)/2, ftello(Minput[0]));exit(1);}
fclose(Minput[0]);

//now can read
#pragma omp parallel for private(thread,threadstart,threadend,i,i2) schedule (static, 1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=pow((double)thread/maxthreads,.5)*ns;
threadend=pow((double)(thread+1)/maxthreads,.5)*ns;

if((Minput[thread]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

fseeko(Minput[thread], 0, SEEK_SET);
Mcurrent[thread]=0;

//read one "row" at a time - this works because indexer is monotonic
Mdatatemp[thread]=malloc(sizeof(float)*count);
for(i=threadstart;i<threadend;i++)
{
if((type==0||type==1)&&count>40000&&thread==0&&i%20000==0&&i*maxthreads<ns)
{printf("Reading kinships for Sample %d out of %d\n", i*maxthreads+1, ns);}

if(indexer[i]!=Mcurrent[thread])	//get to start of Row indexer[i]
{fseeko(Minput[thread], (off_t)sizeof(float)*indexer[i]*(indexer[i]+1)/2, SEEK_SET);}
if(fread(Mdatatemp[thread], sizeof(float), indexer[i]+1, Minput[thread])!=indexer[i]+1)
{printf("Error reading Row %d of %s\n\n", indexer[i]+1, filename);exit(1);};
Mcurrent[thread]=indexer[i]+1;

if(kins!=NULL)	//saving as doubles
{
if(type==0||type==2)	//just read kins
{
for(i2=0;i2<i;i2++)
{
kins[(size_t)indexer2[i2]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i2]]*scale;
kins[(size_t)indexer2[i]*ns+indexer2[i2]]=Mdatatemp[thread][indexer[i2]]*scale;
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i]]*scale;
}

if(type==1||type==3)	//add kins
{
for(i2=0;i2<i;i2++)
{
kins[(size_t)indexer2[i2]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i2]]*scale;
kins[(size_t)indexer2[i]*ns+indexer2[i2]]+=Mdatatemp[thread][indexer[i2]]*scale;
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i]]*scale;
}

if(type==4)	//read to upper triangle
{
if(flag==0)	//indexer2 in order
{
for(i2=0;i2<i;i2++)
{
kins[(size_t)indexer2[i]*ns+indexer2[i2]]=Mdatatemp[thread][indexer[i2]]*scale;
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i]]*scale;
}
else	//indexer2 not in order
{
for(i2=0;i2<i;i2++)
{
if(indexer2[i2]<indexer2[i]){kins[(size_t)indexer2[i]*ns+indexer2[i2]]=Mdatatemp[thread][indexer[i2]];}
else{kins[(size_t)indexer2[i2]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i2]];}
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i]]*scale;
}
}	//end of type=4

if(type==5)	//add to upper triangle
{
if(flag==0)	//indexer2 in order
{
for(i2=0;i2<i;i2++)
{
kins[(size_t)indexer2[i]*ns+indexer2[i2]]+=Mdatatemp[thread][indexer[i2]]*scale;
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i]]*scale;
}
else	//indexer2 not in order
{
for(i2=0;i2<i;i2++)
{
if(indexer2[i2]<indexer2[i]){kins[(size_t)indexer2[i]*ns+indexer2[i2]]+=Mdatatemp[thread][indexer[i2]];}
else{kins[(size_t)indexer2[i2]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i2]];}
}
kins[(size_t)indexer2[i]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i]]*scale;
}
}	//end of type=5
}
else	//saving as floats
{
if(type==0||type==2)	//just read kins
{
for(i2=0;i2<i;i2++)
{
kins_single[(size_t)indexer2[i2]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i2]]*scale;
kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]=Mdatatemp[thread][indexer[i2]]*scale;
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i]]*scale;
}

if(type==1||type==3)	//add kins
{
for(i2=0;i2<i;i2++)
{
kins_single[(size_t)indexer2[i2]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i2]]*scale;
kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]+=Mdatatemp[thread][indexer[i2]]*scale;
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i]]*scale;
}

if(type==4)	//read to upper triangle
{
if(flag==0)	//indexer2 in order
{
for(i2=0;i2<i;i2++)
{
kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]=Mdatatemp[thread][indexer[i2]]*scale;
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i]]*scale;
}
else	//indexer2 not in order
{
for(i2=0;i2<i;i2++)
{
if(indexer2[i2]<indexer2[i]){kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]=Mdatatemp[thread][indexer[i2]];}
else{kins_single[(size_t)indexer2[i2]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i2]];}
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]=Mdatatemp[thread][indexer[i]]*scale;
}
}	//end of type=4

if(type==5)	//add to upper triangle
{
if(flag==0)	//indexer2 in order
{
for(i2=0;i2<i;i2++)
{
kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]+=Mdatatemp[thread][indexer[i2]]*scale;
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i]]*scale;
}
else	//indexer2 not in order
{
for(i2=0;i2<i;i2++)
{
if(indexer2[i2]<indexer2[i]){kins_single[(size_t)indexer2[i]*ns+indexer2[i2]]+=Mdatatemp[thread][indexer[i2]];}
else{kins_single[(size_t)indexer2[i2]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i2]];}
}
kins_single[(size_t)indexer2[i]*ns+indexer2[i]]+=Mdatatemp[thread][indexer[i]]*scale;
}
}	//end of type=5
}	//end of single

}	//end of i loop

fclose(Minput[thread]);
free(Mdatatemp[thread]);
}	//end of thread loop

if(type==0||type==1){printf("\n");}

free(Mcurrent);free(Minput);free(Mdatatemp);
for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);

return(0);
}	//end of read_kins

////////

double read_kin_trace(char *kinstem, int ns, char **ids3, double *kin_diags, int type, int diagonal, int maxthreads)
//type=0 - normal, type=1 - quiet
{
int i, count, count2;
int thread, threadstart, threadend;
int *indexer, *indexer2;
double sum, *Msum, *datatemp;

char **wantids;
float readfloat;

char filename[500];
FILE **Minput;


if(type==0){printf("Reading trace of kinship matrix with stem %s\n",kinstem);}

//first get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count2!=ns){printf("Error find indexes read_kin_trace, please tell Doug\n\n");exit(1);}

if(diagonal==0)	//normal case, have full kinship matrix
{
Minput=malloc(sizeof(FILE *)*maxthreads);
Msum=malloc(sizeof(double*)*maxthreads);

//check size
sprintf(filename, "%s.grm.bin", kinstem);
if((Minput[0]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

fseeko(Minput[0], 0, SEEK_END);
if(ftello(Minput[0])!=(off_t)sizeof(float)*count*(count+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*count*(count+1)/2, ftello(Minput[0]));exit(1);}
fclose(Minput[0]);

//now tally trace
#pragma omp parallel for private(thread,threadstart,threadend,i,readfloat) schedule (static, 1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=thread*((ns-1)/maxthreads+1);
threadend=(thread+1)*((ns-1)/maxthreads+1);
if(threadend>ns){threadend=ns;}

if((Minput[thread]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

Msum[thread]=0;
for(i=threadstart;i<threadend;i++)
{
//get to ith diagonal element
fseeko(Minput[thread], (off_t)sizeof(float)*indexer[i]*(indexer[i]+3)/2, SEEK_SET);
if(fread(&readfloat, sizeof(float), 1, Minput[thread])!=1)
{printf("Error reading Row %d of %s\n\n", i+1, filename);exit(1);};
Msum[thread]+=readfloat;
if(kin_diags!=NULL){kin_diags[indexer2[i]]=readfloat;}
}

fclose(Minput[thread]);
}

sum=0;for(thread=0;thread<maxthreads;thread++){sum+=Msum[thread];}

free(Minput);free(Msum);
}
else	//read diagonal values
{
datatemp=malloc(sizeof(double)*count);
sprintf(filename, "%s.grm.diag", kinstem);
read_values(filename, datatemp, count, NULL, 1, 0, 0);

sum=0;for(i=0;i<ns;i++){sum+=datatemp[indexer[i]];}
if(kin_diags!=NULL)
{
for(i=0;i<ns;i++){kin_diags[indexer2[i]]=datatemp[indexer[i]];}
}

free(datatemp);
}

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);

return(sum/ns);
}

///////////////////////////

int read_details(char *filename, char **kpreds, int *kindex, double *kcentres, double *kmults, double *kweights, char *kal1, char *kal2, double *kexps, int length)
{
int j;

int readint;
float readfloat;
char readchar, *rs;

FILE *input;

rs=malloc(sizeof(char)*10000000);


if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}

for(j=0;j<length;j++)
{
if(fscanf(input, "%s %d %lf %lf %lf %c %c %lf %f ", rs, &readint, kcentres+j, kmults+j, kweights+j, kal1+j, kal2+j, kexps+j, &readfloat)!=9)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation\n\n", j+2, filename);exit(1);}
copy_string(kpreds,j,rs);
kindex[j]=readint-1;
}
fclose(input);

free(rs);
return(0);
}	//end of read_details

///////////////////////////

int read_kinsgz(char *kinstem, float *kins_single, int ns, char **ids3, int partial)
{
size_t scount;
int i, i2, readi, readi2, count, count2, found, size, offset;
int *indexer, *indexer2, *rec;

double readdouble, *datatemp;
char **wantids, *gzbuffer;

char filename[500];
gzFile inputgz;


//size of buffer
size=1000;

printf("Reading %s.grm.id and %s.grm.gz\n\n", kinstem, kinstem);

//get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, filename, NULL, NULL, NULL, 3);
if(count2<ns){printf("Error 83E, please tell Doug %d %d\n", count2, ns);exit(1);}

//get ready to read
sprintf(filename, "%s.grm.gz", kinstem);
if((inputgz=gzopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

gzbuffer=malloc(sizeof(char)*size);

if(partial==0)	//expect full kinship matrix
{
datatemp=malloc(sizeof(double)*count);

//read set of i rows at a time
found=0;scount=0;
for(i=0;i<count;i++)
{
for(i2=0;i2<=i;i2++)
{
if(gzgets(inputgz,gzbuffer,size)==NULL)
{printf("Error reading %s; there is no Row %jd\n\n", filename, scount);exit(1);}

if(strlen(gzbuffer)==size-1)
{printf("Error reading %s; Row %jd is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", filename, scount+1, (int)strlen(gzbuffer), size);exit(1);}

//check starting elements
if(sscanf(gzbuffer, "%d %d %n", &readi, &readi2, &offset)!=2)
{printf("Error reading start of Row %jd of %s\n\n", scount+1, filename);exit(1);}
if(readi!=i+1||readi2!=i2+1)
{printf("Error reading %s; Row %jd starts %d %d (not %d %d)\nIf this is a partial kinship matrix (only contains values for a subset of samples), please add \"--partial YES\"\n\n", filename, scount+1, readi, readi2, i+1, i2+1);exit(1);}

//then save kinships
if(sscanf(gzbuffer+offset, "%lf %lf ", &readdouble, datatemp+i2)!=2)
{printf("Error reading end of Row %jd of %s\n\n", scount+1, filename);exit(1);}
scount++;
}

if(i==indexer[found])	//using i, so save
{
for(i2=0;i2<found;i2++)
{
kins_single[(size_t)indexer2[i2]*ns+indexer2[found]]=datatemp[indexer[i2]];
kins_single[(size_t)indexer2[found]*ns+indexer2[i2]]=datatemp[indexer[i2]];
}
kins_single[(size_t)indexer2[found]*ns+indexer2[found]]=datatemp[indexer[found]];
found++;
if(found==ns){break;}
}
}

free(datatemp);
}	//end of partial=0
else	//not sure how many kinships will have
{
//rec indicates which i to save
rec=malloc(sizeof(int)*ns);
for(i=0;i<ns;i++){rec[i]=-1;}
for(i=0;i<count;i++){rec[indexer[i]]=i;}

//first set kins_single to identify matrix
for(i2=0;i2<ns;i2++)
{
for(i=0;i<ns;i++){kins_single[(size_t)i2*ns+i]=0;}
kins_single[(size_t)i2*ns+i2]=1;
}

//read one row at time until reach end
scount=0;
while(1)
{
if(gzgets(inputgz,gzbuffer,size)==NULL){printf("Finished reading %s after %jd rows\n\n", filename, scount);break;}

if(strlen(gzbuffer)==size-1)
{printf("Error reading %s; Row %jd is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", filename, scount+1, (int)strlen(gzbuffer), size);exit(1);}

//get elements and kinships
if(sscanf(gzbuffer, "%d %d %lf %lf ", &readi, &readi2, &readdouble, &readdouble)!=4)
{printf("Error reading Row %jd of %s\n\n", scount+1, filename);exit(1);}

i=readi-1;
i2=readi2-1;

if(rec[i]!=-1&&rec[i2]!=-1)	//save this element
{
kins_single[(size_t)indexer2[rec[i]]*ns+indexer2[rec[i2]]]=readdouble;
kins_single[(size_t)indexer2[rec[i2]]*ns+indexer2[rec[i]]]=readdouble;
}

scount++;
}

free(rec);
}

gzclose(inputgz);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(gzbuffer);

return(0);
}	//end of read_kinfilegz

////////

int read_kinsraw(char *kinstem, float *kins_single, int ns, char **ids3)
{
int i, i2, count, count2, found;
int *indexer, *indexer2;

double *datatemp; 
char **wantids;

char filename[500];
FILE *input;


printf("Reading %s.grm.id and %s.grm.raw\n\n", kinstem, kinstem);

//get indexes of individuals we want
sprintf(filename, "%s.grm.id", kinstem);
count=countrows(filename);
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
count2=find_strings(wantids, count, ids3, ns, indexer, indexer2, filename, NULL, NULL, NULL, 3);
if(count2<ns){printf("Error 83D, please tell Doug %d %d\n", count2, ns);exit(1);}

//get ready to read
sprintf(filename, "%s.grm.raw", kinstem);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

datatemp=malloc(sizeof(double)*count);
found=0;
for(i=0;i<count;i++)
{
for(i2=0;i2<count;i2++)
{
if(fscanf(input, "%lf\n", datatemp+i2)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", i2+1, i+1, filename);exit(1);}
}

if(i==indexer[found])	//using i, so save
{
for(i2=0;i2<ns;i2++){kins_single[(size_t)indexer2[found]*ns+indexer2[i2]]=datatemp[indexer[i2]];}
found++;
if(found==ns){break;}
}
}
fclose(input);

for(i=0;i<count;i++){free(wantids[i]);}free(wantids);
free(indexer);free(indexer2);free(datatemp);

return(0);
}	//end of read_kinfileraw

///////////////////////////

int write_kins(char *outfile, double *kins, float *kins_single, int ns, char **ids1, char **ids2, int kindetails, char **preds, int *keeppreds_use, double *centres, double *mults, double *weights, char *al1, char *al2, double *exps, int length, char *datafile, double power, int kingz, int kinraw, int type)
//type=0 - quiet, type=1 - normal, type=2 - null for association testing, type=3 - trun/pca/square/gxemm kinships, type=4 - convert, type=5 - diagonal
{
int i, i2, j, count;
float *datatemp;
double sum, sumsq, value, value2, denom;

char filename[500], filename2[500], filename3[500], filename4[500], filename5[500], cmd[500];
FILE *output, *output2, *output3, *output4, *output5;


if(type==1||type==2)	//get denom
{
sum=0;
for(i=0;i<ns;i++)
{
if(kins!=NULL){sum+=kins[(size_t)i*ns+i];}
else{sum+=kins_single[(size_t)i*ns+i];}
}
denom=sum/ns;
}
else{denom=1;}

//write kins to grm.bin
sprintf(filename,"%s.grm.bin", outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

datatemp=malloc(sizeof(float)*ns);
sum=0;sumsq=0;
for(i2=0;i2<ns;i2++)
{
if(type!=5)	//have full matrix
{
for(i=0;i<=i2;i++)
{
if(kins!=NULL){datatemp[i]=kins[(size_t)i2*ns+i]/denom;}
else{datatemp[i]=kins_single[(size_t)i2*ns+i]/denom;}
}
}
else	//only have diagonal
{
for(i=0;i<i2;i++){datatemp[i]=0;}
if(kins!=NULL){datatemp[i2]=kins[i2];}
else{datatemp[i2]=kins_single[i2];}
}
fwrite(datatemp, sizeof(float), i2+1, output);
for(i=0;i<i2;i++){sum+=datatemp[i];sumsq+=pow(datatemp[i],2);}
}
fclose(output);
free(datatemp);

if(type==5)	//will also write just diagonal
{
sprintf(filename,"%s.grm.diag", outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
for(i=0;i<ns;i++)
{
if(kins!=NULL){fprintf(output, "%.6f\n", kins[i]);}
else{fprintf(output, "%.6f\n", kins_single[i]);}
}
fclose(output);
}

if(kindetails==1)	//write details and adjustments
{
//get sum of exps
value=0;
for(j=0;j<length;j++)
{
if(mults[j]!=-9999){value+=exps[j];}
}

sprintf(filename2,"%s.grm.details", outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor Index Centre Scaling Weight A1 A2 Exp_Heritability Share\n");
count=0;value2=0;
for(j=0;j<length;j++)
{
if(mults[j]!=-9999)
{
fprintf(output2, "%s %d %.6f %.6f %.6f %c %c %.6f %.6e\n", preds[j], keeppreds_use[j]+1, centres[j], mults[j], weights[j], al1[j], al2[j], exps[j], exps[j]/value);
count++; 
value2+=weights[j];
}
}
fclose(output2);

sprintf(filename3,"%s.grm.adjust", outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

if(datafile!=NULL){fprintf(output3, "Datafile %s\n", datafile);}
else{fprintf(output3, "Multiple Datafiles\n");}
if(power!=-9999){fprintf(output3, "Power %.6f\n", power);}
else{fprintf(output3, "Multiple Powers\n");}

fprintf(output3, "Num_Preds %d\nSum_Weights %.4f\nDenominator %.4f\nOff_Diag_Variance %e\n", count, value2, denom, 2*sumsq/ns/(ns-1)-pow(2*sum/ns/(ns-1),2));
fclose(output3);
}

//write IDs
sprintf(filename4,"%s.grm.id", outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
for(i=0;i<ns;i++){fprintf(output4, "%s %s\n", ids1[i], ids2[i]);}
fclose(output4);

if(type==1||type==4)
{
printf("Kinship matrix saved in files with prefix %s", outfile);
}

////////

if(kingz==1)
{
sprintf(filename5, "%s.grm", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<ns;i++)
{
for(i2=0;i2<=i;i2++)
{
if(type!=5)	//have full matrix
{
if(kins!=NULL){fprintf(output5, "%d %d %d %.6f\n", i+1, i2+1, count, kins[(size_t)i*ns+i2]/denom);}
else{fprintf(output5, "%d %d %d %.6f\n", i+1, i2+1, count, kins_single[(size_t)i*ns+i2]/denom);}
}
else	//only have diagonal
{
if(i!=i2){fprintf(output5, "%d %d %d 0\n", i+1, i2+1, count);}
else
{
if(kins!=NULL){fprintf(output5, "%d %d %d %.6f\n", i+1, i2+1, count, kins[i]/denom);}
else{fprintf(output5, "%d %d %d %.6f\n", i+1, i2+1, count, kins_single[i]/denom);}
}
}
}
}
fclose(output5);

//try to compress the kinship file
sprintf(cmd, "gzip -f %s", filename5);
system(cmd);
if(type==1||type==4||type==5){printf(", with gzipped version saved in %s.gz", filename5);}
}

if(kinraw==1)
{
sprintf(filename5, "%s.grm.raw", outfile);
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}
for(i=0;i<ns;i++)
{
for(i2=0;i2<ns;i2++)
{
if(type!=5)	//have full matrix
{
if(kins!=NULL){fprintf(output5,"%.6f ", kins[(size_t)i*ns+i2]/denom);}
else{fprintf(output5,"%.6f ", kins_single[(size_t)i*ns+i2]/denom);}
}
else	//only have diagonal
{
if(i!=i2){fprintf(output5,"0 ");}
else
{
if(kins!=NULL){fprintf(output5,"%.6f ", kins[i]/denom);}
else{fprintf(output5,"%.6f ", kins_single[i]/denom);}
}
}
}
fprintf(output5,"\n");
}
fclose(output5);
if(type==1||type==4||type==5)
{
if(kingz==1){printf(" and plain text version saved in %s", filename5);}
else{printf(", with plain text version saved in %s", filename5);}
}
}

if(type==1||type==4){printf("\n\n");}

return(0);
}	//end of write_kins

///////////////////////////

int write_eigen(char *outfile, double *U, double *E, int ns, char **ids1, char **ids2, char *kinstem, int eigenraw)
{
int i, i2;
float *datatemp;

char filename[500], filename2[500], filename3[500], filename4[500];
FILE *output, *output2, *output3, *output4;


sprintf(filename,"%s.eigen.bin",outfile);
if((output=fopen(filename,"wb"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename);exit(1);}

datatemp=malloc(sizeof(float)*ns);
for(i=0;i<ns;i++){datatemp[i]=E[i];}
fwrite(datatemp, sizeof(float), ns, output);
for(i2=0;i2<ns;i2++)
{
for(i=0;i<ns;i++){datatemp[i]=U[(size_t)i2*ns+i];}
fwrite(datatemp, sizeof(float), ns, output);
}
fclose(output);
free(datatemp);

sprintf(filename2,"%s.eigen.id",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename2);exit(1);}
for(i=0;i<ns;i++){fprintf(output2,"%s %s\n", ids1[i], ids2[i]);}
fclose(output2);

sprintf(filename3,"%s.eigen.root",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename3);exit(1);}
fprintf(output3,"Kinship %s\n", kinstem);
fclose(output3);

if(eigenraw==1)
{
sprintf(filename4,"%s.eigen.raw",outfile);
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n", filename4);exit(1);}
for(i=ns-1;i>=0;i--){fprintf(output4, "%.6f ", E[i]);}
fprintf(output4, "\n");
for(i2=ns-1;i2>=0;i2--)
{
for(i=0;i<ns;i++){fprintf(output4, "%.6f ", U[(size_t)i2*ns+i]);}
fprintf(output4, "\n");
}
fclose(output4);
}

printf("Eigen-decomposition saved in files with prefix %s", outfile);
if(eigenraw==1)
{printf(", with plain text version saved in %s", filename4);}
printf("\n\n");

return(0);
}	//end of write_eigen

////////

int read_eigens(char *eigenfile, double *U, double *E, int ns, char **ids3, int axes, int maxthreads)
{
int i, i2, count, count2;
int thread, threadstart, threadend;
int *indexer;
float **Mdatatemp;
char **wantids;

char filename[500];
FILE **Minput;


Minput=malloc(sizeof(FILE *)*maxthreads);
Mdatatemp=malloc(sizeof(float*)*maxthreads);

printf("Reading eigen-decomposition with stem %s\n", eigenfile);

//first get indexes of individuals we want
sprintf(filename, "%s.eigen.id", eigenfile);
count=countrows(filename);
if(count!=ns)
{printf("Error, the decomposition is not suitable (%s contains %d samples, not %d)\nWhen using \"--decompose\" make sure you provide the same sample filterings as now (and if there are missing phenotypes, also use \"--pheno\")\n\n", filename, count, ns);exit(1);}
wantids=malloc(sizeof(char*)*count);
read_ids(filename, NULL, NULL, wantids, count, NULL, 0, 0);

indexer=malloc(sizeof(int)*count);
count2=find_strings(ids3, ns, wantids, count, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count2==0)
{printf("Error, the decomposition is not suitable (%s contains none of the %d samples)\nWhen using \"--decompose\" make sure you provide the same sample filterings as now (and if there are missing phenotypes, also use \"--pheno\")\n\n", filename, ns);exit(1);}
if(count2<ns)
{printf("Error, the decomposition is not suitable (%s contains only %d of the %d samples)\nWhen using \"--decompose\" make sure you provide the same sample filterings as now (and if there are missing phenotypes, also use \"--pheno\")\n\n", filename, count2, ns);exit(1);}

//check size then read eigenvalues
sprintf(filename,"%s.eigen.bin",eigenfile);
if((Minput[0]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

fseeko(Minput[0], 0, SEEK_END);
if(ftello(Minput[0])!=(off_t)sizeof(float)*(count+1)*count)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*(count+1)*count, ftello(Minput[0]));exit(1);}

for(thread=0;thread<maxthreads;thread++){Mdatatemp[thread]=malloc(sizeof(float)*count);}

if(axes==-9999)	//read all eigenvalues
{
fseeko(Minput[0], 0, SEEK_SET);
if(fread(Mdatatemp[0], sizeof(float), count, Minput[0])!=count)
{printf("Error reading eigenvalues from %s\n\n", filename);exit(1);}
for(i=0;i<count;i++){E[i]=Mdatatemp[0][i];}
}
else	//read last axes eigenvalues
{
fseeko(Minput[0], (off_t)sizeof(float)*(count-axes), SEEK_SET);
if(fread(Mdatatemp[0], sizeof(float), axes, Minput[0])!=axes)
{printf("Error reading eigenvalues from %s\n\n", filename);exit(1);}
for(i=0;i<axes;i++){E[i]=Mdatatemp[0][i];}
}
fclose(Minput[0]);

if(axes==-9999)	//read all eigenvectors
{
#pragma omp parallel for private(thread,threadstart,threadend,i,i2) schedule (static, 1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=thread*ns/maxthreads;
threadend=(thread+1)*ns/maxthreads;
if(threadend>ns){threadend=ns;}

if((Minput[thread]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

for(i=threadstart;i<threadend;i++)
{
fseeko(Minput[thread], (off_t)sizeof(float)*(1+i)*count, SEEK_SET);
if(fread(Mdatatemp[thread], sizeof(float), count, Minput[thread])!=count)
{printf("Error reading values for Eigenvector %d from %s\n\n", i+1, filename);exit(1);}
for(i2=0;i2<count;i2++){U[(size_t)i*count+indexer[i2]]=Mdatatemp[thread][i2];}
}

fclose(Minput[thread]);
}
}
else	//read last axes eigenvalues
{
if((Minput[0]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}

fseeko(Minput[0], (off_t)sizeof(float)*(1+count-axes)*count, SEEK_SET);
for(i=0;i<axes;i++)
{
if(fread(Mdatatemp[0], sizeof(float), count, Minput[0])!=count)
{printf("Error reading values for Eigenvector %d from %s\n\n", count-axes+i+1, filename);exit(1);}
for(i2=0;i2<count;i2++){U[(size_t)i*count+indexer[i2]]=Mdatatemp[0][i2];}
}
fclose(Minput[0]);
}

for(thread=0;thread<maxthreads;thread++){free(Mdatatemp[thread]);}
free(Minput);free(Mdatatemp);
for(i=0;i<count;i++){free(wantids[i]);};free(wantids);
free(indexer);

printf("\n");
return(0);
}	//end of read_eigenfile

///////////////////////////


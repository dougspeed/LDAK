/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for reading in data

///////////////////////////

void read_bed_full(char *bedfile, unsigned char **data_char, int num_samples_use, int *keepsamps, int length, int *keeppreds, int num_samples, int num_preds, int maxthreads)
{
int i, j,count, count2, rowlength;
int thread, threadstart, threadend, threadlength, *Mcurrent;
unsigned char startchars[3], **Mrowchars, twobits;

FILE **Minput;


rowlength=(num_samples-1)/4+1;
threadlength=(length-1)/maxthreads+1;

Mcurrent=malloc(sizeof(int)*maxthreads);
Minput=malloc(sizeof(FILE *)*maxthreads);
Mrowchars=malloc(sizeof(unsigned char*)*maxthreads);
for(thread=0;thread<maxthreads;thread++){Mrowchars[thread]=malloc(sizeof(unsigned char)*rowlength);}

//check the length and that the first three digits ok
if((Minput[0]=fopen(bedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",bedfile);exit(1);}

fseeko(Minput[0], 0, SEEK_END);
if(ftello(Minput[0])!=(off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3)
{printf("Error reading %s; should have size %jd (%d samples x %d predictors), but instead has size %jd\n\n", bedfile, (off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3, num_samples, num_preds, ftello(Minput[0]));exit(1);}

fseeko(Minput[0], 0, SEEK_SET);
if(fread(startchars, sizeof(unsigned char), 3, Minput[0])!=3)
if(startchars[0]!=108||startchars[1]!=27)
{printf("Error reading %s; does not appear to be in binary PLINK format\n\n", bedfile);exit(1);}
if(startchars[2]!=1)
{printf("Error reading %s; can only read in SNP-major mode\n\n", bedfile);exit(1);}
fclose(Minput[0]);

//now can read
#pragma omp parallel for private(thread,threadstart,threadend,j,i,count,count2,twobits) schedule (static, 1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=thread*threadlength;
threadend=(thread+1)*threadlength;
if(threadend>length){threadend=length;}

if((Minput[thread]=fopen(bedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",bedfile);exit(1);}
if(fseeko(Minput[thread], sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error finding the first predictor in %s\n\n", bedfile);exit(1);}
Mcurrent[thread]=0;

for(j=threadstart;j<threadend;j++)
{
if(thread==0&&j%20000==0&&j*maxthreads<length)
{printf("Reading predictor %d out of %d\n", j*maxthreads+1, length);}

if(keeppreds[j]!=Mcurrent[thread])
{
if(fseeko(Minput[thread], (off_t)sizeof(unsigned char)*rowlength*keeppreds[j]+sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d\n\n", bedfile, keeppreds[j]+1);exit(1);}
}
if(fread(Mrowchars[thread], sizeof(unsigned char), rowlength, Minput[thread])!=rowlength)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, bedfile);exit(1);}
Mcurrent[thread]=keeppreds[j]+1;

for(i=0;i<num_samples_use;i++)
{
if(i%4==0){data_char[j][i/4]=0;}
count=keepsamps[i]/4;
count2=keepsamps[i]%4;
twobits=(Mrowchars[thread][count]>>(2*count2)) & 3;
switch(twobits)
{
case 3: break;
case 2: data_char[j][i/4]+=(1<<(2*(i%4)));break;
case 0: data_char[j][i/4]+=(2<<(2*(i%4)));break;
case 1: data_char[j][i/4]+=(3<<(2*(i%4)));break;
}
}
}	//end of j loop

fclose(Minput[thread]);
}
printf("\n");

free(Mcurrent);free(Minput);
for(thread=0;thread<maxthreads;thread++){free(Mrowchars[thread]);}free(Mrowchars);
}	//end of read_bed_full

////////

int read_bed_fly(char *bedfile, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, double *bedbytes)
{
int i, j, current, rowlength;
unsigned char startchars[3], *rowchars;
double *readdoubles;

FILE *input;


rowlength=(num_samples-1)/4+1;
rowchars=malloc(sizeof(unsigned char)*rowlength);
readdoubles=malloc(sizeof(double)*(num_samples+4));

if((input=fopen(bedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",bedfile);exit(1);}

//check the length and that the first three digits ok
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3)
{printf("Error reading %s; should have size %jd (%d samples x %d predictors), but instead has size %jd\n\n", bedfile, (off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3, num_samples, num_preds, ftello(input));exit(1);}

fseeko(input, 0, SEEK_SET);
if(fread(startchars, sizeof(unsigned char), 3, input)!=3)
{printf("Error reading first three values of %s\n\n", bedfile);exit(1);}
if(startchars[0]!=108||startchars[1]!=27)
{printf("Error reading %s; does not appear to be in binary PLINK format\n\n", bedfile);exit(1);}
if(startchars[2]!=1)
{printf("Error reading %s; can only read in SNP-major mode\n\n", bedfile);exit(1);}
current=0;

for(j=start;j<end;j++)
{
if(keeppreds[j]!=current)
{
if(fseeko(input, (off_t)sizeof(unsigned char)*rowlength*keeppreds[j]+sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d %d\n\n", bedfile, j, keeppreds[j]+1);exit(1);}
}
if(fread(rowchars, sizeof(unsigned char), rowlength, input)!=rowlength)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, bedfile);exit(1);}
current=keeppreds[j]+1;

for(i=0;i<num_samples;i+=4){memcpy(readdoubles+i,bedbytes+4*rowchars[i/4],sizeof(double)*4);}
for(i=0;i<num_samples_use;i++){data[(size_t)(j-start)*num_samples_use+i]=readdoubles[keepsamps[i]];}

/*
for(i=0;i<num_samples_use;i++)
{
count=keepsamps[i]/4;
count2=keepsamps[i]%4;
twobits=(rowchars[count]>>(2*count2)) & 3;

switch (twobits)
{
case 3: data[(size_t)(j-start)*num_samples_use+i]=0;break;
case 2: data[(size_t)(j-start)*num_samples_use+i]=1;break;
case 0: data[(size_t)(j-start)*num_samples_use+i]=2;break;
case 1: data[(size_t)(j-start)*num_samples_use+i]=missingvalue;
}
}
*/
}	//end of j loop

fclose(input);
free(rowchars);
free(readdoubles);

return(0);
}	//end of read_bed_fly

/////////////////////

int read_sped_fly(char *spedfile, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, double missingvalue, double threshold, int nonsnp)
{
int i, j, current;
float *rowfloats;
double value;

FILE *input;


rowfloats=malloc(sizeof(float)*num_samples);

if((input=fopen(spedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",spedfile);exit(1);}

//check the length
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*num_samples*num_preds)
{printf("Error reading %s; should have size %jd (%d samples x %d predictors), but instead has size %jd\n\n", spedfile, (off_t)sizeof(float)*num_samples*num_preds, num_samples, num_preds, ftello(input));exit(1);}
fseeko(input, 0, SEEK_SET);
current=0;

for(j=start;j<end;j++)
{
if(keeppreds[j]!=current)
{
if(fseeko(input, (off_t)sizeof(float)*num_samples*keeppreds[j], SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d\n\n", spedfile, keeppreds[j]+1);exit(1);}
}
if(fread(rowfloats, sizeof(float), num_samples, input)!=num_samples)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, spedfile);exit(1);}
current=keeppreds[j]+1;

if(nonsnp==0)	//check values between 0 and 2
{
for(i=0;i<num_samples;i++)
{
if(rowfloats[i]!=missingvalue&&(rowfloats[i]<0||rowfloats[i]>2))
{printf("Error reading %s; the value (%.6f) for Predictor %d, Individual %d is outside [0,2]\nIf the data do not represent biallelic SNP allele counts you should use \"--SNP-data NO\"\n\n", spedfile, rowfloats[i], current+1, i+1);exit(1);}
}}

for(i=0;i<num_samples_use;i++){data[(size_t)(j-start)*num_samples_use+i]=rowfloats[keepsamps[i]];}

if(threshold!=-9999)	//threshold values
{
for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)(j-start)*num_samples_use+i];
if(value!=missingvalue)
{
data[(size_t)(j-start)*num_samples_use+i]=missingvalue;
if(value<=1-threshold){data[(size_t)(j-start)*num_samples_use+i]=0.0;}
if(value>=1+threshold){data[(size_t)(j-start)*num_samples_use+i]=2.0;}
if(value>=threshold&&value<=2-threshold){data[(size_t)(j-start)*num_samples_use+i]=1.0;}
}
}
}
}	//end of j loop

fclose(input);
free(rowfloats);

return(0);
}	//end of read_sped_fly

////////

float read_speed_size(char *speedfile)
{
float size;

FILE *input;


if((input=fopen(speedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",speedfile);exit(1);}

fseeko(input, 0, SEEK_SET);
if(fread(&size, sizeof(float), 1, input)!=1)
{printf("Error reading header of first predictor in %s\n\n", speedfile);exit(1);}
if(size!=1&&size!=2)
{printf("Error reading %s; values should be stored using either 1 or 2 bytes (not %.1f)\n\n", speedfile, size);exit(1);}
fclose(input);

return(size);
}

////////

void read_speed_full(char *speedfile, float *speedstarts, float *speedscales, unsigned char **data_char, unsigned short **data_short, int num_samples_use, int *keepsamps, int length, int *keeppreds, int num_samples, int num_preds, int nonsnp, int maxthreads)
{
int i, j;
int thread, threadstart, threadend, threadlength, *Mcurrent;
float size, startfloats[16], minfloat, maxfloat;
unsigned char **Mrowchars;
unsigned short **Mrowshorts;

FILE **Minput;


threadlength=(length-1)/maxthreads+1;

Mcurrent=malloc(sizeof(int)*maxthreads);
Minput=malloc(sizeof(FILE *)*maxthreads);

size=read_speed_size(speedfile);
if(size==1)
{
Mrowchars=malloc(sizeof(unsigned char*)*maxthreads);
for(thread=0;thread<maxthreads;thread++){Mrowchars[thread]=malloc(sizeof(unsigned char)*num_samples);}
}
else
{
Mrowshorts=malloc(sizeof(unsigned short*)*maxthreads);
for(thread=0;thread<maxthreads;thread++){Mrowshorts[thread]=malloc(sizeof(unsigned short)*num_samples);}
}

//check the length
if((Minput[0]=fopen(speedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",speedfile);exit(1);}

fseeko(Minput[0], 0, SEEK_END);
if(ftello(Minput[0])!=(sizeof(float)*16+size*num_samples)*num_preds)
{printf("Error reading %s; should have size %jd (%d samples and %d predictors), but instead has size %jd\n\n", speedfile, (off_t)(sizeof(float)*16+size*num_samples)*num_preds, num_samples, num_preds, ftello(Minput[0]));exit(1);}
fclose(Minput[0]);

//now can read
#pragma omp parallel for private(thread,threadstart,threadend,j,i,minfloat,maxfloat) schedule (static, 1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=thread*threadlength;
threadend=(thread+1)*threadlength;
if(threadend>length){threadend=length;}

if((Minput[thread]=fopen(speedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",speedfile);exit(1);}
Mcurrent[thread]=0;

for(j=threadstart;j<threadend;j++)
{
if(thread==0&&j%20000==0&&j*maxthreads<length)
{printf("Reading predictor %d out of %d\n", j*maxthreads+1, length);}

if(keeppreds[j]!=Mcurrent[thread])
{
if(fseeko(Minput[thread], (off_t)(sizeof(float)*16+size*num_samples)*keeppreds[j], SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d\n\n", speedfile, j+1);exit(1);}
}

if(fread(startfloats, sizeof(float), 16, Minput[thread])!=16)
{printf("Error reading headers of Predictor %d in %s\n\n", keeppreds[j]+1, speedfile);exit(1);}
if(startfloats[0]!=size)
{printf("Error reading %s; the first predictor used %d bytes per value, while Predictor %d uses %.1f\n\n", speedfile, (int)size, keeppreds[j]+1, startfloats[0]);exit(1);}

minfloat=startfloats[1];maxfloat=startfloats[2];
if(nonsnp==0&&(minfloat!=0||maxfloat!=2))
{printf("Error, %s does not contain biallelic SNP allele counts; you should use \"--SNP-data NO\"\n\n", speedfile);exit(1);}

speedstarts[j]=minfloat;
if(size==1){speedscales[j]=(maxfloat-minfloat)/254;}
else{speedscales[j]=(maxfloat-minfloat)/65534;}

if(size==1)
{
if(fread(Mrowchars[thread], sizeof(unsigned char), num_samples, Minput[thread])!=num_samples)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, speedfile);exit(1);}
for(i=0;i<num_samples_use;i++){data_char[j][i]=Mrowchars[thread][keepsamps[i]];}
}
else
{
if(fread(Mrowshorts[thread], sizeof(unsigned short), num_samples, Minput[thread])!=num_samples)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, speedfile);exit(1);}
for(i=0;i<num_samples_use;i++){data_short[j][i]=Mrowshorts[thread][keepsamps[i]];}
}

Mcurrent[thread]=keeppreds[j]+1;
}	//end of j loop

fclose(Minput[thread]);
}
printf("\n");

free(Mcurrent);free(Minput);
if(size==1){for(thread=0;thread<maxthreads;thread++){free(Mrowchars[thread]);}free(Mrowchars);}
else{for(thread=0;thread<maxthreads;thread++){free(Mrowshorts[thread]);}free(Mrowshorts);}
}	//end of read_speed_full

////////

int read_speed_fly(char *speedfile, double *data, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, double missingvalue, double threshold, int nonsnp)
{
int i, j, current;
float size, startfloats[16], minfloat, maxfloat;
double value;
unsigned char *rowchars, onechar;
unsigned short *rowshorts, oneshort;

FILE *input;


if((input=fopen(speedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",speedfile);exit(1);}

size=read_speed_size(speedfile);
if(size==1){rowchars=malloc(sizeof(unsigned char)*num_samples);}
else{rowshorts=malloc(sizeof(unsigned short)*num_samples);}

//check the length
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(sizeof(float)*16+size*num_samples)*num_preds)
{printf("Error reading %s; should have size %jd (%d samples and %d predictors), but instead has size %jd\n\n", speedfile, (off_t)(sizeof(float)*16+size*num_samples)*num_preds, num_samples, num_preds, ftello(input));exit(1);}
fseeko(input, 0, SEEK_SET);
current=0;

for(j=start;j<end;j++)
{
if(keeppreds[j]!=current)
{
if(fseeko(input, (off_t)(sizeof(float)*16+size*num_samples)*keeppreds[j], SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d\n\n", speedfile, j+1);exit(1);}
}

if(fread(startfloats, sizeof(float), 16, input)!=16)
{printf("Error reading headers of Predictor %d from %s\n\n", keeppreds[j]+1, speedfile);exit(1);}
if(startfloats[0]!=size)
{printf("Error reading %s; the first predictor used %d bytes per value, while Predictor %d uses %.1f\n\n", speedfile, (int)size, keeppreds[j]+1, startfloats[0]);exit(1);}

minfloat=startfloats[1];maxfloat=startfloats[2];
if(nonsnp==0&&(minfloat!=0||maxfloat!=2))
{printf("Error, %s does not contain biallelic SNP allele counts; you should use \"--SNP-data NO\"\n\n", speedfile);exit(1);}

if(size==1)
{
if(fread(rowchars, sizeof(unsigned char), num_samples, input)!=num_samples)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, speedfile);exit(1);}
}
else
{
if(fread(rowshorts, sizeof(unsigned short), num_samples, input)!=num_samples)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, speedfile);exit(1);}
}
current=keeppreds[j]+1;

//convert to doubles
for(i=0;i<num_samples_use;i++)
{
if(size==1)	//predictors can have 256 values
{
onechar=rowchars[keepsamps[i]];
if(onechar==255){data[(size_t)(j-start)*num_samples_use+i]=missingvalue;}
else{data[(size_t)(j-start)*num_samples_use+i]=(double)onechar/254*(maxfloat-minfloat)+minfloat;}
}
else	//each value is two bytes - can have 65536 values
{
oneshort=rowshorts[keepsamps[i]];
if(oneshort==65535){data[(size_t)(j-start)*num_samples_use+i]=missingvalue;}
else{data[(size_t)(j-start)*num_samples_use+i]=(double)oneshort/65534*(maxfloat-minfloat)+minfloat;}
}
}

if(threshold!=-9999)	//threshold values
{
for(i=0;i<num_samples_use;i++)
{
value=data[(size_t)(j-start)*num_samples_use+i];
if(value!=missingvalue)
{
data[(size_t)(j-start)*num_samples_use+i]=missingvalue;
if(value<1-threshold){data[(size_t)(j-start)*num_samples_use+i]=0.0;}
if(value>1+threshold){data[(size_t)(j-start)*num_samples_use+i]=2.0;}
if(value>=threshold&&value<=2-threshold){data[(size_t)(j-start)*num_samples_use+i]=1.0;}
}
}
}
}	//end of j loop

fclose(input);

if(size==1){free(rowchars);}
else{free(rowshorts);}

return(0);
}	//end of read_speed_fly

/////////////////////

int read_gen_fly(char *genfile, double *data, float **probs, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, gzFile inputgz, int current, int num_samples, int num_preds, int genskip, int genheaders, int genprobs, double missingvalue, double threshold, double minprob, int nonsnp)
{
int i, k, found, wcount, offset, mark, size, flag;
double prob0, prob1, prob2, prob3, proball, unifrand, value;
double *datatemp, *datap0, *datap1;

char readstring[500], *rs;
char *gzbuffer, *buffptr, *buffptr2;

rs=malloc(sizeof(char)*10000000);


datatemp=malloc(sizeof(double)*num_samples);
datap0=malloc(sizeof(double)*num_samples);
datap1=malloc(sizeof(double)*num_samples);

size=30000000+num_samples*(genprobs*20+(genprobs==0)*4);
gzbuffer=malloc(sizeof(char)*size);

found=0;wcount=0;mark=-1;
while(1)
{
while(current<keeppreds[start+found])	//skip rows
{
if(gzgets(inputgz,gzbuffer,size)==NULL)
{printf("Error reading %s; there is no Row %d\n\n", genfile, genskip+current+1);exit(1);}
if(strlen(gzbuffer)>=size-1)
{printf("Error reading %s; Row %d is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", genfile, genskip+current+1, (int)strlen(gzbuffer), size);exit(1);}
current++;
}

if(keeppreds[start+found]==mark)	//then copy in previous row
{
for(i=0;i<num_samples_use;i++)
{
data[(size_t)found*num_samples_use+i]=data[(size_t)(found-1)*num_samples_use+i];
if(probs!=NULL)
{
probs[0][(size_t)found*num_samples_use+i]=probs[0][(size_t)(found-1)*num_samples_use+i];
probs[1][(size_t)found*num_samples_use+i]=probs[1][(size_t)(found-1)*num_samples_use+i];
}
}
found++;
}
else	//so new predictor
{
if(gzgets(inputgz,gzbuffer,size)==NULL)
{printf("Error reading %s; there is no Row %d\n\n", genfile, genskip+current+1);exit(1);}
if(strlen(gzbuffer)==size-1)
{printf("Error reading %s; Row %d is longer (%d) than expected/allowed (%d)\nPlease tell Doug\n\n", genfile, genskip+current+1, (int)strlen(gzbuffer), size);exit(1);}

buffptr=gzbuffer;
for(k=0;k<genheaders;k++)	//skip genheaders elements
{
if(sscanf(buffptr,"%s %n", rs, &offset)!=1)
{printf("Error reading header %d on Row %d of %s\n\n", k+1, genskip+current+1, genfile);exit(1);}
buffptr+=offset;
}

//now start reading probabilities, and standardize while here
for(i=0;i<num_samples;i++)
{
if(genprobs==0)	//have haps
{
if(buffptr[2]=='-'){buffptr[2]=buffptr[0];}
if(buffptr[0]=='0'&&buffptr[1]==' '&&buffptr[2]=='0')
{datatemp[i]=2.0;flag=1;}
if(buffptr[0]=='0'&&buffptr[1]==' '&&buffptr[2]=='1')
{datatemp[i]=1.0;flag=1;}
if(buffptr[0]=='1'&&buffptr[1]==' '&&buffptr[2]=='0')
{datatemp[i]=1.0;flag=1;}
if(buffptr[0]=='1'&&buffptr[1]==' '&&buffptr[2]=='1')
{datatemp[i]=0.0;flag=1;}
if(flag==0)
{
printf("Error reading haplotypes for Predictor %d, Individual %d in %s (%c%c%c%c)\n\n", genskip+current+1, i+1, genfile, buffptr[0], buffptr[1], buffptr[2], buffptr[3]);exit(1);
}
buffptr+=4;
}

if(genprobs==1)	//have dosages
{
if((buffptr[0]=='0'||buffptr[0]=='1'||buffptr[0]=='2')&&buffptr[1]==' ')
{datatemp[i]=(int)buffptr[0]-48;buffptr+=2;}
else
{
if(sscanf(buffptr, "%[0-9eE.+-] %n", readstring, &offset)==1)
{datatemp[i]=atof(readstring);buffptr+=offset;}
else
{
if(sscanf(buffptr, "%s %n", readstring, &offset)!=1)
{printf("Error reading value for Predictor %d, Individual %d in %s\n\n", genskip+current+1, i+1, genfile);exit(1);}
buffptr+=offset;
if(strcmp(readstring,"NA")==0){datatemp[i]=missingvalue;}
else	//it is not NA
{printf("Error reading %s; unknown value (%s) for Predictor %d, Individual %d\n", genfile, readstring, genskip+current+1, i+1);exit(1);}
}
}

if(nonsnp==0&&datatemp[i]!=missingvalue&&(datatemp[i]<0||datatemp[i]>2))
{printf("Error reading %s; the value (%.6f) for Predictor %d, Individual %d is outside [0,2]\nIf the data do not represent biallelic SNP allele counts you should use \"--SNP-data NO\"\n\n", genfile, datatemp[i], current+1, i+1);exit(1);}
}

////////

if(genprobs>1)	//have probabilities
{
//prob0
if((buffptr[0]=='0'||buffptr[0]=='1')&&buffptr[1]==' ')
{prob0=(int)buffptr[0]-48;offset+=2;buffptr+=2;}
else
{prob0=strtod(buffptr, &buffptr2);buffptr=buffptr2+1;}
datap0[i]=prob0;

//prob1
if((buffptr[0]=='0'||buffptr[0]=='1')&&buffptr[1]==' ')
{prob1=(int)buffptr[0]-48;offset+=2;buffptr+=2;}
else
{prob1=strtod(buffptr, &buffptr2);buffptr=buffptr2+1;}
datap1[i]=prob1;

if(genprobs==2)	//can only check 2 probs not too high
{
if(prob0+prob1>1.02)
{printf("Error reading %s; the two probabilities for Predictor %d, Individual %d (Columns %d and %d) sum to %.6f (greater than 1)\n\n", genfile, genskip+current+1, i+1, genheaders+i*genprobs+1, genheaders+i*genprobs+2, prob0+prob1);exit(1);}
prob2=1.0-prob0-prob1;
if(prob0+prob1>1.0)	//difference treated as rounding error
{
value=prob0+prob1;
prob0=prob0/value;prob1=prob1/value;
}
prob2=1.0-prob0-prob1;
proball=1.0;
prob3=0.0;
}

if(genprobs>=3)
{
//prob2
if((buffptr[0]=='0'||buffptr[0]=='1')&&buffptr[1]==' ')
{prob2=(int)buffptr[0]-48;offset+=2;buffptr+=2;}
else
{prob2=strtod(buffptr, &buffptr2);buffptr=buffptr2+1;}
}

if(genprobs==3)	//check 3 probs not too high (will later warn if too low)
{
if(prob0+prob1+prob2>1.02)
{printf("Error reading %s; the three probabilities for Predictor %d, Individual %d (Columns %d to %d) sum to %.6f (greater than 1)\n\n", genfile, genskip+current+1, i+1, genheaders+i*genprobs+1, genheaders+i*genprobs+3, prob0+prob1+prob2);exit(1);}
if(prob0+prob1+prob2>=.98)	//difference treated as rounding error
{
value=prob0+prob1+prob2;
prob0=prob0/value;prob1=prob1/value;prob2=prob2/value;prob3=0;
}
proball=prob0+prob1+prob2;
prob3=1.0-proball;
}

if(genprobs==4)	//check 4 probs not too high nor too low (will later warn if 3 probs too low)
{
//prob3
if((buffptr[0]=='0'||buffptr[0]=='1')&&buffptr[1]==' ')
{prob3=(int)buffptr[0]-48;offset+=2;buffptr+=2;}
else
{prob3=strtod(buffptr, &buffptr2);buffptr=buffptr2+1;}

if(prob0+prob1+prob2+prob3>1.02)
{printf("Error reading %s; the four probabilities for Predictor %d, Individual %d (Columns%d to %d) sum to %.6f (greater than 1)\n\n", genfile, genskip+current+1, i+1, genheaders+i*genprobs+1, genheaders+i*genprobs+4, prob0+prob1+prob2+prob3);exit(1);}
if(prob0+prob1+prob2+prob3<0.98)
{printf("Error reading %s; the four probabilities for Predictor %d, Individual %d (Columns%d to %d) sum to %.6f (less than 1)\n\n", genfile, genskip+current+1, i+1, genheaders+i*genprobs+1, genheaders+i*genprobs+4, prob0+prob1+prob2+prob3);exit(1);}
if(prob0+prob1+prob2>=.98)	//difference treated as rounding error
{
value=prob0+prob1+prob2;
prob0=prob0/value;prob1=prob1/value;prob2=prob2/value;prob3=0;
}
proball=prob0+prob1+prob2;
}

if(minprob==-9999)	//saving dosage
{
if(proball<0.98)
{
if(wcount<5)
{printf("Warning reading %s; state probabilities for Predictor %d, Individual %d (Columns %d to %d) sum to only %.6f (less than 0.98), so value will be treated as missing. This is uncommon using imputed data\n", genfile, genskip+current+1, i+1, genheaders+i*genprobs+1, genheaders+(i+1)*genprobs, proball);}
datatemp[i]=missingvalue;
wcount++;
}
else	//so must be above 0.98 (in which case has been rounded to 1)
{
datatemp[i]=prob1+2*prob0;

if(threshold!=-9999)	//convert this to hard call
{
value=missingvalue;
if(datatemp[i]<1-threshold){value=0.0;}
if(datatemp[i]>1+threshold){value=2.0;}
if(datatemp[i]>=threshold&&datatemp[i]<=2-threshold){value=1.0;}
datatemp[i]=value;
}
}
}
else	//will sample based on probabilities or select most likely
{
if(proball<0.98)
{
if(wcount<5)
{printf("Warning reading %s; state probabilities for Predictor %d, Individual %d (Columns %d to %d) sum to only %.6f (less than 0.98). This is uncommon using imputed data\n", genfile, genskip+current+1, i+1, genheaders+i*genprobs+1, genheaders+(i+1)*genprobs, proball);}
datatemp[i]=missingvalue;
wcount++;
}
if(minprob==0)	//sample from probabilities
{
//unifrand=(double)rand()/RAND_MAX*(prob0+prob1+prob2+prob3);
unifrand=genrand_real1()*(prob0+prob1+prob2+prob3);
datatemp[i]=2.0;
if(unifrand>prob0){datatemp[i]=1.0;}
if(unifrand>prob0+prob1){datatemp[i]=0.0;}
if(unifrand>prob0+prob1+prob2){datatemp[i]=missingvalue;}
}
else	//select most likely
{
datatemp[i]=missingvalue;
if(prob0>minprob){datatemp[i]=2.0;}
if(prob1>minprob){datatemp[i]=1.0;}
if(prob2>minprob){datatemp[i]=0.0;}
}
}
}	//end of have probabilities
}	//end of i loop

//load up data
for(i=0;i<num_samples_use;i++)
{
data[(size_t)found*num_samples_use+i]=datatemp[keepsamps[i]];
if(probs!=NULL)
{
probs[0][(size_t)found*num_samples_use+i]=datap0[keepsamps[i]];
probs[1][(size_t)found*num_samples_use+i]=datap1[keepsamps[i]];
}
}
mark=keeppreds[start+found];found++;current++;
if(found==end-start){break;}
}
}	//end of while loop
if(wcount>5){printf("In total %d sets of probabilities summed to less than 0.98\n", wcount);}
if(wcount>0){printf("\n");}

free(datatemp);free(datap0);free(datap1);free(rs);free(gzbuffer);

return(current);
}	//end of read_gen_fly

/////////////////////


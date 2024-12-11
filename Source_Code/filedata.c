/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Functions for reading in data

///////////////////////////

void subset_bits(unsigned char *rowchars, unsigned char *rowchars2, int num_samples_use, int *keepsamps)
{
for(int i=0;i<num_samples_use;i++)
{
int i2=keepsamps[i];
//want to put position i of rowchars into position i2 of rowchars2

const int i_index=i/4;
const int i_pos=2*(i%4);
const int i2_index=i2/4;
const int i2_pos=2*(i2%4);
const int newbits=(rowchars[i2_index] >> i2_pos) & 3;

rowchars2[i_index] &= ~(3 << i_pos);
rowchars2[i_index] |= newbits << i_pos;
}
}

////////

void subset_bits_sparse(const int * __restrict rowchars, int * __restrict rowchars2, unsigned int * __restrict keep16s, int num_samples_use)
{
int cur_block, new_word, new_shift, new_count;
const int final_block=num_samples_use/16;
const int final_shift=num_samples_use%16;


//process blocks of 16 samples
cur_block=0;
new_word=0;
new_shift=0;
new_count=0;
while(1)
{
int cur_ind=keep16s[cur_block];

if(cur_ind)	//will use some samples from this block
{
int cur_word=rowchars[cur_block];
do
{
//find the position of the next sample to use
int next_bit=__builtin_ctz(cur_ind);

//add the sample's value to new_word
new_word |= ((cur_word >> (next_bit*2)) & 3) << (new_shift*2);
new_shift++;

if(new_shift==16)	//new_word is full, so must save and update / reset values
{
rowchars2[new_count]=new_word;
new_count++;
new_shift=0;
new_word=0;
}

//this removes the last one from the sample indicator
cur_ind &= cur_ind - 1;
} while(cur_ind);

if(new_count==final_block)	//are on final output block
{
if(new_shift==final_shift)	//have finished
{
if(new_shift)	//save the final word
{rowchars2[new_count]=new_word;}
break;
}
}
}	//end of using some samples from this block

cur_block++;
}
}

////////

void subset_bits_block(const int * __restrict rowchars, int * __restrict rowchars2, unsigned int * __restrict keep16s, int num_samples_use)
{
int cur_block, new_word, new_shift, new_count;
const int final_block=num_samples_use/16;
const int final_shift=num_samples_use%16;


//process blocks of 16 samples
cur_block=0;
new_word=0;
new_shift=0;
new_count=0;
while(1)
{
int cur_ind=keep16s[cur_block];

if(cur_ind)	//will use some samples from this block
{
int cur_word=rowchars[cur_block];
do
{
//find the position of the next sample to use
int next_bit=__builtin_ctz(cur_ind);

//get remainder of block	
int cur_word_rest= (cur_word >> (next_bit*2));

//see how many we can read at once
int cur_ind_inv = (~cur_ind) >> next_bit;
int num_bits=__builtin_ctz(cur_ind_inv);

//how much space in new word
int new_limit=16-new_shift;

//add all remaining values to new_word (might be adding too many)
new_word |= cur_word_rest << (new_shift*2);

if(num_bits<new_limit)	//update new_shift, then blank values we did not want to add
{
new_shift+=num_bits;
new_word &= ( 1 << (new_shift*2) ) -1;
}
else	//new word is full, so must save and update / reset values
{
rowchars2[new_count]=new_word;
new_count++;
new_shift=num_bits-new_limit;

if(new_shift)	//put extra bits at start of new word
{new_word = (cur_word_rest >> (new_limit*2)) & ((1 << (2 * new_shift)) - 1);}
else	//no extra bits, so new word is empty
{new_word=0;}
}

//this removes next_bits + num_bits ones from the sample indicator
cur_ind &= (~(1 << (next_bit + num_bits))) + 1;

} while(cur_ind);

if(new_count==final_block)	//are on final output block
{
if(new_shift==final_shift)	//have finished
{
if(new_shift)	//save the final word
{rowchars2[new_count]=new_word;}
break;
}
}
}	//end of using some samples from this block

cur_block++;
}
}

////////

void convert_doubles(double * __restrict data, int num_samples_use, unsigned char * __restrict rowchars, const double bedbytes[32])
{
int i, i2;

i=0;

//process blocks of 4 samples
while(i+4<=num_samples_use)
{
const int fourbits1=rowchars[i/4] & 15;
memcpy(data, bedbytes+2*fourbits1, sizeof(double)*2); 
data+=2;

const int fourbits2=(rowchars[i/4]>>4) & 15;
memcpy(data, bedbytes+2*fourbits2, sizeof(double)*2);
data+=2;
i+=4;
}

//load any remaining samples
for(i2=i;i2<num_samples_use;i2++)
{
const int twobits=(rowchars[i2/4]>>(2*(i2%4))) & 3;
memcpy(data, bedbytes+2*twobits, sizeof(double));    
data++;
}
}

///////////////////////////

void read_bed_fast(char *bedfile, double *data, double *centres, double *mults, double *sqdevs, double *rates, double *infos, int num_samples_use, int *keepsamps, int length, int *keeppreds, int num_samples, int num_preds, double missingvalue, int *bedzeros, int *bedones, int *bedtwos, int type)
{
//type=0 - got stats already, type=1 - get stats without standardizing, type=2 - get stats with standardizing
int i, i2, j, current, indcount, flag;
int c0, c1, c2;
double mean, var, value;

int num_16s, rowlength, rowlength2, rowlength3;
unsigned int *keep16s;
unsigned char startchars[3], *allchars, *newchars, *rowchars;
double *nulldoubles;

FILE *input;


//see whether necessary to subset samples
flag=0;
for(i=0;i<num_samples_use;i++)
{
if(keepsamps[i]!=i){flag=1;break;}
}

if(flag==1)	//see whether strictly monotonic
{
flag=4;
for(i=1;i<num_samples_use;i++)
{
if(keepsamps[i]<=keepsamps[i-1]){flag=1;break;}
}

if(flag==4)	//monotonic, so will subset using 16-sample indicators
{
num_16s=(num_samples-1)/16+1;
keep16s=malloc(sizeof(unsigned int)*num_16s);

for(i=0;i<num_16s;i++){keep16s[i]=0;}
for(i=0;i<num_samples_use;i++){keep16s[keepsamps[i]/16]|=(1 << (keepsamps[i]%16));}

//decide whether to use sparse subset (flag=2) or block subset (flag=3)
if(num_samples_use*3<num_samples*2){flag=2;}
else{flag=3;}
}
}

if(flag==1){printf("Warning, subset is not monotonic\n\n");}

//set rowlengths (either numbers of char for each predictor or number rounded to the next multiple of four)
rowlength=(num_samples-1)/4+1;
rowlength2=4*((num_samples-1)/16+1);
rowlength3=4*((num_samples_use-1)/16+1);

//allocate (allchars and newchars are a bit longer than necessary, to ensure no overflow when subsetting)
allchars=malloc(sizeof(unsigned char)*rowlength2);
newchars=malloc(sizeof(unsigned char)*rowlength3);
nulldoubles=malloc(sizeof(double)*num_samples_use);

//set nulldoubles to zero
for(i=0;i<num_samples_use;i++){nulldoubles[i]=0;}

//open bed file, check the length and that the first three digits ok

if((input=fopen(bedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",bedfile);exit(1);}

fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3)
{printf("Error reading %s; should have size %jd (%d samples x %d predictors), but instead has size %jd\n\n", bedfile, (off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3, num_samples, num_preds, ftello(input));exit(1);}

fseeko(input, 0, SEEK_SET);
if(fread(startchars, sizeof(unsigned char), 3, input)!=3)
{printf("Error reading first three values of %s\n\n", bedfile);exit(1);}
if(startchars[0]!=108||startchars[1]!=27)
{printf("Error reading %s; does not appear to be in binary PLINK format\n\n", bedfile);exit(1);}
if(startchars[2]!=1)
{printf("Error reading %s; can only read binary PLINK files stored in SNP-major mode\n\n", bedfile);exit(1);}
current=0;

//read and process one predictor at a time
for(j=0;j<length;j++)
{
if(keeppreds[j]!=current)
{
if(fseeko(input, (off_t)sizeof(unsigned char)*rowlength*keeppreds[j]+sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d %d\n\n", bedfile, j, keeppreds[j]+1);exit(1);}
}

if(fread(allchars, sizeof(unsigned char), rowlength, input)!=rowlength)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, bedfile);exit(1);}
current=keeppreds[j]+1;

if(flag==1)	//general subset
{subset_bits(allchars, newchars, num_samples_use, keepsamps);}
if(flag==2)	//sparse subset
{subset_bits_sparse((int *)allchars, (int *)newchars, keep16s, num_samples_use);}
if(flag==3)	//block subset
{subset_bits_block((int *)allchars, (int *)newchars, keep16s, num_samples_use);}

if(type!=0)	//need to get stats
{
//tally numbers of zeros, ones and twos
i=0;
c0=0;c1=0;c2=0;
if(flag==0){rowchars=allchars;}
else{rowchars=newchars;}

//process blocks of four samples
while(i+4<=num_samples_use)
{
c0+=bedzeros[rowchars[0]];
c1+=bedones[rowchars[0]];
c2+=bedtwos[rowchars[0]];
rowchars++;
i+=4;
}

//now any remaining samples (will be at most three, so no need to update rowchars)
for(i2=i;i2<num_samples_use;i2++)
{
switch ((rowchars[0]>>(2*(i2%4))) & 3)
{
case 3: c0++;break;
case 2: c1++;break;
case 0: c2++;
}
}

//now compute stats
indcount=c0+c1+c2;
if(indcount>0){mean=(double)(c1+2*c2)/indcount;var=(double)(c1+4*c2)/indcount-pow(mean,2);}
else{mean=0;var=0;}

centres[j]=mean;
if(var>0)
{
if(type==1){mults[j]=1;}
else{mults[j]=pow(var*indcount/num_samples_use,-0.5);}
}
else{mults[j]=-9999;}
sqdevs[j]=var*indcount/num_samples_use;
rates[j]=(double)indcount/num_samples_use;
value=centres[j]*(1-centres[j]/2);
if(value>0){infos[j]=sqdevs[j]/value;}
else{infos[j]=0;}
}	//end of type!=0

//convert to doubles
if(mults[j]!=-9999)
{
const double v0=(2-centres[j])*mults[j];
const double v1=0;
const double v2=(1-centres[j])*mults[j];
const double v3=(0-centres[j])*mults[j];
const double bedbytes[32]={v0,v0, v1,v0, v2,v0, v3,v0, v0,v1, v1,v1, v2,v1, v3,v1, v0,v2, v1,v2, v2,v2, v3,v2, v0,v3, v1,v3, v2,v3, v3,v3};

if(flag==0){convert_doubles(data+(size_t)j*num_samples_use, num_samples_use, allchars, bedbytes);}
else{convert_doubles(data+(size_t)j*num_samples_use, num_samples_use, newchars, bedbytes);}
}
else{memcpy(data+(size_t)j*num_samples_use, nulldoubles, sizeof(double)*num_samples_use);}
}	//end of j loop

fclose(input);

if(flag==2||flag==3){free(keep16s);}
free(allchars);free(newchars);free(nulldoubles);
}	//end of read_bed_fast

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
{printf("Error reading %s; can only read binary PLINK files stored in SNP-major mode\n\n", bedfile);exit(1);}
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

///////////////////////////

int read_bed_fly(char *bedfile, double *data, int num_samples_use, int *keepsamps, int length, int *keeppreds, int num_samples, int num_preds, double missingvalue)
{
int i, j, current, flag;

int num_16s, rowlength, rowlength2, rowlength3;
unsigned int *keep16s;
unsigned char startchars[3], *allchars, *newchars;

const double bedbytes[32]={2.0,2.0, missingvalue,2.0, 1.0,2.0, 0.0,2.0,
2.0,missingvalue, missingvalue,missingvalue, 1.0,missingvalue, 0.0,missingvalue,
2.0,1.0, missingvalue,1.0, 1.0,1.0, 0.0,1.0,
2.0,0.0, missingvalue,0.0, 1.0,0.0, 0.0,0.0};

FILE *input;


//see whether necessary to subset samples
flag=0;
for(i=0;i<num_samples_use;i++)
{
if(keepsamps[i]!=i){flag=1;break;}
}

if(flag==1)	//see whether strictly monotonic
{
flag=4;
for(i=1;i<num_samples_use;i++)
{
if(keepsamps[i]<=keepsamps[i-1]){flag=1;break;}
}

if(flag==4)	//monotonic, so will subset using 16-sample indicators
{
num_16s=(num_samples-1)/16+1;
keep16s=malloc(sizeof(unsigned int)*num_16s);

for(i=0;i<num_16s;i++){keep16s[i]=0;}
for(i=0;i<num_samples_use;i++){keep16s[keepsamps[i]/16]|=(1 << (keepsamps[i]%16));}

//decide whether to use sparse subset (flag=2) or block subset (flag=3)
if(num_samples_use*3<num_samples*2){flag=2;}
else{flag=3;}
}
}

//set rowlengths (either numbers of char for each predictor or number rounded to the next multiple of four)
rowlength=(num_samples-1)/4+1;
rowlength2=4*((num_samples-1)/16+1);
rowlength3=4*((num_samples_use-1)/16+1);

//allocate (allchars and newchars are a bit longer than necessary, to ensure no overflow when subsetting)
allchars=malloc(sizeof(unsigned char)*rowlength2);
newchars=malloc(sizeof(unsigned char)*rowlength3);

//open bed file, check the length and that the first three digits ok

if((input=fopen(bedfile,"rb"))==NULL)
{printf("Error opening %s\n\n",bedfile);exit(1);}

fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3)
{printf("Error reading %s; should have size %jd (%d samples x %d predictors), but instead has size %jd\n\n", bedfile, (off_t)sizeof(unsigned char)*rowlength*num_preds+sizeof(unsigned char)*3, num_samples, num_preds, ftello(input));exit(1);}

fseeko(input, 0, SEEK_SET);
if(fread(startchars, sizeof(unsigned char), 3, input)!=3)
{printf("Error reading first three values of %s\n\n", bedfile);exit(1);}
if(startchars[0]!=108||startchars[1]!=27)
{printf("Error reading %s; does not appear to be in binary PLINK format\n\n", bedfile);exit(1);}
if(startchars[2]!=1)
{printf("Error reading %s; can only read binary PLINK files stored in SNP-major mode\n\n", bedfile);exit(1);}
current=0;

//read and process one predictor at a time
for(j=0;j<length;j++)
{
if(keeppreds[j]!=current)
{
if(fseeko(input, (off_t)sizeof(unsigned char)*rowlength*keeppreds[j]+sizeof(unsigned char)*3, SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d %d\n\n", bedfile, j, keeppreds[j]+1);exit(1);}
}

if(fread(allchars, sizeof(unsigned char), rowlength, input)!=rowlength)
{printf("Error reading values for Predictor %d from %s\n\n", keeppreds[j]+1, bedfile);exit(1);}
current=keeppreds[j]+1;

if(flag==1)	//general subset
{subset_bits(allchars, newchars, num_samples_use, keepsamps);}
if(flag==2)	//sparse subset
{subset_bits_sparse((int *)allchars, (int *)newchars, keep16s, num_samples_use);}
if(flag==3)	//block subset
{subset_bits_block((int *)allchars, (int *)newchars, keep16s, num_samples_use);}

//convert to doubles
if(flag==0){convert_doubles(data+(size_t)j*num_samples_use, num_samples_use, allchars, bedbytes);}
else{convert_doubles(data+(size_t)j*num_samples_use, num_samples_use, newchars, bedbytes);}
}	//end of j loop

fclose(input);

if(flag==2||flag==3){free(keep16s);}
free(allchars);free(newchars);

return(0);
}	//end of read_bed_fly

/////////////////////

void convert_probs(Bytef *uncompdata2, int num_samples, double *datap0, double *datap1, char bits)
{
int i;


if(bits==8)
{
for(i=0;i<num_samples;i++)
{
datap0[i]=((double)(*(unsigned char *)(uncompdata2+2*i)))/255;
datap1[i]=((double)(*(unsigned char *)(uncompdata2+2*i+1)))/255;
}
}
else	//must be 16 bits
{
for(i=0;i<num_samples;i++)
{
datap0[i]=((double)(*(unsigned short *)(uncompdata2+4*i)))/65535;
datap1[i]=((double)(*(unsigned short *)(uncompdata2+4*i+2)))/65535;
}
}
}

////////

int read_bgen_fly(char *bgenfile, double *data, float **probs, int num_samples_use, int *keepsamps, int start, int end, int *keeppreds, int num_samples, int num_preds, size_t *bgen_indexes, double missingvalue, double threshold, double minprob)
{
int i, j, count, wcount, vcount, flag;
double prob0, prob1, prob2, value, unifrand;
double *datatemp, *datap0, *datap1;

uLongf maxLen1, maxLen2a, maxLen2b, sourceLen, destLen;
Bytef *compdata1, *uncompdata1, *compdata2, *uncompdata2;

unsigned char flags[4];
short readshort;
int hblock, bgen_comp, bgen_layout, readint, readint2;
char ploid1, ploid2, ploid3, phased, bits;

FILE *input;


datatemp=malloc(sizeof(double)*num_samples);
datap0=malloc(sizeof(double)*num_samples);
datap1=malloc(sizeof(double)*num_samples);

//for layout 1, can predict max compressed size (while uncompressed size is fixed)
maxLen1=(uLongf)(13*num_samples+12);
compdata1=malloc(sizeof(Bytef)*maxLen1);
uncompdata1=malloc(sizeof(Bytef)*6*num_samples);

//for layout 2, must guess max size for both compressed and uncompressed
maxLen2a=100*num_samples;
maxLen2b=100*num_samples;
compdata2=malloc(sizeof(Bytef)*maxLen2a);
uncompdata2=malloc(sizeof(Bytef)*maxLen2b);

//open bgenfile
if((input=fopen(bgenfile,"rb"))==NULL)
{printf("Error opening %s\n\n", bgenfile);exit(1);}
fseeko(input, 0, SEEK_SET);

//read header length
fseeko(input, 4, SEEK_SET);
if(fread(&hblock, 4, 1, input)!=1)
{printf("Error reading second value of %s\n\n", bgenfile);exit(1);}

//read and check number of predictors
fseeko(input, 8, SEEK_SET);
if(fread(&count, 4, 1, input)!=1)
{printf("Error reading third value of %s\n\n", bgenfile);exit(1);}
if(count!=num_preds){printf("Error N3DK, please tell Doug %d %d\n\n", count, num_preds);exit(1);}

//read and check number of samples
if(fread(&count, 4, 1, input)!=1)
{printf("Error reading fourth value of %s\n\n", bgenfile);exit(1);}
if(count!=num_samples){printf("Error N3DL, please tell Doug %d %d\n\n", count, num_samples);exit(1);}

//read the flags
fseeko(input, 4+hblock-4, SEEK_SET);
if(fread(flags, 1, 4, input)!=4)
{printf("Error reading final header value from %s\n\n", bgenfile);exit(1);}

bgen_comp=flags[0] & 3;
bgen_layout=(flags[0]>> 2) & 3;

if(bgen_comp==2){printf("Error N3DM, please tell Doug\n\n");exit(1);}

////////

wcount=0;vcount=0;
for(j=start;j<end;j++)
{
if(fseeko(input, bgen_indexes[keeppreds[j]], SEEK_SET)!=0)
{printf("Error reading %s; unable to find Predictor %d %d\n\n", bgenfile, j, keeppreds[j]+1);exit(1);}

if(bgen_layout==1)	//read 3n two-byte probabilities, and convert to dosages
{
if(bgen_comp==0)	//uncompressed
{
if(fread(uncompdata1, 1, 6*num_samples, input)!=6*num_samples)
{printf("Error reading probabilities for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
}
else	//compressed
{
if(fread(&readint, 4, 1, input)!=1)
{printf("Error reading storage size for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
if(readint>(int)maxLen1)
{printf("Error, size of compressed probabilities for Predictor %d (%d) is larger than the maximum allowed (%d), please tell Doug\n\n", j+1, readint, (int)maxLen1);exit(1);}

if(fread(compdata1, 1, readint, input)!=readint)
{printf("Error reading compressed probabilities for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}

sourceLen=(uLongf)(readint);
destLen=(uLongf)(6*num_samples);
if(uncompress(uncompdata1, &destLen, compdata1, sourceLen)!= Z_OK)
{printf("Error decompressing probabilities for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
}

for(i=0;i<num_samples;i++)
{
prob0=(double)(*(unsigned short *)(uncompdata1+6*i));
prob1=(double)(*(unsigned short *)(uncompdata1+6*i+2));
prob2=(double)(*(unsigned short *)(uncompdata1+6*i+4));
if(prob0+prob1+prob2==0){datatemp[i]=missingvalue;datap0[i]=0;datap1[i]=0;}
else
{value=prob0+prob1+prob2;datatemp[i]=(2*prob0+prob1)/value;datap0[i]=prob0/value;datap1[i]=prob1/value;}
}
}
else	//layout 2
{
if(bgen_comp==0)	//uncompressed - get size of uncompressed data, then read data
{
if(fread(&readint2, 4, 1, input)!=1)
{printf("Error reading storage size for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
if(readint2>(int)maxLen2b)
{
free(uncompdata2);
maxLen2b=(uLongf)readint2;
uncompdata2=malloc(sizeof(Bytef)*maxLen2b);
}
if(fread(uncompdata2, 1, readint2, input)!=readint2)
{printf("Error reading probabilities for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
}
else	//compressed - get size of uncompressed and compressed data, then read data
{
if(fread(&readint, 4, 1, input)!=1)
{printf("Error reading storage size for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
if(readint-4>(int)maxLen2a)
{
free(compdata2);
maxLen2a=(uLongf)(readint-4);
compdata2=malloc(sizeof(Bytef)*maxLen2a);
}

if(fread(&readint2, 4, 1, input)!=1)
{printf("Error reading uncompressed size for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
if(readint2>(int)maxLen2b)
{
free(uncompdata2);
maxLen2b=(uLongf)readint2;
uncompdata2=malloc(sizeof(Bytef)*maxLen2b);
}

if(fread(compdata2, 1, readint-4, input)!=readint-4)
{printf("Error reading compressed probabilities for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}

sourceLen=(uLongf)(readint-4);
destLen=(uLongf)readint2;
if(uncompress(uncompdata2, &destLen, compdata2, sourceLen)!= Z_OK)
{printf("Error decompressing probabilities for Predictor %d from %s\n\n", keeppreds[j]+1, bgenfile);exit(1);}
}

//data starts n (4), nalleles (2), min ploid (1), max ploid (1), ploids (n), flag (1), bits(1)
count=*(int *)uncompdata2;
readshort=*(short *)(uncompdata2+4);
ploid1=*(unsigned char *)(uncompdata2+6);
ploid2=*(unsigned char *)(uncompdata2+7);
phased=*(unsigned char *)(uncompdata2+8+num_samples);
bits=*(unsigned char *)(uncompdata2+8+num_samples+1);

if(count!=num_samples)
{printf("Error N3DP, please tell Doug %d %d\n\n", count, num_samples_use);exit(1);}
if(readshort!=2)
{printf("Error N3DO, please tell Doug %hd\n\n", readshort);exit(1);}

if(ploid1!=2||ploid2!=2)
{
if(wcount<5){printf("Warning, predictor %d is not haploid for all samples (the ploidy ranges from %d to %d), so its values will be set to missing\n\n", keeppreds[j]+1, (int)ploid1, (int)ploid2);}

for(i=0;i<num_samples;i++){datatemp[i]=missingvalue;datap0[i]=0;datap1[i]=0;}
wcount++;
}
else
{
if(bits!=8&&bits!=16)
{
if(vcount<5){printf("Warning, predictor %d is coded using %d bits, so its values will be set to missing (LDAK requires either 8 or 16 bits)\n\n", keeppreds[j]+1, (int)bits);}

for(i=0;i<num_samples;i++){datatemp[i]=missingvalue;datap0[i]=0;datap1[i]=0;}
vcount++;
}
else
{
//because we restrict to haploids, we always have 2n probabilities
//correspond to genotypes A1A1 and A1A2 (phase=0) or A1 allele for each haplotype (phase=1)

count=10+num_samples+2*(int)bits/8*num_samples;
if(readint2!=count)
{printf("Error, uncompressed probabilities for predictor %d should have size %d (not %d)\n\n", keeppreds[j]+1, count, readint2);exit(1);}

//convert data to probabilities (those corresponding to missing are redundant)
convert_probs(uncompdata2+8+num_samples+2, num_samples, datap0, datap1, bits);

if(phased==1)	//convert haplotype probs to genotype probs
{
for(i=0;i<num_samples;i++)
{
prob0=datap0[i];
prob1=datap1[i];
datap0[i]=prob0*prob1;
datap1[i]=prob0+prob1-2*prob0*prob1;
}
}

//check for missing and convert to dosages
for(i=0;i<num_samples;i++)
{
ploid3=*(unsigned char *)(uncompdata2+8+i);
flag=(int)((ploid3>> 7) & 1);
if(flag==1){datatemp[i]=missingvalue;}
else{datatemp[i]=2*datap0[i]+datap1[i];}
}
}	//end of correct bits
}	//end of haploid
}	//end of layout2

///////

if(threshold!=-9999)	//convert dosage to hard call
{
for(i=0;i<num_samples;i++)
{
value=missingvalue;
if(datatemp[i]<1-threshold){value=0.0;}
if(datatemp[i]>1+threshold){value=2.0;}
if(datatemp[i]>=threshold&&datatemp[i]<=2-threshold){value=1.0;}
datatemp[i]=value;
}
}

if(minprob==0)	//sample genotypes from probabilities (overwrite dosages)
{
for(i=0;i<num_samples;i++)
{
if(datatemp[i]!=missingvalue)
{
unifrand=genrand_real1();
datatemp[i]=2.0;
if(unifrand>datap0[i]){datatemp[i]=1.0;}
if(unifrand>datap0[i]+datap1[i]){datatemp[i]=0.0;}
}
}
}

if(minprob>0)	//select most likely genotype (minprob must be above 0.5, so only one of below is true)
{
for(i=0;i<num_samples;i++)
{
if(datatemp[i]!=missingvalue)
{
datatemp[i]=missingvalue;
if(datap0[i]>minprob){datatemp[i]=2;}
if(datap1[i]>minprob){datatemp[i]=1;}
if(1-datap0[i]-datap1[i]>minprob){datatemp[i]=0;}
}
}
}

//load up data
for(i=0;i<num_samples_use;i++)
{
data[(size_t)(j-start)*num_samples_use+i]=datatemp[keepsamps[i]];
if(probs!=NULL)
{
probs[0][(size_t)(j-start)*num_samples_use+i]=datap0[keepsamps[i]];
probs[1][(size_t)(j-start)*num_samples_use+i]=datap1[keepsamps[i]];
}
}
}	//end of j loop
if(wcount>5){printf("In total %d predictors are not haploid\n", wcount);}
if(vcount>5){printf("In total %d predictors have unsuitable bit codings\n", vcount);}
if(wcount>0||vcount>0){printf("\n");}

fclose(input);

free(datatemp);free(datap0);free(datap1);
free(compdata1);free(uncompdata1);
free(compdata2);free(uncompdata2);

return(0);
}	//end of read_bgen_fly

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
//read prob0
if((buffptr[0]=='0'||buffptr[0]=='1')&&buffptr[1]==' ')
{prob0=(int)buffptr[0]-48;offset+=2;buffptr+=2;}
else
{prob0=strtod(buffptr, &buffptr2);buffptr=buffptr2+1;}

//read prob1
if((buffptr[0]=='0'||buffptr[0]=='1')&&buffptr[1]==' ')
{prob1=(int)buffptr[0]-48;offset+=2;buffptr+=2;}
else
{prob1=strtod(buffptr, &buffptr2);buffptr=buffptr2+1;}

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
//read prob2
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
//read prob3
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

//save first two (normalized) probs
datap0[i]=prob0;
datap1[i]=prob1;

if(minprob==-9999)	//get dosage genotype
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
else	//sample genotypes based on probabilities or select most likely
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
else	//select most likely (minprob must be above 0.5, so only one of below is true)
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

free(datatemp);free(datap0);free(datap1);free(rs);free(gzbuffer);

return(current);
}	//end of read_gen_fly

/////////////////////


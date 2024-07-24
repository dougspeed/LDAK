/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Sorting functions

///////////////////////////

//for use with qsort

int compare_int (const void *a, const void *b)
{
return(*(int*)a-*(int*)b);
}

int compare_float (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}

int compare_float_rev (const void * a, const void * b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fb > fa) - (fb < fa);
}

int compare_double (const void * a, const void * b)
{
  double fa = *(const double*) a;
  double fb = *(const double*) b;
  return (fa > fb) - (fa < fb);
}

int compare_double_rev (const void * a, const void * b)
{
  double fa = *(const double*) a;
  double fb = *(const double*) b;
  return (fb > fa) - (fb < fa);
}

int compare_string (const void *a, const void *b)
{
  const char *fa = *(const char **) a;
  const char *fb = *(const char **) b;
  return(strcmp(fa,fb));
}

int compare_sorting_double (const void *a, const void *b)
{
  const struct sorting_double fa = *(const struct sorting_double *) a;
  const struct sorting_double fb = *(const struct sorting_double *) b;
  return (fa.value > fb.value) - (fa.value < fb.value);
}

int compare_sorting_double_rev (const void *a, const void *b)
{
  const struct sorting_double fa = *(const struct sorting_double *) a;
  const struct sorting_double fb = *(const struct sorting_double *) b;
  return (fb.value > fa.value) - (fb.value < fa.value);
}

int compare_sorting_string (const void *a, const void *b)
{
  const struct sorting_string fa = *(const struct sorting_string *) a;
  const struct sorting_string fb = *(const struct sorting_string *) b;
  return(strcmp(fa.ptr,fb.ptr));
}

///////////////////////////

void check_dups(char **str, int length, char *filename, int *order, int type)
//type=0 - warn, type=1 - error, type=2 - only sort
{
int j;

struct sorting_string *sptrs;


sptrs=malloc(sizeof(struct sorting_string)*length);
for(j=0;j<length;j++){sptrs[j].ptr=str[j];sptrs[j].index=j;}
qsort(sptrs, length, sizeof(struct sorting_string), compare_sorting_string);

if(order!=NULL)	//save sorting order
{
for(j=0;j<length;j++){order[j]=sptrs[j].index;}
}

if(type!=2)
{
for(j=1;j<length;j++)
{
if(strcmp(sptrs[j-1].ptr,sptrs[j].ptr)==0)
{
if(type==1){printf("Error, %s appears twice in %s\n\n", sptrs[j].ptr, filename);exit(1);}
printf("Warning, %s appears twice in %s\n\n", sptrs[j].ptr, filename);break;
}
}
}

free(sptrs);
}

////////

int find_strings(char **stra, int lengtha, char **strb, int lengthb, int *indexa, int *indexb, char *filenamea, char *filenameb, int *predordera, int *predorderb, int type)
//sees which of stra are in strb - gets indexes of stra used and to which these match in strb
//if filenamea/filenameb provided, should check for duplicates
//type states which list needs sorting (type=3, must sort both)
{
int a, b, count;
int *ordera, *orderb, *founds;

struct sorting_string *sptrsa, *sptrsb;


//get order for stra
ordera=malloc(sizeof(int)*lengtha);
if(type==1||type==3)	//stra not sorted
{
if(predordera!=NULL)	//order provided
{
for(a=0;a<lengtha;a++){ordera[a]=predordera[a];}
}
else
{
sptrsa=malloc(sizeof(struct sorting_string)*lengtha);
for(a=0;a<lengtha;a++){sptrsa[a].ptr=stra[a];sptrsa[a].index=a;
}
qsort(sptrsa, lengtha, sizeof(struct sorting_string), compare_sorting_string);
for(a=0;a<lengtha;a++){ordera[a]=sptrsa[a].index;}
free(sptrsa);
}
}
else	//srta already sorted
{
for(a=0;a<lengtha;a++){ordera[a]=a;}
}

//get order for strb
orderb=malloc(sizeof(int)*lengthb);
if(type==2||type==3)	//strb not sorted
{
if(predorderb!=NULL)	//order provided
{
for(b=0;b<lengthb;b++){orderb[b]=predorderb[b];}
}
else
{
sptrsb=malloc(sizeof(struct sorting_string)*lengthb);
for(b=0;b<lengthb;b++){sptrsb[b].ptr=strb[b];sptrsb[b].index=b;}
qsort(sptrsb, lengthb, sizeof(struct sorting_string), compare_sorting_string);
for(b=0;b<lengthb;b++){orderb[b]=sptrsb[b].index;}
free(sptrsb);
}
}
else	//strb already sorted
{
for(b=0;b<lengthb;b++){orderb[b]=b;}
}

////////

if(filenamea!=NULL)	//check for duplicates
{
for(a=1;a<lengtha;a++)
{
if(strcmp(stra[ordera[a-1]],stra[ordera[a]])==0)
{printf("Error, %s appears twice in %s\n\n", stra[ordera[a]], filenamea);exit(1);}
}
}

if(filenameb!=NULL)	//check for duplicates
{
for(b=1;b<lengthb;b++)
{
if(strcmp(strb[orderb[b-1]],strb[orderb[b]])==0)
{printf("Error, %s appears twice in %s\n\n", strb[orderb[b]], filenameb);exit(1);}
}
}

////////

//founds[a] will point to match of stra[a] in strb (or -1 if not present)
founds=malloc(sizeof(int)*lengtha);
b=0;
for(a=0;a<lengtha;a++)
{
founds[ordera[a]]=-1;
while(b<lengthb)
{
if(strcmp(stra[ordera[a]],strb[orderb[b]])==0)
{founds[ordera[a]]=orderb[b];break;}
if(strcmp(stra[ordera[a]],strb[orderb[b]])<0){break;}
b++;
}
}

//now fill up indexa and indexb
count=0;
for(a=0;a<lengtha;a++)
{
if(founds[a]!=-1)	//stra[a]=strb[founds[a]]
{
if(indexa!=NULL){indexa[count]=a;}
if(indexb!=NULL){indexb[count]=founds[a];}
count++;
}
}

free(ordera);free(orderb);free(founds);

return(count);
}

///////////////////////////

int uni_ids(char **idsa, int na, char **idsb, int nb, char **allids, int type)
//type states which list needs sorting (type=3, must sort both)
{
int a, b, counta, countb, found;
int *indexa, *indexb;

struct sorting_string *sptrsa, *sptrsb;

indexa=malloc(sizeof(int)*na);
indexb=malloc(sizeof(int)*nb);
sptrsa=malloc(sizeof(struct sorting_string)*na);
sptrsb=malloc(sizeof(struct sorting_string)*nb);

for(a=0;a<na;a++){sptrsa[a].ptr=idsa[a];sptrsa[a].index=a;}
if(type==1||type==3)
{qsort(sptrsa, na, sizeof(struct sorting_string), compare_sorting_string);}
for(a=0;a<na;a++){indexa[a]=sptrsa[a].index;}

for(b=0;b<nb;b++){sptrsb[b].ptr=idsb[b];sptrsb[b].index=b;}
if(type==2||type==3)
{qsort(sptrsb, nb, sizeof(struct sorting_string), compare_sorting_string);}
for(b=0;b<nb;b++){indexb[b]=sptrsb[b].index;}

counta=0;countb=0;
if(strcmp(idsa[indexa[0]],idsb[indexb[0]])<0)	//start with a
{copy_string(allids,0,idsa[indexa[0]]);counta=1;found=1;}
else	//start with b
{copy_string(allids,0,idsb[indexb[0]]);countb=1;found=1;}

while(counta<na&&countb<nb)
{
if(strcmp(idsa[indexa[counta]],idsb[indexb[countb]])<0)	//a comes next
{
if(strcmp(allids[found-1],idsa[indexa[counta]])!=0)	//and its new
{copy_string(allids,found,idsa[indexa[counta]]);found++;}
counta++;
}
else	//try b
{
if(strcmp(allids[found-1],idsb[indexb[countb]])!=0)
{copy_string(allids,found,idsb[indexb[countb]]);found++;}
countb++;
}
}

while(counta<na)	//deal with left over a
{
if(strcmp(allids[found-1],idsa[indexa[counta]])!=0)	//its new
{copy_string(allids,found,idsa[indexa[counta]]);found++;}
counta++;
}
while(countb<nb)	//then left over b
{
if(strcmp(allids[found-1],idsb[indexb[countb]])!=0)	//its new
{copy_string(allids,found,idsb[indexb[countb]]);found++;}
countb++;
}

free(indexa);free(indexb);free(sptrsa);free(sptrsb);

return(found);
}

////////

int inter_ids(char **idsa, int na, char **idsb, int nb, char **allids, int type)
//type states which list needs sorting (type=3, must sort both)
{
int a, b, counta, countb, found;
int *indexa, *indexb;

struct sorting_string *sptrsa, *sptrsb;

indexa=malloc(sizeof(int)*na);
indexb=malloc(sizeof(int)*nb);
sptrsa=malloc(sizeof(struct sorting_string)*na);
sptrsb=malloc(sizeof(struct sorting_string)*nb);

for(a=0;a<na;a++){sptrsa[a].ptr=idsa[a];sptrsa[a].index=a;}
if(type==1||type==3)
{qsort(sptrsa, na, sizeof(struct sorting_string), compare_sorting_string);}
for(a=0;a<na;a++){indexa[a]=sptrsa[a].index;}

for(b=0;b<nb;b++){sptrsb[b].ptr=idsb[b];sptrsb[b].index=b;}
if(type==2||type==3)
{qsort(sptrsb, nb, sizeof(struct sorting_string), compare_sorting_string);}
for(b=0;b<nb;b++){indexb[b]=sptrsb[b].index;}

counta=0;countb=0;found=0;
while(counta<na&&countb<nb)
{
if(strcmp(idsa[indexa[counta]],idsb[indexb[countb]])==0)	//a match
{copy_string(allids,found,idsa[indexa[counta]]);counta++;countb++;found++;}
else	//increase a or b as appropriate
{
if(strcmp(idsa[indexa[counta]],idsb[indexb[countb]])<0){counta++;}
else{countb++;}
}
}

free(indexa);free(indexb);free(sptrsa);free(sptrsb);

return(found);
}

///////////////////////////

int uni_preds(int *chra, double *bpa, char **predsa, char *al1a, char *al2a, int Na, int *chrb, double *bpb, char **predsb, char *al1b, char *al2b, int Nb, int *allchr, double *allbp, char **allnames, char *allal1, char *allal2, int k)
//a is existing, so always takes priority - can only change alleles in b
{
int counta, countb, found, dcount;


//get first element - know each list must have at least one non-negative
counta=0;countb=0;
while(bpa[counta]<=0){counta++;}
while(bpb[countb]<=0){countb++;}

if(chra[counta]<chrb[countb]||(chra[counta]==chrb[countb]&&bpa[counta]<=bpb[countb]))	//start with a
{
allchr[0]=chra[counta];allbp[0]=bpa[counta];
copy_string(allnames,0,predsa[counta]);
allal1[0]=al1a[counta];allal2[0]=al2a[counta];
counta++;found=1;
}
else	//start with b
{
allchr[0]=chrb[countb];allbp[0]=bpb[countb];
copy_string(allnames,0,predsb[countb]);
allal1[0]=al1b[countb];allal2[0]=al2b[countb];
countb++;found=1;
}

while(counta<Na){if(bpa[counta]>0){break;}counta++;}
while(countb<Nb){if(bpb[countb]>0){break;}countb++;}

//now rest
dcount=0;
while(counta<Na&&countb<Nb)
{
if(chra[counta]<chrb[countb]||(chra[counta]==chrb[countb]&&bpa[counta]<=bpb[countb]))	//a is next
{
//and must be new
allchr[found]=chra[counta];allbp[found]=bpa[counta];
copy_string(allnames,found,predsa[counta]);
allal1[found]=al1a[counta];allal2[found]=al2a[counta];
found++;counta++;
}
else	//b next
{
if(chrb[countb]!=allchr[found-1]||bpb[countb]!=allbp[found-1])	//and its new
{
allchr[found]=chrb[countb];allbp[found]=bpb[countb];
copy_string(allnames,found,predsb[countb]);
allal1[found]=al1b[countb];allal2[found]=al2b[countb];
found++;
}
else	//not new, check compatible
{dcount+=check_comp(allnames, allal1, allal2, found-1, predsb, al1b, al2b, countb, allchr, allbp, dcount);}
countb++;
}

while(counta<Na){if(bpa[counta]>0){break;}counta++;}
while(countb<Nb){if(bpb[countb]>0){break;}countb++;}
}

//process remainder (either remaining a or remaining b)
while(counta<Na)
{
allchr[found]=chra[counta];allbp[found]=bpa[counta];
copy_string(allnames,found,predsa[counta]);
allal1[found]=al1a[counta];allal2[found]=al2a[counta];found++;
counta++;

while(counta<Na){if(bpa[counta]>0){break;}counta++;}
}

while(countb<Nb)
{
if(chrb[countb]!=allchr[found-1]||bpb[countb]!=allbp[found-1])	//and its new
{
allchr[found]=chrb[countb];allbp[found]=bpb[countb];
copy_string(allnames,found,predsb[countb]);
allal1[found]=al1b[countb];allal2[found]=al2b[countb];
found++;
}
else	//not new, check compatible
{dcount+=check_comp(allnames, allal1, allal2, found-1, predsb, al1b, al2b, countb, allchr, allbp, dcount);}
countb++;

while(countb<Nb){if(bpb[countb]>0){break;}countb++;}
}

if(dcount>10){printf("In total, %d predictors had different names\n", dcount);}
if(dcount>0){printf("\n");}

return(found);
}

////////

int inter_preds(int *chra, double *bpa, char **predsa, char *al1a, char *al2a, int Na, int *chrb, double *bpb, char **predsb, char *al1b, char *al2b, int Nb, int *allchr, double *allbp, char **allnames, char *allal1, char *allal2, int k)
{
int counta, countb, found, dcount;


//get first element - know each list must have at least one non-negative
counta=0;countb=0;
while(bpa[counta]<=0){counta++;}
while(bpb[countb]<=0){countb++;}

found=0;
dcount=0;
while(counta<Na&&countb<Nb)
{
if(chra[counta]==chrb[countb]&&bpa[counta]==bpb[countb])	//add a then check b compatible
{
allchr[found]=chra[counta];allbp[found]=bpa[counta];
copy_string(allnames,found,predsa[counta]);
allal1[found]=al1a[counta];allal2[found]=al2a[counta];
dcount+=check_comp(allnames, allal1, allal2, found, predsb, al1b, al2b, countb, allchr, allbp, dcount);
counta++;countb++;found++;
}
else	//increase either a or b
{
if(chra[counta]<chrb[countb]||(chra[counta]==chrb[countb]&&bpa[counta]<bpb[countb])){counta++;}
else{countb++;}
}

while(counta<Na){if(bpa[counta]>0){break;}counta++;}
while(countb<Nb){if(bpb[countb]>0){break;}countb++;}
}

if(dcount>10){printf("In total, %d predictors had different names\n", dcount);}
if(dcount>0){printf("\n");}

return(found);
}

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Very simple functions (mainly)

///////////////////////////

int countrows(char *filename)
{
int count, flag;
char readchar;
FILE *input;


if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

if(fscanf(input,"%c",&readchar)!=1)
{printf("Error, %s appears to be empty\n\n", filename);exit(1);}
if(readchar==10)
{printf("Error, %s appears to start with an empty line\n\n", filename);exit(1);}
if(readchar==9||readchar==32)
{printf("Error, %s appears to start with a space\n\n", filename);exit(1);}

count=0;flag=0;
while(fscanf(input,"%c",&readchar)==1)
{
if(flag==1){printf("Error, Row %d of %s appears to be empty\n\n", count+1, filename);exit(1);}
if(readchar==10)	//new row, so add one - allow for blank last row
{
count++;
if(fscanf(input,"%c",&readchar)!=1){break;}
if(readchar==10){flag=1;}
if(readchar==9||readchar==32)
{printf("Error, Row %d of %s appears to start with a space\n\n", count+1, filename);exit(1);}
}
}
fclose(input);

if(readchar!=10&&flag==0)	//did not end with a newline
{count++;}

return(count);
}	//end of countrows

////////

int countrows_min(char *filename, int min)
{
int count, flag;
char readchar;
FILE *input;


if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

if(fscanf(input,"%c",&readchar)!=1)
{printf("Error, %s appears to be empty\n\n", filename);exit(1);}
if(readchar==10)
{printf("Error, %s appears to start with an empty line\n\n", filename);exit(1);}
if(readchar==9||readchar==32)
{printf("Error, %s appears to start with a space\n\n", filename);exit(1);}

count=0;flag=0;
while(fscanf(input,"%c",&readchar)==1)
{
if(flag==1){printf("Error, Row %d of %s appears to be empty\n\n", count+1, filename);exit(1);}
if(readchar==10)	//new row, so add one - allow for blank last row
{
count++;
if(fscanf(input,"%c",&readchar)!=1){break;}
if(readchar==10){flag=1;}
if(readchar==9||readchar==32)
{printf("Error, Row %d of %s appears to start with a space\n\n", count+1, filename);exit(1);}
}
if(count==min){break;}
}
fclose(input);

if(readchar!=10&&flag==0)	//did not end with a newline
{count++;}

return(count);
}	//end of countrows_min

////////

int countels(char *filename)
{
int count;
char readchar, readchar2;

FILE *input;


if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

if(fscanf(input,"%c",&readchar)!=1)
{printf("Error, %s appears to be empty\n\n", filename);exit(1);}
if(readchar==10)
{printf("Error, %s appears to start with an empty line\n\n", filename);exit(1);}
if(readchar==9||readchar==32)
{printf("Error, %s appears to start with a space\n\n", filename);exit(1);}

count=1;readchar2=readchar;
while(fscanf(input,"%c", &readchar)==1)
{
if((readchar2==9||readchar2==10||readchar2==32)&&readchar!=9&&readchar!=10&&readchar!=32){count++;}
readchar2=readchar;
}
fclose(input);

return(count);
}	//end of countels

////////

int countcols(char *filename)
{
int count;

char readchar, readchar2;
FILE *input;


if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n",filename);exit(1);}

if(fscanf(input,"%c",&readchar)!=1)
{printf("Error, %s appears to be empty\n\n", filename);exit(1);}
if(readchar==10)
{printf("Error, %s appears to start with an empty line\n\n", filename);exit(1);}
if(readchar==9||readchar==32)
{printf("Error, %s appears to start with a space\n\n", filename);exit(1);}

count=1;readchar2=readchar;
while(fscanf(input,"%c", &readchar)==1)
{
if(readchar==10){break;}
if((readchar2==9||readchar2==32)&&readchar!=9&&readchar!=32){count++;}
readchar2=readchar;
}
fclose(input);

return(count);
}	//end of countcols

///////////////////////////

void copy_string(char **stra, int a, char *strb)
{
stra[a]=malloc(sizeof(char)*(strlen(strb)+1));
strcpy(stra[a],strb);
}	//end of copy_string

////////

void read_strings(char *strfile, char **str, int length, int *indexer, int col, int head)
{
int j, k, count, found, want;

char readchar, *rs;
FILE *input;

rs=malloc(sizeof(char)*10000000);


count=countrows(strfile)-head;
if(count<0)
{printf("Error reading %s; has %d rows but should have at least %d\n\n", strfile, count+head, head+1);exit(1);}
if(indexer!=NULL)
{
if(indexer[length-1]>=count)
{printf("Error reading %s; has %d rows but should have at least %d\n\n", strfile, count+head, head+indexer[length-1]+1);exit(1);}
}

if(countcols(strfile)<col)
{printf("Error reading %s; has %d columns but should have at least %d\n\n", strfile, countcols(strfile), col);exit(1);}

//open and skip headers
if((input=fopen(strfile,"r"))==NULL)
{printf("Error opening %s\n\n", strfile);exit(1);}
j=0;
while(j<head)	//skip rows
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}j++;
}

found=0;
for(j=0;j<count;j++)	//read in first col-1 columns, then string, then skip to end
{
want=j;
if(indexer!=NULL){want=indexer[found];}
for(k=0;k<col-1;k++)
{
if(fscanf(input, "%s ", rs)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", k+1, head+j+1, strfile);exit(1);}
}
if(fscanf(input, "%s%c", rs, &readchar)!=2)
{printf("Error reading Element %d of Row %d of %s\n\n", col, head+j+1, strfile);exit(1);}
if(j==want)
{
str[found]=malloc(sizeof(char)*(strlen(rs)+1));
strcpy(str[found],rs);
found++;
if(found==length){break;}
}
while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
fclose(input);

free(rs);
}	//end of read_strings

////////

void read_values(char *valfile, double *values, int length, int *indexer, int col, int head, int miss)
{
//miss=0 - missing not allowed, miss=1 - missing allowed
int j, k, count, found, want;

char readchar, *rs;
FILE *input;

rs=malloc(sizeof(char)*10000000);


count=countrows(valfile)-head;
if(count<0)
{printf("Error, %s has %d rows but should have at least %d\n\n", valfile, count+head, head+1);exit(1);}
if(indexer!=NULL)
{
if(indexer[length-1]>=count)
{printf("Error, %s has %d rows but should have at least %d\n\n", valfile, count+head, head+indexer[length-1]+1);exit(1);}
}

if(countcols(valfile)<col)
{printf("Error, %s has %d columns but should have at least %d\n\n", valfile, countcols(valfile), col);exit(1);}

//open and skip headers
if((input=fopen(valfile,"r"))==NULL)
{printf("Error opening %s\n\n", valfile);exit(1);}
j=0;
while(j<head)	//skip rows
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}j++;
}

found=0;
for(j=0;j<count;j++)	//read in first col-1 columns, then value, then skip to end
{
want=j;
if(indexer!=NULL){want=indexer[found];}
for(k=0;k<col-1;k++)
{
if(fscanf(input, "%s ", rs)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", k+1, head+j+1, valfile);exit(1);}
}
if(fscanf(input, "%[0-9eE.+-]%c", rs, &readchar)!=2)
{
if(fscanf(input, "%s%c", rs, &readchar)!=2)
{printf("Error reading Element %d of Row %d of %s\n\n", col, head+j+1, valfile);exit(1);}
if(strcmp(rs,"NA")!=0)
{printf("Error reading %s; Element %d of Row %d is unrecognisable (%s)\n\n", valfile, col, head+j+1, rs);exit(1);}
}
if(j==want)
{
values[found]=atof(rs);
if(strcmp(rs,"NA")==0)
{
if(miss==0)
{printf("Error reading %s; Element %d of Row %d is NA\n\n", valfile, col, head+j+1);exit(1);}
values[found]=-9999;
}
found++;
if(found==length){break;}
}
while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
fclose(input);

free(rs);
}	//end of read_values
////////

void read_integers(char *valfile, int *values, int length, int *indexer, int col, int head, int miss)
{
//miss=0 - missing not allowed, miss=1 - missing allowed
int j, k, count, found, want;

char readchar, *rs;
FILE *input;

rs=malloc(sizeof(char)*10000000);


count=countrows(valfile)-head;
if(count<0)
{printf("Error, %s has %d rows but should have at least %d\n\n", valfile, count+head, head+1);exit(1);}
if(indexer!=NULL)
{
if(indexer[length-1]>=count)
{printf("Error, %s has %d rows but should have at least %d\n\n", valfile, count+head, head+indexer[length-1]+1);exit(1);}
}

if(countcols(valfile)<col)
{printf("Error, %s has %d columns but should have at least %d\n\n", valfile, countcols(valfile), col);exit(1);}

//open and skip headers
if((input=fopen(valfile,"r"))==NULL)
{printf("Error opening %s\n\n", valfile);exit(1);}
j=0;
while(j<head)	//skip rows
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}j++;
}

found=0;
for(j=0;j<count;j++)	//read in first col-1 columns, then value, then skip to end
{
want=j;
if(indexer!=NULL){want=indexer[found];}
for(k=0;k<col-1;k++)
{
if(fscanf(input, "%s ", rs)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", k+1, head+j+1, valfile);exit(1);}
}
if(fscanf(input, "%[0-9eE.+-]%c", rs, &readchar)!=2)
{
if(fscanf(input, "%s%c", rs, &readchar)!=2)
{printf("Error reading Element %d of Row %d of %s\n\n", col, head+j+1, valfile);exit(1);}
if(strcmp(rs,"NA")!=0)
{printf("Error reading %s; Element %d of Row %d is unrecognisable (%s)\n\n", valfile, col, head+j+1, rs);exit(1);}
}
if(j==want)
{
values[found]=atoi(rs);
if(strcmp(rs,"NA")==0)
{
if(miss==0)
{printf("Error reading %s; Element %d of Row %d is NA\n\n", valfile, col, head+j+1);exit(1);}
values[found]=-9999;
}
found++;
if(found==length){break;}
}
while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
fclose(input);

free(rs);
}	//end of read_integers

///////////////////////////

void read_ids(char *idsfile,  char **ids1, char **ids2, char **ids3, int length, int *indexer, int head, int shift)
{
//read ids from columns shift+1 and shift+2
int i, j, count, found, want;

char readchar, *rs1, *rs2;
FILE *input;

rs1=malloc(sizeof(char)*10000000);rs2=malloc(sizeof(char)*10000000);


count=countrows(idsfile)-head;
if(count==0){printf("Error, %s contains no samples\n\n", idsfile);exit(1);}
if(indexer!=NULL)
{
if(indexer[length-1]>=count)
{printf("Error reading %s; has %d rows but should have at least %d\n\n", idsfile, count+head, head+indexer[length-1]+1);exit(1);}
}

if(countcols(idsfile)<shift+2)
{printf("Error reading %s; has %d columns but should have at least %d\n\n", idsfile, countcols(idsfile),shift+2);exit(1);}

if((input=fopen(idsfile,"r"))==NULL)
{printf("Error opening %s\n",idsfile);exit(1);}

i=0;
while(i<head)	//skip rows
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}i++;
}

found=0;
for(i=0;i<count;i++)	//read skip shift elements, read two elements, then skip to end of line
{
want=i;
if(indexer!=NULL){want=indexer[found];}

for(j=0;j<shift;j++)
{
if(fscanf(input, "%s ", rs1)!=1){printf("Error reading Element %d from Row %d of %s\n", j+1, head+i+1, idsfile);exit(1);}
}

if(fscanf(input, "%s %s%c", rs1, rs2, &readchar)!=3)
{printf("Error reading Elements %d and %d from Row %d of %s\n", shift+1, shift+2, head+i+1, idsfile);exit(1);}

if(i==want)
{
if(ids1!=NULL){copy_string(ids1,found,rs1);}
if(ids2!=NULL){copy_string(ids2,found,rs2);}
if(ids3!=NULL)
{
ids3[found]=malloc(sizeof(char)*(strlen(rs1)+strlen(rs2)+4));
sprintf(ids3[found],"%s___%s",rs1,rs2);
}
found++;
if(found==length){break;}
}

while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
fclose(input);

free(rs1);free(rs2);
}	//end of read_ids

///////////////////////////

void print_empty(char *readstring, const char *arg)
{
printf("Error, you can only use one main argument (you have used both %s and %s)\n\n", readstring, arg);
}

////////

void print_scaling(double power, int hwestand)
{
if(power==0){printf("Predictors will be centred, but not scaled (option \"--power\")\n\n");}
else
{
printf("Predictors will be centred then scaled by V^(%.4f/2) (option \"--power\"), ", power);
if(hwestand==0)
{printf("where V is the observed variance of the predictor; to instead scale based on expected variance use \"--hwe-stand YES\"\n\n");}
else
{printf("where V=2*MAF*(1-MAF), the expected variance assuming Hardy-Weinberg Equilibrium; to instead scale based on observed variance use \"--hwe-stand NO\"\n\n");}
}
}

///////////////////////////

int just_check(char *filename)
{
char filename2[500];
struct stat statstruct;

if(strcmp(filename,"blank")!=0)
{
sprintf(filename2,"%s/", filename);
if(stat(filename2, &statstruct)==0){return(2);}
//if(access(filename, R_OK)!=0){return(3);}
if(stat(filename, &statstruct)!=0){return(3);}
return(statstruct.st_size==0);
}

return(0);
}

////////

int append_check(char *filename, char *filename2, char *workdir)
{
char filename3[500];
struct stat statstruct;

if(strcmp(filename2,"blank")!=0)	//prefix with workdir (unless absolute), then see if file/folder already exists
{
if(filename2[0]!='/'){sprintf(filename,"%s%s", workdir, filename2);}
else{strcpy(filename,filename2);}

sprintf(filename3,"%s/", filename);
if(stat(filename3, &statstruct)==0){return(2);}
//if(access(filename, R_OK)!=0){return(3);}
if(stat(filename, &statstruct)!=0){return(3);}
return(statstruct.st_size==0);
}
else{strcpy(filename,"blank");}

return(0);
}	//end of append_check

///////////////////////////

int fill_names(char **datastems, char **bimstems, char **famstems, int k, char *udatafile, int dtype)
{
if(dtype==1)
{sprintf(datastems[k],"%s.bed",udatafile);
sprintf(bimstems[k],"%s.bim",udatafile);
sprintf(famstems[k],"%s.fam",udatafile);}

if(dtype==3)
{sprintf(datastems[k],"%s.sped",udatafile);
sprintf(bimstems[k],"%s.bim",udatafile);
sprintf(famstems[k],"%s.fam",udatafile);}
if(dtype==4)
{sprintf(datastems[k],"%s.speed",udatafile);
sprintf(bimstems[k],"%s.bim",udatafile);
sprintf(famstems[k],"%s.fam",udatafile);}

return(0);
}

///////////////////////////

int check_head_ids(char *filename, int type)
//type=0 - keep quiet, type=1 - advise
{
int head;
char readstring[500], readstring2[500];

FILE *input;


head=0;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of %s\n\n", filename);exit(1);}
if((strcmp(readstring,"FID")==0&&strcmp(readstring2,"IID")==0)||(strcmp(readstring,"ID1")==0&&strcmp(readstring2,"ID2")==0))
{
if(type==1){printf("%s has a header row\n", filename);}
head=1;
}
fclose(input);

return(head);
}

////////

int check_head(char *filename, char *findstring, char *findstring2, int type)
//type=0 - keep quiet, type=1 - advise
{
int head;
char readstring[500];

FILE *input;


head=0;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading first element of %s\n\n", filename);exit(1);}
if(strcmp(readstring,findstring)==0||strcmp(readstring,findstring2)==0)
{
if(type==1){printf("%s has a header row\n", filename);}
head=1;
}
fclose(input);

return(head);
}

////////

int find_head(char *findstring, char *filename, int length)
{
int j, k, count;
char **readlist, *rs;

FILE *input;

rs=malloc(sizeof(char)*10000000);


readlist=malloc(sizeof(char*)*length);

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
for(j=0;j<length;j++)
{
if(fscanf(input, "%s ", rs)!=1)
{printf("Error reading name of Column %d of %s\n\n", j+1, filename);exit(1);}
copy_string(readlist, j, rs);
for(k=0;k<j;k++)
{
if(strcmp(readlist[k],readlist[j])==0)
{printf("Error reading %s; Columns %d and %d have the same name (%s)\n\n", filename, k+1, j+1, readlist[k]);exit(1);}
}
}
fclose(input);

count=-1;
for(j=0;j<length;j++)
{
if(strcmp(findstring,readlist[j])==0){count=j;break;}
}

for(j=0;j<length;j++){free(readlist[j]);}free(readlist);

free(rs);
return(count);
}	//end of find_head

///////////////////////////

int reduce_values(double *values, int length, int *keeppreds_use, int *chr, char **preds, double *cm, double *bp, double *cmbp, char *al1, char *al2, double *centres, double *pvalues, int num_regs, int **regindex)
{
int j, r, count, count2;
int *usedpreds;


//will only use usedpreds if have regions
usedpreds=malloc(sizeof(int)*length);

count=0;
for(j=0;j<length;j++)
{
if(values[j]>0)
{
if(count!=j)
{
values[count]=values[j];
keeppreds_use[count]=keeppreds_use[j];
if(chr!=NULL){chr[count]=chr[j];}
if(cm!=NULL){cm[count]=cm[j];}
if(bp!=NULL){bp[count]=bp[j];}
if(cmbp!=NULL){cmbp[count]=cmbp[j];}
al1[count]=al1[j];al2[count]=al2[j];
free(preds[count]);copy_string(preds,count,preds[j]);
if(centres!=NULL){centres[count]=centres[j];}
if(pvalues!=NULL){pvalues[count]=pvalues[j];}
}
usedpreds[j]=count;
count++;
}
else{usedpreds[j]=-1;}
}
for(j=count;j<length;j++){free(preds[j]);}

for(r=0;r<num_regs;r++)
{
count2=0;
for(j=0;j<regindex[r][0];j++)
{
if(usedpreds[regindex[r][1+j]]!=-1)	//still using this predictor, but it may have been moved
{regindex[r][1+count2]=usedpreds[regindex[r][1+j]];count2++;}
}
if(count2==0){printf("Error, all %d predictors in Region %d have zero weight\n\n", regindex[r][0], r+1);exit(1);}
regindex[r][0]=count2;
}

free(usedpreds);

return(count);
}	//end of reduce_values

///////////////////////////

int try_flip(char *ch1a, char *ch2a, int a, char *ch1b, char *ch2b, int b)
{
//want to break unless all of A, C, G and T are present
if((int)ch1a[a]+(int)ch2a[a]+(int)ch1b[b]+(int)ch2b[b]!=287){return(0);}
if(ch1a[a]!='A'&&ch1a[a]!='C'&&ch1a[a]!='G'&&ch1a[a]!='T'){return(0);}
if(ch2a[a]!='A'&&ch2a[a]!='C'&&ch2a[a]!='G'&&ch2a[a]!='T'){return(0);}
if(ch1b[b]!='A'&&ch1b[b]!='C'&&ch1b[b]!='G'&&ch1b[b]!='T'){return(0);}
if(ch2b[b]!='A'&&ch2b[b]!='C'&&ch2b[b]!='G'&&ch2b[b]!='T'){return(0);}

//check that ch1a and ch2a are not ambiguous (A+T or C+G), in which case, nor are ch1b and ch2b
if((int)ch1a[a]+(int)ch2a[a]==138||(int)ch1a[a]+(int)ch2a[a]==149){return(0);}

//then flip ch1a and ch2a (the first pair)
if(ch1a[a]=='A'||ch1a[a]=='T'){ch1a[a]=149-(int)ch1a[a];}
else{ch1a[a]=138-(int)ch1a[a];}
if(ch2a[a]=='A'||ch2a[a]=='T'){ch2a[a]=149-(int)ch2a[a];}
else{ch2a[a]=138-(int)ch2a[a];}
return(1);
}

////////

int check_comp(char **preds, char *al1, char *al2, int a, char **newnames, char *newal1, char *newal2, int b, int *chr, double *bp, int dcount)
{
//previously used to allow flipping
//if(flip==1){try_flip(newal1,newal2,b,al1,al2,a);}

if((al1[a]!=newal1[b]&&al1[a]!=newal2[b])||(al2[a]!=newal1[b]&&al2[a]!=newal2[b]))
{
printf("Error, %s and %s have the same position (Chr%d:%.0f) but different alleles (%c %c and %c %c)\n", preds[a], newnames[b], chr[a], bp[a], al1[a], al2[a], newal1[b], newal2[b]);
//if(try_flip(newal1,newal2,b,al1,al2,a)==1)
//{printf("You can avoid this error by using \"--allow-flips YES\"\n");}
printf("\n");exit(1);
}

if(strcmp(preds[a],newnames[b])!=0)	//names mismatch, use first unless only second rs
{
if(dcount<10){printf("Chr%d:%.0f is referred to as both %s and %s; ", chr[a], bp[a], preds[a], newnames[b]);}
if((preds[a][0]!='r'||preds[a][1]!='s')&&(newnames[a][0]=='r'&&newnames[a][1]=='s'))
{
if(dcount<10){printf("will use latter\n");free(preds[a]);copy_string(preds,a,newnames[b]);}
}
else
{
if(dcount<10){printf("will use former\n");}
}
dcount=dcount+1;
}	//end of names mismatch

return(dcount);
}	//end of check_comp

///////////////////////////

double compute_info(double *data, float *p0, float *p1, int *indexer, int length, double missingvalue)
{
//store info = cov(data,truth)^2/var(data)var(truth)
//to allow for (slight) rounding error, use 1-p0-p1 in place of p2
int i, want, indcount;
double sumd, sumt, meand, meant, vard, vart, vardt;

sumd=0;sumt=0;indcount=0;
for(i=0;i<length;i++)
{
want=i;
if(indexer!=NULL){want=indexer[i];}
if(data[want]!=missingvalue)
{sumd+=data[want];sumt+=p1[want]+2*p0[want];indcount++;}
}
if(indcount<2){return(0);}
meand=sumd/indcount;meant=sumt/indcount;

vard=0;vart=0;vardt=0;
for(i=0;i<length;i++)
{
want=i;
if(indexer!=NULL){want=indexer[i];}
if(data[want]!=missingvalue)
{
vard+=pow(data[want]-meand,2);
vart+=(1-p0[want]-p1[want])*pow(meant,2)+p1[want]*pow(1-meant,2)+p0[want]*pow(2-meant,2);
vardt+=(data[want]-meand)*(p1[want]+2*p0[want]-meant);
}
else	//fill up probabilities based on hwe
{
vart+=meant*(1-meant/2);
}
}
if(vard==0||vart==0){return(0);}

return(pow(vardt,2)/vard/vart);
}	//end of compute_info

///////////////////////////

int check_copy(char **str, int i, char *name, char *mvalue, char *id1, char *id2, char *filename, int type)
{
if(strcmp(name,mvalue)!=0)	//suggesting something non-missing
{
if(strcmp(str[i],name)!=0&&strcmp(str[i],mvalue)!=0)
{
if(type==1)
{printf("Warning reading %s; the PID (%s) for Sample %s %s does not match that previously seen (%s)\nAll PIDs will be set to missing (0)\n\n",filename,name,id1,id2,str[i]);}
if(type==2)
{printf("Warning reading %s; the MID (%s) for Sample %s %s does not match that previously seen (%s)\nAll MIDs will be set to missing (0)\n\n",filename,name,id1,id2,str[i]);}
if(type==3)
{printf("Warning reading %s; the sex (%s) for Sample %s %s does not match that previously seen (%s)\nAll sexes will be set to missing (0)\n\n",filename,name,id1,id2,str[i]);}
if(type==4)
{printf("Warning reading %s; the phenotype (%s) for Sample %s %s does not match that previously seen (%s)\nAll phenotypes will be set to missing (NA)\n\n",filename,name,id1,id2,str[i]);}
return(1);
}
else{free(str[i]);copy_string(str,i,name);}
}

return(0);
}

///////////////////////////

void rdata_warn(int length, int ns)
{
double size;

size=(double)length/1024*ns/1024*8/1024;
if(size>1)
{printf("Warning, to store the region predictors requires %.1f Gb\n\n", size);}
}

////////

void kin_warn(int num_kins, int ns, int type, int flag)
//flag=0 double, flag=1 float
{
double size;

if(flag==0){size=(double)num_kins*ns/1024*ns/1024*8/1024;}
else{size=(double)num_kins*ns/1024*ns/1024*4/1024;}
if(size>1)
{
printf("Warning, to store kinships requires %.1f Gb", size);
if(type==1){printf("; to reduce use \"--memory-save YES\"");}
printf("\n\n");
}
}

////////

void eigen_warn(int ns)
{
double size;

size=(double)ns/1024*ns/1024*8/1024;
if(size>1)
{printf("Warning, to store the eigen-decomposition requires %.1f Gb; sorry, this can not be reduced\n\n", size);}
}

////////

void decomp_warn(int length)
{
double size;

size=(double)length/1024/1024*8/1024;
if(size>1)
{printf("Warning, to perform the decomposition requires %.1f Gb; sorry, this can not be reduced\n\n", size);}
}

////////

void data_warn(int length1, int length2)
{
double size;

size=(double)length1/1024*length2/1024*8/1024;
if(size>1)
{printf("Warning, to process the data requires %.1f Gb; sorry this can not be reduced\n\n", size);}
}

////////

void data_warn2(int length1, int length2)
{
double size;

size=(double)length1/1024*length2/1024*8/1024;
if(size>1)
{printf("Warning, to process the data requires %.1f Gb; if this is too high, you can reduce using \"--bit-size\" (currently %d)\n\n", size, length1);}
}

////////

void data_warn3(int length1, int length2)
{
double size;

size=(double)length1/1024*length2/1024*8/1024;
if(size>1)
{printf("Warning, to process the data requires %.1f Gb; if this is too high, you can reduce using \"--bit-size\"\n\n", size);}
}

////////

void anal_warn(int length1, int length2)
{
double size;

size=(double)length1/1024*length2/1024*8/1024;
if(size>1)
{printf("Warning, to perform the analysis requires %.1f Gb\n\n", size);}
}

void anal_warn2(int length1, int length2)
{
double size;

size=(double)length1/1024*length2/1024*8/1024;
if(size>1)
{printf("Warning, to perform the analysis requires %.1f Gb; if this is too high, you can reduce using \"--bit-size\" (currently %d)\n\n", size, length1);}
}

////////

void model_warn(int length1, int length2)
{
double size;

size=(double)length1/1024*length2/1024*8/1024;
if(size>1)
{printf("Warning, to store the heritability model requires %.1f Gb\n\n", size);}
}

///////////////////////////

int binary_search(double fraction, int length, double *cumsum, int flag)
//if flag=1, then must move down values
{
int j, try, lower, upper, count;
double value;


if(fraction<0||fraction>cumsum[length])
{printf("Error, fraction (%f\n) is not between 0 and %f, please tell Doug\n\n", fraction, cumsum[length]);exit(1);}

if(cumsum[0]!=0){printf("Error, first cumsum is not zero, please tell Doug\n\n");exit(1);}

lower=0;
upper=length;
count=0;
while(1)
{
try=lower+(upper-lower)/2;

if(fraction<cumsum[try])	//move left
{
upper=try;
try=lower+(upper-lower)/2;
}

if(fraction>=cumsum[try+1])	//move right
{
lower=try;
try=lower+(upper-lower)/2;
}

if(fraction>=cumsum[try]&&fraction<cumsum[try+1]){break;}

count++;
if(count>100){printf("Error, can not find %f out of %f with %d - last try %d boundaries %f %f\n\n", fraction, cumsum[length], length, try+1, cumsum[try], cumsum[try+1]);exit(1);}
}

if(flag==1)	//move down values
{
value=cumsum[try+1]-cumsum[try];
for(j=try;j<length;j++){cumsum[j+1]-=value;}
}

return(try);
}

///////////////////////////



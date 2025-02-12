/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Lite version of he_reg, used for getting starting heritabilities for REML

///////////////////////////

double he_part(double *DTD, double *DTS, int start, int end, int ns, int *indexer, int *indexer2, int num_kins, int num_regs, double *Yadj, double **Mkins, float **Mkins_single, double *kintraces, double *X, int *Xstarts, int *Xends, double *Xsums, int memsave, char **kinstems, int thread, int maxthreads)
{
size_t scount, stotal, smark;
int i, i2, j, k, k2, r, total, kcount;
float value, *datatemp;

double Ysumsq, *Dtemp;
char filename[500];
FILE **Minputs;


//set total, stotal
total=num_kins+num_regs;
stotal=(size_t)ns*(size_t)(ns-1)/2;

//allocate Dtemp
Dtemp=malloc(sizeof(double)*ns*total);

if(num_kins>0&&memsave==1)	//allocate datatemp, then open and check kinships
{
sprintf(filename, "%s.grm.id", kinstems[0]);
kcount=countrows(filename);
datatemp=malloc(sizeof(float)*kcount);

Minputs=malloc(sizeof(FILE *)*num_kins);
for(k=0;k<num_kins;k++)
{
sprintf(filename, "%s.grm.bin", kinstems[k]);
if((Minputs[k]=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n", filename);exit(1);}
fseeko(Minputs[k], 0, SEEK_END);
if(ftello(Minputs[k])!=(off_t)sizeof(float)*kcount*(kcount+1)/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*kcount*(kcount+1)/2, ftello(Minputs[k]));exit(1);}
}	//end of k loop
}

//set values to zero
Ysumsq=0;
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTD[k+k2*total]=0;}
DTS[k]=0;
}

//loop through, filling a sample's worth of contributions to DTSa, DTSb, DTDa and DTDb at a time
scount=0;
for(i=start;i<end;i++)
{
//load up first i-1 rows of Dtemp with scaled kinships
for(k=0;k<num_kins;k++)
{
if(memsave==1)	//fill up datatemp - kinships for i start at i*(i+1)/2
{
smark=sizeof(float)*indexer[i]*(indexer[i]+1)/2;
fseeko(Minputs[k], smark, SEEK_SET);
fread(datatemp, sizeof(float), indexer[i], Minputs[k]);
for(i2=0;i2<i;i2++){Dtemp[i2+k*ns]=datatemp[indexer[i2]]/kintraces[k];}
}
else	//kinships already stored
{
if(Mkins!=NULL)
{
for(i2=0;i2<i;i2++){Dtemp[i2+k*ns]=Mkins[k][(size_t)indexer2[i]*ns+indexer2[i2]]/kintraces[k];}
}
else
{
for(i2=0;i2<i;i2++){Dtemp[i2+k*ns]=Mkins_single[k][(size_t)indexer2[i]*ns+indexer2[i2]]/kintraces[k];}
}
}
}

for(r=0;r<num_regs;r++)	//will compute kinships manually (probably faster to create a temporary matrix and vector, then use dgemv)
{
for(i2=0;i2<i;i2++)
{
value=0;for(j=Xstarts[r];j<Xends[r];j++){value+=X[indexer2[i]+j*ns]*X[indexer2[i2]+j*ns];}
Dtemp[i2+(num_kins+r)*ns]=value/Xsums[r];
}
}

//now loop through pairs
for(i2=0;i2<i;i2++)
{
if(thread==0&&scount%100000000==0&&scount*maxthreads<stotal)
{printf("Processing Sample Pair %jd of %jd\n", scount*maxthreads+1, stotal);}

Ysumsq+=pow(Yadj[indexer2[i]]*Yadj[indexer2[i2]],2);
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTD[k+k2*total]+=Dtemp[i2+k*ns]*Dtemp[i2+k2*ns];}
DTS[k]+=Dtemp[i2+k*ns]*Yadj[indexer2[i]]*Yadj[indexer2[i2]];
}
scount++;
}}	//end i2 and i loop
if(thread==0){printf("\n");}

if(num_kins>0&&memsave==1)	//can shut kinships
{
for(k=0;k<num_kins;k++){fclose(Minputs[k]);}
}

free(Dtemp);
if(num_kins>0&&memsave==1){free(datatemp);free(Minputs);}

return(Ysumsq);
}

////////

int he_starts(double *hers, int ns, int num_covars, int num_envs, int num_tops, int num_kins, int num_regs, double *Y, double *Z, double **Mkins, float **Mkins_single, double *kintraces, double *X, int Xtotal, int *Xstarts, int *Xends, double *Xsums, int memsave, int maxthreads, char **kinstems, char **ids3, double missingvalue)
{
int i, j, k, k2, r, count, kcount;
int thread, threadstart, threadend;
double sum, sumsq;
char **wantids;

int num_fixed, total, *indexer, *indexer2;
double *thetas, *Yadj;

double Ysumsq, *DTD, *DTS, *DTD2, *MYsumsq, **MDTD, **MDTS;
double *Xsave, *Xsumssave;

char filename[500];


//set num_fixed and total
num_fixed=num_covars+num_envs+num_tops;
total=num_kins+num_regs;

//allocate variables

indexer=malloc(sizeof(int)*ns);
indexer2=malloc(sizeof(int)*ns);

thetas=malloc(sizeof(double)*num_fixed);
Yadj=malloc(sizeof(double)*ns);

DTD=malloc(sizeof(double)*total*total);
DTS=malloc(sizeof(double)*total);
DTD2=malloc(sizeof(double)*total);

MYsumsq=malloc(sizeof(double)*maxthreads);
MDTD=malloc(sizeof(double *)*maxthreads);
MDTS=malloc(sizeof(double *)*maxthreads);
for(thread=0;thread<maxthreads;thread++)
{MDTD[thread]=malloc(sizeof(double)*total*total);MDTS[thread]=malloc(sizeof(double)*total);}

if(num_regs>0)	//allocate Xsave and Xsums save, and use to save X and Xsums
{
Xsave=malloc(sizeof(double)*ns*Xtotal);
Xsumssave=malloc(sizeof(double)*num_regs);

for(j=0;j<Xtotal;j++)
{
for(i=0;i<ns;i++){Xsave[i+j*ns]=X[i+j*ns];}
}

for(r=0;r<num_regs;r++){Xsumssave[r]=Xsums[r];}
}

if(num_kins>0&&memsave==1)	//set indexers based on first id file (if nk>1, have checked all id files match)
{
sprintf(filename, "%s.grm.id", kinstems[0]);
kcount=countrows(filename);
wantids=malloc(sizeof(char*)*kcount);
read_ids(filename, NULL, NULL, wantids, kcount, NULL, 0, 0);
count=find_strings(wantids, kcount, ids3, ns, indexer, indexer2, NULL, NULL, NULL, NULL, 3);
if(count!=ns){printf("Doug Error %d %d\n", count, ns);exit(1);}
for(i=0;i<kcount;i++){free(wantids[i]);}free(wantids);
}
else	//order of kinships matches that of phenotypes and regions
{
for(i=0;i<ns;i++){indexer[i]=i;indexer2[i]=i;}
}

//solve covariates (get covher and topher) and fill Yadj
reg_covar_lin(Y, Z, ns, num_covars+num_envs, num_tops, thetas, NULL, NULL, Yadj, 1, NULL, NULL);

if(num_regs>0)	//regress covariates out of regions and update Xsums
{
reg_covar_matrix(X, Z, ns, Xtotal, num_fixed);

for(r=0;r<num_regs;r++)
{
sumsq=0;
for(j=Xstarts[r];j<Xends[r];j++)
{
for(i=0;i<ns;i++){sumsq+=pow(X[i+j*ns],2);}
}
Xsums[r]=sumsq/ns;
}
}

////////

//ready to start - S = products, D = [K1 K2 etc]

#pragma omp parallel for private(thread,threadstart,threadend) schedule (static, 1)
for(thread=0;thread<maxthreads;thread++)
{
threadstart=pow((double)thread/maxthreads,.5)*ns;
threadend=pow((double)(thread+1)/maxthreads,.5)*ns;

MYsumsq[thread]=he_part(MDTD[thread], MDTS[thread], threadstart, threadend, ns, indexer, indexer2, num_kins, num_regs, Yadj, Mkins, Mkins_single, kintraces, X, Xstarts, Xends, Xsums, memsave, kinstems, thread, maxthreads);
}

//sum up
Ysumsq=0;
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTD[k+k2*total]=0;}
DTS[k]=0;
}

for(thread=0;thread<maxthreads;thread++)
{
Ysumsq+=MYsumsq[thread];
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){DTD[k+k2*total]+=MDTD[thread][k+k2*total];}
DTS[k]+=MDTS[thread][k];
}
}

//solve and save hers (first value corresponds to noise)
(void)eigen_invert(DTD, total, DTD2, 1, DTS, 1);
sum=0;for(k=0;k<total;k++){sum+=DTS[k];}
hers[0]=1-sum;
for(k=0;k<total;k++){hers[1+k]=DTS[k];}

free(indexer);free(indexer2);
free(thetas);free(Yadj);
free(DTD);free(DTS);free(DTD2);

if(num_regs>0)	//restore X and Xsums
{
for(j=0;j<Xtotal;j++)
{
for(i=0;i<ns;i++){X[i+j*ns]=Xsave[i+j*ns];}
}
for(r=0;r<num_regs;r++){Xsumssave[r]=Xsums[r];}
free(Xsave);free(Xsumssave);
}

return(0);
}	//end of he_starts

///////////////////////////


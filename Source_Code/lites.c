/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

//Construct prediction models from summary statistics - code allows for exps[j] but dont think this can ever happen

///////////////////////////

if(dougvar==1){printf("MAMA omega\n");}
if(dougvar==2){printf("extreme weights\n");}
if(dougvar==3){printf("only use focal correlations\n");}
if(dougvar==4){printf("do not update psi and delta\n");}
if(dougvar==5){printf("do not scale\n");}
if(dougvar==6){printf("wrong meta\n");}

if(dougvar==7){printf("always sample res var\n");}
if(dougvar==8){printf("purity filter\n");}
if(dougvar==9){printf("no blank\n");}
if(dougvar==10){printf("dopnt scale second\n");}


if(dougvar2!=-9999){printf("scaling is %f\n", dougvar2);}


//currently the correlations use five powers (and I do not plan to change this)
num_pows=5;
powers=malloc(sizeof(double)*num_pows);
powers[0]=-1.0;powers[1]=-0.75;powers[2]=-0.5;powers[3]=-0.25;powers[4]=0;

//set variance grid for use with finemappin
vargrid=malloc(sizeof(double)*num_grids);
vargrid[0]=0;
value=pow(grid_max/grid_min,1.0/(num_grids-2));
for(k=1;k<num_grids;k++){vargrid[k]=grid_min*pow(value,k-1);}

//set sscale to 1 (will overwrite if estimating)
sscale=1;

//set escale to one (will overwrite if using cv)
escale=1.0;

if(finemap==0)   //will consider multiple models
{
#include "setparamsb.c"
#include "setlambdas.c"
}
else    //fine-mapping uses only one model
{num_try=1;}

if(prsvar==1)	//work out how many samplings will be recorded
{
num_blocks=num_mcmc/num_step*num_chains;
if(num_blocks<3){printf("Error, to compute PRS variances, there must be at least three samplings (not %d); either increase \"--MCMC-iterations\" or reduce \"--MCMC-step\" (currently %d and %d, respectively)\n\n", num_blocks, num_mcmc, num_step);exit(1);}
}

if(cur_focal==0)   //make a new progress file
{
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);
}

if(num_focals>1)    //will be using existing progress file
{
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Performing the analysis for Target %d out of %d\n", cur_focal+1, num_focals);
fclose(output);
}

//allocate and load multitrait summary statistics
value=(double)data_length/1024*num_sums3*4/1024/1024*8;
if(value>1){printf("Warning, to store the summary statistics requires %.1f Gb\n\n", value);}

Mnss=malloc(sizeof(double*)*num_sums3);
Mchis=malloc(sizeof(double*)*num_sums3);
Mrhos=malloc(sizeof(double*)*num_sums3);
Ma1freq=malloc(sizeof(double*)*num_sums3);

for(q=1;q<num_sums3;q++)
{
Mnss[q]=malloc(sizeof(double)*data_length);
Mchis[q]=malloc(sizeof(double)*data_length);
Mrhos[q]=malloc(sizeof(double)*data_length);
Ma1freq[q]=malloc(sizeof(double)*data_length);
}

//can just reassign for focal trait
Mnss[0]=nss;Mchis[0]=chis;Mrhos[0]=rhos;Ma1freq[0]=a1freq;

if(num_sums1>1)  //read in other summary statistics (missing will get ns, rhos, chis and freq zero)
{
printf("Reading the other summary statistics (if using multiple threads, these will be read in a random order)\n");
#pragma omp parallel for private(q) schedule(static)
for(q=1;q<num_sums1;q++)
{read_sumsfile(sumstems[q], Mnss[q], Mchis[q], Mrhos[q], Ma1freq[q], data_length, preds, al1, al2, bimfile, amb, 1.0, 3);}
printf("\n");
}

if(num_sums2>0) //read in bulk summary statistics
{
printf("Reading summary statistics for %d traits from %s\n", num_sums2, bulkfile);
read_bulkfile(bulkfile, num_sums2, Mnss+num_sums1, Mchis+num_sums1, Mrhos+num_sums1, Ma1freq+num_sums1, data_length, preds, al1, al2, bimfile);
printf("\n");
}

if(num_sums3>1&&dougvar==2)
{
printf("Checking the non-focal summary statistics for extreme values (predictors that explain over 1%% of phenotypic variation)\n");
for(q=1;q<num_sums3;q++)
{
count=0;
for(j=0;j<data_length;j++)
{
if(Mchis[q][j]>0.01*Mnss[q][j])
{
Mchis[q][j]=0.01*Mnss[q][j];
if(Mrhos[q][j]>0){Mrhos[q][j]=pow(Mchis[q][j]/(Mchis[q][j]+Mnss[q][j]),.5);}
else{Mrhos[q][j]=-pow(Mchis[q][j]/(Mchis[q][j]+Mnss[q][j]),.5);}
count++;
}
}
if(count>0){printf("Warning, %s contains %d extreme values (these have been reduced)\n", sumstems[q], count);}
}
printf("\n");
}

////////

//get an upper limit for bittotal and allocate blocks

sprintf(filename2,"%s.cors.windows", corstems[0]);
total=countrows(filename2)-1;
blockstarts=malloc(sizeof(int)*total);
blockends=malloc(sizeof(int)*total);

//allocate multi-cor variables (and load some)

Dns=malloc(sizeof(int)*num_cors);
Dnp=malloc(sizeof(int)*num_cors);
Dnp2=malloc(sizeof(int)*num_cors);
Dkp=malloc(sizeof(int*)*num_cors);
Dkp2=malloc(sizeof(int*)*num_cors);
Dindexes=malloc(sizeof(off_t*)*num_cors);
Dsizes=malloc(sizeof(int**)*num_cors);
Duse=malloc(sizeof(int**)*num_cors);
Duse2=malloc(sizeof(int**)*num_cors);
Dsigns=malloc(sizeof(int**)*num_cors);
Dchr=malloc(sizeof(int*)*num_cors);
Dpreds=malloc(sizeof(char**)*num_cors);
Dbp=malloc(sizeof(double)*num_cors);
Dal1=malloc(sizeof(char*)*num_cors);
Dal2=malloc(sizeof(char*)*num_cors);

if(num_cors>1){printf("Reading details of the predictor-predictor correlations\n");}

flag=0;
for(s=0;s<num_cors;s++)
{
sprintf(filename2, "%s.cors.root", corstems[s]);
read_integers(filename2, Dns+s, 1, NULL, 2, 3, 0);
read_integers(filename2, Dnp+s, 1, NULL, 2, 4, 0);

Dkp[s]=malloc(sizeof(int)*Dnp[s]);
Dkp2[s]=malloc(sizeof(int)*Dnp[s]);

Dindexes[s]=malloc(sizeof(off_t)*total);
Dsizes[s]=malloc(sizeof(int*)*total);
Duse[s]=malloc(sizeof(int*)*total);
Duse2[s]=malloc(sizeof(int*)*total);
Dsigns[s]=malloc(sizeof(int*)*total);

Dchr[s]=malloc(sizeof(int)*Dnp[s]);
Dpreds[s]=malloc(sizeof(char*)*Dnp[s]);
Dbp[s]=malloc(sizeof(double)*Dnp[s]);
Dal1[s]=malloc(sizeof(char)*Dnp[s]);
Dal2[s]=malloc(sizeof(char)*Dnp[s]);

sprintf(filename2, "%s.cors.bim", corstems[s]);
flag+=read_bimfile_lite(filename2, Dchr[s], Dpreds[s], Dbp[s], Dal1[s], Dal2[s], Dnp[s]);
}

if(flag==1){printf("\n");}
if(flag>1){printf("Warning, the order of predictors with matching basepairs should be consistent across predictor-predictor correlations (else some of the predictors with matching basepairs will be excluded)\n\n");}

//get Dnp2, Dkp and Dkp2 (could set directly for s=sumpops[0]) - also check positions and alleles
for(s=0;s<num_cors;s++)	//find which of Dpreds[s] are in preds
{
j=0;j2=0;Dnp2[s]=0;
while(j<data_length&&j2<Dnp[s])
{
if(strcmp(preds[j],Dpreds[s][j2])==0)
{
if(chr[j]!=Dchr[s][j2]||bp[j]!=Dbp[s][j2])
{printf("Error, the position for %s in %s (%d:%.2f) does not match that in %s (%d:%.2f)\n\n", preds[j], corstems[sumpops[0]], chr[j], bp[j], corstems[s], Dchr[s][j2], Dbp[s][j2]);exit(1);}

if((al1[j]!=Dal1[s][j2]&&al1[j]!=Dal2[s][j2])||(al2[j]!=Dal1[s][j2]&&al2[j]!=Dal2[s][j2]))
{printf("Error, the alleles for %s in %s (%c & %c) are not consistent with those in %s (%c & %c)\n\n", preds[j], corstems[sumpops[0]], al1[j], al2[j], corstems[s], Dal1[s][j2], Dal2[s][j2]);exit(1);}

Dkp[s][Dnp2[s]]=j2;Dkp2[s][Dnp2[s]]=j;
j++;j2++;Dnp2[s]++;
}
else
{
if(chr[j]<Dchr[s][j2]||(chr[j]==Dchr[s][j2]&&bp[j]<Dbp[s][j2])){j++;}
else{j2++;}
}
}
if(Dnp2[s]<data_length)
{printf("Warning, only %d of the %d predictors are in the correlations with stem %s\n", Dnp2[s], data_length, corstems[s]);}
}

if(num_cors>1){printf("\n");}

////////

//cohers2 stores ancestries of focal - only used with metasum=1
if(metasum==1){cohers2=malloc(sizeof(double)*num_cors);}

//compare frequencies of sum stats and correlations (maybe infer ancestry for bulk and meta), then blank missing predictors
#include "freqcheck.c"

//get the heritability model
model_warn(data_length*3/2, num_parts+1);
pindexes=malloc(sizeof(int *)*(num_parts+1));
pweights=malloc(sizeof(double *)*(num_parts+1));

if(useann==0)	//annotations provided manually - if using partions, will always use backpart=1 (ensures no predictors excluded)
{
addpart=get_her_model(num_parts, partpref, pindexes, pweights, data_length, keeppreds_use, num_preds, allpreds, predorder, parttype, 1, allone, -9999);
}
else	//read annotations from the focal correlations (to which num_preds and keeppreds_use correspond)
{
addpart=1;
for(q=0;q<num_parts+addpart;q++)
{pindexes[q]=malloc(sizeof(int)*data_length);pweights[q]=malloc(sizeof(double)*data_length);}

//open annotations file and skip header
sprintf(filename2,"%s.cors.annotations", corstems[sumpops[0]]);
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input2, "%c", &readchar);}

found=0;
for(j=0;j<num_preds;j++)
{
if(j<keeppreds_use[found])	//skip row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input2, "%c", &readchar);}
}
else	//read row
{
if(fscanf(input2, "%s ",readstring)!=1)
{printf("Error reading first element of Row %d of %s\n\n", j+2, filename2);exit(1);}
if(strcmp(readstring,preds[found])!=0){printf("Doug error 89ABC %s %s %d\n\n", readstring, preds[found], found);exit(1);}

for(q=0;q<num_parts+addpart;q++)
{
if(fscanf(input2, "%f ", &readfloat)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+2, j+2, filename2);exit(1);}
pindexes[q][found]=(readfloat!=0);
pweights[q][found]=readfloat;
}
found++;
}
if(found==data_length){break;}
}

fclose(input2);
}

////////

//work out which predictors to remove
usedpreds=malloc(sizeof(int)*data_length);

if(impsums==0)    //remove predictors with missing summary statistics (for focal trait)
{
for(j=0;j<data_length;j++){usedpreds[j]=(nss[j]>0);}
}
else    //only remove predictors in windows missing more than three-quarters of summary statistics
{
for(j=0;j<data_length;j++){usedpreds[j]=1;}

sprintf(filename2,"%s.cors.windows", corstems[sumpops[0]]);
total=countrows(filename2)-1;
bittotal=get_windows(filename2, data_length, keeppreds_use, blockstarts, blockends, total);

wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

count=0;for(j=bitstart;j<bitend;j++){count+=(nss[j]>0);}
if(4*count<bitlength)
{
for(j=bitstart;j<bitend;j++){usedpreds[j]=(nss[j]>0);}

if(wcount<5)
{
if(count==0){printf("Warning, Window %d will be excluded because none of its predictors have summary statistics\n", bit+1);}
else{printf("Warning, Window %d will not be imputed because fewer than a quarter of its predictors (%d out of %d) have summary statistics\n", bit+1, count, bitlength);}
}
wcount++;
}
}
if(wcount>5){printf("In total, %d windows have summary statistics for fewer than a quarter their predictors\n", wcount);}
if(wcount>0){printf("\n");}
}

if(finemap==1&&num_cors>1)  //can only use predictors present in all (active) correlations
{
count=0;for(j=0;j<data_length;j++){count+=usedpreds[j];}

for(q=0;q<num_sums3;q++)    //code a bit redundant if there are multiple summary statistics for non-focal correlations
{
s=sumpops[q];
if(s!=sumpops[0])   //active correlations that are not focal
{
//subtract one from all usedpreds, then add on one for overlapping predictors
for(j=0;j<data_length;j++){usedpreds[j]--;}
for(j=0;j<Dnp2[s];j++){usedpreds[Dkp2[s][j]]++;}
}
}

count2=0;for(j=0;j<data_length;j++){count2+=(usedpreds[j]==1);}
if(count2<count){printf("Warning, %d predictors are excluded because they are not present in all correlations\n\n", count-count2);}
}

//squeeze down predictors
count=0;
for(j=0;j<data_length;j++)
{
if(usedpreds[j]>0)
{
if(count!=j)
{
chr[count]=chr[j];cm[count]=cm[j];bp[count]=bp[j];cmbp[count]=cmbp[j];al1[count]=al1[j];al2[count]=al2[j];
free(preds[count]);copy_string(preds,count,preds[j]);
free(along1[count]);copy_string(along1,count,along1[j]);free(along2[count]);copy_string(along2,count,along2[j]);
keeppreds_use[count]=keeppreds_use[j];

nss[count]=nss[j];chis[count]=chis[j];rhos[count]=rhos[j];a1freq[count]=a1freq[j];
for(q=1;q<num_sums3;q++){Mnss[q][count]=Mnss[q][j];Mchis[q][count]=Mchis[q][j];Mrhos[q][count]=Mrhos[q][j];Ma1freq[q][count]=Ma1freq[q][j];}
for(q=0;q<num_parts;q++){pindexes[q][count]=pindexes[q][j];pweights[q][count]=pweights[q][j];}
}
count++;
}}
for(j=count;j<data_length;j++){free(preds[j]);free(along1[j]);free(along2[j]);}

if(count<data_length)   //lost some predictors
{
data_length=count;

//get Dnp2, Dkp and Dkp2 again
for(s=0;s<num_cors;s++)
{
j=0;j2=0;Dnp2[s]=0;
while(j<data_length&&j2<Dnp[s])
{
if(strcmp(preds[j],Dpreds[s][j2])==0)
{
Dkp[s][Dnp2[s]]=j2;Dkp2[s][Dnp2[s]]=j;
j++;j2++;Dnp2[s]++;
}
else
{
if(chr[j]<Dchr[s][j2]||(chr[j]==Dchr[s][j2]&&bp[j]<Dbp[s][j2])){j++;}
else{j2++;}
}
}
}
}
if(data_length!=Dnp2[sumpops[0]]){printf("Error 127XY, please tell Doug %d %d\n\n", data_length, Dnp2[sumpops[0]]);exit(1);}

free(usedpreds);

//can save some space by freeing Dpreds
for(s=0;s<num_cors;s++)
{
for(j=0;j<Dnp[s];j++){free(Dpreds[s][j]);}
}

printf("The analysis will be performed using %d predictors\n\n", data_length);

////////

//read centres, means, variances, sets of rjksums, datarands and highlds for focal correlations (have already allocated first three)
rjksums=malloc(sizeof(double)*data_length*num_pows);
datarands=malloc(sizeof(double)*data_length);
highlds=malloc(sizeof(int)*data_length);

sprintf(filename2,"%s.cors.bin", corstems[sumpops[0]]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read centres, allowing for filtering 
for(j=0;j<data_length;j++)
{
fseeko(input2, sizeof(double)*keeppreds_use[j], SEEK_SET);
if(fread(centres+j, sizeof(double), 1, input2)!=1)
{printf("Error reading mean for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
}

//read scalings, allowing for filtering
for(j=0;j<data_length;j++)
{
fseeko(input2, sizeof(double)*num_preds+sizeof(double)*keeppreds_use[j], SEEK_SET);
if(fread(mults+j, sizeof(double), 1, input2)!=1)
{printf("Error reading scaling for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
}

//read variances, allowing for filtering
for(j=0;j<data_length;j++)
{
fseeko(input2, sizeof(double)*num_preds*2+sizeof(double)*keeppreds_use[j], SEEK_SET);
if(fread(sqdevs+j, sizeof(double), 1, input2)!=1)
{printf("Error reading variance for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
}

//read sets of rjksums, allowing for filtering and possible missing
for(k=0;k<num_pows;k++)
{
for(j=0;j<data_length;j++)
{
fseeko(input2, sizeof(double)*num_preds*(3+k)+sizeof(double)*keeppreds_use[j], SEEK_SET);
if(fread(rjksums+j+k*data_length, sizeof(double), 1, input2)!=1)
{printf("Error reading Tagging %d for Predictor %d from %s\n\n", k+1, j+1, filename2);exit(1);}
}
}

//read datarands, allowing for filtering
for(j=0;j<data_length;j++)
{
fseeko(input2, sizeof(double)*num_preds*8+sizeof(double)*keeppreds_use[j], SEEK_SET);
if(fread(datarands+j, sizeof(double), 1, input2)!=1)
{printf("Error reading noise for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
}

//read highlds, allowing for filtering
for(j=0;j<data_length;j++)
{
fseeko(input2, sizeof(double)*num_preds*9+sizeof(int)*keeppreds_use[j], SEEK_SET);
if(fread(highlds+j, sizeof(int), 1, input2)!=1)
{printf("Error reading LD status for Predictor %d from %s\n\n", j+1, filename2);exit(1);}
}

fclose(input2);

////////

//work out number of blocks and their boundaries

sprintf(filename2,"%s.cors.windows", corstems[sumpops[0]]);
total=countrows(filename2)-1;
bittotal=get_windows(filename2, data_length, keeppreds_use, blockstarts, blockends, total);

//work out window indexes for each set of correlations

for(s=0;s<num_cors;s++)
{
sprintf(filename2,"%s.cors.windows", corstems[s]);
count=countrows(filename2)-1;

indexer=malloc(sizeof(int)*count);
indexer2=malloc(sizeof(int)*count);
read_integers(filename2, indexer, count, NULL, 2, 1, 0);
read_integers(filename2, indexer2, count, NULL, 3, 1, 0);
for(j=0;j<count;j++){indexer[j]--;}

//get indexes of all windows
Dtemp=malloc(sizeof(off_t)*count);
Dtemp[0]=(off_t)sizeof(double)*Dnp[s]*9+sizeof(int)*Dnp[s];
for(j=1;j<count;j++)
{
count2=indexer2[j-1]-indexer[j-1];
Dtemp[j]=Dtemp[j-1]+sizeof(unsigned short)*count2*(count2-1)/2;
}

start=0;
found=0;
for(bit=0;bit<bittotal;bit++)
{
if(start<Dnp2[s])	//there are some predictors remaining
{
while(Dkp2[s][start]<blockstarts[bit]&&start<Dnp2[s]){start++;}
}

if(start==Dnp2[s])	//must be empty
{end=start;}
else	//might contain predictors
{
if(Dkp2[s][start]>=blockends[bit])	//empty
{end=start;}
else
{
end=start;
while(Dkp2[s][end]<blockends[bit])
{
end++;
if(end==Dnp2[s]){break;}
}
}
}

if(end>start)	//window contains predictors
{
Dsizes[s][bit]=malloc(sizeof(int)*2);
Duse[s][bit]=malloc(sizeof(int)*(end-start));
Duse2[s][bit]=malloc(sizeof(int)*(end-start));
Dsigns[s][bit]=malloc(sizeof(int)*(end-start));

//work out which window we are in
while(Dkp[s][start]>=indexer2[found]){found++;}

if(Dkp[s][end-1]>=indexer2[found]){printf("Error BB33, please tell Doug\n\n");exit(1);}
Dindexes[s][bit]=Dtemp[found];

//first element of Dsizes is original number of predictors in window
Dsizes[s][bit][0]=indexer2[found]-indexer[found];

//second element of Dsizes is number of predictors in window being used
Dsizes[s][bit][1]=end-start;

//Duse indexes predictors being used (relative to start of block)
for(j=start;j<end;j++){Duse[s][bit][j-start]=Dkp[s][j]-indexer[found];}

//Duse2 indexes where predictors map (relative to start of block)
for(j=start;j<end;j++){Duse2[s][bit][j-start]=Dkp2[s][j]-blockstarts[bit];}

//Dsigns says whether we should flip correlation (first check consistent)
for(j=start;j<end;j++)
{
if(al1[Dkp2[s][j]]==Dal1[s][Dkp[s][j]]){Dsigns[s][bit][j-start]=1;}
else{Dsigns[s][bit][j-start]=-1;}
}
}
else	//contains no predictors
{
Dindexes[s][bit]=0;

Dsizes[s][bit]=malloc(sizeof(int)*2);
Duse[s][bit]=malloc(sizeof(int));
Duse2[s][bit]=malloc(sizeof(int));
Dsigns[s][bit]=malloc(sizeof(int));

Dsizes[s][bit][0]=0;
Dsizes[s][bit][1]=0;
Duse[s][bit][0]=0;
Duse2[s][bit][0]=0;
Dsigns[s][bit][0]=0;
}
}	//end of bit loop

free(indexer);free(indexer2);
free(Dtemp);
}	//end of s loop

//get bitmax (max original size, across all correlations)
bitmax=Duse[0][0][0];
for(s=0;s<num_cors;s++)
{
for(bit=0;bit<bittotal;bit++)
{
if(Dsizes[s][bit][0]>bitmax){bitmax=Dsizes[s][bit][0];}
}
}

///////////////////////

//allocate some general variables

value=(double)bitmax/1024*bitmax/1024/1024*10;
if(value>1){printf("Warning, to store the correlations requires %.1f Gb\n\n", value);}

exps=malloc(sizeof(double)*data_length);
cors_short=malloc(sizeof(unsigned short)*bitmax*bitmax);
cors=malloc(sizeof(double)*bitmax*bitmax);

if(num_sums3>1)
{
nss4=malloc(sizeof(double)*data_length);
chis4=malloc(sizeof(double)*data_length);
rhos4=malloc(sizeof(double)*data_length);
}
if(num_sums3>1||metasum==1){cohers=malloc(sizeof(double)*num_cors*bittotal);}

////////

//get best power for focal trait (from the preset options) using sumher (without correction for bias) - will fudge if power provided

rjksums2=malloc(sizeof(double)*data_length);
snss=malloc(sizeof(double)*data_length);
schis=malloc(sizeof(double)*data_length);
stats=malloc(sizeof(double)*12*num_pows);
likes=malloc(sizeof(double)*11*num_pows);
svars=malloc(sizeof(double*));
ssums=malloc(sizeof(double*));
svars[0]=malloc(sizeof(double)*data_length);
ssums[0]=malloc(sizeof(double)*3);

if(noscale==0){printf("Testing five values for the power parameter (-1, -0.75, -0.5, -0.25 and 0), and estimating the inflation of test statistics\n");}
else{printf("Testing five values for the power parameter (-1, -0.75, -0.5, -0.25 and 0)\n");}

for(k=0;k<num_pows;k++)	//solve for kth power and extract likelihood - have not yet imputed, so might be missing values
{
ssums[0][0]=0;ssums[0][1]=data_length;ssums[0][2]=1;
count=0;
for(j=0;j<data_length;j++)
{
if(nss[j]>0)
{
rjksums2[count]=rjksums[j];
snss[count]=nss[j];
schis[count]=chis[j];
if(schis[count]<1e-6){schis[count]=1e-6;}
if(schis[count]>100){schis[count]=100;}
svars[0][count]=rjksums[j+k*data_length];
if(hwestand==1){ssums[0][0]+=weights[j]*pow(centres[j]*(1-centres[j]/2),1+powers[k]);}
else{ssums[0][0]+=weights[j]*pow(sqdevs[j],1+powers[k]);}
count++;
}
}

if(noscale==0){solve_sums(stats+k*12, likes+k*11, NULL, NULL, NULL, 1, 1, 0, -9999, count, 0, NULL, NULL, rjksums2, svars, ssums, snss, schis, 0.001, 100, 1, 7, NULL);}
else{solve_sums(stats+k*12, likes+k*11, NULL, NULL, NULL, 1, 0, 0, -9999, count, 0, NULL, NULL, rjksums2, svars, ssums, snss, schis, 0.001, 100, 1, 7, NULL);}
}

best=0;
for(k=1;k<num_pows;k++)
{
if(likes[1+k*11]>likes[1+best*11]){best=k;}
}

if(power!=-9999)    //compare with provided value, and reset best
{
if(powers[best]!=power){printf("Warning, the best-fitting power is %.2f, which does not match the value provided (%.2f)\n\n", powers[best], power);}
for(k=0;k<num_pows;k++)
{
if(powers[k]==power){best=k;break;}
}
}
power=powers[best];

if(noscale==1){printf("The best-fitting power is %.2f\n\n", power);}
else    //scale focal, then secondary statistics
{
sscale=stats[1+best*12];
printf("The best-fitting value is %.2f, while the estimated inflation is %f\n\n", power, sscale);
if(sscale<0.8){printf("Warning, the scaling is very low, so has been increased to 0.8\n\n");sscale=0.8;}
if(sscale>8){printf("Warning, the scaling is high low, so has been reduced to 8\n\n");sscale=8;}

if(dougvar2!=-9999){sscale=dougvar2;}

value=pow(sscale,-1);
for(j=0;j<data_length;j++)
{
if(nss[j]>0)
{
chis[j]*=value;
if(rhos[j]>0){rhos[j]=pow(chis[j]/(chis[j]+nss[j]),.5);}
else{rhos[j]=-pow(chis[j]/(chis[j]+nss[j]),.5);}
}
}

if(num_sums3>1&&dougvar!=10)
{
//read rjksums for best power for all correlations
rjksums3=malloc(sizeof(double)*data_length*num_cors);

for(s=0;s<num_cors;s++)
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

for(j=0;j<data_length;j++){rjksums3[j+s*data_length]=0;}
for(j=0;j<Dnp2[s];j++)
{
fseeko(input2, sizeof(double)*Dnp[s]*(3+best)+sizeof(double)*Dkp[s][j], SEEK_SET);
if(fread(rjksums3+Dkp2[s][j]+s*data_length, sizeof(double), 1, input2)!=1)
{printf("Error reading Tagging %d for Predictor %d from %s\n\n", best+1, j+1, filename2);exit(1);}
}

fclose(input2);
}

//estimate scaling and scale
for(q=1;q<num_sums3;q++)
{
//values of ssums do not matter as only care about scaling
ssums[0][0]=1;ssums[0][1]=data_length;ssums[0][2]=1;
count=0;
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0)
{
rjksums2[count]=rjksums[j];
snss[count]=Mnss[q][j];
schis[count]=Mchis[q][j];
if(schis[count]<1e-6){schis[count]=1e-6;}
if(schis[count]>100){schis[count]=100;}
svars[0][count]=rjksums3[j+sumpops[q]*data_length];
count++;
}
}

solve_sums(stats, likes, NULL, NULL, NULL, 1, 1, 0, -9999, count, 0, NULL, NULL, rjksums2, svars, ssums, snss, schis, 0.001, 100, 1, 7, NULL);
sscale2=stats[1];

if(sscale2<0.8){printf("Warning, the scaling is very low, so has been increased to 0.8\n\n");sscale2=0.8;}
if(sscale2>8){printf("Warning, the scaling is high low, so has been reduced to 8\n\n");sscale2=8;}

value=pow(sscale2,-1);
for(j=0;j<data_length;j++)
{
if(Mnss[q][j]>0)
{
Mchis[q][j]*=value;
if(Mrhos[q][j]>0){Mrhos[q][j]=pow(Mchis[q][j]/(Mchis[q][j]+Mnss[q][j]),.5);}
else{Mrhos[q][j]=-pow(Mchis[q][j]/(Mchis[q][j]+Mnss[q][j]),.5);}
}
}
}

free(rjksums3);
}
}

free(rjksums2);free(snss);free(schis);free(stats);free(likes);
free(svars[0]);free(svars);free(ssums[0]);free(ssums);

///////////////////////

if(her!=-9999)	//have her (and power), so only need to fill exps
{
if(strcmp(indhers,"blank")!=0)	//exps are stored in weights (and they sum to one)
{
for(j=0;j<data_length;j++){exps[j]=weights[j];}
}
else	//must have simple case (no partitions)
{
for(j=0;j<data_length;j++)
{
if(hwestand==1){exps[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps[j]=weights[j]*pow(sqdevs[j],1+power);}
}

sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}
}
}
else    //estimate heritability for focal trait using sumher - have not yet imputed, so might be missing values
{
printf("Estimating heritability using an approximate version of SumHer\n");
total=num_parts+addpart;

//taggings require dl x (2 + total) while sum_hers require dl x (3 + 2 x total)
value=(double)data_length/1024*(5+3*total)/1024/1024*8;
if(value>1){printf("Warning, to store the taggings requires %.1f Gb\n", value);}
printf("\n");

exps2=malloc(sizeof(double)*data_length);
pweights2=malloc(sizeof(double)*bitmax*total);
cors3=malloc(sizeof(double)*bitmax*total);

snss=malloc(sizeof(double)*data_length);
schis=malloc(sizeof(double)*data_length);
stats=malloc(sizeof(double)*(total+2+total)*3);
stats2=malloc(sizeof(double)*12);

svars=malloc(sizeof(double*)*total);
ssums=malloc(sizeof(double*)*total);
for(q=0;q<total;q++)
{
svars[q]=malloc(sizeof(double)*data_length);
ssums[q]=malloc(sizeof(double)*(total+2));
}

svars2=malloc(sizeof(double*));
ssums2=malloc(sizeof(double*));
svars2[0]=malloc(sizeof(double)*data_length);
ssums2[0]=malloc(sizeof(double)*3);

//load up exps2, making sure it sums to one
for(j=0;j<data_length;j++)
{
if(hwestand==1){exps2[j]=weights[j]*pow(centres[j]*(1-centres[j]/2),1+power);}
else{exps2[j]=weights[j]*pow(sqdevs[j],1+power);}
}

sum=0;for(j=0;j<data_length;j++){sum+=exps2[j];}
value=pow(sum,-1);
for(j=0;j<data_length;j++){exps2[j]*=value;}

//re-open focal correlations
s=sumpops[0];

sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

count=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(bit%100==0)
{
printf("Calculating taggings for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating taggings for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read lower triangle of focal correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//now extract squared correlations, correcting off-diagonals for sample size (no need for signs, and not using shrink)
value=(double)(Dns[s]-1)/(Dns[s]-2);
value2=pow(Dns[s]-2,-1);
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<bitlength;j++)
{cors[(size_t)k*bitlength+j]=pow(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1,2)*value-value2;}
}

//load up annotations for this bit
for(q=0;q<total;q++)
{
for(j=0;j<bitlength;j++){pweights2[j+q*bitlength]=pweights[q][bitstart+j]*exps2[bitstart+j];}
}

//get cors x pweights2
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &total, &alpha, cors, &bitlength, pweights2, &bitlength, &beta, cors3, &bitlength);

//copy cors3 into svars allowing for missing values
for(j=0;j<bitlength;j++)
{
if(nss[bitstart+j]>0)
{
for(q=0;q<total;q++){svars[q][count]=cors3[j+q*bitlength];}
count++;
}
}
}	//end of bit loop
printf("\n");

fclose(input2);

//get ssums
for(q=0;q<total;q++)
{
//ssums[q][q2]indicates how much q2 contributes to q
for(q2=0;q2<total;q2++){ssums[q][q2]=0;}
//ssums[q][total] is the number of snps in category q
ssums[q][total]=0;

for(j=0;j<data_length;j++)
{
if(pindexes[q][j]==1)	//add on contribution from q2, then increase tally
{
for(q2=0;q2<total;q2++){ssums[q][q2]+=pweights[q2][j]*exps2[j];}
ssums[q][total]++;
}
}

//ssums[q][total+1] is the proportion of snps in category q
ssums[q][total+1]=ssums[q][total]/data_length;
}

//load up stats, allowing for missing values and ensuring not too big or small
count2=0;
for(j=0;j<data_length;j++)
{
if(nss[j]>0)
{
snss[count2]=nss[j];
schis[count2]=chis[j];
if(schis[count2]<1e-6){schis[count2]=1e-6;}
if(schis[count2]>100){schis[count2]=100;}
count2++;
}
}

if(count!=count2){printf("Doug Error 722BB %d and %d\n\n", count, count2);exit(1);}

if(total==1||parttype==1)	//solve directly
{solve_sums(stats, NULL, NULL, NULL, NULL, 1, 0, 0, -9999, count, 0, NULL, NULL, rjksums, svars, ssums, snss, schis, 0.001, 100, 1, 7, NULL);}
else
{
//solve first for base
for(j=0;j<data_length;j++){svars2[0][j]=svars[total-1][j];}
ssums2[0][0]=ssums[total-1][total-1];
ssums2[0][1]=ssums[total-1][total];
ssums2[0][2]=ssums[total-1][total+1];

solve_sums(stats2, NULL, NULL, NULL, NULL, 1, 0, 0, -9999, count, 0, NULL, NULL, rjksums, svars2, ssums2, snss, schis, 0.001, 100, 1, 7, NULL);

//load coefficients into stats and solve for full model
for(q=0;q<total;q++){stats[q]=0;}
stats[total-1]=stats2[0];
solve_sums(stats, NULL, NULL, NULL, NULL, total, 0, 0, -9999, count, 0, NULL, NULL, rjksums, svars, ssums, snss, schis, 0.001, 100, 1, 8, NULL);
}

//save enrichments
if(num_focals==1){sprintf(filename2,"%s.enrich",outfile);}
else{sprintf(filename2,"%s.focal%d.enrich", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Component Share Expected Enrichment\n");

for(q=0;q<total;q++)
{
if(parttype==0&&q<total-1){fprintf(output2,"Enrich_A%d ",q+1);}
if(parttype==0&&q==total-1){fprintf(output2,"Enrich_Base ");}
if(parttype==1){fprintf(output2,"Enrich_P%d ",q+1);}

fprintf(output2,"%.6f %.6f %.6f\n", stats[total+1+q]/stats[total], ssums[q][total+1], stats[total+1+q]/stats[total]/ssums[q][total+1]);
}
fclose(output2);

//set exps to zero, then add on contributions from each category
for(j=0;j<data_length;j++){exps[j]=0;}
for(q=0;q<total;q++)	//contribution of category q must sum to its heritability (stored in stats[q])
{
if(stats[q]!=0)
{
sum=0;for(j=0;j<data_length;j++){sum+=pweights[q][j]*exps2[j];}
for(j=0;j<data_length;j++){exps[j]+=pweights[q][j]*exps2[j]/sum*stats[q];}
}
}

//get heritability
her=0;for(j=0;j<data_length;j++){her+=exps[j];}
printf("Estimated heritability is %.4f\n", her);
if(her<0.01){printf("Warning, this is very low, so has been increased to 0.01\n");her=0.01;}
//if(her>maxher){printf("Warning, this is very high, so has been reduced to %.4f\n", maxher);her=maxher;}
printf("\n");

//force negative exps to be something small
for(j=0;j<data_length;j++)
{
if(exps[j]<0){exps[j]=1e-16;}
}

//make exps sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}

//ensure no tiny values
value=minher/data_length;
count=0;
for(j=0;j<data_length;j++)
{
if(exps[j]<value){exps[j]=value;count++;}
}
if(count>0){printf("Warning, %d of the per-predictor heritabilities were negative or very small (so have been set to %.4f times the average value; change this threshold using \"--min-her\")\n\n", count, minher);}

//ensure exps still sum to one
sum=0;for(j=0;j<data_length;j++){sum+=exps[j];}
for(j=0;j<data_length;j++){exps[j]=exps[j]/sum;}

free(exps2);free(pweights2);free(cors3);
free(snss);free(schis);free(stats);free(stats2);
for(q=0;q<total;q++){free(svars[q]);free(ssums[q]);}free(svars);free(ssums);
free(svars2[0]);free(svars2);free(ssums2[0]);free(ssums2);
}	//end of estimating heritability

//save her, power and inflation
if(num_focals==1){sprintf(filename2,"%s.arch",outfile);}
else{sprintf(filename2,"%s.focal%d.arch", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Heritability %.4f\nPower %.4f\nInflation %.4f\n", her, power, sscale);
fclose(output2);

//save per-predictor heritabilities
if(num_focals==1){sprintf(filename2,"%s.ind.hers",outfile);}
else{sprintf(filename2,"%s.focal%d.ind.hers", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor Heritability\n");
for(j=0;j<data_length;j++){fprintf(output2, "%s %.6e\n", preds[j], exps[j]*her);}
fclose(output2);

////////

if(impsums==1)  //impute missing summary statistic for focal trait (if necessary)
{
count=0;for(j=0;j<data_length;j++){count+=(nss[j]==0);}
if(count>0)
{
if(num_sums3==1){printf("Imputing summary statistics for (up to) %d predictors\n", count);}
else{printf("Imputing summary statistics for (up to) %d predictors for focal trait (%s)\n", count, sumstems[0]);}

if(metasum==0){impute_sums(0, Mnss, Mrhos, Mchis, data_length, bittotal, blockstarts, blockends, sumpops[0], NULL, corstems, -9999, Dindexes, Dsizes, Duse, Duse2, Dsigns, shrink, outfile, cors_short, cors, preds);}
else{impute_sums(0, Mnss, Mrhos, Mchis, data_length, bittotal, blockstarts, blockends, sumpops[0], cohers2, corstems, num_cors, Dindexes, Dsizes, Duse, Duse2, Dsigns, shrink, outfile, cors_short, cors, preds);}
}
}

if(num_focals==1){sprintf(filename2,"%s.imputed",outfile);}
else{sprintf(filename2,"%s.focal%d.imputed", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor\tA1\tA2\tZ\tn\n");

for(j=0;j<data_length;j++)
{
if(nss[j]>0)
{
if(rhos[j]>0){value=pow(Mchis[0][j],.5);}
else{value=-pow(Mchis[0][j],.5);}
fprintf(output2,"%s\t%c\t%c\t%.4f\t%.0f\n", preds[j], al1[j], al2[j], value, Mnss[0][j]);
}
}
fclose(output2);

///////////////////////

if(skipcv==0)   //deal with pseudos - only need nss2, rhos2 and rhos3 (not nss3 or chis2 / chis3)
{
nss2=malloc(sizeof(double)*data_length);
rhos2=malloc(sizeof(double)*data_length);
rhos3=malloc(sizeof(double)*data_length);

if(cvprop!=-9999)	//make them from datarands
{
//nss2 is (1-cvprop)*nss
for(j=0;j<data_length;j++){nss2[j]=(1-cvprop)*nss[j];}

//rhos2 is rhos plus datarands*root(cvprop/(1-cvprop)/nss)
value=pow(cvprop/(1-cvprop),.5);
for(j=0;j<data_length;j++){rhos2[j]=rhos[j]+datarands[j]*value*pow(nss[j],-.5);}

//rhos3 is complement (also equals rhos - datarands*root((1-cvprop)/cvprop/nss)
for(j=0;j<data_length;j++){rhos3[j]=(rhos[j]-(1-cvprop)*rhos2[j])/cvprop;}
}
else	//read from pseudo files
{
//will also read nss3, chis2 and chis3 (then delete)
nss3=malloc(sizeof(double)*data_length);
chis2=malloc(sizeof(double)*data_length);
chis3=malloc(sizeof(double)*data_length);

sprintf(filename2,"%s.train.summaries", pseudostem);
printf("Reading training summary statistics from %s\n", filename2);
read_sumsfile(filename2, nss2, chis2, rhos2, NULL, data_length, preds, al1, al2, bimfile, amb, 1.0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis2[0], nss2[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis2[j], nss2[j]);}}
printf("\n\n");

sprintf(filename3,"%s.test.summaries", pseudostem);
printf("Reading test summary statistics from %s\n", filename3);
read_sumsfile(filename3, nss3, chis3, rhos3, NULL, data_length, preds, al1, al2, bimfile, amb, 1.0, -9999);
printf("First few stats and ns are: %s %.3f %.1f", preds[0], chis3[0], nss3[0]);
for(j=1;j<data_length;j++){if(j<3){printf(" | %s %.3f %.1f", preds[j], chis3[j], nss3[j]);}}
printf("\n\n");

free(nss3);free(chis2);free(chis3);
}
}

if(checkld==1)	//using highlds (have already set based on cors file)
{
if(strcmp(ldfile,"blank")==0)	
{
count=0;for(j=0;j<data_length;j++){count+=highlds[j];}
printf("There are %d high-LD predictors (these will be excluded when performing pseudo cross-validation)\n\n", count);
if(count==data_length){printf("Error, there are no non-high-LD predictors; please tell Doug\n\n");exit(1);}
}
else	//reset based on file provided by user
{
count=countrows(ldfile);
printf("Reading list of %d high-LD predictors from %s\n", count, ldfile);
wantpreds=malloc(sizeof(char*)*count);
indexer=malloc(sizeof(int)*data_length);
read_strings(ldfile, wantpreds, count, NULL, 1, 0);

count2=find_strings(preds, data_length, wantpreds, count, indexer, NULL, NULL, NULL, NULL, NULL, 3);
if(count2==0){printf("Warning, none of the predictors are in the data\n\n");}
if(count2>0&&count2<count){printf("Warning, only %d of these are in the data\n\n", count2);}
if(count2==count){printf("All of these are in the data\n\n");}
if(count2==data_length){printf("Error, there are no non-high-LD predictors; please tell Doug\n\n");exit(1);}

for(j=0;j<data_length;j++){highlds[j]=0;}
for(j=0;j<count2;j++){highlds[indexer[j]]=1;}

for(j=0;j<count;j++){free(wantpreds[j]);}free(wantpreds);free(indexer);
}
}
else	//not using highlds (set highlds to 0, but probably unnecesary)
{
for(j=0;j<data_length;j++){highlds[j]=0;}
}

////////

if(num_sums3>1)  //deal with extra traits, and impute if necessary (listens to dougvar 1 and 3)
{
#include "secondaries.c"
}

if(num_sums3==1&&metasum==1)    //put cohers2 into cohers (if num_sums3>1 have already set cohers, while if num_sums3=1 & metasum=0 do not use cohers) 
{
for(bit=0;bit<bittotal;bit++)
{
for(s=0;s<num_cors;s++){cohers[s+bit*num_cors]=cohers2[s];}
}
}

///////////////////////

//allocate variables used to estimate effect sizes

if(mcmcfinal==0){total=num_try;}
else{total=num_try*num_chains;}

if(finemap==0){total2=total;}
else{total2=Lnum;}

value=(double)data_length/1024*(total+num_try)/1024/1024*8;
if(value>1){printf("Warning, to store the effect sizes requires %.1f Gb\n\n", value);}

if(prsvar==1)
{
value=(double)data_length/1024*num_blocks/1024/1024*8;
if(value>1){printf("Warning, to store the MCMC samplings requires %.1f Gb\n\n", value);}
}

YTdata=malloc(sizeof(double)*data_length);
cors2=malloc(sizeof(double)*bitmax*total);
effs=malloc(sizeof(double)*data_length*total);
effs2=malloc(sizeof(double)*bitmax*total);
effs3=malloc(sizeof(double)*bitmax*total2);
probs=malloc(sizeof(double)*data_length*num_try);
probs2=malloc(sizeof(double)*bitmax*total2);
probs3=malloc(sizeof(double)*bitmax*2);
if(finemap==0){residuals=malloc(sizeof(double)*bitsize*total);}
else{residuals=malloc(sizeof(double)*bitmax*total);}
ess=malloc(sizeof(double)*total);
ess2=malloc(sizeof(double)*total);
variances=malloc(sizeof(double)*total);
pars=malloc(sizeof(double)*num_try);

pmeans=malloc(sizeof(double)*bitmax*num_grids);
pvars=malloc(sizeof(double)*bitmax*num_grids);
bfs=malloc(sizeof(double)*bitmax*num_grids);
gridmaxes=malloc(sizeof(double)*num_grids);
gridsums=malloc(sizeof(double)*num_grids);

predvars=malloc(sizeof(double)*total*total);
predtops=malloc(sizeof(double)*total);
predcors=malloc(sizeof(double)*total);
predweights=malloc(sizeof(double)*num_try);

dptrs=malloc(sizeof(struct sorting_double)*bitmax);

if(finemap==1){cors4=malloc(sizeof(double)*bitmax*total2);}

if(prsvar==1){effs4=malloc(sizeof(double)*data_length*num_blocks);}

if(mcmcfinal==1){manyrands=malloc(sizeof(double)*128*total);}

///////////////////////

//STARTOFLOOP1

if(skipcv==0)	//solve using training summaries, then test using test summaries
{
//screen and file print
printf("Estimating training effect sizes for %d models\n", num_try);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating training effect sizes for %d models\n", num_try);
fclose(output);

//get XTY for each predictor
for(j=0;j<data_length;j++){YTdata[j]=rhos2[j]*nss2[j];}

//set starting effects based on type of model
#pragma omp parallel for private(p,j,value,value2) schedule(static)
for(p=0;p<num_try;p++)
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0&&highlds[j]==0)
{
//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*her;
value2=pow(exps[j]*her,-.5);

if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
if(trytypes[p]==2)	//lasso - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, -9999, -9999, -9999, nss2[j], 1.0, -9999, -9999, -9999, -9999, NULL, 2, NULL, NULL)/rjksums[j];}
if(trytypes[p]==3)	//ridge - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, -9999, -9999, -9999, nss2[j], 1.0, -9999, -9999, -9999, -9999, NULL, 3, NULL, NULL)/rjksums[j];}
if(trytypes[p]==4)	//bolt - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss2[j], 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL, NULL)/rjksums[j];}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, nss2[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL, NULL)/rjksums[j];}
if(trytypes[p]==7)	//elastic - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, nss2[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, NULL, NULL)/rjksums[j];}
if(trytypes[p]==8)	//ldpred - posterior mean, weighted by tagging (use bolt function)
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss2[j], 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL, NULL)/rjksums[j];}
}
else	//set to zero
{effs[(size_t)p*data_length+j]=0;}
}
}

//set predvars to zero
for(p2=0;p2<num_try;p2++)
{
for(p=0;p<num_try;p++){predvars[p+p2*num_try]=0;}
}

ecount=0;
wcount=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(bit%100==0)
{
printf("Analyzing Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Analyzing Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//read correlations
if(metasum==0||dougvar==6) //easy case - only using focal correlations
{
s=sumpops[0];

//re-open cors
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read lower triangle of focal correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//now extract correlations for predictors we are using (remember shrink, but no need for signs)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<bitlength;j++)
{cors[(size_t)k*bitlength+j]=(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*shrink;}
}

fclose(input2);
}
else    //hard case - might use multiple correlations
{
//set cors to identity
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1;
for(j=k+1;j<bitlength;j++){cors[(size_t)k*bitlength+j]=0;}
}

for(s=0;s<num_cors;s++)
{
if(cohers2[s]>0&&Dsizes[s][bit][1]>0)	//this population contributes
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

fclose(input2);

//now add on off-diagonal correlations for predictors we are using (remember cohers2, signs and shrink)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<Dsizes[s][bit][1];k++)
{
for(j=k+1;j<Dsizes[s][bit][1];j++)
{cors[(size_t)Duse2[s][bit][k]*bitlength+Duse2[s][bit][j]]+=cohers2[s]*(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*Dsigns[s][bit][j]*Dsigns[s][bit][k]*shrink;}
}
}}
}

//get average sample size
sum=0;for(j=bitstart;j<bitend;j++){sum+=nss[j];}
neff=sum/bitlength;

if(resfix==0)   //compute rhosinvrhos
{rhosinvrhos=cg_solve_lower_norm(cors, bitlength, rhos2+bitstart, 1e-6);}

//put starting effects into effs2
for(p=0;p<num_try;p++)
{
for(j=0;j<bitlength;j++){effs2[j+p*bitlength]=effs[(size_t)p*data_length+bitstart+j];}
}

//set variances to one
for(p=0;p<num_try;p++){variances[p]=1.0;}

//calculate ess for all models - equal to (2YTXbeta - t(beta) XTXbeta)/n
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_try, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

for(p=0;p<num_try;p++){ess[p]=2*ddot_(&bitlength, effs2+p*bitlength, &one, rhos2+bitstart, &one)-ddot_(&bitlength, effs2+p*bitlength, &one, cors2+p*bitlength, &one);}

for(count=0;count<maxiter;count++)
{
//set pars to zero
for(p=0;p<num_try;p++){pars[p]=0;}

for(bitstart2=bitstart;bitstart2<bitend;bitstart2+=bitsize)
{
bitend2=bitstart2+bitsize;
if(bitend2>bitend){bitend2=bitend;}
bitlength2=bitend2-bitstart2;

//get XjXTbeta for all snps in (small) block - this code would work if matrix not triangular
//alpha=1.0;beta=0.0;
//dgemm_("T", "N", &bitlength2, &num_try, &bitlength, &alpha, cors+(size_t)(bitstart2-bitstart)*bitlength, &bitlength, effs2, &bitlength, &beta, residuals, &bitsize);

//get XjXTbeta for all snps in (small) block
token=bitstart2-bitstart;
token2=bitend-bitend2;

//middle square - from bitstart2 to bitend2
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength2, &num_try, &alpha, cors+(size_t)token*bitlength+token, &bitlength, effs2+token, &bitlength, &beta, residuals, &bitsize);

alpha=1.0;beta=1.0;
if(token>0)	//left rectangle - from zero to bitstart2
{dgemm_("N", "N", &bitlength2, &num_try, &token, &alpha, cors+token, &bitlength, effs2, &bitlength, &beta, residuals, &bitsize);}
if(token2>0)	//right rectangle - from bitend2 to bitend
{dgemm_("T", "N", &bitlength2, &num_try, &token2, &alpha, cors+(size_t)token*bitlength+(bitend2-bitstart), &bitlength, effs2+bitend2-bitstart, &bitlength, &beta, residuals, &bitsize);}

#pragma omp parallel for private(p,j,j2,sum,value,value2,postmean,value3) schedule(static)
for(p=0;p<num_try;p++)
{
for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0&&highlds[j]==0)
{
//get XjT residuals
sum=YTdata[j]-tryscales[p]*nss2[j]*(residuals[j-bitstart2+p*bitsize]-effs2[j-bitstart+p*bitlength]);

//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*her;
value2=pow(exps[j]*her,-.5);

//get posterior mean
if(trytypes[p]==1)	//lasso-sparse
{postmean=get_postmean(sum, lambdas[p]*value2*nss2[j], -9999, -9999, -9999, nss2[j], variances[p], -9999, -9999, -9999, -9999, NULL, 1, NULL, pars+p);}
if(trytypes[p]==2)	//lasso
{postmean=get_postmean(sum, lambdas[p]*value2, -9999, -9999, -9999, nss2[j], variances[p], -9999, -9999, -9999, -9999, NULL, 2, NULL, pars+p);}
if(trytypes[p]==3)	//ridge
{postmean=get_postmean(sum, lambdas[p]*value, -9999, -9999, -9999, nss2[j], variances[p], -9999, -9999, -9999, -9999, NULL, 3, NULL, pars+p);}
if(trytypes[p]==4)	//bolt
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss2[j], variances[p], tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL, pars+p);}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, nss2[j], variances[p], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL, pars+p);}
if(trytypes[p]==7)	//elastic
{postmean=get_postmean(sum, lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, nss2[j], variances[p], tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, NULL, pars+p);}
if(trytypes[p]==8)	//ldpred (use bolt function)
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss2[j], variances[p], tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL, pars+p);}

//update residuals for remainder of block
value3=postmean-effs2[j-bitstart+p*bitlength];
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

effs2[j-bitstart+p*bitlength]=postmean;
}}	//end of j loop
}	//end of p loop
}	//end of bitstart2 loop

//save old ess
for(p=0;p<num_try;p++){ess2[p]=ess[p];}

//get new ess
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_try, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

for(p=0;p<num_try;p++){ess[p]=2*ddot_(&bitlength, effs2+p*bitlength, &one, rhos2+bitstart, &one)-ddot_(&bitlength, effs2+p*bitlength, &one, cors2+p*bitlength, &one);}

//see whether converged
cflag=0;for(p=0;p<num_try;p++){cflag+=(fabs(ess[p]-ess2[p])<tol);}
if(cflag==num_try){break;}

if(resfix==0)   //update variance based on expected distribution of rho, assuming prior dist prop to variance^-1.5 (do not allow to go below one)
{
for(p=0;p<num_try;p++) 
{
//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta + variances (terms 2 and 3 equal -ess, computed just above)
sum=rhosinvrhos-ess[p]+pars[p];
variances[p]=sum/bitlength*neff;
if(variances[p]<1.0){variances[p]=1.0;}
}
}
}	//end of count loop

if(count==maxiter){printf("Warning, Window %d did not converge within %d iterations\n", bit+1, maxiter);}

cflag=0;for(p=0;p<num_try;p++){cflag+=(ess[p]>1);}
if(cflag==0)	//models not terrible, so update effect sizes
{
for(p=0;p<num_try;p++)
{
for(j=bitstart;j<bitend;j++){effs[(size_t)p*data_length+j]=effs2[j-bitstart+p*bitlength];}
}
}
else	//one or more models is suspect - leave effect sizes at starting estimates
{
printf("Warning, Window %d went wild, will revert to starting estimates\n", bit+1);
ecount++;
wcount+=bitlength;
}

////////

//recompute cors2, first setting small correlations to zero (unsure whether I should reverse shrinkage)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=k+1;j<bitlength;j++)
{
if(fabs(cors[(size_t)k*bitlength+j])<0.01){cors[(size_t)k*bitlength+j]=0;}
}
}
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_try, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

//add on contribution to predvars
alpha=1.0;beta=1.0;
dgemm_("T","N", &num_try, &num_try, &bitlength, &alpha, effs+bitstart, &data_length, cors2, &bitlength, &beta, predvars, &num_try);
}	//end of bit loop

if(ecount>0){printf("%d windows went wild (in total, %d predictors)\n", ecount, wcount);}
printf("\n");

////////

//test the models, and record the best (note that if exps[j]=0 or highld[j]=1, the predictor will have effect zero)
printf("Testing the models using pseudo test summary statistics\n");

//predtops=rhos3 beta
for(p=0;p<num_try;p++){predtops[p]=ddot_(&data_length, effs+(size_t)p*data_length, &one, rhos3, &one);}

//predcors = predtops / root (predvars)
for(p=0;p<num_try;p++)
{
if(predvars[p+p*num_try]>0){predcors[p]=predtops[p]*pow(predvars[p+p*num_try],-.5);}
else{predcors[p]=-9999;}
}

//find best
best=-1;
for(p=0;p<num_try;p++)
{
if(predcors[p]!=-9999)
{
if(best==-1){best=p;value=predcors[p];}
if(predcors[p]>value){best=p;value=predcors[p];}
}
}

if(best==-1)	//not possible to compute a correlation for any models
{printf("Error, it was not possible to compute a correlation for any of the models (suggesting they all have zero effect sizes)\n\n");exit(1);}

//will print out scaled effect sizes
escale=pow(predvars[best+best*num_try],-.5);
if(dougvar==5){escale=1;}

//save accuracies
if(num_focals==1){sprintf(filename2,"%s.cors",outfile);}
else{sprintf(filename2,"%s.focal%d.cors", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2,"Model\tCorrelation\n");

for(p=0;p<num_try;p++)
{
if(predcors[p]!=-9999){fprintf(output2,"%d\t%.6f\n", p+1, predcors[p]);}
else{fprintf(output2,"%d\tNA\n", p+1);}
}
fclose(output2);

//save the best model
if(num_focals==1){sprintf(filename3,"%s.best",outfile);}
else{sprintf(filename3,"%s.focal%d.best", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3, "Model\tType\tlambda\tscale\tp\tf2\tp1\tp2\tp3\tp4\n");

if(trytypes[best]==1){fprintf(output3, "%d\tlasso-sparse\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1, trylams[best], tryscales[best]);}
if(trytypes[best]==2){fprintf(output3, "%d\tlasso\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1);}
if(trytypes[best]==3){fprintf(output3, "%d\tridge\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", best+1);}
if(trytypes[best]==4){fprintf(output3, "%d\tbolt\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\n", best+1, tryps[best], tryf2s[best]);}
if(trytypes[best]==5){fprintf(output3, "%d\tbayesr\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", best+1, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(trytypes[best]==6){fprintf(output3, "%d\tbayesr-shrink\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\n", best+1, tryps[best], tryp2s[best], tryp3s[best], tryp4s[best]);}
if(trytypes[best]==7){fprintf(output3, "%d\telastic\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\n", best+1, 1-tryp3s[best], tryf2s[best]);}
if(trytypes[best]==8){fprintf(output3, "%d\tLDpred\tNA\tNA\t%.e\tNA\tNA\tNA\tNA\tNA\n", best+1, tryps[best]);}

fclose(output3);

printf("The estimated accuracies are saved in %s, while the parameters corresponding to the best model (which has approximate correlation %.2f) are saved in %s\n\n", filename2, predcors[best], filename3);

if(ensemble==0)	//only best model gets weight
{
for(p=0;p<num_try;p++){predweights[p]=0;}
predweights[best]=1;
}
else	//estimate weights using nnls
{
//set coeffs and var explained to zero
for(p=0;p<num_try;p++){predweights[p]=0;}
sumsq=0;

//will loop from best to worst, ignoring dodgy
dptrs=malloc(sizeof(struct sorting_double)*num_try);
count2=0;
for(p=0;p<num_try;p++)
{
if(predvars[p+p*num_try]>0){dptrs[count2].value=predcors[p];dptrs[count2].index=p;count2++;}
}
qsort(dptrs, count2, sizeof(struct sorting_double), compare_sorting_double_rev);

for(count=0;count<100;count++)
{
for(q=0;q<count2;q++)
{
p=dptrs[q].index;

//get t(X_p)(Y-Xtheta)
sum=predtops[p]+predvars[p+p*num_try]*predweights[p];
for(q2=0;q2<count2;q2++)
{
p2=dptrs[q2].index;
sum-=predvars[p+p2*num_try]*predweights[p2];
}

//update coeff p, ensuring not negative
value=sum/predvars[p+p*num_try];
if(value>0){predweights[p]=value;}
else{predweights[p]=0;}
}

//see whether var explained has converged
sumsq2=sumsq;
sumsq=0;
for(p=0;p<num_try;p++)
{
sumsq+=2*predtops[p]*predweights[p];
for(p2=0;p2<num_try;p2++){sumsq-=predweights[p]*predweights[p2]*predvars[p+p2*num_try];}
}

if(fabs(sumsq2-sumsq)<0.000001){break;}
}

free(dptrs);

//scale weights so sum to one
sum=0;for(p=0;p<num_try;p++){sum+=predweights[p];}
value=pow(sum,-1);
for(p=0;p<num_try;p++){predweights[p]*=value;
if(predweights[p]>0){printf("Keep %d weight %f cor %f\n", p+1, predweights[p], predcors[p]);}
}
}
}	//end of using training and test summary statistics

///////////////////////

//STARTOFLOOP2 - solve using main statistic

if(num_sums3>1)  //replace focal summary statistics with mtag version
{
for(j=0;j<data_length;j++){nss[j]=nss4[j];chis[j]=chis4[j];rhos[j]=rhos4[j];}
}

//pflag indicates whether we might be using non-focal ancestries
if(metasum==0)
{
pflag=0;
for(q=1;q<num_sums3;q++)
{
if(pritraits[q]==1&&sumpops[q]!=sumpops[0]){pflag=1;break;}
}
}
else    //assume we are (might be wrong)
{pflag=1;}

//work out how many models we are constructing
if(finemap==0) //can be many
{
if(skipcv==0&&megasave==0)	//will reduce to models with positive weight
{
num_left=0;
for(p=0;p<num_try;p++)
{
if(predweights[p]>0)
{
if(num_left<p)
{
trytypes[num_left]=trytypes[p];
trylams[num_left]=trylams[p];tryscales[num_left]=tryscales[p];
tryps[num_left]=tryps[p];tryp2s[num_left]=tryp2s[p];
tryp3s[num_left]=tryp3s[p];tryp4s[num_left]=tryp4s[p];
tryf2s[num_left]=tryf2s[p];
lambdas[num_left]=lambdas[p];lambdas2[num_left]=lambdas2[p];
lambdas3[num_left]=lambdas3[p];lambdas4[num_left]=lambdas4[p];
predweights[num_left]=predweights[p];
}
num_left++;
}
}

printf("Estimating final effect sizes for the best model\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating final effect sizes for the best model\n");
fclose(output);
}
else	//will use all models
{
num_left=num_try;

printf("Estimating final effect sizes for %d models\n", num_left);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating final effect sizes for %d models\n", num_left);
fclose(output);
}
}
else    //finemapping - only have one model
{
num_left=1;

printf("Computing posterior inclusion probabilities\n");

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Estimating effect sizes and posterior inclusion probabilities\n");
fclose(output);
}

//get XTY for each predictor
for(j=0;j<data_length;j++){YTdata[j]=rhos[j]*nss[j];}

if(finemap==0)  //set starting effects (perhaps better to adjust rjksums for populations, but should not make a big difference)
{
#pragma omp parallel for private(p,j,value,value2) schedule(static)
for(p=0;p<num_left;p++)
{
for(j=0;j<data_length;j++)
{
if(exps[j]>0)
{
//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*her;
value2=pow(exps[j]*her,-.5);

if(trytypes[p]==1)	//lasso-sparse - set to zero
{effs[(size_t)p*data_length+j]=0;}
if(trytypes[p]==2)	//lasso - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, -9999, -9999, -9999, nss[j], 1.0, -9999, -9999, -9999, -9999, NULL, 2, NULL, NULL)/rjksums[j];}
if(trytypes[p]==3)	//ridge - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, -9999, -9999, -9999, nss[j], 1.0, -9999, -9999, -9999, -9999, NULL, 3, NULL, NULL)/rjksums[j];}
if(trytypes[p]==4)	//bolt - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss[j], 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL, NULL)/rjksums[j];}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, nss[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], NULL, NULL)/rjksums[j];}
if(trytypes[p]==7)	//elastic - posterior mean, weighted by tagging
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, nss[j], 1.0, tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, NULL, NULL)/rjksums[j];}
if(trytypes[p]==8)	//ldpred - posterior mean, weighted by tagging (use bolt function)
{effs[(size_t)p*data_length+j]=get_postmean(YTdata[j], lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss[j], 1.0, tryps[p], tryp2s[p], -9999, -9999, NULL, 4, NULL, NULL)/rjksums[j];}
}
else	//set to zero
{effs[(size_t)p*data_length+j]=0;}
}
}
}
else    //set starting effects to zero (num_left=1)
{
for(j=0;j<data_length;j++){effs[j]=0;}
}

//set probabilities to zero (probably unnecessary)
for(p=0;p<num_left;p++)
{
for(j=0;j<data_length;j++){probs[(size_t)p*data_length+j]=0;}
}

if(condition==1)    //open clumping file
{
if(num_focals==1){sprintf(filename3,"%s.conditional",outfile);}
else{sprintf(filename3,"%s.focal%d.conditional", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Chromosome Basepair Conditional_P Original_P\n");
}

if(finemap==1)  //open credible set files
{
if(num_focals==1){sprintf(filename4,"%s.credibles",outfile);}
else{sprintf(filename4,"%s.focal%d.credibles", outfile, cur_focal+1);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}
fprintf(output4,"Set PIP Size Min_P Purity Predictors\n");
}

ecount=0;
wcount=0;
count4=0;   //counts number of credible sets found
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if((mcmcfinal==0&&bit%100==0)||(mcmcfinal==1&&bit%20==0))
{
printf("Analyzing Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Analyzing Window %d of %d\n", bit+1, bittotal);
fclose(output);

if(condition==1)
{
fclose(output3);
if((output3=fopen(filename3,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename3);exit(1);}
}

if(finemap==1)
{
fclose(output4);
if((output4=fopen(filename4,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename4);exit(1);}
}
}

//read correlations
if(pflag==0) //easy case - only using focal correlations
{
s=sumpops[0];

//re-open cors
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

//read lower triangle of focal correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

 //now extract correlations for predictors we are using (remember shrink, but no need for signs)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<Dsizes[s][bit][1];j++)
{cors[(size_t)k*bitlength+j]=(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*shrink;}
}

fclose(input2);
}
else    //hard case - might use multiple correlations
{
//set cors to identity
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1;
for(j=k+1;j<bitlength;j++){cors[(size_t)k*bitlength+j]=0;}
}

for(s=0;s<num_cors;s++)
{
if(cohers[s+bit*num_cors]>0&&Dsizes[s][bit][1]>0)	//this population contributes
{
sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

fclose(input2);

//now add on off-diagonal correlations for predictors we are using (remember cohers, signs and shrink)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<Dsizes[s][bit][1];k++)
{
for(j=k+1;j<Dsizes[s][bit][1];j++)
{cors[(size_t)Duse2[s][bit][k]*bitlength+Duse2[s][bit][j]]+=cohers[s+bit*num_cors]*(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*Dsigns[s][bit][j]*Dsigns[s][bit][k]*shrink;}
}
}}
}

//get average sample size
sum=0;for(j=bitstart;j<bitend;j++){sum+=nss[j];}
neff=sum/bitlength;

//get sum of effs
sumexps=0;for(j=bitstart;j<bitend;j++){sumexps+=exps[j];}

if(resfix==0)   //compute rhosinvrhos
{rhosinvrhos=cg_solve_lower_norm(cors, bitlength, rhos+bitstart, 1e-6);}

if(condition==1)    //perform conditional analysis
{
#include "clump.c"
}

//flag=0 means that we will do vb
flag=mcmcfinal;

if(mcmcfinal==1)	//first try fixed number of iterations, then save means (can not have trytype=7
{
total=num_left*num_chains;

//use effs2 and probs2 to store sums of effect sizes and probs for each model / chain
for(p=0;p<total;p++)
{
for(j=0;j<bitlength;j++){effs2[j+p*bitlength]=0;probs2[j+p*bitlength]=0;}
}

//put starting effects into effs3
for(p2=0;p2<total;p2++)
{
p=p2/num_chains;
for(j=0;j<bitlength;j++){effs3[j+p2*bitlength]=effs[(size_t)p*data_length+bitstart+j];}
}

//set variances to one
for(p=0;p<total;p++){variances[p]=1.0;}

count2=0;
for(count=0;count<num_burn+num_mcmc;count++)
{
for(bitstart2=bitstart;bitstart2<bitend;bitstart2+=bitsize)
{
bitend2=bitstart2+bitsize;
if(bitend2>bitend){bitend2=bitend;}
bitlength2=bitend2-bitstart2;

//get XjXTbeta for all snps in (small) block
token=bitstart2-bitstart;
token2=bitend-bitend2;

//middle square - from bitstart2 to bitend2
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength2, &total, &alpha, cors+(size_t)token*bitlength+token, &bitlength, effs3+token, &bitlength, &beta, residuals, &bitsize);

alpha=1.0;beta=1.0;
if(token>0)	//left rectangle - from zero to bitstart2
{dgemm_("N", "N", &bitlength2, &total, &token, &alpha, cors+token, &bitlength, effs3, &bitlength, &beta, residuals, &bitsize);}
if(token2>0)	//right rectangle - from bitend2 to bitend
{dgemm_("T", "N", &bitlength2, &total, &token2, &alpha, cors+(size_t)token*bitlength+(bitend2-bitstart), &bitlength, effs3+bitend2-bitstart, &bitlength, &beta, residuals, &bitsize);}

//sample new unifrands (might not use)
for(j=0;j<128*total;j++){manyrands[j]=genrand_real1();}

//first do gaussian-only models - can use pragma
#pragma omp parallel for private(p2,p,start,j,j2,sum,value,postmean,postprob,postsamp,value3) schedule(dynamic)
for(p2=0;p2<total;p2++)
{
p=p2/num_chains;
if(trytypes[p]==3||trytypes[p]==4||trytypes[p]==5||trytypes[p]==6||trytypes[p]==8)
{
start=0;
for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-tryscales[p]*nss[j]*(residuals[j-bitstart2+p2*bitsize]-effs3[j-bitstart+p2*bitlength]);

//value is how much to scale gaussian parameters
value=exps[j]*her;

//get posterior mean and sampling
if(trytypes[p]==3)	//ridge
{postmean=get_postsamp_pragma(sum, lambdas[p]*value, -9999, -9999, -9999, nss[j], variances[p2], -9999, -9999, -9999, -9999, 3, &postprob, &postsamp, manyrands+p2*128, &start, NULL);}
if(trytypes[p]==4)	//bolt
{postmean=get_postsamp_pragma(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss[j], variances[p2], tryps[p], tryp2s[p], -9999, -9999, 4, &postprob, &postsamp, manyrands+p2*128, &start, NULL);}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink
{postmean=get_postsamp_pragma(sum, lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, nss[j], variances[p2], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], trytypes[p], &postprob, &postsamp, manyrands+p2*128, &start, NULL);}
if(trytypes[p]==8)	//ldpred (use bolt function)
{postmean=get_postsamp_pragma(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss[j], variances[p2], tryps[p], tryp2s[p], -9999, -9999, 4, &postprob, &postsamp, manyrands+p2*128, &start, NULL);}

//update residuals for remainder of block - only necessary if estimate has changed
value3=postsamp-effs3[j-bitstart+p2*bitlength];
if(value3!=0)
{
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p2*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}
}

if(count>=num_burn){effs2[j-bitstart+p2*bitlength]+=postmean;probs2[j-bitstart+p2*bitlength]+=postprob;}
effs3[j-bitstart+p2*bitlength]=postsamp;
}}	//end of j loop
}	//end of gaussian models
}	//end of p2 loop

//now do lasso models - must do in serial
for(p2=0;p2<total;p2++)
{
p=p2/num_chains;
if(trytypes[p]==1||trytypes[p]==2)
{
for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-tryscales[p]*nss[j]*(residuals[j-bitstart2+p2*bitsize]-effs3[j-bitstart+p2*bitlength]);

//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*her;
value2=pow(exps[j]*her,-.5);

//get posterior mean and sampling
if(trytypes[p]==1)	//lasso-sparse
{postmean=get_postsamp(sum, lambdas[p]*value2*nss[j], -9999, -9999, -9999, nss[j], variances[p2], -9999, -9999, -9999, -9999, NULL, 1, &postprob, &postsamp, NULL);}
if(trytypes[p]==2)	//lasso
{postmean=get_postsamp(sum, lambdas[p]*value2, -9999, -9999, -9999, nss[j], variances[p2], -9999, -9999, -9999, -9999, NULL, 2, &postprob, &postsamp, NULL);}
if(trytypes[p]==7)	//elastic
{postmean=get_postsamp(sum, lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, nss[j], variances[p2], tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, &postprob, &postsamp, NULL);}

//update residuals for remainder of block (effects will always change)
value3=postsamp-effs3[j-bitstart+p2*bitlength];
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p2*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

if(count>=num_burn){effs2[j-bitstart+p2*bitlength]+=postmean;probs2[j-bitstart+p2*bitlength]+=postprob;}
effs3[j-bitstart+p2*bitlength]=postsamp;
}}	//end of j loop
}	//end of lasso and elastic models
}	//end of p2 loop
}	//end of bitstart2 loop

if(prsvar==1)	//see if saving this iteration - will have num_left=1
{
if(count>=num_burn&&(count-num_burn)%num_step==num_step-1)
{
for(p=0;p<num_chains;p++)
{
for(j=bitstart;j<bitend;j++){effs4[(size_t)count2*data_length+j]=effs3[j-bitstart+p*bitlength];}
count2++;
}
}
}

if(resfix==0)   //sample variance based on conditional distribution of rho, assuming prior dist prop to variance^-.5
{
//get Cbeta
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &total, &alpha, cors, &bitlength, effs3, &bitlength, &beta, cors2, &bitlength);

for(p=0;p<total;p++)
{
if(dougvar!=7)  //use sbayesrc criterion to decide how to move
{
//compute t(beta) beta and t(beta) C beta - two (crude) estimates of variance explained for window
value=ddot_(&bitlength, effs3+p*bitlength, &one, effs3+p*bitlength, &one);
value2=ddot_(&bitlength, effs3+p*bitlength, &one, cors2+p*bitlength, &one);

if((value>1e-8&&value2>1e-8&&value>1.1*value2))    //resample variance from inverse gamma (but do not allow to go below 0.7)
{
//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta
sum=rhosinvrhos-2*ddot_(&bitlength, effs3+p*bitlength, &one, rhos+bitstart, &one)+value2;
variances[p]=0.5*sum*neff/rgamma(.5*bitlength);
if(variances[p]<0.7){variances[p]=0.7;}
}
else    //reset
{variances[p]=1.0;}
}
else    //always move (but do not allow to go below one)
{
//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta
sum=rhosinvrhos-2*ddot_(&bitlength, effs3+p*bitlength, &one, rhos+bitstart, &one)+ddot_(&bitlength, effs3+p*bitlength, &one, cors2+p*bitlength, &one);
variances[p]=0.5*sum*neff/rgamma(.5*bitlength);
if(variances[p]<1.0){variances[p]=1.0;}
}
}
}
}	//end of count loop

//get mean effects and probs for each model - store in first num_left columns of effs2 and probs2 (must be careful if use pragma)
for(j=0;j<bitlength;j++)
{
for(p=0;p<num_left;p++)
{
sum=0;sum2=0;
for(p2=0;p2<num_chains;p2++){sum+=effs2[j+(p*num_chains+p2)*bitlength];sum2+=probs2[j+(p*num_chains+p2)*bitlength];}
effs2[j+p*bitlength]=sum/num_mcmc/num_chains;
probs2[j+p*bitlength]=sum2/num_mcmc/num_chains;
}
}

//calculate ess for mean models - equal to 2YTXbeta - t(beta) XTXbeta (for multiple traits, this is ess of first trait only)
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_left, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

for(p=0;p<num_left;p++){ess[p]=2*ddot_(&bitlength, effs2+p*bitlength, &one, rhos+bitstart, &one)-ddot_(&bitlength, effs2+p*bitlength, &one, cors2+p*bitlength, &one);}

cflag=0;for(p=0;p<num_left;p++){cflag+=(ess[p]>1);}
if(cflag==0)	//model(s) not terrible, so update effect sizes and probs
{
for(p=0;p<num_left;p++)
{
for(j=bitstart;j<bitend;j++)
{effs[(size_t)p*data_length+j]=effs2[j-bitstart+p*bitlength];probs[(size_t)p*data_length+j]=probs2[j-bitstart+p*bitlength];}
}
}
else	//one or more models is suspect - leave effect sizes at starting estimates (probs will be zero), and switch to VB
{
printf("Warning, MCMC failed for Window %d, will revert to variational Bayes\n", bit+1);
flag=0;
ecount++;
wcount+=bitlength;
}
}

////////

if(flag==0)	//either did not try mcmc, or it failed, so use vb
{
if(finemap==0)  //put starting effects into effs2
{
for(p=0;p<num_left;p++)
{
for(j=0;j<bitlength;j++){effs2[j+p*bitlength]=effs[(size_t)p*data_length+bitstart+j];}
}
}
else    //set eff2 and effs3 to zero
{
for(p=0;p<Lnum;p++)
{
for(j=0;j<bitlength;j++){effs3[j+p*bitlength]=0;}
}

for(j=0;j<bitlength;j++){effs2[j]=0;}
}

//set variances to one
for(p=0;p<num_left;p++){variances[p]=1.0;}

//calculate ess for all models - equal to (2YTXbeta - t(beta) XTXbeta)/n
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_left, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

for(p=0;p<num_left;p++){ess[p]=2*ddot_(&bitlength, effs2+p*bitlength, &one, rhos+bitstart, &one)-ddot_(&bitlength, effs2+p*bitlength, &one, cors2+p*bitlength, &one);}

for(count=0;count<maxiter;count++)
{
//set pars to zero
for(p=0;p<num_left;p++){pars[p]=0;}

if(finemap==0)  //can sub-divide bit
{
//set probs2 to zero
for(p=0;p<num_left;p++)
{
for(j=0;j<bitlength;j++){probs2[j+p*bitlength]=0;}
}

for(bitstart2=bitstart;bitstart2<bitend;bitstart2+=bitsize)
{
bitend2=bitstart2+bitsize;
if(bitend2>bitend){bitend2=bitend;}
bitlength2=bitend2-bitstart2;

//get XjXTbeta for all snps in (small) block
token=bitstart2-bitstart;
token2=bitend-bitend2;

//middle square - from bitstart2 to bitend2
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength2, &num_left, &alpha, cors+(size_t)token*bitlength+token, &bitlength, effs2+token, &bitlength, &beta, residuals, &bitsize);

alpha=1.0;beta=1.0;
if(token>0)	//left rectangle - from zero to bitstart2
{dgemm_("N", "N", &bitlength2, &num_left, &token, &alpha, cors+token, &bitlength, effs2, &bitlength, &beta, residuals, &bitsize);}
if(token2>0)	//right rectangle - from bitend2 to bitend
{dgemm_("T", "N", &bitlength2, &num_left, &token2, &alpha, cors+(size_t)token*bitlength+(bitend2-bitstart), &bitlength, effs2+bitend2-bitstart, &bitlength, &beta, residuals, &bitsize);}

#pragma omp parallel for private(p,j,j2,sum,value,value2,postmean,postprob,value3) schedule(static)
for(p=0;p<num_left;p++)
{
for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-tryscales[p]*nss[j]*(residuals[j-bitstart2+p*bitsize]-effs2[j-bitstart+p*bitlength]);

//value and value2 are how much to scale gaussian and laplace parameters
value=exps[j]*her;
value2=pow(exps[j]*her,-.5);

//get posterior mean
if(trytypes[p]==1)	//lasso-sparse
{postmean=get_postmean(sum, lambdas[p]*value2*nss[j], -9999, -9999, -9999, nss[j], variances[p], -9999, -9999, -9999, -9999, NULL, 1, &postprob, pars+p);}
if(trytypes[p]==2)	//lasso
{postmean=get_postmean(sum, lambdas[p]*value2, -9999, -9999, -9999, nss[j], variances[p], -9999, -9999, -9999, -9999, NULL, 2, &postprob, pars+p);}
if(trytypes[p]==3)	//ridge
{postmean=get_postmean(sum, lambdas[p]*value, -9999, -9999, -9999, nss[j], variances[p], -9999, -9999, -9999, -9999, NULL, 3, &postprob, pars+p);}
if(trytypes[p]==4)	//bolt
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss[j], variances[p], tryps[p], tryp2s[p], -9999, -9999, NULL, 4, &postprob, pars+p);}
if(trytypes[p]==5||trytypes[p]==6)	//bayesr or bayesr-shrink
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, lambdas3[p]*value, lambdas4[p]*value, nss[j], variances[p], tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], NULL, trytypes[p], &postprob, pars+p);}
if(trytypes[p]==7)	//elastic
{postmean=get_postmean(sum, lambdas[p]*value2, lambdas2[p]*value2, lambdas3[p]*value, -9999, nss[j], variances[p], tryps[p], tryp2s[p], tryp3s[p], -9999, NULL, 7, &postprob, pars+p);}
if(trytypes[p]==8)	//ldpred (use bolt function)
{postmean=get_postmean(sum, lambdas[p]*value, lambdas2[p]*value, -9999, -9999, nss[j], variances[p], tryps[p], tryp2s[p], -9999, -9999, NULL, 4, &postprob, pars+p);}

//update residuals for remainder of block
value3=postmean-effs2[j-bitstart+p*bitlength];
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

effs2[j-bitstart+p*bitlength]=postmean;
probs2[j-bitstart+p*bitlength]=postprob;
}}	//end of j loop
}	//end of p loop
}	//end of bitstart2 loop
}
else    //finemap=1 - can not (easily) sub-divide bit
{
//set probs2 to zero
for(p=0;p<Lnum;p++)
{
for(j=0;j<bitlength;j++){probs2[j+p*bitlength]=0;}
}

//set start of probs3 to zero
for(j=0;j<bitlength;j++){probs3[j]=0;}

if(count==0)    //set end of probs3 to zero, so can check convergence
{
for(j=0;j<bitlength;j++){probs3[j+bitlength]=0;}
}

for(p=0;p<Lnum;p++) //update values for each effect
{
//subtract pth column of effs3 from effs2
for(j=0;j<bitlength;j++){effs2[j]-=effs3[j+p*bitlength];}

//get XjXTbeta based on total effects minus effects p
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &one, &alpha, cors, &bitlength, effs2, &bitlength, &beta, residuals, &bitlength);

//get statistics corresponding to vargrid, keeping track of max bf
#pragma omp parallel for private(k,j,sum) schedule(static)
for(k=0;k<num_grids;k++)
{
for(j=0;j<bitlength;j++)
{
//get XjT residuals and new effect
sum=YTdata[bitstart+j]-nss[bitstart+j]*residuals[j];

//get mean, var and log bf for each predictor
if(vargrid[k]==0)    //null model
{pmeans[j+k*bitlength]=0;pvars[j+k*bitlength]=0;bfs[j+k*bitlength]=0;gridmaxes[k]=0;}
else
{
pmeans[j+k*bitlength]=sum/(nss[bitstart+j]+variances[0]/vargrid[k]);
pvars[j+k*bitlength]=variances[0]/(nss[bitstart+j]+variances[0]/vargrid[k]);
bfs[j+k*bitlength]=.5*log(pvars[j+k*bitlength]/vargrid[k])+.5*pow(pmeans[j+k*bitlength],2)/pvars[j+k*bitlength];

if(j==0){gridmaxes[k]=bfs[j+k*bitlength];}
if(bfs[j+k*bitlength]>gridmaxes[k]){gridmaxes[k]=bfs[j+k*bitlength];}
}
}
}

//get max bf across all variances
max=gridmaxes[0];
for(k=1;k<num_grids;k++)
{
if(gridmaxes[k]>max){max=gridmaxes[k];}
}

//sum weighted bfs
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<num_grids;k++)
{
gridsums[k]=0;
for(j=0;j<bitlength;j++){gridsums[k]+=exps[bitstart+j]*exp(bfs[j+k*bitlength]-max);}
}

//which sum is best 
best=0;
for(k=1;k<num_grids;k++)
{
if(gridsums[k]>gridsums[best]){best=k;}
}

//update probs2, (approx) probs3, effs2, effs3 and variances (which actually adds on squared expectation)
if(vargrid[best]==0)    //null model
{
for(j=0;j<bitlength;j++){probs2[j+p*bitlength]=0;effs3[j+p*bitlength]=0;}
}
else
{
for(j=0;j<bitlength;j++)
{
value=exps[bitstart+j]*exp(bfs[j+best*bitlength]-max)/gridsums[best];
probs2[j+p*bitlength]=value;
probs3[j]+=value;
effs3[j+p*bitlength]=value*pmeans[j+best*bitlength];
effs2[j]+=value*pmeans[j+best*bitlength];
pars[0]+=value*(pvars[j+best*bitlength]+pow(pmeans[j+best*bitlength],2));
}
}
}
}

//get Cbeta
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &num_left, &alpha, cors, &bitlength, effs2, &bitlength, &beta, cors2, &bitlength);

//save old ess
for(p=0;p<num_left;p++){ess2[p]=ess[p];}

//compute new ess
for(p=0;p<num_left;p++){ess[p]=2*ddot_(&bitlength, effs2+p*bitlength, &one, rhos+bitstart, &one)-ddot_(&bitlength, effs2+p*bitlength, &one, cors2+p*bitlength, &one);}

//see whether converged
if(finemap==0)  //convergence based on ess
{
cflag=0;for(p=0;p<num_left;p++){cflag+=(fabs(ess[p]-ess2[p])<tol);}
if(cflag==num_left){break;}
}
else    //convergence based on max difference in (approx) probs
{
max=0;
for(j=0;j<bitlength;j++)
{
if(fabs(probs3[j]-probs3[j+bitlength])>max){max=fabs(probs3[j]-probs3[j+bitlength]);
}
}
if(max<tol){break;}

for(j=0;j<bitlength;j++){probs3[j+bitlength]=probs3[j];}
}

if(resfix==0)   //estimate variance based on expected distribution of rho, assuming prior dist prop to variance^-.5
{
if(finemap==0)  //always move (but do not allow to go below one)
{
for(p=0;p<num_left;p++) 
{
//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta + variances (note that terms 2 and 3 equal -ess, computed just above)
sum=rhosinvrhos-ess[p]+pars[p];
variances[p]=sum/bitlength*neff;
if(variances[p]<1.0){variances[p]=1.0;}
}
}
else    //only one variance - again, always move (but do not allow to go below one)
{
//get Cbeta for each effect
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &Lnum, &alpha, cors, &bitlength, effs3, &bitlength, &beta, cors4, &bitlength);

//need t(rhos) inv(cors) rhos - 2 t(rhos) beta + t(beta) cors beta + variances - sum of t(beta) R beta (note that terms 2 and 3 equal -ess, computed just above)
sum=rhosinvrhos-ess[0]+pars[0];
for(p=0;p<Lnum;p++){sum-=ddot_(&bitlength, effs3+p*bitlength, &one, cors4+p*bitlength, &one);}
variances[0]=sum/bitlength*neff;
if(variances[0]<1.0){variances[0]=1.0;}
}
}
}	//end of count loop

if(count==maxiter){printf("Warning, Window %d did not converge within %d iterations\n", bit+1, maxiter);}

cflag=0;for(p=0;p<num_left;p++){cflag+=(ess[p]>1);}
if(cflag==0)	//model(s) not terrible, so update effect sizes and probs, and maybe credible sets
{
if(finemap==0)
{
for(p=0;p<num_left;p++)
{
for(j=bitstart;j<bitend;j++)
{effs[(size_t)p*data_length+j]=effs2[j-bitstart+p*bitlength];probs[(size_t)p*data_length+j]=probs2[j-bitstart+p*bitlength];}
}
}
else    //save credible sets, effects and probs
{
for(p=0;p<Lnum;p++)
{
if(dougvar!=9)  //blank predictors whose post probs are not higher than prior probs
{
for(j=0;j<bitlength;j++)
{
if(probs2[j+p*bitlength]<=exps[bitstart+j]/sumexps){probs2[j+p*bitlength]=0;}
}
}

//sum probabilities
sum=0;for(j=0;j<bitlength;j++){sum+=probs2[j+p*bitlength];}

if(sum>=csalpha) //will be able to construct a credible set
{
//sort probabilities
sum=0;for(j=0;j<bitlength;j++){dptrs[j].value=probs2[j+p*bitlength];dptrs[j].index=j;}
qsort(dptrs, bitlength, sizeof(struct sorting_double), compare_sorting_double_rev);

//how many do we require to exceed csalpha?, and what is max test statistic
sum=dptrs[0].value;
count=1;
value=chis[bitstart+dptrs[0].index];
for(j=1;j<bitlength;j++)
{
if(sum>=csalpha){break;}
sum+=dptrs[j].value;count++;
if(chis[bitstart+dptrs[j].index]>value){value=chis[bitstart+dptrs[j].index];}
}

//compute purity between first 100 predictors
value2=1;
for(j=0;j<count;j++)
{
if(j==100){break;}
for(j2=0;j2<j;j2++)
{
if(dptrs[j2].index>=dptrs[j].index){value3=pow(cors[(size_t)dptrs[j].index*bitlength+dptrs[j2].index],2);}
else{value3=pow(cors[(size_t)dptrs[j2].index*bitlength+dptrs[j].index],2);}
if(value3<value2){value2=value3;}
}
}

if(value2<0.5&&dougvar==8)  //blank probs
{
for(j=0;j<bitlength;j++){probs2[j+p*bitlength]=0;}
printf("blank %d %d purity %f\n", bit, p, value2);
}

fprintf(output4,"%d %.4f %d %.4e %.4f %s", count4+1, sum, count, erfc(pow(value,.5)*M_SQRT1_2), value2, preds[bitstart+dptrs[0].index]);
for(j=1;j<count;j++){fprintf(output4,",%s", preds[bitstart+dptrs[j].index]);}
fprintf(output4,"\n");
count4++;
}
}

//load up effects amd probs (1 - chance all zero)
#pragma omp parallel for private(j,p,sum) schedule(static)
for(j=0;j<bitlength;j++)
{
effs[bitstart+j]=effs2[j];
sum=1;for(p=0;p<Lnum;p++){sum*=(1-probs2[j+p*bitlength]);}
probs[bitstart+j]=1-sum;
}
}
}
else	//one or more models is suspect - leave effect sizes at starting estimates (probs will be zero)
{
if(mcmcfinal==0)
{
printf("Warning, Variational Bayes failed for Window %d, will revert to starting estimates\n", bit+1);
ecount++;
wcount+=bitlength;
}
else{printf("Warning, Variational Bayes failed for Window %d\n", bit+1);}
}

if(prsvar==1)	//only here if MCMC failed, so will set effects to the final estimate (will have num_left=1)
{
for(count2=0;count2<num_blocks;count2++)
{
for(j=bitstart;j<bitend;j++){effs4[(size_t)count2*data_length+j]=effs[j];}
}
}
}	//end of vb
}	//end of bit loop
printf("\n");
if(ecount>0){printf("%d windows failed (in total, %d predictors)\n\n", ecount, wcount);}

if(condition==1){fclose(output3);}
if(finemap==1){fclose(output4);}

////////

if(finemap==0)
{
if(skipcv==0)	//save effects and probs for best model (maybe ensemble)
{
if(num_focals==1){sprintf(filename2,"%s.effects",outfile);}
else{sprintf(filename2,"%s.focal%d.effects", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(ensemble==0){fprintf(output2,"Predictor A1 A2 Centre Model%d", best+1);}
else{fprintf(output2,"Predictor A1 A2 Centre Ensemble");}
if(prsvar==1)
{
for(p=0;p<num_blocks;p++){fprintf(output2," Sampling%d", p+1);}
}
fprintf(output2,"\n");

for(j=0;j<data_length;j++)
{
sum=0;for(p=0;p<num_left;p++){sum+=predweights[p]*effs[(size_t)p*data_length+j];}
fprintf(output2, "%s %c %c %.6f %.4e", preds[j], al1[j], al2[j], centres[j], sum*mults[j]*escale);

if(prsvar==1)
{
for(p=0;p<num_blocks;p++){fprintf(output2," %.4e", effs4[(size_t)p*data_length+j]*mults[j]*escale);}
}
fprintf(output2,"\n");
}
fclose(output2);

if(num_focals==1){sprintf(filename3,"%s.probs",outfile);}
else{sprintf(filename3,"%s.focal%d.probs", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Posterior_Probability\n");

for(j=0;j<data_length;j++)
{
sum=0;for(p=0;p<num_left;p++){sum+=predweights[p]*probs[(size_t)p*data_length+j];}
fprintf(output3, "%s %.4f\n", preds[j], sum);
}
fclose(output3);

printf("Effect sizes saved in %s, with posterior probabilities in %s\n\n", filename2, filename3);
}

if(skipcv==1||megasave==1)	//save effects and probs for all models
{
if(num_focals==1){sprintf(filename4,"%s.full.effects",outfile);}
else{sprintf(filename4,"%s.focal%d.full.effects", outfile, cur_focal+1);}
if((output4=fopen(filename4,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename4);exit(1);}

fprintf(output4,"Predictor A1 A2 Centre");
for(p=0;p<num_left;p++){fprintf(output4," Model%d", p+1);}
fprintf(output4,"\n");

for(j=0;j<data_length;j++)
{
fprintf(output4, "%s %c %c %.6f", preds[j], al1[j], al2[j], centres[j]);
for(p=0;p<num_left;p++){fprintf(output4," %.4e", effs[(size_t)p*data_length+j]*mults[j]*escale);}
fprintf(output4,"\n");
}
fclose(output4);

if(num_focals==1){sprintf(filename5,"%s.full.probs",outfile);}
else{sprintf(filename5,"%s.focal%d.full.probs", outfile, cur_focal+1);}
if((output5=fopen(filename5,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename5);exit(1);}

fprintf(output5,"Predictor");
for(p=0;p<num_left;p++){fprintf(output5," Model%d", p+1);}
fprintf(output5,"\n");

for(j=0;j<data_length;j++)
{
fprintf(output5, "%s", preds[j]);
for(p=0;p<num_left;p++){fprintf(output5, " %.4f", probs[(size_t)p*data_length+j]);}
fprintf(output5, "\n");
}
fclose(output5);

printf("Effect sizes for all models saved in %s, with posterior probabilities in %s\n\n", filename4, filename5);
}
}
else    //save effects and probs for finemapping model
{
if(num_focals==1){sprintf(filename2,"%s.effects",outfile);}
else{sprintf(filename2,"%s.focal%d.effects", outfile, cur_focal+1);}
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre Effect\n");
for(j=0;j<data_length;j++){fprintf(output2, "%s %c %c %.6f %.4e\n", preds[j], al1[j], al2[j], centres[j], effs[j]*mults[j]);}
fclose(output2);

if(num_focals==1){sprintf(filename3,"%s.probs",outfile);}
else{sprintf(filename3,"%s.focal%d.probs", outfile, cur_focal+1);}
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor Posterior_Probability\n");

for(j=0;j<data_length;j++){fprintf(output3, "%s %.4f\n", preds[j], probs[j]);}
fclose(output3);

printf("Effect sizes and posterior probabilities saved in %s and %s, ", filename2, filename3);
if(num_focals==1){printf("with credible sets saved in %s.credibles\n\n", outfile);}
else{printf("with credible sets saved in %s.focal%d.credibles\n\n", outfile, cur_focal+1);}
}

if(condition==1)
{
if(num_focals==1){printf("Conditional analysis results are saved in %s.conditional\n\n", outfile);}
else{printf("Conditional analysis results are saved in %s.focal%d.conditional\n\n", outfile, cur_focal+1);}
}

////////

free(powers);
free(vargrid);
if(finemap==0)
{
free(trytypes);free(trylams);free(tryscales);free(tryps);free(tryp2s);free(tryp3s);free(tryp4s);free(tryf2s);
free(lambdas);free(lambdas2);free(lambdas3);free(lambdas4);
}
for(q=1;q<num_sums3;q++){free(Mnss[q]);free(Mchis[q]);free(Mrhos[q]);free(Ma1freq[q]);}
free(Mnss);free(Mchis);free(Mrhos);free(Ma1freq);
free(blockstarts);free(blockends);
for(s=0;s<num_cors;s++)
{
for(bit=0;bit<bittotal;bit++){free(Dsizes[s][bit]);free(Duse[s][bit]);free(Duse2[s][bit]);free(Dsigns[s][bit]);}
free(Dkp[s]);free(Dkp2[s]);free(Dindexes[s]);free(Dsizes[s]);free(Duse[s]);free(Duse2[s]);free(Dsigns[s]);
free(Dchr[s]);free(Dpreds[s]);free(Dbp[s]);free(Dal1[s]);free(Dal2[s]);
}
free(Dns);free(Dnp);free(Dnp2);free(Dkp);free(Dkp2);free(Dindexes);free(Dsizes);free(Duse);free(Duse2);free(Dsigns);free(Dchr);free(Dpreds);free(Dbp);free(Dal1);free(Dal2);
if(metasum==1){free(cohers2);}
for(q=0;q<num_parts+1;q++){free(pindexes[q]);free(pweights[q]);}free(pindexes);free(pweights);
free(rjksums);free(datarands);free(highlds);
if(skipcv==0){free(nss2);free(rhos2);free(rhos3);}
free(exps);free(cors_short);free(cors);
if(num_sums3>1){free(nss4);free(chis4);free(rhos4);}
if(num_sums3>1||metasum==1){free(cohers);}
free(YTdata);free(cors2);free(effs);free(effs2);free(effs3);free(probs);free(probs2);free(probs3);free(residuals);
free(ess);free(ess2);free(variances);free(pars);
free(pmeans);free(pvars);free(bfs);free(gridmaxes);free(gridsums);
free(predvars);free(predtops);free(predcors);free(predweights);
free(dptrs);
if(finemap==1){free(cors4);}
if(prsvar==1){free(effs4);}
if(mcmcfinal==1){free(manyrands);}

///////////////////////////

/*
//old fine mapping codes


if(finemap==1&&dougvar!=8)  //remove (near) duplicates - retain records the tagging predictor, retain2 records number tags, shares indicates share of posterior
{
for(j=0;j<bitlength;j++){retain[j]=j;retain2[j]=1;shares[j]=exps[bitstart+j];}

for(j=0;j<bitlength;j++)
{
if(retain[j]==j)    //j has not been replaced
{
for(j2=j+1;j2<bitlength;j2++)
{
if(retain[j2]==j2)  //remove j2 if highly correlated with j (remember correlations have been shrunk)
{
if(cors[(size_t)j*bitlength+j2]/shrink>dupthresh)    //positively correlated - adjust exps and effs
{
exps[bitstart+j]+=exps[bitstart+j2];
exps[bitstart+j2]=0;
for(p=0;p<num_left;p++){effs[(size_t)p*data_length+bitstart+j]+=effs[(size_t)p*data_length+bitstart+j2];effs[(size_t)p*data_length+bitstart+j2]=0;}
retain[j2]=j;
retain2[j]++;
}
if(cors[(size_t)j*bitlength+j2]/shrink<-dupthresh)    //negatively correlated - adjust exps and effs
{
exps[bitstart+j]+=exps[bitstart+j2];
exps[bitstart+j2]=0;
for(p=0;p<num_left;p++){effs[(size_t)p*data_length+bitstart+j]-=effs[(size_t)p*data_length+bitstart+j2];effs[(size_t)p*data_length+bitstart+j2]=0;}
retain[j2]=j;
retain2[j]++;
}
}}
}}

for(j=0;j<bitlength;j++){shares[j]=shares[j]/exps[bitstart+retain[j]];}
}
else    //need to set retain and retain2 so that clump.c works
{
for(j=0;j<bitlength;j++){retain[j]=j;retain2[j]=1;shares[j]=1.0;}
}


//big count is number of windows with a significant predictor
bigcount=0.0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
for(j=bitstart;j<bitend;j++)
{
if(chis[j]>value){bigcount++;break;}
}
}

if(dougvar==10){bigcount*=2;}
if(dougvar>1000){bigcount=dougvar;}

printf("%d of the %d windows contain a predictor with P<%.4e\n", (int)bigcount, bittotal, cutoff);

//big count can not be zero nor exceed 0.9 her/bigvar
if(bigcount==0){bigcount=1.0;}
if(bigcount>0.9*her/bigvar){bigcount=0.9*her/bigvar;}

//small variance set so expected variance explained equals her
smallvar=(her-bigcount*bigvar)/(data_length-bigcount);

//small variance must be at most 1% of big variance
if(smallvar>bigvar/100){smallvar=bigvar/100;}

printf("The prior distribution expects there are %.1f predictors with big effects (expected variance explained %.2e), while the remainder have small effects (expected variance explained %.2e)\n\n", bigcount, bigvar, smallvar);


else    //finemapping - must have bolt model and num_left=1
{
#pragma omp parallel for private(p,start,j,j2,j3,sum,value,postmean,postprob,postsamp,value3) schedule(dynamic)
for(p=0;p<num_chains;p++)
{
start=0;
for(j3=0;j3<bitlength2;j3++)   //will alternate between moving left to right and right to left
{
if(p%2==0){j=bitstart2+j3;}
else{j=bitend2-1-j3;}

if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-nss[j]*(residuals[j-bitstart2+p*bitsize]-effs3[j-bitstart+p*bitlength]);

//value is probability of large effect
value=exps[j]*hits[p];
if(value>0.99){value=0.99;}

if(dougvar!=7){postmean=get_postsamp_pragma(sum, bigvar, smallvar, -9999, -9999, nss[j], variances[p], value, 1.0-value, -9999, -9999, 4, &postprob, &postsamp, manyrands+p*128, &start);}
else{postmean=get_postsamp_pragma(sum, bigvar, 0.0, -9999, -9999, nss[j], variances[p], value, 1.0-value, -9999, -9999, 4, &postprob, &postsamp, manyrands+p*128, &start);}

//update residuals for remainder of block - only necessary if estimate has changed
value3=postsamp-effs3[j-bitstart+p*bitlength];
if(value3!=0)
{
if(p%2==0) //moving to right
{
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}
}
else    //moving to left
{
for(j2=j-1;j2>=bitstart2;j2--){residuals[j2-bitstart2+p*bitsize]+=value3*cors[(size_t)(j2-bitstart)*bitlength+j-bitstart];}
}
}

if(count>=num_burn){effs2[j-bitstart+p*bitlength]+=postmean;probs2[j-bitstart+p*bitlength]+=postprob;}
effs3[j-bitstart+p*bitlength]=postsamp;
hits2[p]+=postprob;
}}	//end of j loop
}	//end of p loop
}


if(finemap==1)  //sample new hits from beta(1+hits2,data_length/bigcount+data_length-hits2), then multiply by data_length 
{
for(p=0;p<num_chains;p++)
{
//if X~gamma(alpha,1) and Y ~gamma(beta,1), then X/(X+Y) is beta(alpha,beta)
value=rgamma(1+hits2[p]);
value2=rgamma((double)data_length/bigcount+bitlength-hits2[p]);
if(dougvar!=12){hits[p]=value/(value+value2)*data_length;}
if(count>=num_burn){hits3+=hits[p];}
}
}

else    //finemapping - must have bayesr model and num_left=1
{
p=0;
for(j=bitstart2;j<bitend2;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-nss[j]*(residuals[j-bitstart2+p*bitsize]-effs2[j-bitstart+p*bitlength]);

//value is probability of large effect
value=exps[j]*bigcount;
if(value>0.99){value=0.99;}

postmean=get_postmean(sum, bigvar, smallvar, -9999, -9999, nss[j], variances[p], value, 1.0-value, -9999, -9999, NULL, 4, &postprob, &postvar);

//update residuals for remainder of block
value3=postmean-effs2[j-bitstart+p*bitlength];
for(j2=j+1;j2<bitend2;j2++){residuals[j2-bitstart2+p*bitsize]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

effs2[j-bitstart+p*bitlength]=postmean;
probs2[j-bitstart+p*bitlength]=postprob;
pars[j-bitstart+p*bitlength]=postvar;
}}	//end of j loop
}

//inf susie model

if(dougvar!=12&&dougvar!=11&&dougvar!=10) //update final effects (ridge model)
{
//get XjXTbeta based on total effects 
alpha=1.0;beta=0.0;
dsymm_("L", "L", &bitlength, &one, &alpha, cors, &bitlength, effs2, &bitlength, &beta, residuals, &bitlength);

for(j=bitstart;j<bitend;j++)
{
if(exps[j]>0)
{
//get XjT residuals
sum=YTdata[j]-nss[j]*(residuals[j-bitstart]-effs3[j-bitstart+Lbit*bitlength]);
postmean=get_postmean(sum, exps[j]*smallvar, -9999, -9999, -9999, nss[j], variances[0], -9999, -9999, -9999, -9999, NULL, 3, NULL, &postvar);

//update residuals for remainder of block
value3=postmean-effs3[j-bitstart+Lbit*bitlength];
#pragma omp parallel for private(j2) schedule(static)
for(j2=j+1;j2<bitend;j2++){residuals[j2-bitstart]+=value3*cors[(size_t)(j-bitstart)*bitlength+j2-bitstart];}

//update effs2, effs3 and variance
effs2[j-bitstart]+=postmean-effs3[j-bitstart+Lbit*bitlength];
effs3[j-bitstart+Lbit*bitlength]=postmean;
pars[j-bitstart]+=postvar;
}
}	//end of j loop
}

//code for resampling noise


if(dougvar2!=-9999)
{
value=(double)bitmax/1024*bitmax/1024/1024*8;
if(value>1){printf("Warning, to store the decomps requires %.1f Gb\n\n", value);}

U=malloc(sizeof(double)*bitmax*bitmax);
E=malloc(sizeof(double)*bitmax);
}

if(dougvar2!=-9999) //decomp correlations
{
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
for(j=k;j<bitlength;j++){U[(size_t)k*bitlength+j]=cors[(size_t)k*bitlength+j];}
for(j=k;j<bitlength;j++){U[(size_t)j*bitlength+k]=cors[(size_t)k*bitlength+j];}
}

lwork=-1;
dsyev_("V", "L", &bitlength, U, &bitlength, E, &wkopt, &lwork, &info);
if(info!=0){printf("Decomp error 1; please tell Doug\n\n");exit(1);}
lwork=(int)wkopt;

work=malloc(sizeof(double)*lwork);
dsyev_("V", "L", &bitlength, U, &bitlength, E, work, &lwork, &info);
if(info!=0){printf("Decomp error 2; please tell Doug\n\n");exit(1);}
}

if(dougvar2!=-9999)
{
rhos2=malloc(sizeof(double)*data_length);
rhos3=malloc(sizeof(double)*data_length);

for(j=0;j<data_length;j++){rhos2[j]=rhos[j];}
}

if(dougvar2!=-9999) //resample rhos
{
//new values are d x U diag( E/ (E+d) ) UT (rho/d + beta)
for(j=bitstart;j<bitend;j++){rhos[j]=rhos2[j]+dougvar2*effs2[j-bitstart];}
alpha=1.0;beta=0.0;
dgemv_("T", &bitlength, &bitlength, &alpha, U, &bitlength, rhos+bitstart, &one, &beta, rhos3+bitstart, &one);

for(j=bitstart;j<bitend;j++){rhos3[j]*=E[j-bitstart]/(E[j-bitstart]+dougvar2);}
alpha=1.0;beta=0.0;
dgemv_("N", &bitlength, &bitlength, &alpha, U, &bitlength, rhos3+bitstart, &one, &beta, rhos+bitstart, &one);

for(j=bitstart;j<bitend;j++){YTdata[j]=rhos[j]*nss[j];}
}
*/

//tried to infer parameters from summary statistics
/*

if(finemap==0)  //test different models
{
//get expected heritability tagged by each predictor (cant be negative)
for(j=0;j<data_length;j++)
{
exps[j]=0;
for(q=0;q<total;q++){exps[j]+=stats[q]*svars[q][j]/ssums[q][q];}
if(exps[j]<=0){exps[j]=1e-6;}
}

//reload stats, only ensuring not too small
for(j=0;j<data_length;j++)
{
snss[j]=nss[j];
schis[j]=chis[j];
if(schis[j]<1e-6){schis[j]=1e-6;}
}

//get weighted variance of chisq statistics
sum=0;sumsq=0;value=0;
for(j=0;j<data_length;j++){sum+=schis[j]/rjksums[j];sumsq+=pow(schis[j],2)/rjksums[j];value+=pow(rjksums[j],-1);}
mean=sum/value;
var=sumsq/value-pow(mean,2);
printf("observed mean %f var is %f\n", mean, var);

sum=0;for(j=0;j<data_length;j++){sum+=rjksums[j];}
printf("mean ld is %f\n", sum/data_length);

for(p=0;p<num_try;p++)  //get expected variance for each set of parameters
{
sum=0;sumsq=0;value=0;
for(j=0;j<data_length;j++)
{
sum+=(1+snss[j]*exps[j])/rjksums[j];
if(trytypes[p]==5||trytypes[p]==6)  //bayesr models
{sumsq+=(tryps[p]*pow(1+snss[j]*exps[j]*lambdas[p],2)+tryp2s[p]*pow(1+snss[j]*exps[j]*lambdas2[p],2)+tryp3s[p]*pow(1+snss[j]*exps[j]*lambdas3[p],2)+tryp4s[p]*pow(1+snss[j]*exps[j]*lambdas4[p],2))/rjksums[j];}
if(trytypes[p]==8)  //ldpred models
{sumsq+=(tryps[p]*pow(1+snss[j]*exps[j]*lambdas[p],2)+tryp2s[p]*pow(1+snss[j]*exps[j]*lambdas2[p],2))/rjksums[j];}
value+=pow(rjksums[j],-1);
}
mean=sum/value;
value=3*sumsq/value-pow(mean,2);

printf("probs %f %f %f %f mean %f var %f\n", tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], mean, value);

//how about likelihood?
like=0;
for(j=0;j<data_length;j++)
{
if(trytypes[p]==5||trytypes[p]==6)  //bayesr models
{
value=1+snss[j]*exps[j]*lambdas[p];
value2=1+snss[j]*exps[j]*lambdas2[p];
value3=1+snss[j]*exps[j]*lambdas3[p];
value4=1+snss[j]*exps[j]*lambdas4[p];
like+=log(tryps[p]*pow(schis[j]*value,-0.5)*exp(-0.5*schis[j]/value)+tryp2s[p]*pow(schis[j]*value2,-0.5)*exp(-0.5*schis[j]/value2)+tryp3s[p]*pow(schis[j]*value3,-0.5)*exp(-0.5*schis[j]/value3)+tryp4s[p]*pow(schis[j]*value4,-0.5)*exp(-0.5*schis[j]/value4))/rjksums[j];
}
if(trytypes[p]==8)  //ldpred models
{
value=1+snss[j]*exps[j]*lambdas[p];
value2=1+snss[j]*exps[j]*lambdas2[p];
like+=log(tryps[p]*pow(schis[j]*value,-0.5)*exp(-0.5*schis[j]/value)+tryp2s[p]*pow(schis[j]*value2,-0.5)*exp(-0.5*schis[j]/value2))/rjksums[j];
}
if(like!=like||isinf(like)){printf("j %d of %d values %f %f %f %f , chisq %f vales %f %f %f %f, size %e exp %e centres %f mults %f tag %f\n", j+1, data_length ,  tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], schis[j], value, value2, value3, value4, exp(-schis[j]), exps[j], centres[j], mults[j], rjksums[j]);exit(1);}
}

if(p==0){best=0;likenull=like;}
if(like>likenull){best=p;likenull=like;}
printf("probs %f %f %f %f like %f\n", tryps[p], tryp2s[p], tryp3s[p], tryp4s[p], like-likenull);
}
printf("best %f %f %f %f %d\n", tryps[best], tryp2s[best], tryp3s[best], tryp4s[best], best+1);
}
*/

//tried to estimate scaling based on rhosinvrhos
/*
if(noscale==2)  //scale based on rhosinvrhos
{
//re-open focal correlations
s=sumpops[0];

sprintf(filename2,"%s.cors.bin", corstems[s]);
if((input2=fopen(filename2,"rb"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}

sum2=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

if(bit%100==0)
{
printf("Calculating scaling for Window %d of %d\n", bit+1, bittotal);

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Calculating scaling for Window %d of %d\n", bit+1, bittotal);
fclose(output);
}

//get average sample size
sum=0;for(j=bitstart;j<bitend;j++){sum+=nss[j];}
neff=sum/bitlength;

//read lower triangle of focal correlations
fseeko(input2, Dindexes[s][bit], SEEK_SET);
for(k=0;k<Dsizes[s][bit][0]-1;k++)
{
count2=fread(cors_short+(size_t)k*Dsizes[s][bit][0]+k+1, sizeof(unsigned short), Dsizes[s][bit][0]-k-1, input2);
if(count2!=Dsizes[s][bit][0]-k-1)
{printf("Error reading correlations for Window %d from %s\n\n", bit+1, filename2);exit(1);}
}

//now extract correlations for all predictors (remember shrink, but no need for signs)
#pragma omp parallel for private(k,j) schedule(static)
for(k=0;k<bitlength;k++)
{
cors[(size_t)k*bitlength+k]=1.0;
for(j=k+1;j<bitlength;j++)
{cors[(size_t)k*bitlength+j]=(0.00005*cors_short[(size_t)Duse[s][bit][k]*Dsizes[s][bit][0]+Duse[s][bit][j]]-1)*shrink;}
}

//compute rhosinvrhos
rhosinvrhos=cg_solve_lower_norm(cors, bitlength, rhos+bitstart, 1e-6);
sum2+=rhosinvrhos*neff;
}
fclose(input2);

value=sum2/data_length;
printf("2Estimated inflation of test statistics is %f\n", value);

value2=pow(value,-1);
for(j=0;j<data_length;j++)
{
chis[j]*=value2;
if(rhos[j]>0){rhos[j]=pow(chis[j]/(chis[j]+nss[j]),.5);}
else{rhos[j]=-pow(chis[j]/(chis[j]+nss[j]),.5);}
}
}
*/

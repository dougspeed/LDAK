/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.
*/

///////////////////////////

//Impute summary statistics - most code copied from lites.c (and could simplify because must have num_cors=1)

///////////////////////////

sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

//need multi-version pointers

Mnss=malloc(sizeof(double*)*num_sums3);
Mchis=malloc(sizeof(double*)*num_sums3);
Mrhos=malloc(sizeof(double*)*num_sums3);
Ma1freq=malloc(sizeof(double*)*num_sums3);
Mnss[0]=nss;Mchis[0]=chis;Mrhos[0]=rhos;Ma1freq[0]=a1freq;

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
}

////////

//compare frequencies of sum stats and correlations (maybe infer ancestry for bulk), then blank missing predictors
#include "freqcheck.c"

//can save some space by freeing Dpreds
for(s=0;s<num_cors;s++)
{
for(j=0;j<Dnp[s];j++){free(Dpreds[s][j]);}
}

//work out blocks (might be redundant)

sprintf(filename2,"%s.cors.windows", corstems[sumpops[0]]);
total=countrows(filename2)-1;
bittotal=get_windows(filename2, data_length, keeppreds_use, blockstarts, blockends, total);

printf("The analysis will be performed using %d predictors across %d windows\n\n", data_length, bittotal);

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

if(end>start)	//contains predictors
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

////////

//ready to do imputation

value=(double)bitmax/1024*bitmax/1024/1024*10;
if(value>1){printf("Warning, to store the correlations requires %.1f Gb\n\n", value);}

cors_short=malloc(sizeof(unsigned short)*bitmax*bitmax);
cors=malloc(sizeof(double)*bitmax*bitmax);

count=0;
count2=0;
count3=0;
for(bit=0;bit<bittotal;bit++)
{
bitstart=blockstarts[bit];
bitend=blockends[bit];
bitlength=blockends[bit]-blockstarts[bit];

count4=0;for(j=bitstart;j<bitend;j++){count4+=(nss[j]>0);}
if(count4>0&&count4<bitlength)  //will impute window
{
count+=bitlength-count4;
count2=(4*count4<bitlength);
}
count3+=(count4==0);
}
if(count3>0){printf("Warning, %d windows are missing summary statistics for all predictors (so can not be imputed)\n\n", count3);}

if(count>0)
{
printf("Imputing summary statistics for %d predictors\n", count);
if(count2>0){printf("Warning, %d windows have summary statistics for fewer than a quarter of the predictors (their imputation quality may be lower)\n\n", count2);}

impute_sums(0, Mnss, Mrhos, Mchis, data_length, bittotal, blockstarts, blockends, sumpops[0], NULL, corstems, -9999, Dindexes, Dsizes, Duse, Duse2, Dsigns, shrink, outfile, cors_short, cors, preds);
}
else
{
printf("Warning, no summary statistics will be imputed\n\n");
}

sprintf(filename,"%s.imputed",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

if(gotfreq==0){fprintf(output, "Predictor\tA1\tA2\tZ\tn\tP\n");}
else{fprintf(output, "Predictor\tA1\tA2\tZ\tn\tA1Freq\tP\n");}

for(j=0;j<data_length;j++)
{
if(nss[j]>0)
{
if(rhos[j]>0){value=pow(Mchis[0][j],.5);}
else{value=-pow(Mchis[0][j],.5);}
if(gotfreq==0){fprintf(output,"%s\t%c\t%c\t%.4f\t%.0f\t%.4e\n", preds[j], al1[j], al2[j], value, Mnss[0][j], erfc(fabs(value)*M_SQRT1_2));}
else{fprintf(output,"%s\t%c\t%c\t%.4f\t%.0f\t%.4f\t%.4e\n", preds[j], al1[j], al2[j], value, Mnss[0][j], Ma1freq[0][j], erfc(fabs(value)*M_SQRT1_2));}
}
}
fclose(output);

////////

free(Mnss);free(Mchis);free(Mrhos);free(Ma1freq);
free(blockstarts);free(blockends);
for(s=0;s<num_cors;s++)
{
for(bit=0;bit<bittotal;bit++){free(Dsizes[s][bit]);free(Duse[s][bit]);free(Duse2[s][bit]);free(Dsigns[s][bit]);}
free(Dkp[s]);free(Dkp2[s]);free(Dindexes[s]);free(Dsizes[s]);free(Duse[s]);free(Duse2[s]);free(Dsigns[s]);
free(Dchr[s]);free(Dpreds[s]);free(Dbp[s]);free(Dal1[s]);free(Dal2[s]);
}
free(Dns);free(Dnp);free(Dnp2);free(Dkp);free(Dkp2);free(Dindexes);free(Dsizes);free(Duse);free(Duse2);free(Dsigns);free(Dchr);free(Dpreds);free(Dbp);free(Dal1);free(Dal2);
free(cors_short);free(cors);

///////////////////////////


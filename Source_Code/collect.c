/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//combine multi-ancestry effect sizes

///////////////////////////

printf("Merging the ancestry-specific scorefiles\n\n");
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output,"Merging the ancestry-specific scorefiles\n");
fclose(output);

//read bim files for all correlations corresponding
Dns=malloc(sizeof(int)*num_cors);
Dnp=malloc(sizeof(int)*num_cors);
Dchr=malloc(sizeof(int*)*num_cors);
Dpreds=malloc(sizeof(char**)*num_cors);
Dbp=malloc(sizeof(double)*num_cors);
Dal1=malloc(sizeof(char*)*num_cors);
Dal2=malloc(sizeof(char*)*num_cors);

if(num_cors>1){printf("Reading details of the predictor-predictor correlations\n");}

flag=0;
for(s=0;s<num_cors;s++)
{
count=0;for(q=0;q<num_sums1;q++){count+=(pritraits2[q]==2&&sumpops2[q]==s);}

if(count==0){Dnp[s]=0;}
else
{
sprintf(filename2, "%s.cors.root", corstems[s]);
read_integers(filename2, Dnp+s, 1, NULL, 2, 4, 0);

Dchr[s]=malloc(sizeof(int)*Dnp[s]);
Dpreds[s]=malloc(sizeof(char*)*Dnp[s]);
Dbp[s]=malloc(sizeof(double)*Dnp[s]);
Dal1[s]=malloc(sizeof(char)*Dnp[s]);
Dal2[s]=malloc(sizeof(char)*Dnp[s]);

sprintf(filename2, "%s.cors.bim", corstems[s]);
flag+=read_bimfile_lite(filename2, Dchr[s], Dpreds[s], Dbp[s], Dal1[s], Dal2[s], Dnp[s]);
}
}

//work out union of predictors (use Dns as counter)
count=0;for(s=0;s<num_cors;s++){count+=Dnp[s];}
preds=malloc(sizeof(char*)*count);
chr=malloc(sizeof(int)*count);
bp=malloc(sizeof(double)*count);
al1=malloc(sizeof(char)*count);
al2=malloc(sizeof(char)*count);

num_preds_use=0;
for(s=0;s<num_cors;s++){Dns[s]=0;}
while(1)
{
//get position of next predictor
readint=-1;
for(s=0;s<num_cors;s++)
{
if(Dns[s]<Dnp[s])
{
if(readint==-1){s2=s;readint=Dchr[s][Dns[s]];readdouble=Dbp[s][Dns[s]];}
if(Dchr[s][Dns[s]]<readint||(Dchr[s][Dns[s]]==readint&&Dbp[s][Dns[s]]<readdouble)){s2=s;readint=Dchr[s][Dns[s]];readdouble=Dbp[s][Dns[s]];}
}
}
if(readint==-1){break;}

//add this predictor
copy_string(preds,num_preds_use,Dpreds[s2][Dns[s2]]);
chr[num_preds_use]=readint;
bp[num_preds_use]=readdouble;
al1[num_preds_use]=Dal1[s2][Dns[s2]];
al2[num_preds_use]=Dal2[s2][Dns[s2]];
num_preds_use++;

//jump over matching predictors
for(s=0;s<num_cors;s++)
{
if(Dns[s]<Dnp[s])
{
if(Dchr[s][Dns[s]]==readint&&Dbp[s][Dns[s]]==readdouble&&strcmp(Dpreds[s][Dns[s]],preds[num_preds_use-1])==0){Dns[s]++;}
}
}
}

for(s=0;s<num_cors;s++)
{
if(Dnp[s]>0)
{
for(j=0;j<Dnp[s];j++){free(Dpreds[s][j]);}
free(Dchr[s]);free(Dpreds[s]);free(Dbp[s]);free(Dal1[s]);free(Dal2[s]);
}
}
free(Dns);free(Dnp);free(Dchr);free(Dpreds);free(Dbp);free(Dal1);free(Dal2);

//allocate centres and effects, and set to missing and zero
blupcentres=malloc(sizeof(double*)*num_focals);
blupfactors=malloc(sizeof(double*)*num_focals);
for(cur_focal=0;cur_focal<num_focals;cur_focal++)
{
blupcentres[cur_focal]=malloc(sizeof(char*)*num_preds_use);
blupfactors[cur_focal]=malloc(sizeof(char*)*num_preds_use);

for(j=0;j<num_preds_use;j++){blupcentres[cur_focal][j]=-9999;blupfactors[cur_focal][j]=-9999;}
}

//read predictors and effects for each profile, then get indexes and check alleles consistent
rs=malloc(sizeof(char)*1000000);

for(cur_focal=0;cur_focal<num_focals;cur_focal++)
{
sprintf(filename2,"%s.focal%d.effects", outfile, cur_focal+1);
count=countrows(filename2)-1;
printf("Reading effects for %d predictors from %s\n", count, filename2);

allpreds=malloc(sizeof(char*)*count);
allal1=malloc(sizeof(char)*count);
allal2=malloc(sizeof(char)*count);
effs=malloc(sizeof(double)*count);
effs2=malloc(sizeof(double)*count);
indexer=malloc(sizeof(int)*count);

if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input2, "%c", &readchar);}

for(j=0;j<count;j++)
{
if(fscanf(input2, "%s %c %c %lf %lf ", rs, allal1+j, allal2+j, effs+j, effs2+j)!=5)
{printf("Error reading Row %d of %s\n\n", j+2, filename2);exit(1);}
copy_string(allpreds,j,rs);
}
fclose(input2);

found=0;
if(flag<2)  //try fast sort (small risk of problems when multiple correlations have predictors with matching basepairs)
{
k=0;
for(j=0;j<count;j++)
{
while(k<num_preds_use)
{
if(strcmp(allpreds[j],preds[k])==0){indexer[j]=k;found++;break;}
k++;
}
}
}

if(found<count) //use proper sort
{
count2=find_strings(allpreds, count, preds, num_preds_use, NULL, indexer, NULL, NULL, NULL, NULL, 3);
if(count2<count){printf("Error 3NXW, please tell doug (%d and %d)\n\n", count, count2);exit(1);}
}

for(j=0;j<count;j++)
{
j2=indexer[j];
if((allal1[j]!=al1[j2]&&allal1[j]!=al2[j2])||(allal2[j]!=al1[j2]&&allal2[j]!=al2[j2]))	//inconsistent
{printf("Error 4NXW, please tell doug (%c %c and %c %c)\n\n", allal1[j], allal2[j], al1[j], al2[j]);exit(1);}

if(allal1[j]==al1[j2]){blupcentres[cur_focal][j2]=effs[j];blupfactors[cur_focal][j2]=effs2[j];}
else{blupcentres[cur_focal][j2]=2-effs[j];blupfactors[cur_focal][j2]=-effs2[j];}
}

for(j=0;j<count;j++){free(allpreds[j]);}
free(allpreds);free(allal1);free(allal2);free(effs);free(effs2);free(indexer);
}

//save combined effects

sprintf(filename2,"%s.combined.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2");
for(cur_focal=0;cur_focal<num_focals;cur_focal++){fprintf(output2," Centre%d Effect%d", cur_focal+1, cur_focal+1);}
fprintf(output2,"\n");

for(j=0;j<num_preds_use;j++)
{
count=0;for(cur_focal=0;cur_focal<num_focals;cur_focal++){count+=(blupcentres[cur_focal][j]!=-9999);}
if(count>0)
{
fprintf(output2, "%s %c %c", preds[j], al1[j], al2[j]);
for(cur_focal=0;cur_focal<num_focals;cur_focal++)
{
if(blupcentres[cur_focal][j]!=-9999){fprintf(output2, " %.6f %.4e", blupcentres[cur_focal][j], blupfactors[cur_focal][j]);}
else{fprintf(output2, " NA NA");}
}
fprintf(output2,"\n");
}
}
fclose(output2);

printf("Multi-ancestry model saved in %s\n\n", filename2);

if(dougvar==9)  //make meta-analysis versions
{
readdoubles=malloc(sizeof(double)*num_focals);

for(cur_focal=0;cur_focal<num_focals;cur_focal++)   //weights are estimated correlations
{
sprintf(filename2,"%s.focal%d.best", outfile, cur_focal+1);
read_integers(filename2, &best, 1, NULL, 1, 1, 0);
printf("model %s best %d\n", filename2, best);

sprintf(filename2,"%s.focal%d.cors", outfile, cur_focal+1);
read_values(filename2, readdoubles+cur_focal, 1, NULL, 2, best, 0);
printf("correlation is %f\n", readdoubles[cur_focal]);
if(readdoubles[cur_focal]<0.01){readdoubles[cur_focal]=0.01;}
}

sprintf(filename2,"%s.meta.effects",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

fprintf(output2,"Predictor A1 A2 Centre Meta1 Meta2\n");
for(j=0;j<num_preds_use;j++)
{
fprintf(output2, "%s %c %c NA", preds[j], al1[j], al2[j]);

sum=0;sum2=0;sum3=0;
for(cur_focal=0;cur_focal<num_focals;cur_focal++)
{
if(blupcentres[cur_focal][j]!=-9999)
{
sum+=readdoubles[cur_focal]*blupcentres[cur_focal][j];
sum2+=readdoubles[cur_focal]*blupfactors[cur_focal][j];
sum3+=readdoubles[cur_focal];
}
}
if(sum3>0){fprintf(output2, " %.6f", sum2/sum3);}
else{fprintf(output2, " NA");}

sum=0;sum2=0;sum3=0;
for(cur_focal=0;cur_focal<num_focals;cur_focal++)
{
if(blupcentres[cur_focal][j]!=-9999)
{
sum+=blupcentres[cur_focal][j];
sum2+=blupfactors[cur_focal][j];
sum3++;
}
}
if(sum3>0){fprintf(output2, " %.6f\n", sum2/sum3);}
else{fprintf(output2, " NA");}
}
fclose(output2);

printf("Meta-analysis versions saved in %s\n\n", filename2);
}

for(j=0;j<num_preds_use;j++){free(preds[j]);}free(preds);
free(chr);free(bp);free(al1);free(al2);
for(cur_focal=0;cur_focal<num_focals;cur_focal++){free(blupcentres[cur_focal]);free(blupfactors[cur_focal]);}
free(blupcentres);free(blupfactors);
free(rs);


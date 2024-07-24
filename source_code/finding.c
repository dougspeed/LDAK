/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Finding or removing tags
//stypes=1 - normal predictor (potential tag), stypes=2 - top/target pred, stypes=3 - both tag and top/target

///////////////////////////

//get gene details and get genemax
gchr=malloc(sizeof(int)*num_seek);
gbp1=malloc(sizeof(double)*num_seek);
gbp2=malloc(sizeof(double)*num_seek);
gnames=malloc(sizeof(char*)*num_seek);
gstrand=malloc(sizeof(int)*num_seek);
gstarts=malloc(sizeof(int)*num_seek);
gends=malloc(sizeof(int)*num_seek);

for(g=0;g<num_seek;g++)
{
gchr[g]=chr[sindex[g]];
gbp1[g]=cmbp[sindex[g]]-1000*window_kb;
if(gbp1[g]<0){gbp1[g]=0;}
gbp2[g]=cmbp[sindex[g]]+1000*window_kb;
gstrand[g]=1;
copy_string(gnames,g,preds[sindex[g]]);
gstarts[g]=-9999;gends[g]=-9999;
}

cut_genes(gchr, gbp1, gbp2, gnames, gstrand, gstarts, gends, data_length, chr, preds, cmbp, weights, "blank", -9999, -9999, 0, 0, 0, 1, pvafile, pvalues, -9999, -9999, 2, NULL, NULL, NULL, -9999, num_seek);

genemax=1;
for(g=0;g<num_seek;g++)
{
if(gends[g]-gstarts[g]>genemax){genemax=gends[g]-gstarts[g];}
}
printf("The longest window contains %d predictors\n\n", genemax);

////////

//can now allocate variables
data_warn(num_samples_use, genemax);
data=malloc(sizeof(double)*num_samples_use*genemax);

cors=malloc(sizeof(double)*genemax);

if(mode==108)
{
effects=malloc(sizeof(double*)*(num_scores+1));
for(k=0;k<num_scores+1;k++){effects[k]=malloc(sizeof(double)*data_length);}
}

//fill some

if(mode==108)	//set effects to zero
{
for(k=0;k<num_scores+1;k++)
{
for(j=0;j<data_length;j++){effects[k][j]=0;}
}
}
else	//so mode=109 - set mults to 1 (so that predictors not read are not excluded)
{
for(j=0;j<data_length;j++){mults[j]=1;}
}

//prepare for reading data
if(binary==0){open_datagz(&datainputgz, datafile, num_samples, genskip, genheaders, genprobs);}
current=0;start=0;end=0;

//deal with progress and on-the-fly files
sprintf(filename,"%s.progress",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}
fclose(output);

if(mode==108)
{
sprintf(filename2,"%s.replacements",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
fprintf(output2, "Predictor\tBest_Tag\tCorrelation\n"); 
}

//ready for g loop
for(g=0;g<num_seek;g++)
{
total=gends[g]-gstarts[g];

if(g%100==0)
{
printf("Searching for tags for Predictor %d of %d\n", g+1, num_seek);
if((output=fopen(filename,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename);exit(1);}
fprintf(output, "Searching for tags for Predictor %d of %d\n", g+1, num_seek);
fclose(output);

if(mode==108)
{
fclose(output2);
if((output2=fopen(filename2,"a"))==NULL)
{printf("Error re-opening %s\n\n",filename2);exit(1);}
}
}

if(gstarts[g]!=-9999)
{
shuffle=0;
for(j=0;j<end-gstarts[g];j++)	//using values already in data, so shuffle back
{
for(i=0;i<num_samples_use;i++)
{data[(size_t)shuffle*num_samples_use+i]=data[(size_t)(gstarts[g]-start+j)*num_samples_use+i];}
shuffle++;
}

current=read_data_fly(datafile, dtype, data+(size_t)shuffle*num_samples_use, NULL, num_samples_use, keepsamps, gstarts[g]+shuffle, gends[g], keeppreds_use, datainputgz, current, num_samples, num_preds, genskip, genheaders, genprobs, bedbytes, missingvalue, -9999, -9999, nonsnp, maxthreads);
stand_data(data+(size_t)shuffle*num_samples_use, centres+gstarts[g]+shuffle, mults+gstarts[g]+shuffle, sqdevs+gstarts[g]+shuffle, num_samples_use, total-shuffle, missingvalue, -1, 0, 0, NULL, 1, preds+gstarts[g]+shuffle);

if(mode==108)	//searching for best tag
{
if(stypes[sindex[g]]==3){found=sindex[g];value=1;}	//top predictor within the data
else	//not in data
{
found=-1;value=0;
if(mults[sindex[g]]!=-9999)	//get correlations with top predictor and find best tag
{
alpha=1.0/num_samples_use;beta=0.0;
dgemv_("T", &num_samples_use, &total, &alpha, data, &num_samples_use, data+(size_t)(sindex[g]-gstarts[g])*num_samples_use, &one, &beta, cors, &one);

for(j=gstarts[g];j<gends[g];j++)
{
if((stypes[j]==1||stypes[j]==3)&&mults[j]!=-9999)	//potential tag and not trivial
{
if(fabs(cors[j-gstarts[g]])>fabs(value)){found=j;value=cors[j-gstarts[g]];}
}}
}
}

if(found!=-1&&pow(value,2)>=mincor)	//have found a (valid) tag, so load up effects, and print tag
{
for(k=0;k<num_scores;k++){effects[k][found]+=((value>0)-(value<0))*blupfactors[k][g];}
effects[num_scores][found]=1;
fprintf(output2, "%s\t%s\t%.6f\n", preds[sindex[g]], preds[found], value);
}
else	//print that there is no tag
{fprintf(output2, "%s\tNONE\tNA\n", preds[sindex[g]]);}
}
else	//so mode=109 - remove any predictor with correlation squared >= mincor with top predictor
{
alpha=1.0/num_samples_use;beta=0.0;
dgemv_("T", &num_samples_use, &total, &alpha, data, &num_samples_use, data+(size_t)(sindex[g]-gstarts[g])*num_samples_use, &one, &beta, cors, &one);

for(j=gstarts[g];j<gends[g];j++)	//consider all, including itself, already removed and trivial
{
if(pow(cors[j-gstarts[g]],2)>=mincor){mults[j]=-9999;}
}
}

start=gstarts[g];if(gends[g]>end){end=gends[g];}
}}	//end of using g and g loop
printf("\n");
if(mode==108){fclose(output2);}

if(mode==108)	//save new tags
{
sprintf(filename3,"%s.score.tags",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}
fprintf(output3,"Predictor\tA1\tA2\tCentre");
for(k=0;k<num_scores;k++){fprintf(output3, "\tEffect_%d", k+1);}
fprintf(output3,"\n");
for(j=0;j<data_length;j++)
{
if(effects[num_scores][j]==1)
{
fprintf(output3, "%s\t%c\t%c\t%.6f", preds[j], al1[j], al2[j], centres[j]);
for(k=0;k<num_scores;k++){fprintf(output3, "\t%.6f", effects[k][j]);}
fprintf(output3,"\n");
}}
fclose(output3);

printf("New score file saved in %s, with list of tags in %s\n\n", filename3, filename2);
}
else	//so mode=109 - save kept and lost predictors (ignore those with stype=2)
{
sprintf(filename2,"%s.in",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}
sprintf(filename3,"%s.out",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

count=0;count2=0;
for(j=0;j<data_length;j++)
{
if(stypes[j]!=3)
{
if(mults[j]!=-9999){fprintf(output2,"%s\n", preds[j]);count++;}
else{fprintf(output3,"%s\n", preds[j]);count2++;}
}
}

fclose(output2);
fclose(output3);

printf("Search complete: %d predictors kept (saved in %s), %d lost (%s)\n\n", count, filename2, count2, filename3);
}

//free allocations from setdl.c
free(sindex);free(stypes);
if(mode==108)
{
for(k=0;k<num_scores;k++){free(blupcentres[k]);free(blupfactors[k]);}
free(blupcentres);free(blupfactors);
}

//frees from above
for(g=0;g<num_seek;g++){free(gnames[g]);}free(gnames);
free(gchr);free(gbp1);free(gbp2);free(gstrand);free(gstarts);free(gends);

free(data);
free(cors);
if(mode==108)
{
for(k=0;k<num_scores+1;k++){free(effects[k]);}free(effects);
}
if(binary==0){gzclose(datainputgz);}

///////////////////////////


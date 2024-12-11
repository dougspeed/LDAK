/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Merge pathway tagging files

///////////////////////////

//allocations
pcounts=malloc(sizeof(int)*num_tags);
Minput=malloc(sizeof(FILE *)*num_tags);
readstrings=malloc(sizeof(char **)*11);
writestrings=malloc(sizeof(char **)*11);
readfloats=malloc(sizeof(float)*12);
writefloats=malloc(sizeof(float)*12);
for(q=0;q<11;q++){readstrings[q]=malloc(sizeof(char)*500);writestrings[q]=malloc(sizeof(char)*500);}

//check each pathway.sums file and get numbers of partitions
num_parts=0;
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.sums",tagstems[k]);

if(countcols(filename)!=13){printf("Error, %s should have 13 columns (not %d), suggesting the file has been changed since creation with \"--calc-taggings\"\n\n", filename, countcols(filename));exit(1);}

pcounts[k]=countrows(filename);	
if(pcounts[k]<1)
{printf("Error, %s should have at least one row (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, pcounts[k]);exit(1);}

num_parts+=pcounts[k];
}

//check each pathway.tagging file and get number of predictors
sprintf(filename,"%s.pathway.tagging",tagstems[0]);
count=countcols(filename)-1;
if(count<1)
{printf("Error, %s should have at least two columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", filename, 1+count);exit(1);}

for(k=1;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.tagging",tagstems[k]);
count2=countcols(filename)-1;
if(count2!=count){printf("Error, number of predictors in %s (%d) does not match number in %s.pathway.tagging (%d)\n\n", filename, count2, tagstems[0], count);exit(1);}
}

//check each base.tagging file
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.base.tagging",tagstems[k]);

if(countcols(filename)!=11){printf("Error, %s should have 11 columns (not %d), suggesting the file has been changed since creation with \"--calc-taggings\"\n\n", filename, countcols(filename));exit(1);}

count2=countrows(filename)-1;
if(count2!=count){printf("Error, number of predictors in %s (%d) does not match number in %s.pathway.tagging (%d)\n\n", filename, count2, tagstems[k], count);exit(1);}
}

//deal with labels
catlabels=malloc(sizeof(char *)*num_parts);
count2=0;
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.sums",tagstems[k]);
read_strings(filename, catlabels+count2, pcounts[k], NULL, 1, 0);
count2+=pcounts[k];
}

if(strcmp(catlabels[0],"Pathway_1")==0)	//have generic names
{
printf("It appears the pathways have generic labels (the first one is called Pathway_1)");
for(q=0;q<num_parts;q++)
{
free(catlabels[q]);
catlabels[q]=malloc(sizeof(char)*500);
sprintf(catlabels[q],"Pathway_%d",q+1);
}
}
else	//have bespoke names - will test unique
{
printf("Checking none of the %d pathway labels are duplicates\n", num_parts);
check_dups(catlabels, num_parts, "the pathway tagging files", NULL, 1);
printf("\n");
}

////////

//open all base.taggings files
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.base.tagging",tagstems[k]);
if((Minput[k]=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
}

//check they are identical, except for number of categories, and save
printf("Processing the %d base.tagging files\n", num_parts);

sprintf(filename2,"%s.base.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

if(fscanf(Minput[0], "%s %s %s %s %s %s %s %s %s %s %s ", writestrings[0], writestrings[1], writestrings[2], writestrings[3], writestrings[4], writestrings[5], writestrings[6], writestrings[7], writestrings[8], writestrings[9], writestrings[10])!=11)
{printf("Error reading Row 1 of %s.base.tagging\n\n", tagstems[0]);exit(1);}

for(k=1;k<num_tags;k++)
{
if(fscanf(Minput[k], "%s %s %s %s %s %s %s %s %s %s %s ", readstrings[0], readstrings[1], readstrings[2], readstrings[3], readstrings[4], readstrings[5], readstrings[6], readstrings[7], readstrings[8], readstrings[9], readstrings[10])!=11)
{printf("Error reading Row 1 of %s.base.tagging\n\n", tagstems[k]);exit(1);}

for(q=0;q<11;q++)
{
if(strcmp(writestrings[q], readstrings[q])!=0)
{printf("Error, %s.base.tagging and %s.base.tagging do not match for Element %d of Row 1 (%s and %s)\n\n", tagstems[0], tagstems[k], q+1, writestrings[q], readstrings[q]);exit(1);}
}
}

fprintf(output2, "%s %s %s %s %s %s %s %s %s %s %s\n", writestrings[0], writestrings[1], writestrings[2], writestrings[3], writestrings[4], writestrings[5], writestrings[6], writestrings[7], writestrings[8], writestrings[9], writestrings[10]);

for(j=0;j<count;j++)
{
if(fscanf(Minput[0], "%s %s %s %s %s %s %s %d %s %s %s ", writestrings[0], writestrings[1], writestrings[2], writestrings[3], writestrings[4], writestrings[5], writestrings[6], &writeint, writestrings[8], writestrings[9], writestrings[10])!=11)
{printf("Error reading Row %d of %s.base.tagging\n\n", 2+j, tagstems[0]);exit(1);}

for(k=1;k<num_tags;k++)
{
if(fscanf(Minput[k], "%s %s %s %s %s %s %s %d %s %s %s ", readstrings[0], readstrings[1], readstrings[2], readstrings[3], readstrings[4], readstrings[5], readstrings[6], &readint, readstrings[8], readstrings[9], readstrings[10])!=11)
{printf("Error reading Row %d of %s.base.tagging\n\n", 2+j, tagstems[k]);exit(1);}

for(q=0;q<11;q++)
{
if(q!=7&&strcmp(writestrings[q], readstrings[q])!=0)
{printf("Error, %s.base.tagging and %s.base.tagging do not match for Element %d of Row %d (%s and %s)\n\n", tagstems[0], tagstems[k], q+1, 2+j, writestrings[q], readstrings[q]);exit(1);}
}

writeint+=readint-1;
}

fprintf(output2, "%s %s %s %s %s %s %s %d %s %s %s\n", writestrings[0], writestrings[1], writestrings[2], writestrings[3], writestrings[4], writestrings[5], writestrings[6], writeint, writestrings[8], writestrings[9], writestrings[10]);
}

fclose(output2);
for(k=0;k<num_tags;k++){fclose(Minput[k]);}

////////

//concatenate the pathway.tagging files, updating the pathway names
printf("Processing the %d pathway.tagging files\n", num_parts);

sprintf(filename2,"%s.pathway.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

count2=0;
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.tagging",tagstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

for(q=0;q<pcounts[k];q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element 1 of Row %d of %s\n\n", 1+q, filename);}
fprintf(output2, "%s", catlabels[count2+q]);

for(j=0;j<count;j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 2+j, 1+q, filename);}
fprintf(output2, " %s", readstring);
}
fprintf(output2, "\n");
}

count2+=pcounts[k];
fclose(input);
}
fclose(output2);

//concatenate the pathway.sums files, updating the pathway names
printf("Processing the %d pathway.sums files\n", num_parts);

sprintf(filename2,"%s.pathway.sums",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

count2=0;
for(k=0;k<num_tags;k++)
{
sprintf(filename,"%s.pathway.sums",tagstems[k]);
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

for(q=0;q<pcounts[k];q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element 1 of Row %d of %s\n\n", 1+q, filename);}
fprintf(output2, "%s", catlabels[count2+q]);

for(j=0;j<12;j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 2+j, 1+q, filename);}
fprintf(output2, " %s", readstring);
}
fprintf(output2, "\n");
}

count2+=pcounts[k];
fclose(input);
}
fclose(output2);

printf("\nNew pathway tagging file saved with stem %s\n\n", outfile);

free(pcounts);
free(Minput);
for(q=0;q<11;q++){free(readstrings[q]);free(writestrings[q]);}
free(readstrings);free(writestrings);
free(readfloats);free(writefloats);
for(q=0;q<num_parts;q++){free(catlabels[q]);}free(catlabels);

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Reducing a tagging files

///////////////////////////

count=countrows(tagfile)-num_parts-2;

//open tagging file and output file
if((input=fopen(tagfile,"r"))==NULL)
{printf("Error opening %s\n\n",tagfile);exit(1);}
sprintf(filename2,"%s.tagging",outfile);
if((output2=fopen(filename2,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename2);exit(1);}

//read through header then tagging rows
for(j=0;j<1+count;j++)
{
if(j%100000==0){printf("Processing Row %d of %d of %s\n", j+1, count, tagfile);}

//first print out first 9 elements
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element 1 of Row %d of %s\n\n", j+2, tagfile);exit(1);}
fprintf(output2, "%s", readstring);
for(q=1;q<9;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", q+1, j+2, tagfile);exit(1);}
fprintf(output2, " %s", readstring);
}

//then the required columns
for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 10+q, j+2, tagfile);exit(1);}
if(keepparts[q]==1){fprintf(output2, " %s", readstring);}
}
fprintf(output2,"\n");
}

//now sort sums
count2=0;
for(q2=0;q2<num_parts;q2++)
{
if(keepparts[q2]==1||q2==num_parts)	//read and print (parts of) this row
{
if(fscanf(input, "The %s %s %s %s %s %s %s %s ", readstring, readstring, readstring, readstring, readstring, readstring, readstring, readstring)!=8)
{printf("Error reading Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+q2, tagfile);exit(1);}

if(parttype==0)
{
if(count2==num_reds-1){fprintf(output2, "The relative contribution of the Base to each category");}
else{fprintf(output2, "The relative contribution of Annotation %d to each category", count2+1);}
}
else{fprintf(output2, "The relative contribution of Partition %d to each category", count2+1);}

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%lf ", &readdouble)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", 10+q, count+2+q2, tagfile);exit(1);}
if(keepparts[q]==1){fprintf(output2, " %.4f", readdouble);}
}
fprintf(output2,"\n");
count2++;
}
else	//skip this row
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}
}

//and the final row
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+num_parts, tagfile);exit(1);}
if(strcmp(readstring,"The")==0)
{printf("Error, the tagging file %s was made using an older version of LDAK; please remake using this version\n\n", tagfile);exit(1);}
if(strcmp(readstring,"There")!=0)
{printf("Error reading the first element of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+num_parts, tagfile);exit(1);}

if(fscanf(input, "%s %d %s %d %s %d %s %s ", readstring, &readint, readstring, &readint2, readstring, &readint3, readstring, readstring)!=8)
{printf("Error reading first 9 elements of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", count+2+num_parts, tagfile);exit(1);}
fprintf(output2, "There are %d reference %d regression %d heritability predictors", readint, readint2, readint3);

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%d ", &readint)!=1)
{printf("Error reading Element %d of Row %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", 10+q, count+2+num_parts, tagfile);exit(1);}
if(keepparts[q]==1){fprintf(output2, " %d", readint);}
}
fprintf(output2,"\n");

fclose(input);
fclose(output2);

printf("\n");

////////

if(strcmp(matfile,"blank")!=0)	//reduce matrix
{
count=countrows(matfile)-1;

//open heritability matfile and output file
if((input=fopen(matfile,"r"))==NULL)
{printf("Error opening %s\n\n",matfile);exit(1);}
sprintf(filename3,"%s.matrix",outfile);
if((output3=fopen(filename3,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename3);exit(1);}

//read through header then tagging rows
for(j=0;j<1+count;j++)
{
if(j%100000==0){printf("Processing Row %d of %d of %s\n", j+1, count, matfile);}

//first print out first element
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element 1 of Row %d of %s\n\n", j+2, tagfile);exit(1);}
fprintf(output3, "%s", readstring);

//then the required columns
for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of Row %d of %s\n\n", 10+q, j+2, tagfile);exit(1);}
if(keepparts[q]==1){fprintf(output3, " %s", readstring);}
}
fprintf(output3,"\n");
}

fclose(input);
fclose(output3);

printf("\n");
}

printf("New tagging file saved in %s", filename2);
if(strcmp(matfile,"blank")!=0){printf(", with new heritability matrix saved in %s", filename3);}
printf("\n\n");

///////////////////////////


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Check for multiple phenotypes if mode 140 (this code makes parts of parsefiles.c redundant)

///////////////////////////

sprintf(filename,"%s.root", bocofile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--kvik-step3\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=10)
{printf("Error, %s should have ten rows (not %d), suggesting the file has been changed since creation with \"--kvik-step1\"\n\n", filename, count);exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

for(j=0;j<4;j++)	//skip four rows
{
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"One")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, Row 5 of %s should say \"One Phenotype\" or \"Multiple Phenotypes\" (not %s %s), suggesting the file has been changed since creation with \"--kvik-step1\"\n\n",filename, readstring, readstring2);exit(1);}

if(strcmp(readstring,"Multiple")==0)
{
printf("Error, you can not use \"--kvik-step3\" for multi-phenotype analyses; you should instead use \"--cut-genes\", \"--calc-genes-reml\" and \"--join-genes-reml\" for each phenotype in turn\n\n");exit(1);
}

fclose(input);

///////////////////////////


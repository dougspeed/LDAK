/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Check and read details from the prs root file - here for modes 131, 132, 134 and 140 (134 becomes either 131 or 132)

///////////////////////////

if(checkprs==-9999){checkprs=1;}

//check and open root file

sprintf(filename,"%s.root", prsfile);
if(just_check(filename)!=0)
{
if(mode==131||mode==132){printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\"\n\n", filename);}
else{printf("Error reading %s; this file should have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countrows(filename);
if(count!=11)
{printf("Error, %s should have 11 rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read analysis type (think no need to separate fastgwa 1 and 2)

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Analysis")!=0)
{printf("Error reading %s; Row 1 should begin \"Analysis\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

fastgwa=-1;
if(strcmp(readstring2,"KVIK")==0){fastgwa=0;}
if(strcmp(readstring2,"fastGWA")==0){fastgwa=1;}
if(strcmp(readstring2,"KVIK_Pedigree")==0){fastgwa=2;}

if(fastgwa==-1){printf("Error reading %s; Element 2 of Row 1 should be \"KVIK\", \"KVIK_Pedigree\" or \"fastGWA\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring2);exit(1);}

if(strcmp(sampwfile,"blank")!=0&&fastgwa!=1)
{printf("Error, you can only use \"--sample-weights\" if the Step 1 PRS were constructed using \"--fastGWA YES\"\n\n");exit(1);}

//read regression model

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 2 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Regression")!=0)
{printf("Error reading %s; Row 2 should begin \"Regression\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(mode==131)	//check consistent
{
if(strcmp(readstring2,"Logistic")==0){printf("Error, \"--binary YES\" was used when making the PRS, so you should now use \"--logistic\" (instead of \"--linear\")\n\n");exit(1);}
}
if(mode==132)	//check consistent
{
if(strcmp(readstring2,"Linear")==0){printf("Error, \"--binary YES\" was not used when making the PRS, so you should now use \"--linear\" (instead of \"--logistic\")\n\n");exit(1);}
}
if(mode==134)	//switch to 131 or 132
{
if(strcmp(readstring2,"Linear")==0){mode=131;}
else{mode=132;}
}

//read data, phenotypes, covariates and factors

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0)
{printf("Error reading %s; Row 3 should begin \"Datafile\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring3)!=2)
{printf("Error reading Row 4 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Phenotypes")!=0)
{printf("Error reading %s; Row 4 should begin \"Phenotypes\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring4)!=2)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Covariates")!=0)
{printf("Error reading %s; Row 5 should begin \"Covariates\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring5)!=2)
{printf("Error reading Row 6 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Factors")!=0)
{printf("Error reading %s; Row 6 should begin \"Factors\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(mode==131||mode==132)	//check covariates are consistent
{
if(strcmp(covarfile,"blank")==0&&strcmp(readstring4,"none")!=0)
{printf("Error, \"--covar\" was used when making the PRS, but is not used now\n\n");exit(1);}
if(strcmp(covarfile,"blank")!=0&&strcmp(readstring4,"none")==0)
{printf("Error, \"--covar\" is used now, but was not used when making the PRS\n\n");exit(1);}

if(strcmp(factorfile,"blank")==0&&strcmp(readstring5,"none")!=0)
{printf("Error, \"--factors\" was used when making the PRS, but is not used now\n\n");exit(1);}
if(strcmp(factorfile,"blank")!=0&&strcmp(readstring5,"none")==0)
{printf("Error, \"--factors\" is used now, but was not used when making the PRS\n\n");exit(1);}
}

if(checkroot==1&&(mode==131||mode==132))	//check names consistent
{
flag=0;

if(strcmp(respfile,readstring3)!=0)
{printf("Error, the phenotypes file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring3, respfile);flag=1;}
if(strcmp(covarfile,"blank")!=0&&strcmp(covarfile,readstring4)!=0)
{printf("Error, the covariates file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring4, covarfile);flag=1;}
if(strcmp(factorfile,"blank")!=0&&strcmp(factorfile,readstring5)!=0)
{printf("Error, the factors file used when making the PRS (%s) does not match that used now (%s)\n\n", readstring5, factorfile);flag=1;}

if(flag==1){printf("If you are sure you have provided the correct files, use \"--check-root NO\"\n\n");exit(1);}
}

//read phenotype number

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 7 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Phenotype_Number")!=0&&strcmp(readstring,"Num_Phenotypes")!=0)
{printf("Error reading %s; Row 7 should begin \"Phenotype_Number\" or \"Num_Phenotypes\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n",filename, readstring);exit(1);}

if(strcmp(readstring,"Phenotype_Number")==0)	//set kvikparity to zero and check mpheno
{
kvikparity=0;

if(mpheno==-1)
{printf("Error, \"--mpheno ALL\" is used now, but was not used when making the PRS\n\n");exit(1);}
if(mpheno!=-9999&&mpheno!=atoi(readstring2))
{printf("Error, the phenotype number provided now (%d) does not match that used with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" (%d)\n\n", mpheno, atoi(readstring2));exit(1);}
}
else	//set kvikparity to the number of phenotypes and maybe check mpheno
{
kvikparity=atoi(readstring2);

if(mpheno!=-1)
{printf("Error, \"--mpheno ALL\" was used when making the PRS, but is not used now\n\n");exit(1);}
}

//read num_samples_used (save in num_samples_use2)
if(fscanf(input, "%s %d ", readstring, &num_samples_use2)!=2)
{printf("Error reading Row 8 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Num_Samples_Used")!=0)
{printf("Error reading %s; Row 8 should begin \"Num_Samples_Used\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}

//skip one row
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read num_chr (save in num_chr2)
if(fscanf(input, "%s %d ", readstring, &num_chr2)!=2)
{printf("Error reading Row 10 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Num_Chromosomes")!=0)
{printf("Error reading %s; Row 10 should begin \"Num_Chromosomes\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}

//read multivariate
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 11 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Multivariate")!=0)
{printf("Error reading %s; Row 11 should begin \"Multivariate\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}

if((multi==-9999||multi==0)&&strcmp(readstring2,"YES")==0)
{printf("Error, \"--multivariate YES\" was used when making the PRS, but is not used now\n\n");exit(1);}
if(multi==1&&strcmp(readstring2,"NO")==0)
{printf("Error, \"--multivariate YES\" is used now, but was not used when making the PRS\n\n");exit(1);}

fclose(input);

//check size of prs file or files
if(mpheno!=-1)	//only one file
{
sprintf(filename,"%s.loco.prs", prsfile);
if(just_check(filename)!=0)
{
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999){printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);}
else{printf("Error reading %s; this file should have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countrows(filename)-1;
if(count!=num_samples_use2)
{printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, num_samples_use2+1, count+1);exit(1);}

count=countcols(filename);
if(fastgwa==0&&count!=2+num_chr2){printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, 2+num_chr2, count);exit(1);}
if(fastgwa!=0&&count!=3){printf("Error, %s should have three columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}
}
else	//kvikparity files
{
for(m=0;m<kvikparity;m++)
{
sprintf(filename,"%s.pheno%d.loco.prs", prsfile, m+1);
if(just_check(filename)!=0)
{
if(kvikstep==-9999&&gctastep==-9999&&faststep==-9999){printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename);}
else{printf("Error reading %s; this file should have been created using \"--kvik-step1\", \"--GCTA-LOCO-step1\" or \"--fastGWA-step1\"\n\n", filename);}
exit(1);
}

count=countrows(filename)-1;
if(count!=num_samples_use2)
{printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, num_samples_use2+1, count+1);exit(1);}

count=countcols(filename);
if(fastgwa==0&&count!=2+num_chr2){printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, 2+num_chr2, count);exit(1);}
if(fastgwa!=0&&count!=3){printf("Error, %s should have three columns (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}
}
}

if(multi==1)	//check multivariate files exist (must be mode 131 or 132)
{
sprintf(filename,"%s.multivariate.transformation", prsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--multivariate YES\"\n\n", filename);exit(1);}

count=countrows(filename);
count2=countcols(filename);
if(count!=kvikparity||count2!=kvikparity)
{printf("Error, %s should have %d rows and %d columns (not %d and %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, kvikparity, kvikparity, count, count2);exit(1);}

sprintf(filename,"%s.multivariate.inverse", prsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--multivariate YES\"\n\n", filename);exit(1);}

count=countrows(filename);
count2=countcols(filename);
if(count!=kvikparity||count2!=kvikparity)
{printf("Error, %s should have %d rows and %d columns (not %d and %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, kvikparity, kvikparity, count, count2);exit(1);}
}

////////

if(mode==140)	//get power (stored in Row 5 of details file) and set pvafile and sumsfile (bit complicated when multiple phenotypes)
{
if(kvikparity==0)	//previously analyzed just one phenotype (so only one pvafile and one sumsfile)
{
sprintf(filename,"%s.loco.details", prsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=7)
{printf("Error, %s should have seven rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open and skip five rows
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read power and close file
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 6 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Power")!=0)
{printf("Error reading %s; Row 6 should begin \"Power\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
if(strcmp(readstring2,"NA")==0)
{printf("Error, the Step 1 PRS were constructed using \"--fastGWA YES\" (you can only use \"--kvik-step3\" if you have previously run both \"--kvik-step1\" and \"--kvik-step2\")\n\n");exit(1);}
power=atof(readstring2);
fclose(input);

//set pvafile and sumsfile
sprintf(pvafile,"%s.pvalues",assfile);
sprintf(sumsfile,"%s.summaries",assfile);

//check they exist
if(just_check(pvafile)!=0)
{printf("Error reading %s; this file should have been created using \"--linear\" or \"--logistic\" with \"--PRS\"\n\n", pvafile);exit(1);}
if(just_check(sumsfile)!=0)
{printf("Error reading %s; this file should have been created using \"--linear\" or \"--logistic\" with \"--PRS\"\n\n", sumsfile);exit(1);}
}
else	//previously analyzed multiple phenotypes
{
powers=malloc(sizeof(double)*kvikparity);

for(m=0;m<kvikparity;m++)
{
sprintf(filename,"%s.pheno%d.loco.details", prsfile, m+1);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" with \"--LOCO YES\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=7)
{printf("Error, %s should have seven rows (not %d), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, count);exit(1);}

//open and skip five rows
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}

//read power and close file
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 6 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Power")!=0)
{printf("Error reading %s; Row 6 should begin \"Power\" (not %s), suggesting the file has been changed since creation with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\"\n\n", filename, readstring);exit(1);}
if(strcmp(readstring2,"NA")==0)
{printf("Error, the Step 1 PRS were constructed using \"--fastGWA YES\" (you can only use \"--kvik-step3\" if you have previously run both \"--kvik-step1\" and \"--kvik-step2\")\n\n");exit(1);}
powers[m]=atof(readstring2);
fclose(input);

//check pvafile exists
sprintf(filename,"%s.pheno%d.pvalues",assfile, m+1);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--linear\" or \"--logistic\" with \"--PRS\"\n\n", filename);exit(1);}

//check sumsfile exists
sprintf(filename,"%s.pheno%d.summaries",assfile, m+1);
if(just_check(filename)!=0)
{printf("Error reading %s; this file should have been created using \"--linear\" or \"--logistic\" with \"--PRS\"\n\n", filename);exit(1);}
}

//set power, pvafile and sumsfile corresponding to the first phenotype
power=powers[0];
sprintf(pvafile,"%s.pheno1.pvalues",assfile);
sprintf(sumsfile,"%s.pheno1.summaries",assfile);
}
}

///////////////////////////


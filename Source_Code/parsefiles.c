/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Process and check files

///////////////////////////

//data - convenient to set use_data, which states whether using datafiles
//0 - no data, 1 - "normal", 2 - only through regions, 3 - only through top preds
//4 - have data, but only using ids, 5 - multiple datafiles, 6 - have bim (not using ids)
use_data=0;

if(mode==101||mode==102||mode==103||mode==104||mode==105||mode==106||mode==107||mode==108||mode==109||mode==110||mode==111||mode==112||mode==114||(mode==117&&extract==1)||mode==122||mode==127||mode==128||mode==131||mode==132||mode==136||mode==137||mode==138||(mode==139&&cut1!=-9999)||mode==140||mode==141||mode==145||mode==151||mode==152||mode==153||mode==154||mode==156||mode==158||mode==160||mode==162||mode==171||mode==172||mode==173||mode==175||mode==186||mode==187||mode==188||mode==189||mode==191||mode==192||mode==193||mode==194)
{use_data=1;}

if(use_data==0&&num_regs>0)
{use_data=2;}

if(use_data==0&&strcmp(topfile,"blank")!=0)
{use_data=3;}

if(mode==177&&dtype!=-9999)
{use_data=4;}

if(mode==181||mode==182||mode==183||mode==184||mode==185||mode==190)
{
if(strcmp(datalist,"blank")==0){use_data=1;}
else{use_data=5;}
}

if(mode==159)
{use_data=6;}

if((use_data==1||use_data==2||use_data==3)&&dtype==-9999)
{
if(mode==151||mode==152||mode==153||mode==154)
{printf("Error, you must provide a set of genetic data files using \"--bfile\", \"--bgen\" or \"--speed\"\n\n");exit(1);}

printf("Error, you must provide a set of genetic data files using \"--bfile\", \"--bgen\", \"--sp\", \"--sped\", \"--speed\" or \"--gen\"\n\n");exit(1);
}

if((mode==181||mode==182||mode==183||mode==184||mode==185)&&dtype==-9999)
{printf("Error, you must provide one or more sets of genetic data files using \"--bfile\", \"--bgen\", \"--sp\", \"--sped\", \"--speed\", \"--gen\", \"--mbfile\", \"--msp\", \"--msped\" or \"--mgen\"\n\n");exit(1);}

if(mode==190&&dtype==-9999)
{printf("Error, you must provide two sets of genetic data files using \"--mbfile\", \"--msp\", \"--msped\" or \"--mspeed\"\n\n");exit(1);}

if(use_data==0&&dtype!=-9999)
{
if(strcmp(weightsfile,"blank")!=0)
{printf("Warning, genetic data files and weightings are provided but will not be used\n\n");}
else
{printf("Warning, genetic data files are provided but will not be used\n\n");}
dtype=-9999;strcpy(weightsfile,"blank");
}

if((use_data==3||use_data==4)&&strcmp(weightsfile,"blank")!=0)
{printf("Warning, weightings are provided but will not be used\n\n");strcpy(weightsfile,"blank");}

if(extract==1)
{
if(strcmp(topfile,"blank")!=0)
{printf("Warning, the predictor filterings do not apply to top predictors\n\n");}
if(strcmp(targetfile,"blank")!=0)
{printf("Warning, the predictor filterings do not apply to target predictors\n\n");}
}

////////

if(famhead==-9999)	//do not have a separate famfile
{
if(dtype==1||dtype==3||dtype==4){famhead=0;}

if(dtype==5)	//then must be providing a list (else famhead would be set)
{
count=countcols(datalist);
checkcols(datalist,count);
if(count>1){printf("Warning, it is assumed that Column %d of %s provides sample files; if these are instead fam files, please use \"--fam-files YES\"\n\n", count, datalist);}
famhead=2;
}
}

if(use_data==1||use_data==2||use_data==3||use_data==4||use_data==5)	//get stems - dtype must be 1-5 (dtypes 11-15 will have become dtype 5)
{
if(strcmp(datalist,"blank")==0)	//single dataset
{
num_files=1;
datastems=malloc(sizeof(char *));datastems[0]=malloc(sizeof(char)*500);
bimstems=malloc(sizeof(char *));bimstems[0]=malloc(sizeof(char)*500);
famstems=malloc(sizeof(char *));famstems[0]=malloc(sizeof(char)*500);

if(dtype==1||dtype==3||dtype==4){fill_names(datastems, bimstems, famstems, 0, udatafile, dtype);}
if(dtype==2)
{
strcpy(datastems[0],udatafile);
strcpy(bimstems[0],udatafile);
if(strcmp(ufamfile,"blank")!=0){strcpy(famstems[0],ufamfile);}
else{strcpy(famstems[0],udatafile);}
}
if(dtype==5)
{
strcpy(datastems[0],udatafile);
if(strcmp(ubimfile,"blank")!=0){strcpy(bimstems[0],ubimfile);}
else{strcpy(bimstems[0],udatafile);}
if(strcmp(ufamfile,"blank")!=0){strcpy(famstems[0],ufamfile);}
else{strcpy(famstems[0],udatafile);}
}
}
else	//reading from list (can not be dtype=2)
{
num_files=countrows(datalist);

count=countcols(datalist);
checkcols(datalist,count);

if(dtype==1&&count!=1&&count!=3)
{printf("Error, each row of %s must contain either one or three elements (not %d), providing either the data stem or bed, bim and fam files\n\n",datalist,count);exit(1);}
if(dtype==3&&count!=1&&count!=3)
{printf("Error, each row of %s must contain either one or three elements (not %d), providing either the data stem or sped, bim and fam files\n\n",datalist,count);exit(1);}
if(dtype==4&&count!=1&&count!=3)
{printf("Error, each row of %s must contain either one or three elements (not %d), providing either the data stem or speed, bim and fam files\n\n",datalist,count);exit(1);}

if(dtype==5)
{
if(genheaders<4&&count!=3)
{printf("Error, each row of %s must contain three elements (not %d), providing gen, bim and sample/fam files\n\n", datalist, count);exit(1);}

if(count!=2&&count!=3)
{printf("Error, each row of %s must contain two or three elements (not %d), providing gen and sample/fam files or gen, bim and sample/fam files\n\n", datalist, count);exit(1);}

if(count==2&&oxchr==-9999)
{printf("Error, either %s should include bim files, or you should use \"--oxford-single-chr\"\n\n", datalist);exit(1);}
}

datastems=malloc(sizeof(char *)*num_files);
bimstems=malloc(sizeof(char *)*num_files);
famstems=malloc(sizeof(char *)*num_files);

if((input=fopen(datalist,"r"))==NULL)
{printf("Error opening %s\n\n",datalist);exit(1);}

for(k=0;k<num_files;k++)
{
datastems[k]=malloc(sizeof(char)*500);
bimstems[k]=malloc(sizeof(char)*500);
famstems[k]=malloc(sizeof(char)*500);

if(count==1)	//must be dtypes 1, 3, or 4
{
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading stem for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
fill_names(datastems, bimstems, famstems, k, readstring, dtype);
}
if(count==2)	//must be dtype 5
{
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading gen file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(datastems[k],readstring);
strcpy(bimstems[k],datastems[k]);

if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading sample/fam file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(famstems[k],readstring);
}
if(count==3)	//must be dtypes 1, 3, 4 or 5
{
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading predictor file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(datastems[k],readstring);
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading bim file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(bimstems[k],readstring);
if(fscanf(input, "%s ", readstring2)!=1)
{printf("Error reading fam/sample file for Dataset %d from %s\n\n", k+1, datalist);exit(1);}
(void)append_check(readstring,readstring2,workdir);
strcpy(famstems[k],readstring);
}

for(k2=0;k2<k;k2++)	//check predictor file names different to previous
{
if(strcmp(datastems[k],datastems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", datalist, datastems[k]);exit(1);}
}
}	//end of k loop

fclose(input);
}	//end of using datalist
}	//end of use_data!=0

////////

if(mode==190&&num_files!=2)
{printf("Error, when using \"--calc-sim-data\", you must provide a pair of genetic data files using \"--mbfile\", \"--msped\", \"--mspeed\" or \"--mgen\"\n\n");exit(1);}

if(mode==181||mode==182||mode==183||mode==184||mode==185)	//check names different to outfile
{
if(mode==181){sprintf(filename,"%s.bed",outfile);}
if(mode==182){sprintf(filename,"%s.sp",outfile);}
if(mode==183){sprintf(filename,"%s.sped",outfile);}
if(mode==184){sprintf(filename,"%s.speed",outfile);}
if(mode==185){sprintf(filename,"%s.gen",outfile);}
for(k=0;k<num_files;k++)
{
if(strcmp(datastems[k],filename)==0)
{
if(num_files==1)
{printf("Error, the name of the new predictor file (%s) matches that provided\n\n", filename);exit(1);}
else
{printf("Error, the name of the new predictor file (%s) matches that of Dataset %d\n\n", filename, k+1);exit(1);}
}}
}

////////

//check datafiles exist (and are not folders)

for(k=0;k<num_files;k++)
{
if(just_check(datastems[k])!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", datastems[k]);exit(1);}
if(just_check(bimstems[k])!=0)
{printf("Error reading bim file %s; check file exists and is not empty\n\n", bimstems[k]);exit(1);}
if(just_check(famstems[k])!=0)
{printf("Error reading fam/sample file %s; check file exists and is not empty\n\n", famstems[k]);exit(1);}
}

//set datafile, bimfile and famfile (used when there is one set of data files)
if(use_data==0||use_data==5)
{strcpy(datafile,"blank");strcpy(bimfile,"blank");strcpy(famfile,"blank");}
if(use_data==1||use_data==2||use_data==3)
{strcpy(datafile,datastems[0]);strcpy(bimfile,bimstems[0]);strcpy(famfile,famstems[0]);}
if(use_data==4)
{strcpy(datafile,"blank");strcpy(bimfile,"blank");strcpy(famfile,famstems[0]);}
if(use_data==6)	//do not have fam file - currently must be mode 159
{strcpy(datafile,"blank");sprintf(bimfile,"%s.cors.bim", corname);strcpy(famfile,"blank");}

////////

if(dtype==5)	//give some advice
{
if(num_files==1){printf("Expect the predictor file to have ");}
else{printf("Expect each predictor file to have ");}
if(genskip==0){printf("no header row, ");}
if(genskip==1){printf("one header row, ");}
if(genskip>1){printf("%d header rows, ", genskip);}
if(genheaders==0){printf("then for each row to have no header columns ");}
if(genheaders==1){printf("then for each row to have one header column (which is ignored), ");}
if(genheaders==2){printf("then for each row to have two header columns (which are ignored), ");}
if(genheaders==3){printf("then for each row to have three header columns (Predictor, A1, A2), ");}
if(genheaders==4){printf("then for each row to have four header columns (Predictor, BP, A1, A2), ");}
if(genheaders==5){printf("then for each row to have five header columns (IGNORED, Predictor, BP, A1, A2), ");}
if(genheaders==6){printf("then for each row to have six header columns (IGNORED, IGNORED, Predictor, BP, A1, A2), ");}
if(genprobs==0){printf("and contain two binary values (\"0 0\", \"1 0\", \"0 1\" or \"1 1\") per sample");}
if(genprobs==1&&nonsnp==0){printf("and contain one dosage value (expected allele count) per sample");}
if(genprobs==1&&nonsnp==1){printf("and contain one value per sample");}
if(genprobs==2){printf("and contain two probabilities (for A1A1 and A1A2) for each sample");}
if(genprobs==3){printf("and contain three probabilities (for A1A1, A1A2 and A2A2) for each sample");}
if(genprobs==4){printf("and contain four probabilities (for A1A1, A1A2, A2A2 and missing) for each sample; ");}
printf("; these can be changed using \"--gen-skip\", \"--gen-headers\" and \"--gen-probs\"");
if(num_files==1){printf(" (if a bim file is provided, predictor details will be read from this instead)\n\n");}
else{printf(" (if bim files are provided, predictor details will be read from these instead)\n\n");}
}

if((mode==151||mode==152||mode==153||mode==154)&&dtype==4)
{
if((input=fopen(datafile,"rb"))==NULL)
{printf("Error opening %s\n\n",datafile);exit(1);}
if(read_speed_size(datafile)==2)
{printf("Error, you can not use long SPEED format with \"--ridge\", \"--bolt\", \"--bayesr\" or \"--elastic\" (when you make %s, you should not add \"--speed-long YES\")\n\n", datafile);exit(1);}
fclose(input);
}

///////////////////////////

//data filtering

if(num_subs>0)	//check exist
{
for(s=0;s<num_subs;s++)
{
sprintf(filename,"%s%d",subpref, s+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}
}

///////////////////////////

//data scaling (and pvalues) and coding

if(strcmp(centresfile,"blank")!=0)	//check correct size
{
count=countcols(centresfile);
if(count!=4)
{printf("Error, %s should have four columns, providing the predictor names, the A1 and A2 alleles, then the mean values for the A1 alleles (not %d)\n\n", centresfile, count);exit(1);}
checkcols(centresfile,count);
}
	
if(use_data==1||use_data==2||use_data==6)	//will use the weights vector
{
if(strcmp(weightsfile,"blank")==0&&ignoreweights!=1)	//weights not provided
{
if(mode==105||mode==112||mode==114||mode==127||mode==128||mode==136||mode==138||mode==140||mode==141||mode==145||mode==151||mode==152||mode==153||mode==154||mode==159||mode==173||mode==186||mode==187||mode==188||mode==189)	//weights are used - but fine to set to one
{
//printf("Warning, the predictor weightings have been set to one (equivalent to adding \"--ignore-weights YES\"); if you wish to specify different weightings, use \"--weights\"\n\n");
ignoreweights=1;
}
else	//weightings are not used, so set to one
{ignoreweights=1;}
}
}

if(strcmp(weightsfile,"blank")!=0)	//check correct size
{
count=countcols(weightsfile);
if(count!=2&&count!=6)
{printf("Error, %s should either have two or six columns (not %d), with Columns 1 and 2 providing predictor names then weightings\n\n", weightsfile, count);exit(1);}
checkcols(weightsfile,count);

if(count==6)
{
if((input=fopen(weightsfile,"r"))==NULL)
{printf("Error opening %s\n\n",weightsfile);exit(1);}
if(fscanf(input, "Predictor %s ", readstring)!=1)
{printf("Error reading Row 1 of %s; should begin \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--join-weights\" or \"--calc-weights-all\")\"\n\n", weightsfile, readstring);exit(1);}
fclose(input);
}
}

if(strcmp(pvafile,"blank")!=0)	//check correct size
{
count=countcols(pvafile);
if(count!=2)
{printf("Error, %s should have two columns (not %d), providing Predictor names then p-values\n\n", pvafile, count);exit(1);}
checkcols(pvafile,count);
}

if(strcmp(impfile,"blank")!=0)	//check correct size
{
count=countcols(impfile);
if(count!=2)
{printf("Error, %s should have two columns (not %d), providing Predictor names then p-values\n\n", impfile, count);exit(1);}
checkcols(impfile,count);
}

if(power==-9999)
{
if(mode==112||mode==114||mode==127||mode==128||mode==137||mode==138||mode==141||mode==145||mode==159||mode==173||num_regs>0)
{
if(hwestand==1){printf("Error, you must use \"--power\" to specify the assumed relationship between heritability and allele frequency; for human data we recommend -0.25\n\n");exit(1);}
else{printf("Error, you must use \"--power\" to specify the assumed relationship between heritability and predictor variance; for human data we recommend -0.25\n\n");exit(1);}
}

if(mode==172)	//set to zero, then will warn later
{power=0;}

if(mode==186||mode==187||mode==188||mode==189)
{printf("Error, you must use \"--power\" to specify how to standardize predictors; for no standardization (i.e., to simply sum allele counts) use \"--power 0\"\n\n");exit(1);}
}
else
{
if(mode!=112&&mode!=114&&mode!=127&&mode!=128&&mode!=137&&mode!=138&&mode!=140&&mode!=141&&mode!=145&&mode!=151&&mode!=152&&mode!=153&&mode!=154&&mode!=159&&mode!=172&&mode!=173&&mode!=186&&mode!=187&&mode!=188&&mode!=189&&num_regs==0)
{printf("Warning, \"--power\" was provided but will not be used\n\n");power=-9999;}
}

///////////////////////////

//kinships - set num_kins and kinship stems, then do some checks

if(strcmp(kinname,"blank")==0&&strcmp(kinlist,"blank")==0){num_kins=0;}
else
{
if(strcmp(kinname,"blank")!=0)	//will have already been prefixed by folder
{
num_kins=1;
kinstems=malloc(sizeof(char *));
kinstems[0]=malloc(sizeof(char)*500);
strcpy(kinstems[0], kinname);
}
else	//read from file
{
num_kins=countrows(kinlist);
kinstems=malloc(sizeof(char *)*num_kins);

if((input=fopen(kinlist,"r"))==NULL)
{printf("Error opening %s\n\n",kinlist);exit(1);}
for(k=0;k<num_kins;k++)
{
kinstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, kinlist);exit(1);}
(void)append_check(kinstems[k],filename,workdir);

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(kinstems[k],kinstems[k2])==0)
{printf("Error reading %s; the stem %s appears twice\n\n", kinlist, kinstems[k]);exit(1);}
}
}
fclose(input);
}	//end of reading from file
}

////////

if((mode==115||mode==118||mode==119||mode==161||mode==162||mode==163||mode==164||mode==166||mode==167||mode==168||mode==169||mode==170||mode==201||mode==202||mode==203)&&num_kins>1)
{printf("Error, you can only provide one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==117&&extract==1&&num_kins>1)
{printf("Error, when using \"--extract\" or \"--exclude\" to subtract Predictors from a kinship matrix, you can only provide one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==120&&num_kins==1)
{printf("Error, to use \"--calc-sim-grms\" you must provide two or more kinship matrices\n\n");exit(1);}

if(mode==123&&(minrel!=-9999||maxrel!=-9999)&&num_kins!=1)
{printf("Error, you can only use \"--min-rel\" or \"--max-rel\" with \"--pcgc\" if using one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==124&&(minrel!=-9999||maxrel!=-9999)&&num_kins!=1)
{printf("Error, you can only use \"--min-rel\" or \"--max-rel\" with \"--he\" if using one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==131&&num_kins>1)
{printf("Error, with \"--linear\" you can provide at most one kinship matrix (not %d); for complex mixed-model linear regression you should instead perform the analysis in two steps: first use\"--solve-null\" to construct a combined kinship matrix, then use this matrix with \"--linear\"\n\n", num_kins);exit(1);}

if(mode==133&&num_kins==1)
{printf("Error, to use \"--solve-null\" you must provide multiple kinship matrices (otherwise, use \"--linear\")\n\n");exit(1);}

if(mode==138&&num_kins>1)
{printf("Error, with \"--calc-genes-reml\" you can provide at most one kinship matrix (not %d); for complex gene/chunk-based REML you should instead perform the analysis in two steps: first use\"--solve-null\" to construct a combined kinship matrix, then use this matrix with \"--calc-genes-reml\"\n\n", num_kins);exit(1);}

if(strcmp(eigenfile,"blank")!=0&&num_kins>1)
{printf("Error, you can only use \"--eigen\" if providing one kinship matrix (not %d)\n\n", num_kins);exit(1);}

if(mode==116||mode==117||mode==118||mode==119||mode==133||mode==164||mode==167||mode==168)	//check names different to outfile
{
for(k=0;k<num_kins;k++)
{
if(strcmp(kinstems[k],outfile)==0)
{
if(num_kins==1){printf("Error, the stem of the new kinship matrix (%s) matches that provided\n\n", outfile);exit(1);}
else{printf("Error, the stem of the new kinship matrix (%s) matches that of Kinship %d\n\n", outfile, k+1);exit(1);}
}
}
}

////////

if(num_kins>0)	//check required kinship files exists and set kinnums and kinsums
{
kinnums=malloc(sizeof(int)*num_kins);
kinsums=malloc(sizeof(double)*num_kins);

flag=0;	//records how many low-rank kinship matrices
for(k=0;k<num_kins;k++)
{
sprintf(filename,"%s.grm.id",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when calculating kinships\n\n", filename);exit(1);}

if(mode!=118&&mode!=119)	//not converting
{
sprintf(filename,"%s.grm.bin",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when calculating kinships\n\n", filename);exit(1);}
}

if(mode==118)	//convert raw
{
sprintf(filename,"%s.grm.gz",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}

if(mode==119)	//convert raw
{
sprintf(filename,"%s.grm.raw",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}

if(kindetails==1&&mode!=118&&mode!=119)	//check details and adjust exist, then read things from adjust
{
sprintf(filename,"%s.grm.details",kinstems[k]);
if(just_check(filename)!=0)
{
printf("Error reading %s; this file would have been created when calculating kinships; ", filename);
if(mode==121||mode==126)
{printf("you can ignore this warning by adding \"--kinship-details NO\" (however this file will be required to subsequently calculate BLUPs)");}
if(mode!=113&&mode!=116&&mode!=117&&mode!=121&&mode!=122&&mode!=126&&mode!=162)
{printf("you can ignore this warning by adding \"--kinship-details NO\"");}
printf("\n\n");exit(1);
}

if(countcols(filename)==7)
{
printf("Error, %s has seven columns (not nine), indicating it was created using an old version of LDAK; you should remake the corresponding kinship matrix using the latest version of LDAK", filename);
if(mode==121||mode==126)
{printf(", or ignore this warning by adding \"--kinship-details NO\" (however this file will be required to subsequently calculate BLUPs)");}
if(mode!=113&&mode!=116&&mode!=117&&mode!=121&&mode!=122&&mode!=126&&mode!=162)
{printf(", or ignore this warning by adding \"--kinship-details NO\"");}
printf("\n\n");exit(1);
}

sprintf(filename,"%s.grm.adjust",kinstems[k]);
if(just_check(filename)!=0)
{
printf("Error reading %s; this file would have been created when calculating kinships; ", filename);
if(mode==121||mode==126)
{printf("you can ignore this warning by adding \"--kinship-details NO\" (however this file will be required to subsequently calculate BLUPs)");}
if(mode!=113&&mode!=116&&mode!=117&&mode!=121&&mode!=122&&mode!=126&&mode!=162)
{printf("you can ignore this warning by adding \"--kinship-details NO\"");}
printf("\n\n");exit(1);
}

//read top line
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, %s should begin \"Datafile\" or \"Multiple\" (not %s), suggesting it was created using an old version of LDAK\n\n", filename, readstring);exit(1);}

if((mode==122||mode==162)&&checkroot==1)	//check kinship matrix corresponds to the datafile
{
if(strcmp(readstring,"Multiple")==0)
{printf("Error, the kinship matrix %s uses predictors from multiple datasets; you must ensure all the predictors are contained within the data file %s, then add \"--check-root NO\"\n\n", kinstems[k], datafile);exit(1);}
(void)append_check(readstring3,readstring2,workdir);
if(strcmp(datafile,readstring2)!=0&&strcmp(datafile,readstring3)!=0)
{printf("Error, the kinship matrix %s corresponds to the data file %s, which appears to be different to that provided now (%s); if you are sure the data file is correct, add \"--check-root NO\"\n\n", kinstems[k], readstring2, datafile);exit(1);}
}

//skip power line, read kinnums, skip weights line, read kinsums
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}
if(fscanf(input, "Num_Preds %d ", kinnums+k)!=1)
{printf("Error reading Row 3 of %s\n\n", filename);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input,"%c", &readchar);}
if(fscanf(input, "Denominator %lf ", kinsums+k)!=1)
{printf("Error reading Row 5 of %s\n\n", filename);exit(1);}
fclose(input);

if(kinnums[k]<num_samples_use){flag++;}
}	//end of kindetails=1
else
{kinnums[k]=-9999;kinsums[k]=-9999;}

if(diagonal==1)
{
sprintf(filename,"%s.grm.diag",kinstems[k]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when calculating kinships using \"--gxemm-free\"\n\n", filename);exit(1);}
}
}	//end of k loop

if((mode==121||mode==123||mode==124)&&flag>0)
{
if(num_kins==1){printf("Warning, the kinship matrix was computed from fewer predictors (%d) than the number of samples (%d); it might be more efficient to replace this with a region\n\n", kinnums[0], num_samples_use);}
else{printf("Warning, %d of the %d kinship matrices were computed from fewer predictors than the number of samples (%d); it might be more efficient to replace these with regions\n\n", flag, num_kins, num_samples_use);}
}
}	//end of num_kins>0

///////////////////////////

//regions, responses, summaries and fixed

if(num_regs>0)	//check regions exist and contain elements, and set rprune
{
for(r=0;r<num_regs;r++)
{
sprintf(filename,"%s%d",regpref, r+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}

if(rprune==-9999)
{
if(strcmp(weightsfile,"blank")!=0)	//with weights, no pruning should be required
{rprune=2;}
else	//without weights, pruning depends on whether using summaries or not
{
if(strcmp(sumsfile,"blank")!=0){rprune=0.5;}
else{rprune=0.98;}
}
}
}

////////

if(strcmp(respfile,"blank")!=0)	//have phenotypes - some of these checks may be redundant if using prsfile
{
count=countcols(respfile);
if(count<3)
{printf("Error, %s should have at least three columns (not %d)\n\n", respfile, count);exit(1);}
num_resps=count-2;

if(mpheno>num_resps)
{printf("Error, %s contains only %d phenotypes (so you can not use \"--mpheno %d\")\n\n", respfile, num_resps, mpheno);exit(1);}

if(mpheno2>num_resps)
{printf("Error, %s contains only %d phenotypes (so you can not use \"--mpheno2 %d\")\n\n", respfile, num_resps, mpheno2);exit(1);}

if((mode==229||mode==230)&&(mpheno==-9999||mpheno2==-9999)){printf("Error, you must specify two phenotypes using \"--mpheno\" and \"--mpheno2\"\n\n");exit(1);}

/*
if(strcmp(npheno,"blank")!=0)	//set mpheno
{
if(check_head_ids(respfile,0)==0)
{printf("Error, you can not use \"--pheno-name\" because %s does not have a header row\n\n", respfile);exit(1);}

if((input=fopen(respfile,"r"))==NULL)
{printf("Error opening %s\n\n", respfile);exit(1);}
if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading first two elements of %s\n\n", respfile);exit(1);}

mpheno=-9999;
for(j=0;j<num_resps;j++)
{
if(fscanf(input,"%s ", readstring)!=1)
{printf("Error reading element %d of %s\n\n", 3+j, respfile);exit(1);}
if(strcmp(readstring,npheno)==0)
{
if(mpheno!=-9999){printf("Error, %s contains at least two phenotypes called \"%s\"\n\n", respfile, npheno);exit(1);}
mpheno=j+1;
}
}
if(mpheno==-9999){printf("Error, %s does not contain a phenotype called \"%s\"\n\n", respfile, npheno);exit(1);}
fclose(input);
}
*/

if(mpheno==-9999)
{
if(num_resps>1)
{
printf("Error, %s contains multiple phenotypes, so you must specify one using \"--mpheno\"", respfile);
if(mode==121||mode==123||mode==124||mode==126||mode==129||mode==130||mode==131||mode==132||mode==151||mode==152||mode==153||mode==154||mode==163)
{printf(", or use \"--mpheno ALL\" to analyse all phenotypes");}
printf("\n\n");exit(1);
}
else{mpheno=1;}
}

keepresps=malloc(sizeof(int)*num_resps);
if(mpheno==-1)
{
num_resps_use=num_resps;
for(m=0;m<num_resps;m++){keepresps[m]=m;}
}
else{num_resps_use=1;keepresps[0]=mpheno-1;}

if(mode==151||mode==152||mode==153||mode==154)	//check whether the phenotypes are binary, and that dichot is consistent
{
head=check_head_ids(respfile,0);
count=countrows_plus(respfile,2+num_resps)-head;
readdoubles=malloc(sizeof(double)*count);

//see how many binary
count4=0;
for(m=0;m<num_resps_use;m++)
{
read_values(respfile, readdoubles, count, NULL, 3+keepresps[m], head, 1);

//see if all values are 0/1/missing
count2=0;for(i=0;i<count;i++){count2+=(readdoubles[i]!=-9999&&readdoubles[i]!=0&&readdoubles[i]!=1&&readdoubles[i]!=missingvalue);}

//or if all values are 1/2/missing
count3=0;for(i=0;i<count;i++){count3+=(readdoubles[i]!=-9999&&readdoubles[i]!=1&&readdoubles[i]!=2&&readdoubles[i]!=missingvalue);}

count4+=(count2==0||count3==0);
}

if(num_resps_use==1)
{
if(dichot==-9999&&count4==1)
{printf("Error, the response in %s appears to be binary, so you must use either \"--binary YES\" or \"--binary NO\" to tell LDAK whether or not to use a quasi-logistic model\n\n", respfile);exit(1);}
if(dichot==1&&count4==0)
{printf("Error, the response in %s does not appear to be binary, so you can not use \"--binary YES\"\n\n", respfile);exit(1);}
}
else
{
if(dichot==-9999&&count4>0)
{printf("Error, some of the responses in %s appear to be binary, so you must use either \"--binary YES\" or \"--binary NO\" to tell LDAK whether or not to use a quasi-logistic model\n\n", respfile);exit(1);}
if(dichot==1&&count4<num_resps_use)
{printf("Error, some of the responses in %s appear not to be binary, so you can not use \"--binary YES\"\n\n", respfile);exit(1);}
}

free(readdoubles);
}
}	//end of got respfile
else	//do not have phenotype
{
if(mode==115||mode==121||mode==122||mode==125||mode==138||mode==140||mode==172)	//will fake a response
{num_resps_use=1;keepresps=malloc(sizeof(int));keepresps[0]=0;}
else
{num_resps_use=0;}
}

////////

if(strcmp(sumsfile,"blank")!=0)	//check correct columns and set amb and scaling
{
count=countcols(sumsfile);

//assume it is in LDAK format, unless find otherwise
sformat=0;

if(sformat==0)	//check whether in cojo format
{
if(find_head("SNP", sumsfile, count)!=-1&&find_head("A1", sumsfile, count)!=-1&&find_head("A2", sumsfile, count)!=-1&&find_head("b", sumsfile, count)!=-1&&find_head("se", sumsfile, count)!=-1&&find_head("N", sumsfile, count)!=-1)
{
printf("Have detected that %s is in COJO format (has columns labelled \"SNP\", \"A1\", \"A2\", \"b\", \"se\" and \"N\")\n\n", sumsfile);
sformat=1;
}
}

if(sformat==0)	//check whether in prs-cs format
{
if(find_head("SNP", sumsfile, count)!=-1&&find_head("A1", sumsfile, count)!=-1&&find_head("A2", sumsfile, count)!=-1&&find_head("BETA", sumsfile, count)!=-1&&find_head("SE", sumsfile, count)!=-1)
{
printf("Have detected that %s is in PRS-CS format (has columns labelled \"SNP\", \"A1\", \"A2\", \"BETA\" and \"SE\")\n\n", sumsfile);
if(fixn==-9999){printf("Error, you should use \"--fixed-n\" to provide the sample size corresponding to %s\n\n", sumsfile);exit(1);}
sformat=2;
}
}

if(sformat==0)	//check it is in LDAK format
{
if(find_head("Predictor", sumsfile, count)==-1&&find_head("SNP", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"Predictor\" or \"SNP\"\n\n", sumsfile);exit(1);}
if(find_head("A1", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sumsfile);exit(1);}
if(find_head("A2", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sumsfile);exit(1);}

if(mode==146||(mode==147&&plet==1))	//only need unsigned statistics (Z, Stat or P)
{
if(find_head("Z", sumsfile, count)==-1&&find_head("Stat", sumsfile, count)==-1&&find_head("P", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"Z\", \"Stat\" or \"P\"\n\n", sumsfile);exit(1);}
}
else	//will need signed statistics (either Z, Stat+Direction or P+Direction)
{
if(find_head("Z", sumsfile, count)==-1&&(find_head("Stat", sumsfile, count)==-1||find_head("Direction", sumsfile, count)==-1)&&(find_head("P", sumsfile, count)==-1||find_head("Direction", sumsfile, count)==-1))
{printf("Error, %s should have a column named \"Z\", two columns named \"Stat\" and \"Direction\", or two columns named \"P\" and \"Direction\"\n\n", sumsfile);exit(1);}
}

if(fixn==-9999)
{
if(find_head("n", sumsfile, count)==-1&&find_head("N", sumsfile, count)==-1)
{printf("Error, %s should have a column named \"n\" or \"N\" (or use \"--fixed-n\" if sample size is constant)\n\n", sumsfile);exit(1);}
}
else
{
if(find_head("n", sumsfile, count)==-1&&find_head("N", sumsfile, count)==-1)
{printf("Warning, you have used \"--fixed-n\", so the sample sizes in %s will be ignored\n\n", sumsfile);}
}
}

//see whether there is a column called A1Freq
gotfreq=(find_head("A1Freq", sumsfile, count)!=-1);

if(checkfreq==1&&gotfreq==0)
{printf("Error, %s does not contain a column named \"A1Freq\", so it is not possible to use \"--check-frequencies YES\"\n\n", sumsfile);exit(1);}

if(amb==-9999)
{
if(mode==146||(mode==147&&plet==1)||mode==150){amb=1;}
else{amb=0;}
}
if(scaling==-9999){scaling=1;}
}

////////

if(strcmp(sums2file,"blank")!=0)	//check correct columns and set scaling (amb already set)
{
count=countcols(sums2file);

//assume it is in LDAK format, unless find otherwise
sformat2=0;

if(sformat2==0)	//check whether in cojo format
{
if(find_head("SNP", sums2file, count)!=-1&&find_head("A1", sums2file, count)!=-1&&find_head("A2", sums2file, count)!=-1&&find_head("b", sums2file, count)!=-1&&find_head("se", sums2file, count)!=-1&&find_head("N", sums2file, count)!=-1)
{
printf("Have detected that %s is in COJO format (has columns labelled \"SNP\", \"A1\", \"A2\", \"b\", \"se\" and \"N\")\n\n", sums2file);
sformat2=1;
}
}

if(sformat2==0)	//check whether in prs-cs format
{
if(find_head("SNP", sums2file, count)!=-1&&find_head("A1", sums2file, count)!=-1&&find_head("A2", sums2file, count)!=-1&&find_head("BETA", sums2file, count)!=-1&&find_head("SE", sums2file, count)!=-1)
{
printf("Have detected that %s is in PRS-CS format (has columns labelled \"SNP\", \"A1\", \"A2\", \"BETA\" and \"SE\")\n\n", sums2file);
if(fixn2==-9999){printf("Error, you should use \"--fixed-n2\" to provide the sample size corresponding to %s\n\n", sums2file);exit(1);}
sformat2=2;
}
}

if(sformat2==0)	//check it is in LDAK format
{
if(find_head("Predictor", sums2file, count)==-1&&find_head("SNP", sums2file, count)==-1)
{printf("Error, %s should have a column named \"Predictor\" or \"SNP\"\n\n", sums2file);exit(1);}
if(find_head("A1", sums2file, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sums2file);exit(1);}
if(find_head("A2", sums2file, count)==-1)
{printf("Error, %s should have a column named \"A1\"\n\n", sums2file);exit(1);}

if(mode==146||(mode==147&&plet==1))	//only need unsigned statistics (Z, Stat or P)
{
if(find_head("Z", sums2file, count)==-1&&find_head("Stat", sums2file, count)==-1&&find_head("P", sums2file, count)==-1)
{printf("Error, %s should have a column named \"Z\", \"Stat\" or \"P\"\n\n", sums2file);exit(1);}
}
else	//will need signed statistics (either Z, Stat+Direction or P+Direction)
{
if(find_head("Z", sums2file, count)==-1&&(find_head("Stat", sums2file, count)==-1||find_head("Direction", sums2file, count)==-1)&&(find_head("P", sums2file, count)==-1||find_head("Direction", sums2file, count)==-1))
{printf("Error, %s should have a column named \"Z\", columns named \"Stat\" and \"Direction\", or columns named \"P\" and \"Direction\"\n\n", sums2file);exit(1);}
}

if(fixn2==-9999)
{
if(find_head("n", sums2file, count)==-1&&find_head("N", sums2file, count)==-1)
{printf("Error, %s should have a column named \"n\" or \"N\" (or use \"--fixed-n2\" if sample size is constant)\n\n", sums2file);exit(1);}
}
else
{
if(find_head("n", sums2file, count)==-1&&find_head("N", sums2file, count)==-1)
{printf("Warning, you have used \"--fixed-n2\", so the sample sizes in %s will be ignored\n\n", sums2file);}
}
}

if(scaling2==-9999){scaling2=1;}
}

////////

if(strcmp(covarfile,"blank")!=0)	//have quantitative covariates (will add an intercept)
{
count=countcols(covarfile);
if(count<3){printf("Error, %s should have at least three columns (not %d)\n\n", covarfile, count);exit(1);}
checkcols(covarfile,count);

//get covariate indexes
keepcovars=malloc(sizeof(int)*(count-2));
if(strcmp(covarnums,"blank")==0&&strcmp(covarnames,"blank")==0)	//no filtering
{
num_quants=count-1;
for(j=0;j<count-2;j++){keepcovars[j]=j;}
}
else	//filtering
{
if(strcmp(covarnums,"blank")!=0)	//filtering based on numbers
{num_quants=1+find_covar_numbers(covarnums, keepcovars, count-2, covarfile);}
else	//filtering based on names
{num_quants=1+find_covar_names(covarnames, keepcovars, count-2, covarfile);}
}
}
else	//do not have quantitative covariates
{
if(mode==121||mode==122||mode==123||mode==124||mode==125||mode==126||mode==127||mode==128||mode==129||mode==130||mode==229||mode==230||mode==131||mode==132||mode==133||mode==138||mode==140||mode==151||mode==152||mode==153||mode==154||mode==156||mode==169||mode==170||mode==172||mode==173||mode==175||mode==194)	//will be using covar (these are all cases except adjust-grm, where covariates required)
{num_quants=1;}
else
{num_quants=0;}
}

if(strcmp(envfile,"blank")!=0)	//have environmental variables
{
count=countcols(envfile);
if(count<3){printf("Error, %s should have at least three columns (not %d)\n\n", envfile, count);exit(1);}
checkcols(envfile,count);
num_envs=count-2;

if(discenv==1&&num_envs<2)
{printf("Error, when using \"--subgroups YES\" there must be at least 2 environmental variables\n\n");exit(1);}

if(discenv==1&&(mode==121||mode==123||mode==124))	//estimating variances - check num_kins
{
if(num_kins!=2*num_envs){printf("Error, %s you must provide %d kinship matrices, created using \"--gxemm-free\" with \"--subgroups YES\" (not %d)\n\n", kinlist, 2*num_envs, num_kins);exit(1);}
}

if(mode==131&&num_envs!=1)
{printf("Error, when using \"--linear\", you can only have one environmental variable (not %d)\n\n", num_envs);exit(1);}
}
else
{num_envs=0;}

if(strcmp(topfile,"blank")!=0)	//have top preds
{num_tops=countrows(topfile);}
else
{num_tops=0;}

if(strcmp(factorfile,"blank")!=0)	//have factors
{
count=countcols(factorfile);
if(count<3){printf("Error, %s should have at least three columns (not %d)\n\n", factorfile, count);exit(1);}
checkcols(factorfile,count);
}

if(strcmp(povarfile,"blank")!=0)	//have prs covar
{
count=countcols(povarfile);
if(count<4){printf("Error, %s should have at least four columns (not %d)\n\n", povarfile, count);exit(1);}
checkcols(povarfile,count);

//how many times does profile appear in top row
if((input=fopen(povarfile,"r"))==NULL)
{printf("Error opening %s\n\n",povarfile);exit(1);}

num_prs=0;
for(j=0;j<count;j++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of %s\n\n", j+1, povarfile);exit(1);}
num_prs+=(strcmp(readstring,"Profile")==0);
}
if(num_prs==0){printf("Error, %s does not appear to be correctly formatted\n\n", povarfile);exit(1);}

fclose(input);

//get number of chrs
num_chr3=(count-2)/num_prs-1;
if(count!=2+num_prs*(1+num_chr3)){printf("Error, %s does not appear to be correctly formatted\n\n", povarfile);exit(1);}
printf("%s contains values for %d PRS and %d chromosomes\n\n", povarfile, num_prs, num_chr3);

//read chr numbers
chrindex3=malloc(sizeof(int)*num_chr3);

if((input=fopen(povarfile,"r"))==NULL)
{printf("Error opening %s\n\n",povarfile);exit(1);}

//skip three elements
if(fscanf(input, "%s %s %s ", readstring, readstring, readstring)!=3)
{printf("Error reading first three elements %d of %s\n\n", j+1, povarfile);exit(1);}

for(j=0;j<num_chr3;j++)
{
if(fscanf(input, "Chr%d ", chrindex3+j)!=1)
{printf("Error reading Element %d of %s\n\n", j+4, povarfile);exit(1);}
}

for(k=1;k<num_prs;k++)	//check same chr
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Element %d of %s\n\n", 3+k*(1+num_chr3), povarfile);exit(1);}

for(j=0;j<num_chr3;j++)
{
if(fscanf(input, "Chr%d ", &readint)!=1)
{printf("Error reading Element %d of %s\n\n", 4+j+k*(1+num_chr3), povarfile);exit(1);}

if(readint!=chrindex3[j])
{printf("Error, chromosomes for PRS %d do not match those for PRS 1\n\n", k+1);exit(1);}
}
}

fclose(input);
}
else{num_prs=0;}

////////

if(strcmp(offsetfile,"blank")!=0)	//have offsets
{
count=countcols(offsetfile);
if(count!=3){printf("Error, %s should have three columns (not %d)\n\n", offsetfile, count);exit(1);}
checkcols(offsetfile,count);
}

///////////////////////////

//calculating kinships

if(strcmp(partpref,"blank")!=0)	//check partitions exist and contain elements
{
for(q=0;q<num_parts;q++)
{
sprintf(filename,"%s%d",partpref, q+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}
}

if(strcmp(invsfile,"blank")!=0)	//check inverse files exist and correct size
{
sprintf(filename,"%s.inverse.predictors", invsfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--solve-gre\"\n\n", filename);exit(1);}
count=countrows(filename);

sprintf(filename,"%s.inverse.bin", invsfile);	
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--solve-gre\"\n\n", filename);exit(1);}

if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(float)*(count+1)*count/2)
{printf("Error reading %s; should have size %jd not %jd\n\n", filename, (off_t)sizeof(float)*(count+1)*count/2, ftello(input));exit(1);}
fclose(input);
}

///////////////////////////

//reml, blup and he/pcgc

if(strcmp(hersfile,"blank")!=0)	//check size and that values add to one
{
count=countcols(hersfile);
if(count!=1)
{printf("Error, %s should have one column (not %d)\n\n", hersfile, count);exit(1);}
checkcols(hersfile,count);

count2=countrows(hersfile);
if(count2!=num_kins+num_regs)
{printf("Error, %s should have %d rows (not %d)\n\n", hersfile, num_kins+num_regs, count2);exit(1);}

vstarts=malloc(sizeof(double)*(num_kins+num_regs));
read_values(hersfile, vstarts, num_kins+num_regs, NULL, 1, 0, 0);
sum=0;for(k=0;k<num_kins+num_regs;k++){sum+=vstarts[k];}
if(sum<0||sum>1)
{printf("Warning, the heritabilities in %s sum to %.4f (normally the sum is between zero and one)\n\n", hersfile, sum);}
free(vstarts);
}

if(strcmp(oversfile,"blank")!=0)	//check size and read proportions into ssums
{
count=countcols(oversfile);
if(count!=num_kins){printf("Error, the number of categories in %s (%d) does not match the number of kinship matrices (%d)\n\n", oversfile, count, num_kins);exit(1);}
checkcols(oversfile,count);

count2=countrows(oversfile);
if(count2!=count+3){printf("Error, %s should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-overlaps\"\n\n", oversfile, count+3, count2);exit(1);}

ssums=malloc(sizeof(double*)*num_kins);
for(k=0;k<num_kins;k++)
{
ssums[k]=malloc(sizeof(double)*(num_kins+2));
read_values(oversfile, ssums[k], num_kins+2, NULL, k+1, 1, 0);
}
}

////////

if(strcmp(remlfile,"blank")!=0)	//must be blup or pred - check numbers match those provided (will read covariates in later) 
{
if((input=fopen(remlfile,"r"))==NULL)
{printf("Error opening %s\n\n",remlfile);exit(1);}

printf("Reading details from %s\n", remlfile);

if(fscanf(input, "Num_Kinships %d ", &readint)!=1)
{printf("Error reading Row 1 of %s (should begin \"Num_Kinships\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(readint>0&&num_kins==0)
{
if(readint==1){printf("Error, you must use \"--grm\" to provide the kinship matrix used when performing REML\n\n");exit(1);}
else{printf("Error, you must use \"--mgrm\" to provide the %d kinship matrices used when performing REML\n\n", readint);exit(1);}
}
if(readint!=num_kins)
{printf("Error, the number of kinship matrices provided now (%d) is different to the number in %s (%d)\n\n", num_kins, remlfile, readint);exit(1);}

if(fscanf(input, "Num_Regions %d ", &num_scores)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Regions\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}

if(num_scores>0&&mode==125){printf("Error, it is not possible to use \"--reml-predict\" if regions were used when performing REML\n\n");exit(1);}

if(fscanf(input, "Num_Top_Predictors %d ", &num_tops)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Top_Predictors\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}

if(num_tops>0&&mode==125){printf("Error, it is not possible to use \"--reml-predict\" if top predictors were used when performing REML\n\n");exit(1);}

if(fscanf(input, "Num_Covariates %d ", &readint)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Covariates\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(covarfile,"blank")!=0)	//covariates provided
{
if(readint==1){printf("Error, the covariates provided now were not used when performing REML\n\n");exit(1);}
if(num_quants!=readint){printf("Error, the number of covariates in %s (%d) does not match the number used when performing REML (%d)\n\n", covarfile, num_quants-1, readint-1);exit(1);}
}
else	//not provided
{
if(readint>1){printf("Error, you must use \"--covar\" to provide the covariates used when performing REML\n\n");exit(1);}
}

if(fscanf(input, "Num_Environments %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Environments\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(envfile,"blank")!=0)	//environments provided
{
if(readint==0){printf("Error, the environmental variables provided now were not used when performing REML\n\n");exit(1);}
if(num_envs!=readint){printf("Error, the number of environmental variables in %s (%d) does not match the number in %s (%d)\n\n", envfile, num_envs, remlfile, readint);exit(1);}
}
else	//not provided
{
if(readint>0){printf("Error, you must use \"--enviro\" to provide the environmental variables used when performing REML\n\n");exit(1);}
}

if(fscanf(input, "Blupfile %s ", readstring)!=1)
{printf("Error reading Row 6 of %s (should begin \"Blupfile\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(readstring,"none")!=0)
{
strcpy(blupfile2,readstring);
if(append_check(blupfile,blupfile2,workdir)!=0)
{printf("Error reading blupfile %s; this file would have been created when performing REML\n\n", blupfile);exit(1);}

if(countcols(blupfile)!=2+2*(num_kins+num_scores))
{printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation when performing REML\n\n ", blupfile, 2*(num_kins+num_scores), countcols(blupfile));exit(1);}
}
else{strcpy(blupfile,"blank");}

if(fscanf(input, "Regfile %s ", readstring)!=1)
{printf("Error reading Row 7 of %s (should begin \"Regfile\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(strcmp(readstring,"none")!=0)
{
strcpy(regfile2,readstring);
if(append_check(regfile,regfile2,workdir)!=0)
{printf("Error reading regfile %s; this file would have been created when performing REML\n\n", regfile);exit(1);}

if(countcols(regfile)!=4+num_scores+(num_tops>0))
{printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation when performing REML\n\n", regfile, 4+num_scores+(num_tops>0), countcols(regfile));exit(1);}
}

if(fscanf(input, "Coeffsfile %s ", cofile2)!=1)
{printf("Error reading Row 8 of %s (should begin \"Coeffsfile\"), suggesting the file has been changed since creation when performing REML\n\n",remlfile);exit(1);}
if(append_check(cofile,cofile2,workdir)!=0)
{printf("Error reading coeffsfile %s; this file would have been created when performing REML\n\n", cofile);exit(1);}
fclose(input);

printf("\n");
}	//end of using remlfile

////////

if(mode==123||mode==124)	//he or pcgc
{
//first set adjusted
if(adjusted==-9999)
{
if(num_quants+num_envs+num_tops==1){adjusted=0;}
else{adjusted=1;}
}

if(adjusted==1&&num_kins>0)
{
//first check all kinships match
sprintf(filename,"%s.grm.root", kinstems[0]);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--adjust-grm\"; if you are not using adjusted kinships, use \"--adjusted NO\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
if(fscanf(input, "Covariates %s ", readstring2)!=1)
{printf("Error reading Row 2 of %s (should begin \"Covariates\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
if(fscanf(input, "Top_Predictors %s ", readstring3)!=1)
{printf("Error reading Row 3 of %s (should begin \"Top_Predictors\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
if(fscanf(input, "Environments %s ", readstring4)!=1)
{printf("Error reading Row 4 of %s (should begin \"Environments\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename);exit(1);}
fclose(input);

for(k=1;k<num_kins;k++)
{
sprintf(filename2,"%s.grm.root", kinstems[k]);
if(just_check(filename2)!=0)
{printf("Error reading %s; this file would have been created using \"--adjust-grm\"; if you are not using adjusted kinships, use \"--adjusted NO\"\n\n", filename2);exit(1);}
if((input2=fopen(filename2,"r"))==NULL)
{printf("Error opening %s\n\n",filename2);exit(1);}
if(fscanf(input2, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(fscanf(input2, "Covariates %s ", readstring)!=1)
{printf("Error reading Row 2 of %s (should begin \"Covariates\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(strcmp(readstring,readstring2)!=0)
{printf("Error, the covariate file listed in %s (%s) does not match that listed in %s (%s)\n\n", filename2, readstring, filename, readstring2);exit(1);}
if(fscanf(input2, "Top_Predictors %s ", readstring)!=1)
{printf("Error reading Row 3 of %s (should begin \"Top_Predictors\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(strcmp(readstring,readstring3)!=0)
{printf("Error, the top predictors file listed in %s (%s) does not match that listed in %s (%s)\n\n", filename2, readstring, filename, readstring3);exit(1);}
if(fscanf(input2, "Environments %s ", readstring)!=1)
{printf("Error reading Row 4 of %s (should begin \"Environments\"), suggesting the file has been changed since creation with \"--adjust-grm\"\n\n",filename2);exit(1);}
if(strcmp(readstring,readstring4)!=0)
{printf("Error, the environmental variables file listed in %s (%s) does not match that listed in %s (%s)\n\n", filename2, readstring, filename, readstring4);exit(1);}
fclose(input2);
}

if(strcmp(covarfile,"blank")==0&&strcmp(readstring2,"none")!=0)
{printf("Error, you must use \"--covar\" to provide the covariates used with \"--adjust-grm\"\n\n");exit(1);}
if(strcmp(covarfile,"blank")!=0&&strcmp(readstring2,"none")==0)
{printf("Error, covariates are provided now, but were not used with \"--adjust-grm\"\n\n");exit(1);}

if(strcmp(topfile,"blank")==0&&strcmp(readstring3,"none")!=0)
{printf("Error, you must use \"--top-preds\" to provide the top predictors used with \"--adjust-grm\"\n\n");exit(1);}
if(strcmp(topfile,"blank")!=0&&strcmp(readstring3,"none")==0)
{printf("Error, top predictors are provided now, but were not used with \"--adjust-grm\"\n\n");exit(1);}

if(strcmp(envfile,"blank")==0&&strcmp(readstring4,"none")!=0)
{printf("Error, you must use \"--enviro\" to provide the environmental variables used with \"--adjust-grm\"\n\n");exit(1);}
if(strcmp(envfile,"blank")!=0&&strcmp(readstring4,"none")==0)
{printf("Error, environmental variables are provided now, but were not used with \"--adjust-grm\"\n\n");exit(1);}

if(checkroot==1&&strcmp(readstring2,"none")!=0)	//check covarfiles match
{
(void)append_check(readstring,readstring2,workdir);
if(strcmp(covarfile,readstring)!=0&&strcmp(covarfile,readstring2)!=0)
{printf("Error, the kinship matrix %s corresponds to the covariate file %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", kinstems[0], readstring2, covarfile);exit(1);}
}
if(checkroot==1&&strcmp(readstring3,"none")!=0)	//check top predictor files match
{
(void)append_check(readstring,readstring3,workdir);
if(strcmp(topfile,readstring)!=0&&strcmp(topfile,readstring3)!=0)
{printf("Error, the kinship matrix %s corresponds to the top predictors file %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", kinstems[0], readstring3, topfile);exit(1);}
}
if(checkroot==1&&strcmp(readstring4,"none")!=0)	//check envfiles match
{
(void)append_check(readstring,readstring4,workdir);
if(strcmp(envfile,readstring)!=0&&strcmp(envfile,readstring4)!=0)
{printf("Error, the kinship matrix %s corresponds to the environmental variables file %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", kinstems[0], readstring4, envfile);exit(1);}
}
}
}

////////

if(strcmp(relfile,"blank")!=0)
{
count=countcols(relfile);
if(count!=5&&count!=6){printf("Error, %s should have either five or six columns (not %d), containing pairs of samples, their estimated relatedness, and perhaps also their environmental similarity\n\n", relfile, count);exit(1);}
checkcols(relfile,count);
}

///////////////////////////

//association analysis

if(mode==131&&trios==1)
{
count=countcols(famfile);
if(count<4){printf("Error, %s should have at least four columns (not %d), with the 3rd and 4th columns containing paternal and maternal IDs (PID and MID)\n\n", famfile, count);exit(1);}
checkcols(famfile,count);
}

if(strcmp(sampwfile,"blank")!=0)
{
count=countcols(sampwfile);
if(count!=3)
{printf("Error, %s should have three columns (not %d), containing the two sample IDs then a positive weight\n\n", sampwfile, count);exit(1);}
checkcols(sampwfile,count);
}

if(strcmp(transfile,"blank")!=0)
{
count=countcols(transfile);
if(count!=num_resps)
{printf("Error, %s should have %d columns (not %d), one for each phenotype in %s\n\n", transfile, num_resps, count, respfile);exit(1);}
count=countrows(transfile);
if(count!=num_resps)
{printf("Error, %s should have %d rows (not %d), one for each phenotype in %s\n\n", transfile, num_resps, count, respfile);exit(1);}
}

////////

if(strcmp(genefile,"blank")!=0)	//check correct size
{
count=countcols(genefile);
if(count!=4&&count!=5)
{printf("Error, %s should either four or five columns (not %d), providing the name, chromosome, start and end basepairs, and perhaps the orientation (\"+\" or \"-\") of each genomic region\n\n", genefile, count);exit(1);}
checkcols(genefile,count);

if(count!=5&&(up_buffer>0||down_buffer>0))
{printf("Error, to use \"--up-buffer\" and \"--down-buffer\", %s should have five columns (not %d), providing the name, chromosome, start and end basepairs, and the orientation (\"+\" or \"-\") of each genomic region\n\n", genefile, count);exit(1);}
}

///////////////////////////

//sumher

if(strcmp(annpref,"blank")!=0)	//check annotations exist and contain elements
{
for(q=0;q<num_anns;q++)
{
sprintf(filename,"%s%d",annpref, q+1);
if(just_check(filename)!=0)
{printf("Error reading %s; check file exists and is not empty\n\n", filename);exit(1);}
}
}

if(strcmp(labfile,"blank")!=0&&num_parts>0)
{
count=countcols(labfile);
if(count!=1)
{printf("Error, %s should have one column (not %d)\n\n", labfile, count);exit(1);}
checkcols(labfile,count);

if(countrows(labfile)!=num_parts)
{printf("Error, the number of rows of %s (%d) does not match the number of partitions (%d)\n\n", labfile, countrows(labfile), num_parts);exit(1);}
}

////////

if(strcmp(taglist,"blank")!=0)	//have taglist (might also have matlist - can not have matlist without taglist)
{
num_tags=countrows(taglist);
if(num_tags==1){printf("Error, %s must contain at least two tagging files (not one)\n\n", taglist);exit(1);}

tagstems=malloc(sizeof(char *)*num_tags);
if((input=fopen(taglist,"r"))==NULL)
{printf("Error opening %s\n\n",taglist);exit(1);}
for(k=0;k<num_tags;k++)
{
tagstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, taglist);exit(1);}
if(append_check(tagstems[k],filename,workdir)!=0)
{printf("Error reading tagging file %s; check file exists and is not empty\n\n", tagstems[k]);exit(1);}

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(tagstems[k],tagstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", taglist, tagstems[k]);exit(1);}
}
}
fclose(input);
}

if(strcmp(matlist,"blank")!=0)
{
if(countrows(matlist)!=num_tags)
{printf("Error, the number of rows of %s (%d) does not match the number of rows of %s (%d)\n\n", matlist, countrows(matlist), taglist, num_tags);exit(1);}

matstems=malloc(sizeof(char *)*num_tags);
if((input=fopen(matlist,"r"))==NULL)
{printf("Error opening %s\n\n",matlist);exit(1);}
for(k=0;k<num_tags;k++)
{
matstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, matlist);exit(1);}
if(append_check(matstems[k],filename,workdir)!=0)
{printf("Error reading heritability matrix %s; check file exists and is not empty\n\n", matstems[k]);exit(1);}

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(matstems[k],matstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", matlist, matstems[k]);exit(1);}
}
}
fclose(input);
}

////////

if(strcmp(pathlist,"blank")!=0)	//have pathlist
{
num_tags=countrows(pathlist);
if(num_tags==1){printf("Error, %s must contain at least two tagging files (not one)\n\n", pathlist);exit(1);}

tagstems=malloc(sizeof(char *)*num_tags);
if((input=fopen(pathlist,"r"))==NULL)
{printf("Error opening %s\n\n",pathlist);exit(1);}
for(k=0;k<num_tags;k++)
{
tagstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, pathlist);exit(1);}
sprintf(filename2,"%s.base.tagging",filename);
(void)append_check(tagstems[k],filename,workdir);

for(k2=0;k2<k;k2++)	//check different to previous
{
if(strcmp(tagstems[k],tagstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", pathlist, tagstems[k]);exit(1);}
}
}
fclose(input);
}

////////

if(strcmp(tagfile,"blank")!=0)	//check correct size and set num_parts and parttype (0=annotations, 1=partitions)
{
count=countcols(tagfile);
if(count<10)
{printf("Error, %s should have at least 10 columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, count);exit(1);}
num_parts=count-9;

count=countrows_min(tagfile,num_parts+3);
if(count<num_parts+3)
{printf("Error, %s should have at least %d rows, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, num_parts+3);exit(1);}

if((input=fopen(tagfile,"r"))==NULL)
{printf("Error opening %s\n\n",tagfile);exit(1);}
if(fscanf(input, "%s %s %s %s %s %s %s %s %s ", readstring, readstring2, readstring2, readstring2, readstring2, readstring2, readstring2, readstring2, readstring2)!=9)
{printf("Error reading Row 1 of %s\n\n", tagfile);exit(1);}
if(strcmp(readstring,"Predictor")!=0)
{printf("Error reading %s; first element should be \"Predictor\" (not %s), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", tagfile, readstring);exit(1);}

for(q=0;q<num_parts;q++)
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading element %d of %s, suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", j+1, tagfile);exit(1);}
}
parttype=(strcmp(readstring,"Base")!=0);

fclose(input);
}

////////

if(strcmp(pathfile,"blank")!=0)	//check have files, check sizes and get num_parts
{
sprintf(filename,"%s.base.tagging", pathfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}

count=countcols(filename);
if(count!=11){printf("Error, %s should have 11 columns (not %d), suggesting the file has been changed since creation with \"--calc-taggings\"\n\n", filename, count);exit(1);}

count2=countrows(filename)-1;

sprintf(filename,"%s.pathway.sums", pathfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}

count=countcols(filename);
if(count!=13){printf("Error, %s should have 13 columns (not %d), suggesting the file has been changed since creation with \"--calc-taggings\"\n\n", filename, count);exit(1);}

num_parts=countrows(filename);

sprintf(filename,"%s.pathway.tagging", pathfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-taggings\"\n\n", filename);exit(1);}

count=countcols(filename);
if(count!=count2+1){printf("Error, %s should have %d columns (not %d), suggesting the file has been changed since creation with \"--calc-taggings\"\n\n", filename, count2+1, count);exit(1);}
}

////////

if(strcmp(altfile,"blank")!=0)	//check correct size
{
count=countcols(altfile);
if(count!=2){printf("Error, %s should have two columns, providing predictor names then taggings (not %d)\n\n", altfile, count);exit(1);}
checkcols(altfile,count);
}

//dont check cvsfile

if(strcmp(catfile,"blank")!=0)	//set num_reds, keeppart and keeppart2
{
keepparts=malloc(sizeof(int)*num_parts);
keepparts2=malloc(sizeof(int)*num_parts);

count=countcols(catfile);
if(count!=1){printf("Error, %s should have one column (not %d)\n\n", catfile, count);exit(1);}
checkcols(catfile,count);

count=countrows(catfile);
printf("Reading the %d categories from %s\n", count, catfile);

for(j=0;j<num_parts;j++){keepparts[j]=0;}
if(parttype==0){keepparts[num_parts-1]=1;}

if((input=fopen(catfile,"r"))==NULL)
{printf("Error opening %s\n\n",catfile);exit(1);}

for(j=0;j<count;j++)
{
if(fscanf(input, "%d ", &readint)!=1)
{printf("Error reading Element %d of %s\n\n", j+1, catfile);exit(1);}

if(readint<0||(readint==0&&(parttype==1||count>1)))
{printf("Error reading Element %d of %s; all values must be positive integers (not %d)", j+1, catfile, readint);exit(1);}

if(readint!=0)	//readint=0 is shorthand for only the base
{
if(parttype==0&&readint>num_parts-1){printf("Error reading %s; Row %d specifies Annotation %d, which is greater than the number of annotations (%d)\n\n", catfile, j+1, readint, num_parts-1);exit(1);}
if(parttype==1&&readint>num_parts){printf("Error reading %s; Row %d specifies Partition %d, which is greater than the number of partitions (%d)\n\n", catfile, j+1, readint, num_parts);exit(1);}

if(keepparts[readint-1]==1){printf("Error, %d appears twice in %s\n\n", readint, catfile);exit(1);}
keepparts[readint-1]=1;
}
}
fclose(input);

num_reds=0;
for(j=0;j<num_parts;j++)
{
if(keepparts[j]==1){keepparts2[num_reds]=j;num_reds++;}
}
}

//check taufile in required.c

////////

if(strcmp(powfile,"blank")!=0)	//check contains numbers
{
count=countcols(powfile);
if(count!=1){printf("Error, %s should have one column (not %d)\n\n", powfile, count);exit(1);}
checkcols(powfile,count);

count=countrows(powfile);
printf("Checking the %d values in %s\n", count, powfile);

powers=malloc(sizeof(double)*count);
read_values(powfile,powers,count,NULL,1,0,0);

for(k=0;k<count;k++)
{
if(powers[k]<-1.25||powers[k]>0.25)
{printf("Warning, the power on Row %d is %.4f (usually values are between -1.25 and 0.25)\n", k+1, powers[k]);}
}

/*
if((input=fopen(powfile,"r"))==NULL)
{printf("Error opening %s\n\n",powfile);exit(1);}
for(k=0;k<count;k++)
{
if(fscanf(input, "%[0-9eE.+-] ", readstring)==1)
{
if(atof(readstring)<-1.25||atof(readstring)>0.25)
{printf("Warning, the power on Row %d is %.4f (usually values are between -1.25 and 0.25)\n", k+1, atof(readstring));}
}
else
{
if(fscanf(input, "%s ", readstring)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, powfile);exit(1);}
printf("Error reading %s; the first value on Row %d is unrecognisable (%s)\n\n", powfile, k+1, readstring);exit(1);
}
}
fclose(input);
*/

free(powers);
printf("\n");
}

////////

if(strcmp(matfile,"blank")!=0)	//check correct size (will have tagfile, and have set num_parts)
{
count=countcols(matfile);
if(count<2)
{printf("Error, %s should have at least two columns (not %d), suggesting the file has been changed since creation with \"--calc-tagging\"\n\n", matfile, count);exit(1);}
if(count-1!=num_parts){printf("Error, the number of categories in %s (%d) does not match the number in %s (%d)\n\n", matfile, count-1, tagfile, num_parts);exit(1);}
}

///////////////////////////

//individual-level data prediction, then megaprs

if(strcmp(indhers,"blank")!=0)	//check correct size
{
count=countcols(indhers);
if(count!=2)
{printf("Error, %s should have two columns, providing predictor names then info score (not %d)\n\n", indhers, countcols(indhers));exit(1);}
checkcols(indhers,count);
}

////////

if(strcmp(fastfile,"blank")!=0)	//check samples and overlap
{
//first check ID file
sprintf(filename,"%s.grm.id", fastfile);
if(just_check(filename)!=0)
{printf("Error reading %s\n\n", filename);exit(1);}
count=countrows(filename);

//now check pairs file
sprintf(filename,"%s.grm.sp", fastfile);
if(just_check(filename)!=0)
{printf("Error reading %s\n\n", filename);exit(1);}
count2=countrows(filename);

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

for(j=0;j<count2;j++)
{
if(fscanf(input, "%d %d %lf ", &readint, &readint2, &readdouble)!=3)
{printf("Error reading Row %d of %s\n\n", j+1, filename);exit(1);}
if(readint>=count||readint2>=count)
{printf("Error reading Row %d of %s; one of the indexes (%d and %d) is higher than the total number of individuals (%d)\n\n", j+1, filename, readint, readint2, count);exit(1);}
}

fclose(input);
}

////////

//dont check bvsfile

////////

if(strcmp(blockfile,"blank")!=0)
{
count=countcols(blockfile);
if(count!=2)
{printf("Error, %s should have two columns (not %d), providing the chromosome and basepair of each break point\n\n", blockfile, count);exit(1);}
checkcols(blockfile,count);

count=countrows(blockfile)-1;

if((input=fopen(blockfile,"r"))==NULL)
{printf("Error opening %s\n\n",blockfile);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading the headers of %s\n\n", blockfile);exit(1);}
if(strcmp(readstring,"Chr")!=0&&strcmp(readstring,"Chromosome")!=0)
{printf("Error reading %s; the column names should be Chr and BP (not %s and %s)\n\n", blockfile, readstring, readstring2);exit(1);}
if(strcmp(readstring2,"BP")!=0&&strcmp(readstring2,"Basepair")!=0)
{printf("Error reading %s; the column names should be Chr and BP (not %s and %s)\n\n", blockfile, readstring, readstring2);exit(1);}

if(fscanf(input, "%d %lf ", &readint, &readdouble)!=2)
{printf("Error reading Row 2 of %s\n\n", blockfile);exit(1);}

for(j=1;j<count;j++)
{
if(fscanf(input, "%d %lf ", &readint2, &readdouble2)!=2)
{printf("Error reading Row %d of %s\n\n", j+2, blockfile);exit(1);}

if(readint2<readint){printf("The chromosome on Row %d (%d) is smaller than that on Row %d (%d); the breakpoints must be in order\n\n", j+2, readint2, j+1, readint);exit(1);}

if(readint2==readint&&readdouble2<readdouble){printf("The location on Row %d (Chr %d, BP %.2f) is smaller than that on Row %d (Chr %d, BP %.2f); the breakpoints must be in order\n\n", j+2, readint2, readdouble2, j+1, readint, readdouble);exit(1);}

readint=readint2;
readdouble=readdouble2;
}
fclose(input);
}

////////

if(strcmp(corslist,"blank")!=0)
{
num_cors=countrows(corslist);
corstems=malloc(sizeof(char *)*num_cors);
if((input=fopen(corslist,"r"))==NULL)
{printf("Error opening %s\n\n",corslist);exit(1);}

for(k=0;k<num_cors;k++)
{
corstems[k]=malloc(sizeof(char)*500);
if(fscanf(input, "%s ", filename)!=1)
{printf("Error reading Row %d of %s\n\n", k+1, corslist);exit(1);}
(void)append_check(corstems[k],filename,workdir);

if(k>0)	//check different to previous
{
for(k2=0;k2<k;k2++)
{
if(strcmp(corstems[k],corstems[k2])==0)
{printf("Error reading %s; the file %s appears twice\n\n", corslist, corstems[k]);exit(1);}
}
}
}
fclose(input);
}

////////

if(strcmp(corname,"blank")!=0)	//must be mode=159 - check files then read windows
{
//first read root - need num_preds_used, num_windows, num_pairs, and num_blocks (for checking sizes of some other files)
sprintf(filename,"%s.cors.root", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=8)
{printf("Error reading %s; should have eight rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename, count);exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

if(fscanf(input, "%s %s ", readstring, readstring2)!=2)
{printf("Error reading Row 1 of %s\n\n", filename);exit(1);}
if(strcmp(readstring,"Datafile")!=0&&strcmp(readstring,"Multiple")!=0)
{printf("Error, %s should begin \"Datafile\" or \"Multiple\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename, readstring);exit(1);}

if(fscanf(input, "Num_Samples %s ", readstring)!=1)
{printf("Error reading Row 2 of %s (should begin \"Num_Samples\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Predictors %s ", readstring)!=1)
{printf("Error reading Row 3 of %s (should begin \"Num_Predictors\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}
if(fscanf(input, "Num_Samples_Used %s ", readstring)!=1)
{printf("Error reading Row 4 of %s (should begin \"Num_Samples_Used\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Predictors_Used %d ", &readint)!=1)
{printf("Error reading Row 5 of %s (should begin \"Num_Predictors_Used\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Windows %d ", &readint2)!=1)
{printf("Error reading Row 6 of %s (should begin \"Num_Windows\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Pairs %jd ", &scount)!=1)
{printf("Error reading Row 7 of %s (should begin \"Num_Pairs\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}

if(fscanf(input, "Num_Jackknifes %d ", &num_blocks)!=1)
{printf("Error reading Row 8 of %s (should begin \"Num_Jackknifes\"), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename);exit(1);}

fclose(input);

//now check bin file exists and correct size
sprintf(filename,"%s.cors.bin", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);
if(ftello(input)!=(off_t)sizeof(double)*readint*4+sizeof(float)*scount)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)sizeof(double)*readint*4+sizeof(float)*scount, ftello(input));exit(1);}
fclose(input);

//check bim file exists and correct size
sprintf(filename,"%s.cors.bim", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename, readint, count);exit(1);}

//check windows file exists and correct size
sprintf(filename,"%s.cors.windows", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=readint2+1)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename, readint2, count);exit(1);}

//check highLD file exists
sprintf(filename,"%s.cors.highld", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}

//check noise file exists and correct size
sprintf(filename,"%s.cors.noise", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}

count=countrows(filename);
if(count!=readint)
{printf("Error reading %s; should have %d rows (not %d), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n",filename, readint, count);exit(1);}

if(prsvar==1)	//check jackknife file exists and correct size
{
sprintf(filename,"%s.cors.jackknifes", corname);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}
if((input=fopen(filename,"rb"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
fseeko(input, 0, SEEK_END);

if(ftello(input)!=(off_t)sizeof(float)*num_blocks*readint)
{printf("Error reading %s; should have size %jd, but instead has size %jd\n\n", filename, (off_t)sizeof(float)*num_blocks*readint, ftello(input));exit(1);}
fclose(input);
}

//read windows

sprintf(filename,"%s.cors.windows",corname);
bittotal=countrows(filename)-1;
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}

//read top line
if(fscanf(input, "%s %s %s %s %s", readstring, readstring2, readstring2, readstring2, readstring2)!=5)
{printf("Error reading Row 1 of %s, suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n", filename);exit(1);}
if(strcmp(readstring,"Window")!=0)
{printf("Error, %s should begin \"Window\" (not %s), suggesting the file has been changed since creation with \"--calc-cors\" or \"--join-cors\"\n\n", filename, readstring);exit(1);}

blockstarts=malloc(sizeof(int)*bittotal);
blockends=malloc(sizeof(int)*bittotal);
blocklengths=malloc(sizeof(int)*bittotal);
blockindexes=malloc(sizeof(size_t)*bittotal);

for(j=0;j<bittotal;j++)
{
if(fscanf(input,"%d %d %d %s %s ", &readint, blockstarts+j, blockends+j, readstring, readstring)!=5)
{printf("Error reading Row %d of %s\n\n", 2+j, filename);exit(1);}
blockstarts[j]--;
blocklengths[j]=blockends[j]-blockstarts[j];
if(j==0){blockindexes[0]=0;}
else{blockindexes[j]=blockindexes[j-1]+(size_t)blocklengths[j-1]*(blocklengths[j-1]-1)/2;}
}
fclose(input);
}

////////

if(strcmp(pseudostem,"blank")!=0)	//check training and test summaries exist and look ok
{
sprintf(filename,"%s.train.summaries", pseudostem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pseudo-summaries\"\n\n", filename);exit(1);}
if(countcols(filename)!=6){printf("Error, %s should have six columns (not %d), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n", filename, countcols(filename));exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Predictor %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Predictor\"), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n",filename);exit(1);}
fclose(input);

sprintf(filename,"%s.test.summaries", pseudostem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pseudo-summaries\"\n\n", filename);exit(1);}
if(countcols(filename)!=6){printf("Error, %s should have six columns (not %d), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n", filename, countcols(filename));exit(1);}

if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Predictor %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Predictor\"), suggesting the file has been changed since creation with \"--pseudo-summaries\"\n\n",filename);exit(1);}
fclose(input);
}

////////

if(strcmp(bestfile,"blank")!=0)	//check has correct size
{
if(countrows(bestfile)!=2){printf("Error, %s should have two rows (not %d), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, countrows(bestfile));exit(1);}
if(countcols(bestfile)!=11){printf("Error, %s should have 11 columns (not %d), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, countcols(bestfile));exit(1);}

if((input=fopen(bestfile,"r"))==NULL)
{printf("Error opening %s\n\n",bestfile);exit(1);}
if(fscanf(input, "Model %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Model\"), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
fclose(input);
}

////////

if(strcmp(fracfile,"blank")!=0)	//either bolt, bayesr, elastic or mega with ptype 2, 3, 4, 5 or 6
{
count=countcols(fracfile);
checkcols(fracfile,count);

//if(mode==151)	//ridge - not allowed fracfile
if(mode==152)	//bolt
{
if(count!=2){printf("Error, %s should have two columns (not %d), providing p and f2\n\n", fracfile, count);exit(1);}
}
if(mode==153)	//bayesr
{
if(pointmass==1)
{
if(count!=4){printf("Error, %s should have four columns (not %d), providing  the fractions for the point mass, then the small, medium and large Gaussian distributions\n\n", fracfile, count);exit(1);}
}
else
{
if(count!=4){printf("Error, %s should have four columns (not %d), providing the fractions for the tiny, small, medium and large Gaussian distributions\n\n", fracfile, count);exit(1);}
}
}
if(mode==154)	//elastic
{
if(count!=2){printf("Error, %s should have two columns (not %d), providing p and f2\n\n", fracfile, count);exit(1);}
}

if(ptype==2||ptype==3)	//lasso and ridge
{
if(count!=1){printf("Error, %s should have one column (not %d), providing the heritability\n\n", fracfile, count);exit(1);}
}
if(ptype==4)	//bolt
{
if(count!=3){printf("Error, %s should have three columns (not %d), providing the heritability, p then f2\n\n", fracfile, count);exit(1);}
}
if(ptype==5)	//bayes
{
if(count!=5){printf("Error, %s should have five columns (not %d), providing the heritability, then fractions for the point mass and the small, medium and large Gaussian distributions\n\n", fracfile, count);exit(1);}
}
if(ptype==6)	//bayes-shrink
{
if(count!=5){printf("Error, %s should have five columns (not %d), providing the heritability, then fractions for the tiny, small, medium and large Gaussian distributions\n\n", fracfile, count);exit(1);}
}
}

///////////////////////////

//pca, decompose and adjust-grm

if(strcmp(pcastem,"blank")!=0)
{
sprintf(filename,"%s.vect", pcastem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created when using \"--pca\"\n\n", filename);exit(1);}
count=countcols(filename);
if(count<3)
{printf("Error, %s should have at least three columns (not %d), suggesting the file has been changed since creation with \"--pca\"\n\n",filename, count);exit(1);}
axes=count-2;

sprintf(filename,"%s.values", pcastem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pca\"\n\n", filename);exit(1);}

sprintf(filename,"%s.root", pcastem);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--pca\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--pca\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches kinship
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(kinstems[0],readstring)!=0&&strcmp(kinstems[0],readstring2)!=0)
{printf("Error, the pca %s corresponds to the kinship matrix %s, which appears to be different to that provided now (%s); if you are sure the kinship matrix is correct, use \"--check-root NO\"\n\n", pcastem, readstring, kinstems[0]);exit(1);}
}
fclose(input);
}

////////

if(strcmp(eigenfile,"blank")!=0)
{
sprintf(filename,"%s.eigen.id", eigenfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--decompose\"\n\n", filename);exit(1);}
sprintf(filename,"%s.eigen.bin", eigenfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--decompose\"\n\n", filename);exit(1);}

sprintf(filename,"%s.eigen.root", eigenfile);
if(just_check(filename)!=0)
{printf("Error reading %s; this file would have been created using \"--decompose\"\n\n", filename);exit(1);}
if((input=fopen(filename,"r"))==NULL)
{printf("Error opening %s\n\n",filename);exit(1);}
if(fscanf(input, "Kinship %s ", readstring)!=1)
{printf("Error reading Row 1 of %s (should begin \"Kinship\"), suggesting the file has been changed since creation with \"--decompose\"\n\n",filename);exit(1);}
if(checkroot==1)	//check matches kinship
{
(void)append_check(readstring2,readstring,workdir);
if(strcmp(kinstems[0],readstring)!=0&&strcmp(kinstems[0],readstring2)!=0)
{printf("Error, the eigen-decomposition %s corresponds to the kinship matrix %s, which appears to be different to that provided now (%s); if you are sure the eigen-decomposition is correct, use \"--check-root NO\"\n\n", eigenfile, readstring, kinstems[0]);exit(1);}
}
fclose(input);
}

///////////////////////////

//stats, scores, making phenotypes, snps, jackknifing, folds, find gaussian

if(strcmp(scorefile,"blank")!=0)	//check correct size and get num_scores (and maybe num_blocks and jackprop)
{
count=countcols(scorefile);
if(count<5)
{printf("Error, %s should have least five columns: \"Predictor\", \"A1\", \"A2\" and \"Centre\" (the average allele count; set to NA if unknown), followed by a column of effect sizes for each profile (not %d)\n\n", scorefile, count);exit(1);}
checkcols(scorefile,count);

if((input=fopen(scorefile,"r"))==NULL)
{printf("Error opening %s\n\n",scorefile);exit(1);}
if(fscanf(input, "%s %s %s %s", readstring, readstring2, readstring3, readstring4)!=4)
{printf("Error reading first four elements of %s\n\n", scorefile);exit(1);}

if(strcmp(readstring,"Predictor")!=0)
{printf("Error, the first column of %s should be named \"Predictor\" (not %s)", scorefile, readstring);exit(1);}
if(strcmp(readstring2,"A1")!=0)
{printf("Error, the second column of %s should be named \"A1\" (not %s)", scorefile, readstring2);exit(1);}
if(strcmp(readstring3,"A2")!=0)
{printf("Error, the third column of %s should be named \"A2\" (not %s)", scorefile, readstring3);exit(1);}
if(strcmp(readstring4,"Centre")!=0)
{printf("Error, the fourth column of %s should be named \"Centre\" (not %s)", scorefile, readstring4);exit(1);}
fclose(input);

num_scores=count-4;

if(prsvar==1)
{
if(num_scores<3){printf("Error, %s should contain at least four sets of effect sizes (not %d)\n\n", scorefile, num_scores);exit(1);}
num_blocks=num_scores-1;

//get jackprop
if((input=fopen(scorefile,"r"))==NULL)
{printf("Error opening %s\n\n",scorefile);exit(1);}
if(fscanf(input, "Predictor A1 A2 Centre Real Jackknife_%lf", &jackprop)!=1)
{printf("Error reading first six elements of %s (sixth element should take the form \"Jackknife_X\", where X is the jackknife proportion\n\n", scorefile);exit(1);}
fclose(input);

printf("jackprop is %f\n", jackprop);
}
}

if(strcmp(cofile,"blank")!=0)	//mode must be 122, 125 or 172
{
count=countrows(cofile);
count2=countcols(cofile);
checkcols(cofile,count2);
count3=count-(count2==4);

if(mode==122||mode==125)	//blup or pred - have already checked number of fixed effects 
{
if(count!=1+num_quants+num_envs&&count2!=4)
{printf("Error, %s should have %d rows and 4 columns (not %d and %d), suggesting the file has been changed since creation with \"--reml\"\n\n", cofile, 1+num_quants+num_envs, count, count2);exit(1);}
}
else	//so calc-scores
{
if(count2!=1&&count2!=4)
{printf("Error, %s should have either four columns, with Column 2 providing the coefficients (the format produced by \"--reml\"), or one column (not %d)\n\n", cofile, count2);exit(1);}

if(count3>1&&strcmp(covarfile,"blank")==0&&strcmp(envfile,"blank")==0)
{printf("Error, %s contains %d coefficients, so you must use \"--covars\" and/or \"--enviro\" to provide covariates and/or environmental variables\n\n", cofile, count3);exit(1);}
if(count3!=num_quants+num_envs)
{printf("Error, the number of coefficients in %s (%d) does not match the total number of covariates and environmental variables (%d)\n\n", cofile, count3, num_quants+num_envs);exit(1);}
}
}

if(strcmp(finalfile,"blank")!=0)	//check correct number of scores
{
count=countcols(finalfile);
if(count<5)
{printf("Error, %s should have least five columns: \"Predictor\", \"A1\", \"A2\" and \"Centre\" (the average allele count; set to NA if unknown), followed by a column of effect sizes for each profile (not %d)\n\n", finalfile, count);exit(1);}
if(count!=num_scores+4){printf("Error, %s should contain %d scores (not %d), to match the number in %s\n", finalfile, num_scores, count-4, scorefile);exit(1);} 
}

///////////////////////////


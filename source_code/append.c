/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Append slashes and prefix workdir to folder / file names and checks things exist

///////////////////////////

//if workdir has been specified, make sure it ends in "/", then check it exists

if(strcmp(workdir2,"blank")!=0)
{
count=strlen(workdir2);
if(workdir2[count-1]!='/'){sprintf(workdir,"%s/",workdir2);}
else{strcpy(workdir,workdir2);}

if((dir=opendir(workdir))==NULL)
{printf("Error, working directory %s does not exist; correct this first\n\n",workdir);exit(1);}
closedir(dir);
}
else
{strcpy(workdir,"");}

//make sure folder ends in "/", prefix with workdir (unless absolute reference), then check exists or can be made

if(strcmp(folder2,"blank")!=0)
{
count=strlen(folder2);

if(folder2[0]=='/')
{
if(folder2[count-1]!='/'){sprintf(folder,"%s/",folder2);}
else{strcpy(folder,folder2);}
}
else
{
if(folder2[count-1]!='/'){sprintf(folder,"%s%s/",workdir,folder2);}
else{sprintf(folder,"%s%s", workdir,folder2);}
}

if(stat(folder, &statstruct)!=0)	//folder does not exist (relies on folder ending in /)
{
if(mkdir(folder, S_IRWXU | S_IROTH)!=0)	//but can not make
{printf("Error, failed to create folder %s\nCheck write permissions or that a file of that name does not already exist\n\n",folder);exit(1);}
}
}
else
{strcpy(folder,"blank");}

//prefix outfile with workdir (unless absolute reference)
(void)append_check(outfile,outfile2,workdir);

///////////////////////////
//prefix other files (and for many, check whether exists)

(void)append_check(udatafile,udatafile2,workdir);
if(append_check(ufamfile,ufamfile2,workdir)!=0)
{printf("Error reading fam/sample file %s; check file exists and is not empty\n\n", ufamfile);exit(1);}
if(append_check(ubimfile,ubimfile2,workdir)!=0)
{printf("Error reading bim file %s; check file exists and is not empty\n\n", ubimfile);exit(1);}
(void)append_check(unamefile,unamefile2,workdir);
if(append_check(datalist,datalist2,workdir)!=0)
{printf("Error reading list of data files %s; check file exists and is not empty\n\n", datalist);exit(1);}

////////

if(append_check(bsampfile,bsampfile2,workdir)!=0)
{printf("Error reading keep file %s; check file exists and is not empty\n\n", bsampfile);exit(1);}
if(append_check(csampfile,csampfile2,workdir)!=0)
{printf("Error reading remove file %s; check file exists and is not empty\n\n", csampfile);exit(1);}
(void)append_check(subpref,subpref2,workdir);

if(append_check(bpredfile,bpredfile2,workdir)!=0)
{printf("Error reading extract file %s; check file exists and is not empty\n\n", bpredfile);exit(1);}
if(append_check(cpredfile,cpredfile2,workdir)!=0)
{printf("Error reading exclude file %s; check file exists and is not empty\n\n", cpredfile);exit(1);}

////////

if(append_check(centresfile,centresfile2,workdir)!=0)
{printf("Error reading predictor means file %s; check file exists and is not empty\n\n", centresfile);exit(1);}
if(append_check(weightsfile,weightsfile2,workdir)!=0)
{printf("Error reading weightings file %s; check file exists and is not empty\n\n", weightsfile);exit(1);}

if(append_check(pvafile,pvafile2,workdir)!=0)
{printf("Error reading p-values file %s; check file exists and is not empty\n\n", pvafile);exit(1);}
if(append_check(impfile,impfile2,workdir)!=0)
{printf("Error reading importances file %s; check file exists and is not empty\n\n", impfile);exit(1);}

////////

(void)append_check(kinname,kinname2,workdir);
if(append_check(kinlist,kinlist2,workdir)!=0)
{printf("Error reading list of kinships %s; check file exists and is not empty\n\n", kinlist);exit(1);}

(void)append_check(regpref,regpref2,workdir);

if(append_check(respfile,respfile2,workdir)!=0)
{printf("Error reading phenotype file %s; check file exists and is not empty\n\n", respfile);exit(1);}

if(append_check(sumsfile,sumsfile2,workdir)!=0)
{printf("Error reading summary file %s; check file exists and is not empty\n\n", sumsfile);exit(1);}
if(append_check(sums2file,sums2file2,workdir)!=0)
{printf("Error reading second summary file %s; check file exists and is not empty\n\n", sums2file);exit(1);}

if(append_check(covarfile,covarfile2,workdir)!=0)
{printf("Error reading covariates file %s; check file exists and is not empty\n\n", covarfile);exit(1);}
if(append_check(topfile,topfile2,workdir)!=0)
{printf("Error reading top predictors file %s; check file exists and is not empty\n\n", topfile);exit(1);}
if(append_check(envfile,envfile2,workdir)!=0)
{printf("Error reading environmental variables file %s; check file exists and is not empty\n\n", envfile);exit(1);}

if(append_check(offsetfile,offsetfile2,workdir)!=0)
{printf("Error reading offset file %s; check file exists and is not empty\n\n", offsetfile);exit(1);}

////////

if(append_check(infosfile,infosfile2,workdir)!=0)
{printf("Error reading info scores file %s; check file exists and is not empty\n\n", infosfile);exit(1);}

if(append_check(targetfile,targetfile2,workdir)!=0)
{printf("Error reading target predictors file %s; check file exists and is not empty\n\n", targetfile);exit(1);}

////////

(void)append_check(partpref,partpref2,workdir);

////////

if(append_check(malesfile,malesfile2,workdir)!=0)
{printf("Error reading males file %s; check file exists and is not empty\n\n", malesfile);exit(1);}

(void)append_check(invsfile,invsfile2,workdir);

////////

if(append_check(hersfile,hersfile2,workdir)!=0)
{printf("Error reading starting heritabilities file %s; check file exists and is not empty\n\n", hersfile);exit(1);}

if(append_check(oversfile,oversfile2,workdir)!=0)
{printf("Error reading overlaps file %s; check file exists and is not empty\n\n", oversfile);exit(1);}

////////

if(append_check(remlfile,remlfile2,workdir)!=0)
{printf("Error reading reml file %s; check file exists and is not empty\n\n", remlfile);exit(1);}

if(append_check(relfile,relfile2,workdir)!=0)
{printf("Error reading relationships file %s; check file exists and is not empty\n\n", relfile);exit(1);}

////////

if(append_check(sampwfile,sampwfile2,workdir)!=0)
{printf("Error reading sample weights file %s; check file exists and is not empty\n\n", sampwfile);exit(1);}

(void)append_check(locofile,locofile2,workdir);
(void)append_check(bocofile,bocofile2,workdir);

if(append_check(genefile,genefile2,workdir)!=0)
{printf("Error reading gene annotations file %s; check file exists and is not empty\n\n", genefile);exit(1);}

////////

(void)append_check(annpref,annpref2,workdir);
(void)append_check(labfile,labfile2,workdir);

if(append_check(printfile,printfile2,workdir)!=0)
{printf("Error reading regression predictors file %s; check file exists and is not empty\n\n", printfile);exit(1);}

if(append_check(herfile,herfile2,workdir)!=0)
{printf("Error reading heritability predictors file %s; check file exists and is not empty\n\n", herfile);exit(1);}

if(append_check(taglist,taglist2,workdir)!=0)
{printf("Error reading list of tagging files %s; check file exists and is not empty\n\n", taglist);exit(1);}
if(append_check(matlist,matlist2,workdir)!=0)
{printf("Error reading list of heritability matrices %s; check file exists and is not empty\n\n", matlist);exit(1);}

if(append_check(tagfile,tagfile2,workdir)!=0)
{printf("Error reading tagging file %s; check file exists and is not empty\n\n", tagfile);exit(1);}
if(append_check(catfile,catfile2,workdir)!=0)
{printf("Error reading reduce categories file %s; check file exists and is not empty\n\n", catfile);exit(1);}
if(append_check(altfile,altfile2,workdir)!=0)
{printf("Error reading alternative taggings file %s; check file exists and is not empty\n\n", altfile);exit(1);}
if(append_check(cvsfile,cvsfile2,workdir)!=0)
{printf("Error reading cross-validation predictors file %s; check file exists and is not empty\n\n", cvsfile);exit(1);}

if(append_check(powfile,powfile2,workdir)!=0)
{printf("Error reading powers file %s; check file exists and is not empty\n\n", powfile);exit(1);}

if(append_check(matfile,matfile2,workdir)!=0)
{printf("Error reading heritability matrix file %s; check file exists and is not empty\n\n", matfile);exit(1);}
if(append_check(taufile,taufile2,workdir)!=0)
{printf("Error reading tau file %s; check file exists and is not empty\n\n", taufile);exit(1);}

if(append_check(expfile,expfile2,workdir)!=0)
{printf("Error reading expectations file %s; check file exists and is not empty\n\n", expfile);exit(1);}

////////

if(append_check(indhers,indhers2,workdir)!=0)
{printf("Error reading per-predictor heritabilities file %s; check file exists and is not empty\n\n", indhers);exit(1);}

if(append_check(bvsfile,bvsfile2,workdir)!=0)
{printf("Error reading cross-validation samples file %s; check file exists and is not empty\n\n", bvsfile);exit(1);}

if(append_check(corslist,corslist2,workdir)!=0)
{printf("Error reading list of correlation files %s; check file exists and is not empty\n\n", corslist);exit(1);}

(void)append_check(corname,corname2,workdir);
(void)append_check(pseudostem,pseudostem2,workdir);

if(append_check(bestfile,bestfile2,workdir)!=0)
{printf("Error reading best model file %s; check file exists and is not empty\n\n", bestfile);exit(1);}

if(append_check(ldfile,ldfile2,workdir)!=0)
{printf("Error reading high-LD file %s; check file exists and is not empty\n\n", ldfile);exit(1);}

if(append_check(fracfile,fracfile2,workdir)!=0)
{printf("Error reading fractions file %s; check file exists and is not empty\n\n", fracfile);exit(1);}

////////

(void)append_check(pcastem,pcastem2,workdir);

(void)append_check(eigenfile,eigenfile2,workdir);

////////

if(append_check(scorefile,scorefile2,workdir)!=0)
{printf("Error reading score file %s; check file exists and is not empty\n\n", scorefile);exit(1);}
if(append_check(cofile,cofile2,workdir)!=0)
{printf("Error reading coefficients file %s; check file exists and is not empty\n\n", cofile);exit(1);}
if(append_check(finalfile,finalfile2,workdir)!=0)
{printf("Error reading final effects file %s; check file exists and is not empty\n\n", finalfile);exit(1);}

if(append_check(probsfile,probsfile2,workdir)!=0)
{printf("Error reading probabilties file %s; check file exists and is not empty\n\n", probsfile);exit(1);}
if(append_check(causalsfile,causalsfile2,workdir)!=0)
{printf("Error reading causals file %s; check file exists and is not empty\n\n", causalsfile);exit(1);}
if(append_check(effectsfile,effectsfile2,workdir)!=0)
{printf("Error reading effects file %s; check file exists and is not empty\n\n", effectsfile);exit(1);}

if(append_check(predlista,predlista2,workdir)!=0)
{printf("Error reading predictor list %s; check file exists and is not empty\n\n", predlista);exit(1);}
if(append_check(predlistb,predlistb2,workdir)!=0)
{printf("Error reading predictor list %s; check file exists and is not empty\n\n", predlistb);exit(1);}

if(append_check(jackfile,jackfile2,workdir)!=0)
{printf("Error reading datapoint pairs file %s; check file exists and is not empty\n\n", jackfile);exit(1);}
if(append_check(proffile,proffile2,workdir)!=0)
{printf("Error reading profile file %s; check file exists and is not empty\n\n", proffile);exit(1);}

if(append_check(likefile,likefile2,workdir)!=0)
{printf("Error reading likelihoods file %s; check file exists and is not empty\n\n", likefile);exit(1);}

////////

(void)append_check(greout,greout2,workdir);

///////////////////////////


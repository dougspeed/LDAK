/*
Copyright 2026 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//set parameters for prediction models (megaprs)

///////////////////////////

//set num_try
if(strcmp(bestfile,"blank")!=0)	//considering a single model
{
num_try=1;
}
else
{
if(strcmp(fracfile,"blank")==0)	//defaults depend on mode
{
if(ptype==0)	//mega - do lasso, ridge, bayesr (minus ridge) and elastic
{num_try=45;}
if(ptype==1)	//lasso-sparse
{num_try=80;}
if(ptype==2)	//lasso
{num_try=1;}
if(ptype==3)	//ridge
{num_try=1;}
if(ptype==4)	//bolt
{num_try=18;}
if(ptype==5)	//bayesr
{
if(exparam==0){num_try=35;}
else{num_try=56;}
}
if(ptype==6)	//bayesr-shrink
{
if(exparam==0){num_try=35;}
else{num_try=56;}
}
if(ptype==7)	//elastic
{num_try=9;}
if(ptype==8)	//ldpred
{num_try=9;}
}
else	//will set models based on values in fracfile
{
num_try=countrows(fracfile);
}
}

//allocate parameters
trytypes=malloc(sizeof(int)*num_try);
trylams=malloc(sizeof(double)*num_try);
tryscales=malloc(sizeof(double)*num_try);
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

//set scale to 1, then remaining values to NA (then will change when required)
for(p=0;p<num_try;p++)
{
trytypes[p]=-9999;trylams[p]=-9999;tryscales[p]=1;
tryps[p]=-9999;tryp2s[p]=-9999;tryp3s[p]=-9999;tryp4s[p]=-9999;tryf2s[p]=-9999;
}

if(strcmp(bestfile,"blank")!=0)	//considering a single model
{
//open bestfile, skip the header row and get model
if((input=fopen(bestfile,"r"))==NULL)
{printf("Error opening %s\n\n",bestfile);exit(1);}
readchar=0;while(readchar!=10){readchar=10;(void)fscanf(input, "%c", &readchar);}
if(fscanf(input, "%d %s ", &readint, readstring)!=2)
{printf("Error reading Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
if(strcmp(readstring,"lasso-sparse")!=0&&strcmp(readstring,"lasso")!=0&&strcmp(readstring,"ridge")!=0&&strcmp(readstring,"bolt")!=0 &&strcmp(readstring,"bayesr")!=0&&strcmp(readstring,"bayesr-shrink")!=0&&strcmp(readstring,"elastic")!=0)
{printf("Error, Element 2 of Row 2 of %s should be  \"lasso-sparse\", \"lasso\", \"ridge\", \"bolt\", \"bayesr\", \"bayesr-shrink\" or \"--elastic\" (not %s), suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile, readstring);exit(1);}

if(strcmp(readstring,"lasso-sparse")==0)	//lasso-sparse - read lambda and s
{
trytypes[0]=1;
if(fscanf(input, "%lf %lf ", trylams, tryscales)!=2){printf("Error reading lambda and s from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

//if(strcmp(readstring,"lasso")==0)	nothing to read
//if(strcmp(readstring,"ridge")==0)	nothing to read

if(strcmp(readstring,"bolt")==0)	//bolt - read p and f2 (and set p2/p3)
{
trytypes[0]=4;
if(fscanf(input, "NA NA %lf %lf ", tryps, tryf2s)!=2){printf("Error reading p and f2 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
tryp2s[0]=1-tryps[0];
}

if(strcmp(readstring,"bayesr")==0)	//bayesr-sparse - read p, p2, p3 and p4
{
trytypes[0]=5;
if(fscanf(input, "NA NA NA NA %lf %lf %lf %lf ", tryps, tryp2s, tryp3s, tryp4s)!=4){printf("Error reading p, p2, p3 and p4 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"bayesr-shrink")==0)	//bayesr-shrink - read p, p2, p3 and p4
{
trytypes[0]=6;
if(fscanf(input, "NA NA NA NA %lf %lf %lf %lf ", tryps, tryp2s, tryp3s, tryp4s)!=4){printf("Error reading p, p2, p3 and p4 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

if(strcmp(readstring,"elastic")==0)	//elastic - read p and f2
{
trytypes[0]=7;
if(fscanf(input, "NA NA %lf %lf ", &readdouble, tryf2s)!=2){printf("Error reading p and f2 from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
tryps[0]=readdouble/2;
tryp2s[0]=readdouble/2;
tryp3s[0]=1-readdouble;
}

if(strcmp(readstring,"LDpred")==0)	//LDpred - read p
{
trytypes[0]=8;
if(fscanf(input, "NA NA %lf ", tryps)!=1){printf("Error reading p from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
tryp2s[0]=1-tryps[0];
}

if(strcmp(readstring,"CS")==0)	//CS - read lambda (phi)
{
trytypes[0]=9;
if(fscanf(input, "%lf ", trylams)!=1){printf("Error reading lambda from Row 2 of %s, suggesting the file has been changed since creation with \"--mega-prs\"\n\n", bestfile);exit(1);}
}

fclose(input);
}
else	//not bestfile
{
if(strcmp(fracfile,"blank")==0)	//using defaults
{
count=0;

if(ptype==0)	//mega - using lasso, ridge, bayesr (minus ridge) and elastic
{
//lasso and ridge
trytypes[count]=2;count++;
trytypes[count]=3;count++;

//bayesr (squeezed down)
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;

for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryp4s[count]=loads[j];tryp3s[count]=loads[j2];
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
while(tryp4s[count]==0){tryp4s[count]=tryp3s[count];tryp3s[count]=tryp2s[count];tryp2s[j]=0;}
}
}}}
free(loads);

//elastic
loads=malloc(sizeof(double)*3);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.1;loads[2]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=7;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}	//end of ptype=0

if(ptype==1)	//lasso-sparse
{
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.9;trylams[count]=0.001*pow(100,(double)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.5;trylams[count]=0.001*pow(100,(double)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.2;trylams[count]=0.001*pow(100,(double)p/19);count++;}
for(p=0;p<20;p++){trytypes[count]=1;tryscales[count]=0.1;trylams[count]=0.001*pow(100,(double)p/19);count++;}
}

if(ptype==2)	//lasso
{
trytypes[count]=2;count++;
}

if(ptype==3)	//ridge
{
trytypes[count]=3;count++;
}

if(ptype==4)	//bolt
{
loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
for(j=0;j<6;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=4;tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}	//end of ptype=4

if(ptype==5)	//bayesr (first is ridge)
{
if(exparam==0)
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;

trytypes[0]=5;tryps[0]=0;tryp2s[0]=0;tryp3s[0]=0;tryp4s[0]=1;
count=1;

for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryp4s[count]=loads[j]/100;tryp3s[count]=loads[j2]/10;
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}
else
{
loads=malloc(sizeof(double)*6);
loads[0]=0;loads[1]=.001;loads[2]=.01;loads[3]=.05;loads[4]=.1;loads[5]=.2;

trytypes[0]=5;tryps[0]=0;tryp2s[0]=0;tryp3s[0]=0;tryp4s[0]=1;
count=1;

for(j=0;j<6;j++){
for(j2=j;j2<6;j2++){
for(j3=j2;j3<6;j3++){
if(j+j2+j3>0)
{
trytypes[count]=5;tryp4s[count]=loads[j]/100;tryp3s[count]=loads[j2]/10;
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}

//squeeze down the models (make sure tryps4 non-zero)
for(j=0;j<count;j++)
{
while(tryp4s[j]==0)
{tryp4s[j]=tryp3s[j];tryp3s[j]=tryp2s[j];tryp2s[j]=0;}
}
}	//end of ptype=5

if(ptype==6)	//bayes-shrink (first is ridge)
{
if(exparam==0)
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;

trytypes[0]=6;tryps[0]=0;tryp2s[0]=0;tryp3s[0]=0;tryp4s[0]=1;
count=1;

for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
trytypes[count]=6;tryp4s[count]=loads[j]/100;tryp3s[count]=loads[j2]/10;
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}
else
{
loads=malloc(sizeof(double)*6);
loads[0]=0;loads[1]=.001;loads[2]=.01;loads[3]=.05;loads[4]=.1;loads[5]=.2;

trytypes[0]=6;tryps[0]=0;tryp2s[0]=0;tryp3s[0]=0;tryp4s[0]=1;

count=1;
for(j=0;j<6;j++){
for(j2=j;j2<6;j2++){
for(j3=j2;j3<6;j3++){
if(j+j2+j3>0)
{
trytypes[count]=6;tryp4s[count]=loads[j]/100;tryp3s[count]=loads[j2]/10;
tryp2s[count]=loads[j3];tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);
}

//squeeze down the models (make sure tryps4 non-zero)
for(j=0;j<count;j++)
{
while(tryp4s[j]==0)
{tryp4s[j]=tryp3s[j];tryp3s[j]=tryp2s[j];tryp2s[j]=0;}
}
}	//end of ptype=6

if(ptype==7)	//elastic
{
loads=malloc(sizeof(double)*3);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.1;loads[2]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
trytypes[count]=7;tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}	//end of ptype=7

if(ptype==8)	//ldpred
{
trytypes[0]=8;tryps[0]=1.0;tryp2s[0]=0.0;
trytypes[1]=8;tryps[1]=0.3;tryp2s[1]=0.7;
trytypes[2]=8;tryps[2]=0.1;tryp2s[2]=0.9;
trytypes[3]=8;tryps[3]=0.03;tryp2s[3]=0.97;
trytypes[4]=8;tryps[4]=0.01;tryp2s[4]=0.99;
trytypes[5]=8;tryps[5]=0.003;tryp2s[5]=0.997;
trytypes[6]=8;tryps[6]=0.001;tryp2s[6]=0.999;
trytypes[7]=8;tryps[7]=0.0003;tryp2s[7]=0.9997;
trytypes[8]=8;tryps[8]=0.0001;tryp2s[8]=0.9999;
count=9;
}	//end of ptype=8

if(count!=num_try){printf("Doug error %d and %d\n", num_try, count);exit(1);}
}
else	//read fracfile - must be model 4, 5, 6, 7, 8 or 9
{
for(p=0;p<num_try;p++){trytypes[p]=ptype;}

if(ptype==4)	//bolt - read p and f2 (also set p2)
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<=0||tryps[p]>1){printf("Error, probability in Row %d of %s (%.4f) is not within (0,1])\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryp2s[p]=1-tryps[p];}

read_values(fracfile, tryf2s, num_try, NULL, 2, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]==1&&tryf2s[p]!=0){printf("Error, the probability in Row %d of %s is one, so the fraction must be zero (not %.4f)\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
if(tryf2s[p]<0||tryf2s[p]>1){printf("Error, fraction in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
}
}

if(ptype==5||ptype==6)	//bayesr or bayesr-shrink - read p, p2, p3 and p4
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
read_values(fracfile, tryp2s, num_try, NULL, 2, 0, 0);
read_values(fracfile, tryp3s, num_try, NULL, 3, 0, 0);
read_values(fracfile, tryp4s, num_try, NULL, 4, 0, 0);

for(p=0;p<num_try;p++)
{
if(tryps[p]>=1){printf("Error, the first probability in %s must be less than one (not %f)\n\n", fracfile, tryps[p]);exit(1);}
sum=tryps[p]+tryp2s[p]+tryp3s[p]+tryp4s[p];
if(sum<0.99||sum>1.01)
{printf("Error, probabilities in Row %d of %s sum to %.4f (not one)\n\n", p+1, fracfile, sum);exit(1);}
tryps[p]=tryps[p]/sum;
tryp2s[p]=tryp2s[p]/sum;
tryp3s[p]=tryp3s[p]/sum;
tryp4s[p]=tryp4s[p]/sum;
}
}

if(ptype==7)	//elastic - read 2p and f2 (also set p2 and p3)
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<0||tryps[p]>1){printf("Error, probability in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryps[p]=tryps[p]/2;tryp2s[p]=tryps[p];tryp3s[p]=1-tryps[p]-tryp2s[p];}

read_values(fracfile, tryf2s, num_try, NULL, 3, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryf2s[p]<0||tryf2s[p]>1){printf("Error, fraction in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
if(tryps[p]==1&&tryf2s[p]>0){printf("Error in Row %d of %s; can not have p=1 and f2>0\n\n", p+1, fracfile);exit(1);}
if(tryps[p]==0&&tryf2s[p]<1){printf("Error in Row %d of %s; can not have p=0 and f2<1\n\n", p+1, fracfile);exit(1);}
}
}

if(ptype==8)	//ldpred
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<0||tryps[p]>1){printf("Error, probability in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryp2s[p]=1-tryps[p];}
}
}	//end of reading fracfile
}	//end of not using best model

//save parameters
if(num_focals==1){sprintf(filename,"%s.parameters",outfile);}
else{sprintf(filename,"%s.focal%d.parameters", outfile, cur_focal+1);}
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

fprintf(output, "Model\tType\tLasso_lambda\tLasso_s\tBolt_p\tBolt_f2\tBayesR_p1\tBayesR_p2\tBayesR_p3\tBayesR_p4\tElastic_p\tElastic_f2\n");
for(p=0;p<num_try;p++)
{
if(trytypes[p]==1){fprintf(output, "%d\tlasso-sparse\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, trylams[p], tryscales[p]);}
if(trytypes[p]==2){fprintf(output, "%d\tlasso\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1);}
if(trytypes[p]==3){fprintf(output, "%d\tridge\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1);}
if(trytypes[p]==4){fprintf(output, "%d\tbolt\tNA\tNA\t%.4f\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryps[p], tryf2s[p]);}
if(trytypes[p]==5){fprintf(output, "%d\tbayesr\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\tNA\tNA\n", p+1, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==6){fprintf(output, "%d\tbayesr-shrink\tNA\tNA\tNA\tNA\t%.4f\t%.4f\t%.4f\t%.4f\tNA\tNA\n", p+1, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
if(trytypes[p]==7){fprintf(output, "%d\telastic\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%.4f\t%.4f\n", p+1, 1-tryp3s[p], tryf2s[p]);}
if(trytypes[p]==8){fprintf(output, "%d\tLDpred\tNA\tNA\t%.4f\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", p+1, tryps[p]);}
}
fclose(output);

///////////////////////////


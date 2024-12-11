/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//set parameters for prediction models (for ind-level data) - for LOCO, require first model to be ridge

///////////////////////////

//now parameters

/*
if(highstruct==1)
{
num_try=1;
if(strcmp(fracfile,"blank")!=0){printf("Warning, will ignore the parameters provided in %s (due to switching to ridge regression because of strong structure)\n\n", fracfile);}
}
else
*/
{
if(strcmp(fracfile,"blank")==0)	//defaults depend on mode
{
if(mode==151)	//ridge
{num_try=1;}
if(mode==152)	//bolt
{
if(ldpred==0){num_try=18;}
else{num_try=9;}
}
if(mode==153)	//bayesr
{num_try=35;}
if(mode==154)	//elastic
{num_try=10;}
}
else	//will set models based on values in fracfile
{
num_try=countrows(fracfile);
}
}

//allocate parameters
tryps=malloc(sizeof(double)*num_try);
tryp2s=malloc(sizeof(double)*num_try);
tryp3s=malloc(sizeof(double)*num_try);
tryp4s=malloc(sizeof(double)*num_try);
tryf2s=malloc(sizeof(double)*num_try);

/*
if(highstruct==1)	//just need to save a blank parameters file
{
//set values to NA, then change when required
tryps[0]=-9999;tryp2s[0]=-9999;tryp3s[0]=-9999;tryp4s[0]=-9999;tryf2s[0]=-9999;

//if(mode==151)	//nothing to set (and probably will not be here)

if(mode==152){tryps[0]=0.5;tryp2s[0]=0.5;tryf2s[0]=0.5;}

if(mode==153){tryps[0]=0;tryp2s[0]=0;tryp3s[0]=0;tryp4s[0]=1;}

if(mode==154){tryps[0]=0.0;tryp2s[0]=0.0;tryp3s[0]=1.0;tryf2s[0]=1.0;}
}
else
*/
{
//set values to NA, then change when required
for(p=0;p<num_try;p++){tryps[p]=-9999;tryp2s[p]=-9999;tryp3s[p]=-9999;tryp4s[p]=-9999;tryf2s[p]=-9999;}

if(strcmp(fracfile,"blank")==0)	//using defaults
{
//if(mode==151)	//ridge - default is 1 model - no parameters to set

if(mode==152)	//bolt - default is 18 models (first is ridge)
{
if(ldpred==0)
{
loads=malloc(sizeof(double)*6);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.2;loads[2]=.1;loads[3]=.05;loads[4]=.02;loads[5]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;

count=0;
for(j=0;j<6;j++)
{
for(j2=0;j2<3;j2++)
{
tryps[count]=loads[j];tryp2s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}
else
{
tryps[0]=1.0;tryp2s[0]=0.0;tryf2s[0]=0.0;
tryps[1]=0.3;tryp2s[1]=0.7;tryf2s[1]=0.0;
tryps[2]=0.1;tryp2s[2]=0.9;tryf2s[2]=0.0;
tryps[3]=0.03;tryp2s[3]=0.97;tryf2s[3]=0.0;
tryps[4]=0.01;tryp2s[4]=0.99;tryf2s[4]=0.0;
tryps[5]=0.003;tryp2s[5]=0.997;tryf2s[5]=0.0;
tryps[6]=0.001;tryp2s[6]=0.999;tryf2s[6]=0.0;
tryps[7]=0.0003;tryp2s[7]=0.9997;tryf2s[7]=0.0;
tryps[8]=0.0001;tryp2s[8]=0.9999;tryf2s[8]=0.0;
}
}

if(mode==153)	//bayesr - default is 35 models (first is ridge)
{
loads=malloc(sizeof(double)*5);
loads[0]=0;loads[1]=.01;loads[2]=.05;loads[3]=.1;loads[4]=.2;

tryps[0]=0;tryp2s[0]=0;tryp3s[0]=0;tryp4s[0]=1;

count=1;
for(j=0;j<5;j++){
for(j2=j;j2<5;j2++){
for(j3=j2;j3<5;j3++){
if(j+j2+j3>0)
{
tryp4s[count]=loads[j];tryp3s[count]=loads[j2];tryp2s[count]=loads[j3];
tryps[count]=1-tryp2s[count]-tryp3s[count]-tryp4s[count];count++;
}
}}}
free(loads);

//squeeze down the models (make sure tryps4 non-zero)
for(j=0;j<count;j++)
{
while(tryp4s[j]==0)
{tryp4s[j]=tryp3s[j];tryp3s[j]=tryp2s[j];tryp2s[j]=0;}
}
}

if(mode==154)	//elastic - default is 10 models (first is ridge)
{
loads=malloc(sizeof(double)*3);
loads2=malloc(sizeof(double)*3);
loads[0]=.5;loads[1]=.1;loads[2]=.01;
loads2[0]=.5;loads2[1]=.3;loads2[2]=.1;

tryps[0]=0.0;tryp2s[0]=0.0;tryp3s[0]=1.0;tryf2s[0]=1.0;

count=1;
for(j=0;j<3;j++)
{
for(j2=0;j2<3;j2++)
{
tryps[count]=loads[j]/2;tryp2s[count]=loads[j]/2;
tryp3s[count]=1-loads[j];tryf2s[count]=loads2[j2];count++;
}}
free(loads);free(loads2);
}
}
else	//read fracfile (must have loco=0, so do not need first model ridge)
{
//if(mode==151)	//ridge - not allowed fracfile

if(mode==152)	//bolt - read p and f2 (also set p2)
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryps[p]<=0||tryps[p]>=1){printf("Error, probability in Row %d of %s (%.4f) is not within (0,1)\n\n", p+1, fracfile, tryps[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryp2s[p]=1-tryps[p];}

read_values(fracfile, tryf2s, num_try, NULL, 2, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryf2s[p]<0||tryf2s[p]>1){printf("Error, fraction in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
}
}

if(mode==153)	//bayesr or bayesr-shrink - read p, p2, p3 and p4
{
read_values(fracfile, tryps, num_try, NULL, 1, 0, 0);
read_values(fracfile, tryp2s, num_try, NULL, 2, 0, 0);
read_values(fracfile, tryp3s, num_try, NULL, 3, 0, 0);
read_values(fracfile, tryp4s, num_try, NULL, 4, 0, 0);

for(p=0;p<num_try;p++)
{
sum=tryps[p]+tryp2s[p]+tryp3s[p]+tryp4s[p];
if(sum<0.99||sum>1.01)
{printf("Error, probabilities in Row %d of %s sum to %.4f (not one)\n\n", p+1, fracfile, sum);exit(1);}
tryps[p]=tryps[p]/sum;
tryp2s[p]=tryp2s[p]/sum;
tryp3s[p]=tryp3s[p]/sum;
tryp4s[p]=tryp4s[p]/sum;
}
}

if(mode==154)	//elastic - read 2p and f2 (also set p2 and p3)
{
loads=malloc(sizeof(double)*num_try);
read_values(fracfile, loads, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(loads[p]<0||loads[p]>1){printf("Error, probability in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, loads[p]);exit(1);}
}
for(p=0;p<num_try;p++){tryps[p]=loads[p]/2;tryp2s[p]=loads[p]/2;tryp3s[p]=1-loads[p];}

read_values(fracfile, tryf2s, num_try, NULL, 2, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryf2s[p]<0||tryf2s[p]>1){printf("Error, fraction in Row %d of %s (%.4f) is not within [0,1]\n\n", p+1, fracfile, tryf2s[p]);exit(1);}
if(tryps[p]==1&&tryf2s[p]>0){printf("Error in Row %d of %s; can not have p=1 and f2>0\n\n", p+1, fracfile);exit(1);}
if(tryps[p]==0&&tryf2s[p]<1){printf("Error in Row %d of %s; can not have p=0 and f2<1\n\n", p+1, fracfile);exit(1);}
}
free(loads);
}
}	//end of reading fracfile
}

if(verbose==1)	//save parameters
{
sprintf(filename,"%s.parameters",outfile);
if((output=fopen(filename,"w"))==NULL)
{printf("Error writing to %s; check you have permission to write and that there does not exist a folder with this name\n\n",filename);exit(1);}

if(mode==151)
{
fprintf(output, "Model\n");
for(p=0;p<num_try;p++){fprintf(output, "%d\n", p+1);}
}
if(mode==152)
{
fprintf(output, "Model\tp\tf2\n");
for(p=0;p<num_try;p++){fprintf(output, "%d\t%.4f\t%.4f\n", p+1, tryps[p], tryf2s[p]);}
}
if(mode==153)
{
fprintf(output, "Model\tp1\tp2\tp3\tp4\n");
for(p=0;p<num_try;p++){fprintf(output, "%d\t%.4f\t%.4f\t%.4f\t%.4f\n", p+1, tryps[p], tryp2s[p], tryp3s[p], tryp4s[p]);}
}
if(mode==154)
{
fprintf(output, "Model\tp\tf2\n");
for(p=0;p<num_try;p++){fprintf(output, "%d\t%.4f\t%.4f\n", p+1, 1-tryp3s[p], tryf2s[p]);}
}
fclose(output);
}

///////////////////////////

/*
old code for reading her scaling
read_values(fracfile, tryhers, num_try, NULL, 1, 0, 0);
for(p=0;p<num_try;p++)
{
if(tryhers[p]<=0){printf("Error, heritability scaling in Row %d of %s is not positive(%.4f) \n\n", p+1, fracfile, tryhers[p]);exit(1);}
if(tryhers[p]>1){printf("Warning, heritability scaling in Row %d of %s is greater than one (%.4f)\n\n", p+1, fracfile, tryhers[p]);exit(1);}
}
*/


/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Compute REML derivatives

///////////////////////////

//want AI (minus second derivatives) and BI (first derivatives)
//require YTPKPY, YTPKPKPY and tr(PK) for each kinship or pair of kinships
//sweights must be 1 when shortcut=1

if(Xtotal>0)	//will need PX and XTPY - must have shortcut 0 or 1 (for shortcut=1, actually get UTPX and XTPY)
{
if(shortcut==0)	//get both directly
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &Xtotal, &ns, &alpha, P, &ns, X, &ns, &beta, PX, &ns);
dgemv_("T", &ns, &Xtotal, &alpha, PX, &ns, Y, &one, &beta, XTPY, &one);
}

if(shortcut==1)	//PX = (invD - DUTXF - BUTZH) UTX, can get XTPY directly
{
for(j=0;j<Xtotal;j++)
{
for(i=0;i<ns;i++){PX[i+j*ns]=UTX[i+j*ns]/D[i];}
}
alpha=1.0;beta=0.0;
dgemm_("N", "N", &Xtotal, &Xtotal, &ns, &alpha, F, &Xtotal, UTX, &ns, &beta, FUTX, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &Xtotal, &Xtotal, &alpha, DUTX, &ns, FUTX, &Xtotal, &beta, PX, &ns);
alpha=1.0;beta=0.0;
dgemm_("N", "N", &num_fixed, &Xtotal, &ns, &alpha, H, &num_fixed, UTX, &ns, &beta, HUTX, &num_fixed);
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &Xtotal, &num_fixed, &alpha, BUTZ, &ns, HUTX, &num_fixed, &beta, PX, &ns);
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &Xtotal, &alpha, PX, &ns, UTY, &one, &beta, XTPY, &one);
}
}

if(shortcut==2)	//get PR = invVR - invVZ (ZTVZ) ZTVR - already have VZZTVZ
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_vects, &ns, &alpha, Z, &ns, VR, &ns, &beta, ZTVR, &num_fixed);

for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){PR[i+g*ns]=VR[i+g*ns];}
}
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &num_vects, &num_fixed, &alpha, VZZTVZ, &ns, ZTVR, &num_fixed, &beta, PR, &ns);
}

////////

//get KPYs and trace PKs (for shortcut 2, latter is E(t(PR) KR)

//get KPY and trace PK for noise (i.e., for K = invW)

for(i=0;i<ns;i++){KPY[0][i]=PY[i]/sweights[i];}

if(shortcut==0)
{
traces[0]=0;
for(i=0;i<ns;i++){traces[0]+=P[(size_t)i*ns+i]/sweights[i];}
}

if(shortcut==1)	//trace P = trace (invD - DUTXF - BUTZH)
{
traces[0]=0;
for(i=0;i<ns;i++)
{
traces[0]+=1.0/D[i];
for(j=0;j<Xtotal;j++){traces[0]-=DUTX[i+j*ns]*F[j+i*Xtotal];}
for(j=0;j<num_fixed;j++){traces[0]-=BUTZ[i+j*ns]*H[j+i*num_fixed];}
}
}

if(shortcut==2)	//trace P invW = E(t(PR) * invW R)
{
sum=0;
for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){sum+=PR[i+g*ns]*R_single[i+g*ns]/sweights[i];}
}
traces[0]=sum/num_vects;
}

if(shortcut==3)	//trace P invW = trace (invD - BUTZH)
{
traces[0]=0;
for(i=0;i<ns;i++)
{
traces[0]+=1.0/D[i]/sweights[i];
for(j=0;j<num_fixed;j++){traces[0]-=BUTZ[i+j*ns]*H[j+i*num_fixed]/sweights[i];}
}
}

////////

if(shortcut==2&&memsave==1)	//will be using V_single to read kinships, so save its values in V_full
{
i2=0;j2=0;
for(j=0;j<ns;j++)
{
for(i=0;i<=j;i++)
{
V_full[(size_t)j2*ns+i2]=V_single[(size_t)j*ns+i];
i2++;
if(i2==ns){j2++;i2=j2;}
}
}
}

for(k=0;k<num_kins;k++)	//get KPYs and traces for kinships
{
if(shortcut==0)
{
if(memsave==0)	//kins already stored
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, Mkins[k], &ns, PY, &one, &beta, KPY[1+k], &one);

traces[1+k]=0;
for(scount=0;scount<stotal;scount++){traces[1+k]+=P[scount]*Mkins[k][scount];}
}
else	//must read kins - will save space by saving them in V
{
read_kins(kinstems[k], V, NULL, 1.0, ns, ids3, 2, maxthreads);

alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, V, &ns, PY, &one, &beta, KPY[1+k], &one);

traces[1+k]=0;
for(scount=0;scount<stotal;scount++){traces[1+k]+=P[scount]*V[scount];}
}
}

if(shortcut==1)	//must have k=0 and tr(PK) = tr(U(...)UTUEUT) = tr(E(invD - DUTXF - BUTZH))
{
for(i=0;i<ns;i++){KPY[1][i]=E[i]*PY[i];}

traces[1]=0;
for(i=0;i<ns;i++)
{
traces[1]+=E[i]/D[i];
for(j=0;j<Xtotal;j++){traces[1]-=E[i]*DUTX[i+j*ns]*F[j+i*Xtotal];}
for(j=0;j<num_fixed;j++){traces[1]-=E[i]*BUTZ[i+j*ns]*H[j+i*num_fixed];}
}
}

if(shortcut==2)
{
//fill PY_single
for(i=0;i<ns;i++){PY_single[i]=PY[i];}

//get KPY_single
if(memsave==0)	//kins already stored
{
alpha_single=1.0;beta_single=0.0;
ssymv_("U", &ns, &alpha_single, Mkins_single[k], &ns, PY_single, &one, &beta_single, KPY_single, &one);
}
else	//must read kins
{
read_kins(kinstems[k], NULL, V_single, 1.0, ns, ids3, 4, maxthreads);

alpha_single=1.0;beta_single=0.0;
ssymv_("U", &ns, &alpha_single, V_single, &ns, PY_single, &one, &beta_single, KPY_single, &one);
}

//fill KPY
for(i=0;i<ns;i++){KPY[1+k][i]=KPY_single[i];}

//trace = E(t(PR) * KR)
sum=0;
for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){sum+=PR[i+g*ns]*KR_single[k][i+g*ns];}
}
traces[1+k]=sum/num_vects;
}

if(shortcut==3)	//tr(PK) = tr(diag(K) (invD - BUTZH))
{
for(i=0;i<ns;i++){KPY[1+k][i]=kin_diags[i+k*ns]*PY[i];}

traces[1+k]=0;
for(i=0;i<ns;i++)
{
traces[1+k]+=kin_diags[i+k*ns]/D[i];
for(j=0;j<num_fixed;j++){traces[1+k]-=kin_diags[i+k*ns]*BUTZ[i+j*ns]*H[j+i*num_fixed];}
}
}

}

for(r=0;r<num_regs;r++)	//get KPYs and traces for regions - must have shortcut 0 or 1
{
token=Xends[r]-Xstarts[r];
if(shortcut==0)
{
alpha=1.0/Xsums[r];beta=0.0;
dgemv_("N", &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, XTPY+Xstarts[r], &one, &beta, KPY[1+num_kins+r], &one);

traces[1+num_kins+r]=0;
for(i=0;i<ns;i++)
{
for(j=Xstarts[r];j<Xends[r];j++){traces[1+num_kins+r]+=PX[i+j*ns]*X[i+j*ns]/Xsums[r];}
}
}

if(shortcut==1)
{
alpha=1.0/Xsums[r];beta=0.0;
dgemv_("N", &ns, &token, &alpha, UTX+Xstarts[r]*ns, &ns, XTPY+Xstarts[r], &one, &beta, KPY[1+num_kins+r], &one);

traces[1+num_kins+r]=0;
for(i=0;i<ns;i++)
{
for(j=Xstarts[r];j<Xends[r];j++){traces[1+num_kins+r]+=PX[i+j*ns]*UTX[i+j*ns]/Xsums[r];}
}
}
}	//end of r loop

if(shortcut==2&&memsave==1)	//restore V_single
{
i2=0;j2=0;
for(j=0;j<ns;j++)
{
for(i=0;i<=j;i++)
{
V_single[(size_t)j*ns+i]=V_full[(size_t)j2*ns+i2];
i2++;
if(i2==ns){j2++;i2=j2;}
}
}
}

////////

//get PKPYs

if(shortcut==2)	//first get invV KPYs
{
for(k=0;k<total;k++)
{
for(i=0;i<ns;i++){RHS2_single[i+k*ns]=KPY[k][i];}
}

if(ldlt==0)	//using cholesky
{
spotrs_("U", &ns, &total, V_single, &ns, RHS2_single, &ns, &info);
if(info!=0){printf("Error second Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, ns);exit(1);}
}
else	//using ldlt
{
ssytrs_("U", &ns, &total, V_single, &ns, ipiv, RHS_single, &ns, &info);
if(info!=0){printf("Error first LDLT solve failed, please tell Doug (info %d, length %d)\n\n", info, ns);exit(1);}
}

//fill up VKPY
for(k=0;k<total;k++)
{
for(i=0;i<ns;i++){VKPY[i+k*ns]=RHS2_single[i+k*ns];}
}
}

for(k=0;k<total;k++)
{
if(shortcut==0)
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, KPY[k], &one, &beta, PKPY[k], &one);
}

if(shortcut==1)	//get (invD - DUTXF - BUTZH) KPY
{
//set PKPY = invD KPY
for(i=0;i<ns;i++){PKPY[k][i]=KPY[k][i]/D[i];}
if(Xtotal>0)	//using X, so subtract DUTXF KPY
{
alpha=1.0;beta=0.0;
dgemv_("N", &Xtotal, &ns, &alpha, F, &Xtotal, KPY[k], &one, &beta, FKPY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &Xtotal, &alpha, DUTX, &ns, FKPY, &one, &beta, PKPY[k], &one);
}
//subtract BUTZH KPY
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &ns, &alpha, H, &num_fixed, KPY[k], &one, &beta, HKPY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, BUTZ, &ns, HKPY, &one, &beta, PKPY[k], &one);
}

if(shortcut==2)	//PKPY = invV KPY - invVZ invZTVZ ZT invVPKY
{
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, VKPY+k*ns, &one, &beta, ZTVKPY, &one);

for(i=0;i<ns;i++){PKPY[k][i]=VKPY[i+k*ns];}
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, VZZTVZ, &ns, ZTVKPY, &one, &beta, PKPY[k], &one);
}

if(shortcut==3)	//get (invD - BUTZH) KPY
{
//set PKPY = invD KPY
for(i=0;i<ns;i++){PKPY[k][i]=KPY[k][i]/D[i];}
//subtract BUTZH KPY
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &ns, &alpha, H, &num_fixed, KPY[k], &one, &beta, HKPY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, BUTZ, &ns, HKPY, &one, &beta, PKPY[k], &one);
}
}	//end of k loop

///////////////////////////

//fill AI = -2nd deriv = .5*YTPKPKPY and BI = 1st deriv = .5*(YTPKPKY -trace(PK))

//start assuming all terms fixed (AI diag, BI=0)
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){AI[k+k2*total]=0;}
AI[k+k*total]=1;BI[k]=0;
}

//now fill up
for(k=0;k<total;k++)
{
if(fixed[k]<3)
{
//diagonals and B
gam2=0;for(i=0;i<ns;i++){gam2+=KPY[k][i]*PY[i];}
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[k][i]*KPY[k][i];}
AI[k+k*total]=.5*gam3;
BI[k]=.5*(gam2-traces[k]);

//off-diagonals
for(k2=0;k2<k;k2++)
{
if(fixed[k2]<3)
{
gam3=0;for(i=0;i<ns;i++){gam3+=PKPY[k][i]*KPY[k2][i];}
AI[k+k2*total]=.5*gam3;
AI[k2+k*total]=.5*gam3;
}}	//end of k2 free and k2 loop
}}	//end of k free and k loop

//for stability, scale AI so has trace one
sum=0;for(k=0;k<total;k++){sum+=AI[k+k*total];}
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){AI[k+k2*total]*=pow(sum,-1);}
}

//invert then set fixed to zero
(void)eigen_invert(AI, total, AI2, -1, AI3, 1);
for(k=0;k<total;k++)
{
if(fixed[k]>=3){AI[k+k*total]=0;}
}

//undo scaling of AI
for(k=0;k<total;k++)
{
for(k2=0;k2<total;k2++){AI[k+k2*total]*=pow(sum,-1);}
}

///////////////////////////


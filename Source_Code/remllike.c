/*
Copyright 2024 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

//Compute REML likelihood

///////////////////////////

//V has form var0 invW + var1 K1 + var2 K2 + ... + var(1+nk) X1X1T + var(1+nk+1) X2X2T + ... 
//when shortcut=1, will have nk=0 or nk=1 and write V as var0 I + var1 UEUT + XCXT
//sweights must be 1 when shortcut=1

if(shortcut==0)	//get invV directly
{
//load up V (only use upper triangle)
for(scount=0;scount<stotal;scount++){V[scount]=0;}
for(i=0;i<ns;i++){V[(size_t)i*ns+i]=vars[0]/sweights[i];}

for(k=0;k<num_kins;k++)	//add on varK for kinships
{
if(vars[1+k]!=0)
{
if(memsave==0)	//have kinships saved in Mkins
{
for(scount=0;scount<stotal;scount++){V[scount]+=vars[1+k]*Mkins[k][scount];}
}
else{read_kins(kinstems[k], V, NULL, vars[1+k], ns, ids3, 5, maxthreads);}
}
}

for(r=0;r<num_regs;r++)	//and for regions
{
if(vars[1+num_kins+r]!=0)
{
token=Xends[r]-Xstarts[r];
alpha=vars[1+num_kins+r]/Xsums[r];beta=1.0;
dgemm_("N", "T", &ns, &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, V, &ns);
}
}

//invert V
detV=cholesky_invert(V, ns, -1, NULL, &info, 0);
if(info!=0)	//try ldlt
{
printf("Cholesky decomposition failed, switching to LDLT decomposition\n");

//reload V
for(scount=0;scount<stotal;scount++){V[scount]=0;}
for(i=0;i<ns;i++){V[(size_t)i*ns+i]=vars[0]/sweights[i];}

for(k=0;k<num_kins;k++)	//add on varK for kinships
{
if(vars[1+k]!=0)
{
if(memsave==0)	//have kinships saved in Mkins
{
for(scount=0;scount<stotal;scount++){V[scount]+=vars[1+k]*Mkins[k][scount];}
}
else{read_kins(kinstems[k], V, NULL, vars[1+k], ns, ids3, 5, maxthreads);}
}
}

for(r=0;r<num_regs;r++)	//and for regions
{
if(vars[1+num_kins+r]!=0)
{
token=Xends[r]-Xstarts[r];
alpha=vars[1+num_kins+r]/Xsums[r];beta=1.0;
dgemm_("N", "T", &ns, &ns, &token, &alpha, X+Xstarts[r]*ns, &ns, X+Xstarts[r]*ns, &ns, &beta, V, &ns);
}
}

detV=ldlt_invert(V, ns, -1, NULL, &info, 1);
if(info!=0)	//another fail
{printf("Error, both Cholesky and LDLT decompositions failed, unable to continue; please tell Doug\n\n");exit(1);}

//would previously also try an eigendecomposition (saving space by using P as workspace)
//must first remake V then run detV=eigen_invert(V, ns, V2, -1, P, 1);
}
}	//end of shortcut=0

if(shortcut==1)	//get invV = U (invD - invD UTX F) UT, where F= (t(UTX) invDUTX + invC)^-1 t(invDUTX)
{
//find D and detD where var0 I + var1 K1 = UDUT (D = vars0 I or D = vars0 I + vars1 E)
detD=0;
for(i=0;i<ns;i++)
{
if(num_kins==1){D[i]=vars[0]+vars[1]*E[i];}
else{D[i]=vars[0];}
if(D[i]!=0){detD+=log(fabs(D[i]));}
}

if(Xtotal>0)	//using X, so find F
{
//get DUTX = invD UTX
for(r=0;r<num_regs;r++)
{
for(j=Xstarts[r];j<Xends[r];j++)
{
if(vars[1+num_kins+r]!=0)	//fill column
{
for(i=0;i<ns;i++){DUTX[i+j*ns]=UTX[i+j*ns]/D[i];}
}
else	//blank column
{
for(i=0;i<ns;i++){DUTX[i+j*ns]=0;}
}
}}

//get XTVX = t(UTX) DUTX + invC, where C = vars/Xsums (recall for regions with var=0, DUTX is blank)
//detV = |D + UTX invC XTU| = |D| |C| |XUT invD UTX + invC| = |D| |C| |XTVX|
alpha=1.0;beta=0.0;
dgemm_("T", "N", &Xtotal, &Xtotal, &ns, &alpha, UTX, &ns, DUTX, &ns, &beta, XTVX, &Xtotal);

detC=0;
for(r=0;r<num_regs;r++)
{
for(j=Xstarts[r];j<Xends[r];j++)
{
if(vars[1+num_kins+r]!=0)	//add invC to diagonal and contribution to detC
{
XTVX[j+j*Xtotal]+=Xsums[r]/vars[1+num_kins+r];
detC+=log(fabs(vars[1+num_kins+r]/Xsums[r]));
}
else	//blank row except for 1 on diagonal (no contribution to detC)
{
for(j2=0;j2<Xtotal;j2++){XTVX[j+j2*Xtotal]=0;}
XTVX[j+j*Xtotal]=1;
}
}}	//end of j and r loops

//get F = invXTVX t(DUTX) - save XTVX in case ldlt fails
for(j=0;j<Xtotal;j++)
{
for(j2=0;j2<Xtotal;j2++){XTVX3[j+j2*Xtotal]=XTVX[j+j2*Xtotal];}
}

//first set F = t(DUTX), then solve
for(i=0;i<ns;i++)
{
for(j=0;j<Xtotal;j++){F[j+i*Xtotal]=DUTX[i+j*ns];}
}
detXTVX=ldlt_invert(XTVX, Xtotal, ns, F, &info, 1);
if(info!=0)	//inversion failed, so use eigen (XTVX saved in XTVX3, F will still equal DUTX)
{detXTVX=eigen_invert(XTVX3, Xtotal, XTVX2, ns, F, 1);}
detV=detD+detC+detXTVX;
}	//end of using X
else	//not using X, so F not used and V = UDUT
{detV=detD;}
}	//end of shortcut=1

if(shortcut==2)	//get invV [Z, Y, vect] - will only use invV vect later
{
//load up V_single (only use upper triangle)
for(scount=0;scount<stotal;scount++){V_single[scount]=0;}
for(i=0;i<ns;i++){V_single[(size_t)i*ns+i]=vars[0]/sweights[i];}

for(k=0;k<num_kins;k++)	//add on varK for kinships
{
if(vars[1+k]!=0)
{
if(memsave==0)	//kins already stored
{
for(scount=0;scount<stotal;scount++){V_single[scount]+=vars[1+k]*Mkins_single[k][scount];}
}
else{read_kins(kinstems[k], NULL, V_single, vars[1+k], ns, ids3, 5, maxthreads);}
}
}

//decompose V
if(ldlt==0)	//using cholesky
{
spotrf_("U", &ns, V_single, &ns, &info);
if(info!=0)
{printf("Error, Cholesky decomposition failed; to switch to a (slower) LDLT decomposion, restart adding \"--LDLT YES\"\n\n");exit(1);}
}
else	//using ldlt
{
ipiv=malloc(sizeof(int)*ns);
lwork=-1;
ssytrf_("U", &ns, V_single, &ns, ipiv, &wkopt_single, &lwork, &info);
if(info!=0)
{printf("Error, LDLT priming failed; please tell Doug (info %d, length %d)\n\n", info, ns);exit(1);}
lwork=(int)wkopt_single;
work_single=malloc(sizeof(float)*lwork);
ssytrf_("U", &ns, V_single, &ns, ipiv, work_single, &lwork, &info);
free(work_single);
if(info!=0)
{printf("Error, LDLT decomposition failed, please tell Doug (info %d, length %d)\n\n", info, ns);exit(1);}
}

//get detV - for cholesky, must square, for ldlt, dont
detV=0;
for(i=0;i<ns;i++)
{
if(fabs(V_single[(size_t)i*ns+i])>=0.000001){detV+=(2-ldlt)*log(fabs(V_single[(size_t)i*ns+i]));}
else{detV+=(2-ldlt)*log(0.000001);}
}

//fill RHS_single with Z, Y and R
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){RHS_single[i+j*ns]=Z[i+j*ns];}
}
for(i=0;i<ns;i++){RHS_single[i+num_fixed*ns]=Y[i];}
for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){RHS_single[i+(1+num_fixed+g)*ns]=R_single[i+g*ns];}
}

//get invV RHS_single
if(ldlt==0)	//using cholesky
{
spotrs_("U", &ns, &total2, V_single, &ns, RHS_single, &ns, &info);
if(info!=0){printf("Error first Cholesky solve failed, please tell Doug (info %d, length %d)\n\n", info, ns);exit(1);}
}
else	//using ldlt
{
ssytrs_("U", &ns, &total2, V_single, &ns, ipiv, RHS_single, &ns, &info);
if(info!=0){printf("Error first LDLT solve failed, please tell Doug (info %d, length %d)\n\n", info, ns);exit(1);}
}

//fill VZ, VY and VR from RHS_single
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){VZ[i+j*ns]=RHS_single[i+j*ns];}
}
for(i=0;i<ns;i++){VY[i]=RHS_single[i+num_fixed*ns];}
for(g=0;g<num_vects;g++)
{
for(i=0;i<ns;i++){VR[i+g*ns]=RHS_single[i+(1+num_fixed+g)*ns];}
}
}	//end of shortcut=2

if(shortcut==3)	//get invV = invD, where D = var0 I + var1 diag(K1) + ...
{
//find D and detV
detV=0;
for(i=0;i<ns;i++)
{
D[i]=vars[0]/sweights[i];
for(k=0;k<num_kins;k++){D[i]+=vars[1+k]*kin_diags[i+k*ns];}
if(D[i]!=0){detV+=log(fabs(D[i]));}
}
}

////////

//ultimately want PY, where P = invV - invVZ invZTVZ ZTinvV or U (invD-DUTXF-BUTZH) UT, where BUTZ=UTinvV UTZ and H = invZTVZ t(BUTZ)

//get invZTVZ

if(shortcut==0)	//get invVZ then ZTVZ = ZT invVZ
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_fixed, &ns, &alpha, V, &ns, Z, &ns, &beta, VZ, &ns);
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, VZ, &ns, &beta, ZTVZ, &num_fixed);
}

if(shortcut==1)	//get BUTZ = UTinvV UTZ, then ZTVZ = ZTU BUTZ
{
//set BUTZ = invD UTZ
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){BUTZ[i+j*ns]=UTZ[i+j*ns]/D[i];}
}
if(Xtotal>0)	//using X, so subtract DUTZ F UTZ
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &Xtotal, &num_fixed, &ns, &alpha, F, &Xtotal, UTZ, &ns, &beta, FUTZ, &Xtotal);
alpha=-1.0;beta=1.0;
dgemm_("N", "N", &ns, &num_fixed, &Xtotal, &alpha, DUTX, &ns, FUTZ, &Xtotal, &beta, BUTZ, &ns);
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, UTZ, &ns, BUTZ, &ns, &beta, ZTVZ, &num_fixed);
}

if(shortcut==2)	//already have VZ
{
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, VZ, &ns, &beta, ZTVZ, &num_fixed);
}

if(shortcut==3)	//get BUTZ = invV Z, then ZTVZ = ZT BUTZ
{
//set BUTZ = invD Z
for(j=0;j<num_fixed;j++)
{
for(i=0;i<ns;i++){BUTZ[i+j*ns]=Z[i+j*ns]/D[i];}
}
alpha=1.0;beta=0.0;
dgemm_("T", "N", &num_fixed, &num_fixed, &ns, &alpha, Z, &ns, BUTZ, &ns, &beta, ZTVZ, &num_fixed);
}

//invert ZTVZ
detZTVZ=eigen_invert(ZTVZ, num_fixed, ZTVZ2, -1, ZTVZ3, 1);

if(shortcut==0)	//P = invV - invVZ invZTVZ ZTinvV
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_fixed, &num_fixed, &alpha, VZ, &ns, ZTVZ, &num_fixed, &beta, VZZTVZ, &ns);

for(scount=0;scount<stotal;scount++){P[scount]=V[scount];}
alpha=-1.0;beta=1.0;
dgemm_("N", "T", &ns, &ns, &num_fixed, &alpha, VZZTVZ, &ns, VZ, &ns, &beta, P, &ns);
}

if(shortcut==1)	//P = U (invD-DUTXF-BUTZH) UT, where H = invZTVZ t(BUTZ)
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_fixed, &ns, &num_fixed, &alpha, ZTVZ, &num_fixed, BUTZ, &ns, &beta, H, &num_fixed);
}

if(shortcut==2)	//only get invVZ invZTVZ
{
alpha=1.0;beta=0.0;
dgemm_("N", "N", &ns, &num_fixed, &num_fixed, &alpha, VZ, &ns, ZTVZ, &num_fixed, &beta, VZZTVZ, &ns);
}

if(shortcut==3)	//P = invD-BUTZH, where H = invZTVZ t(BUTZ)
{
alpha=1.0;beta=0.0;
dgemm_("N", "T", &num_fixed, &ns, &num_fixed, &alpha, ZTVZ, &num_fixed, BUTZ, &ns, &beta, H, &num_fixed);
}

////////

//get PY (for shortcut=1, actually get UTPY) and gam = YTPY

if(shortcut==0)	//get direct
{
alpha=1.0;beta=0.0;
dgemv_("N", &ns, &ns, &alpha, P, &ns, Y, &one, &beta, PY, &one);
gam=0;for(i=0;i<ns;i++){gam+=PY[i]*Y[i];}
}

if(shortcut==1) //get UTPY = (invD - DUTXF - BUTZH) UTY
{
//set PY = invD UTY
for(i=0;i<ns;i++){PY[i]=UTY[i]/D[i];}
if(Xtotal>0)	//using X, so subtract DUTXF UTY
{
alpha=1.0;beta=0.0;
dgemv_("N", &Xtotal, &ns, &alpha, F, &Xtotal, UTY, &one, &beta, FUTY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &Xtotal, &alpha, DUTX, &ns, FUTY, &one, &beta, PY, &one);
}
//subtract BUTZH UTY
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &ns, &alpha, H, &num_fixed, UTY, &one, &beta, HUTY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, BUTZ, &ns, HUTY, &one, &beta, PY, &one);

gam=0;for(i=0;i<ns;i++){gam+=PY[i]*UTY[i];}
}

if(shortcut==2)	//PY = invV Y - invVZ invZTVZ ZTinvV Y
{
alpha=1.0;beta=0.0;
dgemv_("T", &ns, &num_fixed, &alpha, Z, &ns, VY, &one, &beta, ZTVY, &one);

for(i=0;i<ns;i++){PY[i]=VY[i];}
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, VZZTVZ, &ns, ZTVY, &one, &beta, PY, &one);

gam=0;for(i=0;i<ns;i++){gam+=PY[i]*Y[i];}
}

if(shortcut==3) //get PY = (invD - BUTZH) Y
{
//set PY = invD Y
for(i=0;i<ns;i++){PY[i]=Y[i]/D[i];}
//subtract BUTZH UTY
alpha=1.0;beta=0.0;
dgemv_("N", &num_fixed, &ns, &alpha, H, &num_fixed, Y, &one, &beta, HUTY, &one);
alpha=-1.0;beta=1.0;
dgemv_("N", &ns, &num_fixed, &alpha, BUTZ, &ns, HUTY, &one, &beta, PY, &one);

gam=0;for(i=0;i<ns;i++){gam+=PY[i]*Y[i];}
}

///////////////////////////


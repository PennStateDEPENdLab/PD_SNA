*code written by Kimberly Nolf in 6/2013
*creates PD scores for the SIDP-IV used in the couples study protocol

MISSING VALUES 
negtv1,negtv2,negtv3,negtv4,negtv5,negtv6,negtv7,szoid1,szoid2,szoid3,szoid4,
szoid5stypl8,szoid6,szoid7,histr1,histr2,histr3,histr4,histr5,histr6,histr7,histr8,narci1,narci2,narci3,narci4,narci5,
narci6,narci7,narci8,narci9,antso1,antso2,antso3,antso4,antso5,antso6,antso7,cantso,bordl1,bordl2,bordl3,bordl4,
bordl5,bordl6,bordl7,bordl8,bordl9,avoid1,avoid2,avoid3,avoid4,avoid5,avoid6,avoid7,depen1,depen2,depen3,depen4,
depen5,depen6,depen7,depen8,obcmp1,obcmp2,obcmp3,obcmp4,obcmp5,obcmp6,obcmp7,obcmp8,parnd1,parnd2,
parnd3,parnd4,parnd5,parnd6,parnd7,stypl1,stypl2,stypl3,stypl4,stypl5,stypl6,stypl7,szoid5stypl8,stypl9,deprs1,
deprs2,deprs3,deprs4,deprs5,deprs6,deprs7(99) .



compute negtv_sidp=sum(negtv1,negtv2,negtv3,negtv4,negtv5,negtv6,negtv7) .
compute szoid_sidp=sum(szoid1,szoid2,szoid3,szoid4,szoid5stypl8,szoid6,szoid7) .
compute histr_sidp=sum(histr1,histr2,histr3,histr4,histr5,histr6,histr7,histr8) .
compute narci_sidp=sum(narci1,narci2,narci3,narci4,narci5,narci6,narci7,narci8,narci9) .
compute antso_sidp=sum(antso1,antso2,antso3,antso4,antso5,antso6,antso7) .
compute bordl_sidp=sum(bordl1,bordl2,bordl3,bordl4,bordl5,bordl6,bordl7,bordl8,bordl9) .
compute avoid_sidp=sum(avoid1,avoid2,avoid3,avoid4,avoid5,avoid6,avoid7) .
compute depen_sidp=sum(depen1,depen2,depen3,depen4,depen5,depen6,depen7,depen8) .
compute obcmp_sidp=sum(obcmp1,obcmp2,obcmp3,obcmp4,obcmp5,obcmp6,obcmp7,obcmp8) .
compute parnd_sidp=sum(parnd1,parnd2,parnd3,parnd4,parnd5,parnd6,parnd7) .
compute stypl_sidp=sum(stypl1,stypl2,stypl3,stypl4,stypl5,stypl6,stypl7,szoid5stypl8,stypl9) .
compute deprs_sidp=sum(deprs1,deprs2,deprs3,deprs4,deprs5,deprs6,deprs7) .
compute OPD_sidp=sum(negtv_sidp,szoid_sidp,histr_sidp,narci_sidp,antso_sidp,avoid_sidp,depen_sidp,obcmp_sidp,parnd_sidp,stypl_sidp,deprs_sidp) .
EXECUTE .

COMPUTE NumMiss = NMISS(negtv1,negtv2,negtv3,negtv4,negtv5,negtv6,negtv7,szoid1,szoid2,szoid3,szoid4,
szoid5stypl8,szoid6,szoid7,histr1,histr2,histr3,histr4,histr5,histr6,histr7,histr8,narci1,narci2,narci3,narci4,narci5,
narci6,narci7,narci8,narci9,antso1,antso2,antso3,antso4,antso5,antso6,antso7,cantso,bordl1,bordl2,bordl3,bordl4,
bordl5,bordl6,bordl7,bordl8,bordl9,avoid1,avoid2,avoid3,avoid4,avoid5,avoid6,avoid7,depen1,depen2,depen3,depen4,
depen5,depen6,depen7,depen8,obcmp1,obcmp2,obcmp3,obcmp4,obcmp5,obcmp6,obcmp7,obcmp8,parnd1,parnd2,
parnd3,parnd4,parnd5,parnd6,parnd7,stypl1,stypl2,stypl3,stypl4,stypl5,stypl6,stypl7,
szoid5stypl8,stypl9,deprs1,deprs2,deprs3,deprs4,deprs5,deprs6,deprs7) .
EXECUTE .

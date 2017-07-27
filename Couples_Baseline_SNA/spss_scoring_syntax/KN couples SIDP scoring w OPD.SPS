*code written by Kimberly Nolf in 6/2013
*creates PD scores for the SIDP-IV used in the couples study protocol

Numeric UsrID (F12.0).
COMPUTE UsrID=(PTNUM * 10) + DyadID.
VARIABLE LABELS  UsrID 'combination of PTNUM/dyad'.
EXECUTE.


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

VARIABLE LABELS
negtv_sidp 'negtv dimensional score'
szoid_sidp 'szoid dimensional score'
histr_sidp 'histr dimensional score'
narci_sidp 'narci dimensional score'
antso_sidp 'antso dimensional score'
bordl_sidp 'bordl dimensional score'
avoid_sidp 'avoid dimensional score'
depen_sidp 'depen dimensional score'
obcmp_sidp 'obcmp dimensional score'
parnd_sidp 'parnd dimensional score'
stypl_sidp 'stypl dimensional score'
deprs_sidp 'deprs dimensional score'
OPD_sidp 'OtherPD or dimensional score for all PDs except BPD' .
EXECUTE.

COMPUTE NumMiss = NMISS(negtv1,negtv2,negtv3,negtv4,negtv5,negtv6,negtv7,szoid1,szoid2,szoid3,szoid4,
szoid5stypl8,szoid6,szoid7,histr1,histr2,histr3,histr4,histr5,histr6,histr7,histr8,narci1,narci2,narci3,narci4,narci5,
narci6,narci7,narci8,narci9,antso1,antso2,antso3,antso4,antso5,antso6,antso7,cantso,bordl1,bordl2,bordl3,bordl4,
bordl5,bordl6,bordl7,bordl8,bordl9,avoid1,avoid2,avoid3,avoid4,avoid5,avoid6,avoid7,depen1,depen2,depen3,depen4,
depen5,depen6,depen7,depen8,obcmp1,obcmp2,obcmp3,obcmp4,obcmp5,obcmp6,obcmp7,obcmp8,parnd1,parnd2,
parnd3,parnd4,parnd5,parnd6,parnd7,stypl1,stypl2,stypl3,stypl4,stypl5,stypl6,stypl7,
szoid5stypl8,stypl9,deprs1,deprs2,deprs3,deprs4,deprs5,deprs6,deprs7) .
EXECUTE .


COUNT negtvCount=negtv1 negtv2 negtv4 negtv7 negtv5 negtv3 negtv6(1,2,3).
VARIABLE LABELS  negtvCount 'negtv symptom count'.
EXECUTE.
COUNT szoidCount=szoid1 szoid2 szoid3 szoid4 szoid5stypl8 szoid6 szoid7(1,2,3).
VARIABLE LABELS szoidCount 'szoid symptom count'.
Execute.
COUNT histrCount=histr1 histr2 histr3 histr4 histr5 histr6 histr7 histr8(1,2,3).
VARIABLE Labels histrCount 'histr symptom count'.
EXECUTE.
COUNT narciCount=narci1 narci2 narci3 narci4 narci5 narci6 narci7 narci8 narci9(1,2,3).
VARIABLE LABELS narciCount 'narci symptom count'.
EXECUTE.
COUNT antsoCount=antso1 antso2 antso3 antso4 antso5 antso6 antso7(1,2,3).
VARIABLE LABELS antsoCount 'antso symptom count'.
EXECUTE.
COUNT bordlCount=bordl1 bordl2 bordl3 bordl4 bordl5 bordl6 bordl7 bordl8 bordl9(1,2,3).
VARIABLE LABELS bordlCount 'bordl symptom count'.
EXECUTE.
COUNT avoidCount=avoid1 avoid2 avoid3 avoid4 avoid5 avoid6 avoid7(1,2,3).
VARIABLE LABELS avoidCount 'avoid symptom count'.
EXECUTE.
COUNT depenCount=depen1 depen2 depen3 depen4 depen5 depen6 depen7 depen8(1,2,3).
VARIABLE LABELS depenCount 'depen symptom count'.
EXECUTE.
COUNT obcmpCount=obcmp1 obcmp2 obcmp3 obcmp4 obcmp5 obcmp6 obcmp7 obcmp8(1,2,3).
VARIABLE LABELS obcmpCount 'obcmp symptom count'.
EXECUTE.
COUNT parndCount=parnd1 parnd2 parnd3 parnd4 parnd5 parnd6 parnd7(1,2,3).
VARIABLE LABELS parndCount 'parnd symptom count'.
EXECUTE.
COUNT styplCount=stypl1 stypl2 stypl3 stypl4 stypl5 stypl6 stypl7 szoid5stypl8,stypl9(1,2,3).
VARIABLE LABELS styplCount 'stypl symptom count'.
EXECUTE.
COUNT deprsCount=deprs1,deprs2,deprs3,deprs4,deprs5,deprs6,deprs7(1,2,3).
VARIABLE LABELS deprsCount 'deprs symptom count'.
EXECUTE.

RECODE slfinj1 slfinj2 slfinj3 slfinj4 sectionQpatient sectionQinterviewer (' '='8888').
EXECUTE.

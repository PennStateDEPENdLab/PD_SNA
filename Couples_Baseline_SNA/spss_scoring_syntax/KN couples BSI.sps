*code written by kimberly nolf for the couples study 

compute bsiGSI=MEAN.18(bsi1 to bsi18).
execute.

compute bsiSOM=MEAN.6(bsi1,bsi4,bsi7,bsi10,bsi13,bsi16).
execute.

compute bsiDEP=MEAN.6(bsi2,bsi5,bsi8,bsi11,bsi14,bsi17).
execute.

compute bsiANX=MEAN.6(bsi3,bsi6,bsi9,bsi12,bsi15,bsi18).
execute.

recode bsiGSI bsiSOM bsiDEP bsiANX (SYSMIS=9999) .
execute .

variable labels bsiGSI 'BSI Global Symptom Index' bsiSOM 'BSI Somatization' bsiANX 'BSI Anxiety' bsiDEP 'BSI Depression'.


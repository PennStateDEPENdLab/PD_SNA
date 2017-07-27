*code originally written by Jennifer Morse and editted by Kimberly Nolf in 3.2014 for use with the Couples study data. 

COMPUTE IIP_PD1=MEAN.10(iip1,iip35,iip36,iip42,iip51,iip55,iip60,iip78,iip79,iip81,iip86) .

COMPUTE IIP_PD2=MEAN.9(iip13,iip14,iip26,iip28,iip32,iip34,iip38,iip40,iip41,iip84) .

COMPUTE IIP_PD3=MEAN.6(iip50,iip53,iip58,iip63,iip77,iip80,iip88) .

COMPUTE IIP_C1=MEAN.8(iip2,iip9,iip16,iip48,iip59,iip66,iip72,iip74,iip75) .

COMPUTE IIP_C2=MEAN.9(iip3,iip7,iip17,iip19,iip22,iip33,iip43,iip49,iip71,iip85) .

COMPUTE IIP_PD=MEAN(IIP_PD1,IIP_PD2,IIP_PD3) .

COMPUTE IIP_C=mean(IIP_C1,IIP_C2) .

COMPUTE iip_pa=mean(iip21,iip40,iip57,iip58,iip65,iip68,iip76,iip80) .

COMPUTE iip_bc=mean(iip1,iip26,iip28,iip38,iip41,iip50,iip73,iip88) .

COMPUTE iip_de=mean(iip11,iip18,iip20,iip24,iip27,iip31,iip46,iip82) .

COMPUTE iip_fg=mean(iip3,iip7,iip17,iip22,iip43,iip45,iip71,iip85) .

COMPUTE iip_hi=mean(iip5,iip6,iip8,iip9,iip12,iip15,iip23,iip49) .

COMPUTE iip_jk=mean(iip2,iip10,iip29,iip44,iip48,iip54,iip69,iip83) .

COMPUTE iip_lm=mean(iip25,iip37,iip47,iip59,iip64,iip67,iip70,iip87) .

COMPUTE iip_no=mean(iip4,iip30,iip39,iip52,iip56,iip61,iip62,iip78) .

COMPUTE iip_bpd=mean(iip51,iip53,iip55,iip66,iip77,iip80,iip89,iip90) .
EXECUTE .

VARIABLE LABELS
 iip_pa 'IIP domineering' iip_bc 'IIP vindictive' iip_de 'IIP cold' iip_fg 'IIP socially avoidant'
 iip_hi 'IIP nonassertive' iip_jk 'IIP exploitable' iip_lm 'IIP overly nurturant' iip_no 'IIP intrusive' iip_bpd 'Clifton BPD scale' .

IF (IIP_PD >=1.1) PD_on_IIP = 1 .
EXECUTE .

DO IF (~ MISSING(IIP_PD)) .
RECODE
  PD_on_IIP  (SYSMIS=0)  .
END IF .
EXECUTE .

COMPUTE IIPSumValue = SUM.80(iip1,iip2,iip3,iip4,iip5,iip6,iip7,iip8,iip9,iip10,iip11,iip12,iip13,iip14,iip15,iip16,iip17,iip18
 ,iip19,iip20,iip21,iip22,iip23,iip24,iip25,iip26,iip27,iip28,iip29,iip30,iip31,iip32,iip33,iip34,iip35,iip36,iip37,iip38
 ,iip39,iip40,iip41,iip42,iip43,iip44,iip45,iip46,iip47,iip49,iip48,iip50,iip51,iip52,iip53,iip54,iip55,iip56,iip57,iip58
 ,iip59,iip60,iip61,iip62,iip63,iip64,iip65,iip66,iip67,iip68,iip69,iip70,iip71,iip72,iip73,iip74,iip75,iip76,iip77,iip78
 ,iip79,iip80,iip81,iip82,iip83,iip84,iip85,iip86,iip87,iip88, iip89, iip90) .
EXECUTE .

COMPUTE NumMiss = NMISS(iip1,iip2,iip3,iip4,iip5,iip6,iip7,iip8,iip9,iip10,iip11,iip12,iip13,iip14,iip15,iip16,iip17,iip18
 ,iip19,iip20,iip21,iip22,iip23,iip24,iip25,iip26,iip27,iip28,iip29,iip30,iip31,iip32,iip33,iip34,iip35,iip36,iip37,iip38
 ,iip39,iip40,iip41,iip42,iip43,iip44,iip45,iip46,iip47,iip49,iip48,iip50,iip51,iip52,iip53,iip54,iip55,iip56,iip57,iip58
 ,iip59,iip60,iip61,iip62,iip63,iip64,iip65,iip66,iip67,iip68,iip69,iip70,iip71,iip72,iip73,iip74,iip75,iip76,iip77,iip78
 ,iip79,iip80,iip81,iip82,iip83,iip84,iip85,iip86,iip87,iip88, iip89, iip90) .
EXECUTE .

variable labels IIPSumValue 'IIP Sum' PD_on_IIP 'PD on IIP'.
EXECUTE.
value labels pd_on_iip 0 'no' 1 'yes'.
EXECUTE.

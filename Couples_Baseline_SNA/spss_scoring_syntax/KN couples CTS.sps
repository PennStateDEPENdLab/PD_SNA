*CTS scoring created by Kimberly Nolf on 3.12.14 to score the CTS used in the couples study.

GET DATA
  /TYPE=ODBC
  /CONNECT='DSN=Pilkonis_Couples;Description=Pilkonis_Couples;UID=;Trusted_Connection=Yes;APP='+
    'IBM SPSS Products: Statistics Common;WSID=OACB124402-W7;DATABASE=Pilkonis_Couples;Network='+
    'DBMSSOCN'
  /SQL='SELECT PTNUM, dyadID, mth, CTS1, CTS2, CTS3, CTS4, CTS5, CTS6, CTS7, CTS8, CTS9, CTS10, '+
    'CTS11, CTS12, CTS13, CTS14, CTS15, CTS16, CTS17, CTS18, CTS19, CTS20, CTS21, CTS22, CTS23, '+
    'CTS24, CTS25, CTS26, CTS27, CTS28, CTS29, CTS30, CTS31, CTS32, CTS33, CTS34, CTS35, CTS36, '+
    'CTS37, CTS38, CTS39, CTS40, UsrID FROM Pilkonis_Couples.Assessments.[CTS]'
  /ASSUMEDSTRWIDTH=255.
CACHE.
EXECUTE.

RECODE CTS1 to CTS40
(0=0) (1=1) (2=2) (3=4) (4=8) (5=15) (6=25) (7=98) .
EXECUTE .

MISSING VALUES CTS1 TO CTS40 (98) .
EXECUTE .


*DEFINITIONS:  v=victim, p=perpetrator, m=minor, s=severe

compute PsychAggV = sum.7(CTS1,CTS11,CTS15,CTS19,CTS27,CTS33,CTS35,CTS37) .
execute .
compute PsychAggP = sum.7(CTS2,CTS12,CTS16,CTS20,CTS28,CTS34,CTS36,CTS38) .
execute .

compute PhysAssV = sum.11(CTS3,CTS5,CTS7,CTS9,CTS13,CTS17,CTS21,CTS23,CTS25,CTS29,CTS31,CTS39) .
execute .
compute PhysAssP = sum.11(CTS4,CTS6,CTS8,CTS10,CTS14,CTS18,CTS22,CTS24,CTS26,CTS30,CTS32,CTS40) .
execute .

compute PsychAggVM = sum.4(CTS1,CTS19,CTS27,CTS35) .
execute .
compute PsychAggPM = sum.4(CTS2,CTS20,CTS28,CTS36) .
execute .

compute PsychAggVS = sum.4(CTS11,CTS15,CTS33,CTS37) .
execute .
compute PsychAggPS = sum.4(CTS12,CTS16,CTS34,CTS38).
execute .

compute PhysAssVM = sum.5(CTS3,CTS5,CTS7,CTS25,CTS29) .
execute .
compute PhysAssPM = sum.5(CTS4,CTS6,CTS8,CTS26,CTS30) .
execute .

compute PhysAssVS = sum.7(CTS9,CTS13,CTS17,CTS21,CTS23,CTS31,CTS39) .
execute .
compute PhysAssPS = sum.7(CTS10,CTS14,CTS18,CTS22,CTS24,CTS32,CTS40) .
execute .

compute Victim=sum(PsychAggV,PhysAssV) .
compute Perp=sum(PsychAggP,PhysAssP) . 
execute .

VARIABLE LABELS 
PsychAggV 'CTS psychological aggression-victim',
PsychAggP 'CTS psychological aggression-perpetrator',
PhysAssV 'CTS physical assault-victim',
PhysAssP 'CTS physical assault-perpetrator',
PsychAggVM 'CTS psychological aggression-victim (minor)',
PsychAggPM 'CTS psychological aggression-perpetrator (minor)',
PsychAggVS 'CTS psychological aggression-victim (severe)',
PsychAggPS 'CTS psychological aggression-perpetrator (severe)',
PhysAssVM 'CTS physical assault-victim (minor)',
PhysAssPM 'CTS physical assault-perpetrator (minor)',
PhysAssVS 'CTS physical assault-victim (severe)',
PhysAssPS 'CTS physical assault-perpetrator (severe)',
Victim 'CTS victim total',
Perp 'CTS perpetrator total'.
Execute.


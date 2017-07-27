GET DATA /TYPE=XLSX 
  /FILE='C:\Users\cushmank\Desktop\PANASPrePost.xlsx' 
  /SHEET=name 'PANASPre' 
  /CELLRANGE=full 
  /READNAMES=on 
  /ASSUMEDSTRWIDTH=32767. 
EXECUTE. 

ALTER TYPE PTNUM DyadID UserID ScrnID mth ratercode raterID PrePost (F2.0).
EXECUTE.
ALTER TYPE PANAS1 to PANAS20 (F2.0).
Execute.

VALUE LABELS prepost 3 'pre-interaction' 4 'post-interaction'.
Execute.
VALUE LABELS panas1 to panas20 1 'very slightly or not at all' 2 'a little' 3 'moderately' 4 'quite a bit' 5 'extremely'.
EXECUTE.
VARIABLE LABELS panas1 'interested' panas2 'distressed' panas3 'excited' panas4 'upset' panas5 'strong' panas6 'guilty' panas7 'scared'
panas8 'hostile' panas9 'enthusiastic' panas10 'proud' panas11 'irritable' panas12 'alert' panas13 'ashamed' panas14 'inspired' panas15 'nervous'
panas16 'determined' panas17 'attentive' panas18 'jittery' panas19 'active' panas20 'afraid'.
EXECUTE.

COMPUTE
NegAffect=sum.9(PANAS2,PANAS4,PANAS6,PANAS7,PANAS8,PANAS11,PANAS13,PANAS15,PANAS18,PANAS20) .
EXECUTE.
COMPUTE 
PosAffect=sum.9(PANAS1,PANAS3,PANAS5,PANAS9,PANAS10,PANAS12,PANAS14,PANAS16,PANAS17,PANAS19).
EXECUTE .
ALTER TYPE NegAffect PosAffect (F2.0).
EXECUTE.

RECODE PANAS1 to PANAS20 (SYSMIS=999).
EXECUTE.
RECODE NegAffect PosAffect (SYSMIS=999).
EXECUTE.
MISSING VALUES 
PANAS1 to PANAS20 (999).
EXECUTE.
MISSING VALUES NegAffect PosAffect (999).
EXECUTE.

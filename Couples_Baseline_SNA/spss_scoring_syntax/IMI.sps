
GET DATA
  /TYPE=ODBC
  /CONNECT='DSN=MS Access Database;DBQ=X:\Couples Study\Couples Study Single '+
    'Entry.accdb;DriverId=25;FIL=MS Access;MaxBufferSize=2048;PageTimeout=5;'
  /SQL='SELECT PTNUM, DyadID, ScrnID, mth, ratercode, raterID, pre, IMI1, IMI2, IMI3, IMI4, '+
    'IMI5, IMI6, IMI7, IMI8, IMI9, IMI10, IMI11, IMI12, IMI13, IMI14, IMI15, IMI16, IMI17, IMI18, '+    
    'IMI19, IMI20, IMI21, IMI22, IMI23, IMI24, IMI25, IMI26, IMI27, IMI28, IMI29, IMI30, IMI31, '+
    'IMI32 FROM IMIC_PRE'
  /ASSUMEDSTRWIDTH=255.

CACHE.
EXECUTE.

GET DATA
  /TYPE=ODBC
  /CONNECT='DSN=MS Access Database;DBQ=X:\Couples Study\Couples Study Single '+
    'Entry.accdb;DriverId=25;FIL=MS Access;MaxBufferSize=2048;PageTimeout=5;'
  /SQL='SELECT PTNUM, DyadID, ScrnID, mth, ratercode, raterID, post, IMI1, IMI2, IMI3, IMI4, '+
    'IMI5, IMI6, IMI7, IMI8, IMI9, IMI10, IMI11, IMI12, IMI13, IMI14, IMI15, IMI16, IMI17, IMI18, '+    
    'IMI19, IMI20, IMI21, IMI22, IMI23, IMI24, IMI25, IMI26, IMI27, IMI28, IMI29, IMI30, IMI31, '+
    'IMI32 FROM IMIC_Post'
  /ASSUMEDSTRWIDTH=255.

CACHE.
EXECUTE.
DATASET NAME DataSet2 WINDOW=FRONT.

VALUE LABELS prepost 3 'pre-interaction' 4 'post-interaction'.
Execute.
VALUE LABELS imi1 to imi32 1 'not at all' 2 'somewhat' 3 'moderately so' 4 'very much so'.
EXECUTE.
VARIABLE LABELS imi1 'hearing my partner makes me feel bossed around' imi2 'hearing my partner makes me feel distant from him/her' imi3 'hearing my partner makes me feel important'
imi4 'hearing my partner makes me feel in charge' imi5 'hearing my partner makes me feel appreciated by him/her' imi6 'hearing my partner makes me feel complimented'
imi7 'hearing my partner makes me feel dominant' imi8 'hearing my partner makes me feel welcome by him/her' imi9 'hearing my partner makes me feel annoyed'
imi10 'hearing my partner makes me feel taken charge of' imi11 'my partner makes me feel like I could lean on him/her for support' imi12 'my partner makes me feel like I am an intruder'
imi13 'my partner makes me feel like I should tell him/her he/she is often inconsiderate' imi14 'my partner makes me feel like I can relax and s/he will take charge' 
imi15 'my partner makes me feel like I could ask him/her to do anything' imi16 'based on interactions w/ partner it appears that s/he doesnt want to get involved w/ me'
imi17 'based on interactions w/ partner it appears that s/he is most comfotable withdrawing' imi18 'based on interactions w/ partner it appears that s/he wants me to put him/her on pedestal'
imi19 'based on interactions w/ partner it appears that s/he thinks he/she cant do anything for self' imi20 'based on interactions w/ partner it appears that s/he thinks its every person for self' 
imi21 'based on interactions w/ partner it appears that s/he would accept whatever I said' imi22 'based on interactions w/ partner it appears that s/he wants to be charming one'
imi23 'based on interactions w/ partner it appears that s/he thinks s/he is always in control' imi24 'based on interactions w/ partner it appears that s/he thinks s/hes inadequate' 
imi25 'based on interactions w/ partner it appears that s/he thinks I have most of the answers' imi26 'based on interactions w/ partner it appears that s/he enjoys being with people' 
imi27 'based on interactions w/ partner it appears that s/he would rather be left alone' imi28 'based on interactions w/ partner it appears that s/he sees me as superior'
imi29 'based on interactions w/ partner it appears that s/he is carrying a grudge' imi30 'based on interactions w/ partner it appears that s/he is nervous around me'
imi31 'based on interactions w/ partner it appears that s/he trusts me' imi32 'based on interactions w/ partner it appears that s/he thinks others find him/her interesting'.
EXECUTE.

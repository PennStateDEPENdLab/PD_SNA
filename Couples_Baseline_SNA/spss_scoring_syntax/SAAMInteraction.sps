GET DATA /TYPE=XLSX
  /FILE='C:\Users\cushmank\Desktop\IntrctnData\SAAM\SAAMPrePost.xlsx'
  /SHEET=name 'SAAMPre'
  /CELLRANGE=full
  /READNAMES=on
  /ASSUMEDSTRWIDTH=32767.
EXECUTE.

ALTER TYPE PTNUM DyadID UserID ScrnID mth ratercode raterID PrePost (F2.0).
EXECUTE.
ALTER TYPE SAAM1 to SAAM21 (F2.0).
EXECUTE.

VALUE LABELS prepost 3 'pre-interaction' 4 'post-interaction'.
EXECUTE.
VALUE Labels SAAM1 to SAAM21 1 'disagree strongly' 2 '...' 3 '...' 4 'neutral/mixed' 5 '...' 6 '...' 7 'agree strongly'.
EXECUTE.

VARIABLE LABELS saam1 'I wish someone would tell me they really love me' saam2 'I would be uncomfortable having a good friend or partner close to me' saam3 'I feel alone and dont feel like getting close to others'
saam4 'I feel loved' saam5 'I wish someone close could see me now' saam6 'If something went wrong I feel like I could depend on someone' saam7 'I feel like others care about me' saam8 'I feel a strong need to be unconditionally loved'
saam9 'Im afraid someone will want to get too close to me' saam10 'If someone tried to get close to me I would try to keep my distance' saam11 'I feel relaxed knowing that close others are there for me right now'
saam12 'I really need to feel loved right now' saam13 'I feel like I have someone to rely on' saam14 'I want to share my feelings with someone' saam15 'I feel like I am loved by others but really dont care'
saam16 'The idea of being emotionally close to someone makes me nervous' saam17 'I want to talk with someone who cares for me about things that are worrying me'
saam18 'I feel secure and close to other people' saam19 'I really need someones emotional support' saam20 'I feel I can trust those close to me' saam21 'I have mixed feeligns about being close to others'.
EXECUTE.

COMPUTE
Anxiety=sum.7(saam1,saam5,saam8,saam12,saam14,saam17,saam19).
EXECUTE.
COMPUTE
Avoidance=sum.7(saam2,saam3,saam9,saam10,saam15,saam16,saam21).
EXECUTE.
COMPUTE
Security=sum.7(saam4,saam6,saam7,saam11,saam13,saam18,saam20).
EXECUTE.
ALTER TYPE anxiety, avoidance, security (F2.0).
EXECUTE.


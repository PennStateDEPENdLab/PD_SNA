* bipolar mood disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (1 thru 2=1) (51 thru 52=1)  (35=1) (85=1)  (ELSE=0)  INTO  bip1  bip2 bip3 bip4 bip5 bip6 bip7 bip8 bip9 bip10 bip11 bip12 .
EXECUTE .
IF (bip1 = 1 OR bip2 = 1 OR bip3 = 1 OR bip4 = 1 OR bip5 = 1 OR bip6 = 1 OR bip7 = 1 OR bip8 = 1 OR bip9 = 1 OR bip10 = 1 OR bip11 = 1 or bip12 = 1) life_bip_dx = 1 .
EXECUTE .
RECODE
  life_bip_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (1 thru 2=1)  (35=1)  (ELSE=0)  INTO  cbip1  cbip2 cbip3 cbip4 cbip5 cbip6 cbip7 cbip8 cbip9 cbip10 cbip11 cbip12.
EXECUTE .
IF (cbip1 = 1 OR cbip2 = 1 OR cbip3 = 1 OR cbip4 = 1 OR cbip5 = 1 OR cbip6 = 1 OR cbip7 = 1 OR cbip8 = 1 OR cbip9 = 1 OR cbip10 = 1 OR cbip11 = 1 OR cbip12 = 1) curr_bip_dx = 1 .
EXECUTE .
RECODE
  curr_bip_dx  (MISSING=0)  .
EXECUTE .

* unipolar mood disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (3 thru 5=1) (53 thru 55=1) (8=1) (58=1) (ELSE=0) INTO mdd1 mdd2 mdd3 mdd4 mdd5 mdd6 mdd7 mdd8 mdd9 mdd10 mdd11 mdd12.
EXECUTE .
IF (mdd1 = 1 OR mdd2 = 1 OR mdd3 = 1 OR mdd4 = 1 OR mdd5 = 1 OR mdd6 = 1 OR mdd7 = 1 OR mdd8 = 1 OR mdd9 = 1 OR mdd10 = 1 OR mdd11 = 1 OR mdd12 = 1) life_dep_dx = 1 .
RECODE
  life_dep_dx  (MISSING=0) .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (3 thru 5=1) (8=1) (ELSE=0) INTO cmdd1 cmdd2 cmdd3 cmdd4 cmdd5 cmdd6 cmdd7 cmdd8 cmdd9 cmdd10 cmdd11 cmdd12.
EXECUTE .
IF (cmdd1 = 1 OR cmdd2 = 1 OR cmdd3 = 1 OR cmdd4 = 1 OR cmdd5 = 1 OR cmdd6 = 1 OR cmdd7 = 1 OR cmdd8 = 1 OR cmdd9 = 1 OR cmdd10 = 1 OR cmdd11 = 1 OR cmdd12 = 1) curr_dep_dx = 1 .
RECODE
  curr_dep_dx  (MISSING=0) .
EXECUTE .

* any mood disorder .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (1 thru 5=1) (51 thru 55=1)  (8=1)  (58=1) (35=1) (85=1)  (ELSE=0)  INTO  mood1  mood2 mood3 mood4 mood5 mood6 mood7 mood8 mood9 mood10 mood11 mood12.
EXECUTE .
IF (mood1 = 1 OR mood2 = 1 OR mood3 = 1 OR mood4 = 1 OR mood5 = 1 OR mood6 = 1 OR mood7 = 1 OR mood8 = 1 OR mood9 = 1 OR mood10 = 1 OR mood11 = 1 OR mood12 = 1) life_mood_dx = 1 .
EXECUTE .
RECODE
  life_mood_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (1 thru 5=1)  (8=1)  (35=1)  (ELSE=0)  INTO  cmood1  cmood2 cmood3 cmood4 cmood5 cmood6 cmood7 cmood8 cmood9 cmood10 cmood11 cmood12.
EXECUTE .
IF (cmood1 = 1 OR cmood2 = 1 OR cmood3 = 1 OR cmood4 = 1 OR cmood5 = 1 OR cmood6 = 1 OR cmood7 = 1 OR cmood8 OR cmood9 = 1 OR cmood10 = 1 OR cmood11 = 1 OR cmood12 = 1) curr_mood_dx = 1 .
EXECUTE .
RECODE curr_mood_dx  (MISSING=0)  .
EXECUTE .

* drug abuse disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (12 thru 20=1) (62 thru 80 = 1) (ELSE=0) INTO drug1 drug2 drug3 drug4 drug5 drug6 drug7 drug8 drug9 drug10 drug11 drug12.
EXECUTE .
IF (drug1 = 1 OR drug2 = 1 OR drug3 = 1 OR drug4 = 1 OR drug5 = 1 OR drug6 = 1 OR drug7 = 1 OR drug8 = 1 OR drug9 = 1 OR drug10 = 1 OR drug11 = 1 OR drug12 =1) life_drug_dx = 1 .
EXECUTE .
RECODE
  Life_drug_dx (MISSING=0) .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (12 thru 20=1) (ELSE=0) INTO cdrug1 cdrug2 cdrug3 cdrug4 cdrug5 cdrug6 cdrug7 cdrug8 cdrug9 cdrug10 cdrug11 cdrug12.
EXECUTE .
IF (cdrug1 = 1 OR cdrug2 = 1 OR cdrug3 = 1 OR cdrug4 = 1 OR cdrug5 = 1 OR cdrug6 = 1 OR cdrug7 = 1 OR cdrug8 = 1 OR cdrug9 = 1 OR cdrug10 = 1 OR cdrug11 = 1 OR cdrug12 = 1) curr_drug_dx = 1 .
EXECUTE .
RECODE
  curr_drug_dx (MISSING=0) .
EXECUTE .

* anxiety disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (40=1) (90=1)  (21 thru 26=1) (71 thru 76=1) (36 thru 37=1) (86 thru 87=1)  (ELSE=0)  INTO  anx1  anx2 anx3 anx4 anx5 anx6 anx7 anx8 anx9 anx10 anx11 anx12.
EXECUTE .
IF (anx1 = 1 OR anx2 = 1 OR anx3 = 1 OR anx4 = 1 OR anx5 = 1 OR anx6 = 1 OR anx7 = 1 OR anx8 = 1 OR anx9 = 1 OR anx10 = 1 OR anx11 = 1 OR anx12 = 1) life_anx_dx = 1 .
EXECUTE .
RECODE
  life_anx_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (40=1) (21 thru 26=1) (36 thru 37=1)  (ELSE=0)  INTO  canx1 canx2 canx3 canx4 canx5 canx6 canx7 canx8 canx9 canx10 canx11 canx12.
EXECUTE .
IF (canx1 = 1 OR canx2 = 1 OR canx3 = 1 OR canx4 = 1 OR canx5 = 1 OR canx6 = 1 OR canx7 = 1 OR canx8 = 1 OR canx9 = 1 OR canx10 = 1 OR canx11 = 1 OR canx12 = 1 ) curr_anx_dx = 1 .
EXECUTE .
RECODE
  curr_anx_dx  (MISSING=0)  .
EXECUTE .

* somatoform disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (27 thru 30=1)  (77 thru 80=1)  (38=1)  (88=1) (ELSE=0)  INTO  som1  som2 som3 som4 som5 som6 som7 som8 som9 som10 som11 som12.
EXECUTE .
IF (som1 = 1 OR som2 = 1 OR som3 = 1 OR som4 = 1 OR som5 = 1 OR som6 = 1 OR som7 = 1 OR som8 = 1 OR som9 = 1 OR som10 = 1 OR som11 = 1 OR som12 = 1) Life_som_dx = 1 .
EXECUTE .
RECODE
  life_som_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (27 thru 30=1)  (38=1) (ELSE=0)  INTO  csom1  csom2 csom3 csom4 csom5 csom6 csom7 csom8 csom9 csom10 csom11 csom12.
EXECUTE .
IF (csom1 = 1 OR csom2 = 1 OR csom3 = 1 OR csom4 = 1 OR csom5 = 1 OR csom6 = 1 OR csom7 = 1 OR csom8 = 1 OR csom9 = 1 OR csom10 = 1 OR csom11 = 1 OR csom12 = 1) curr_som_dx = 1 .
EXECUTE .
RECODE
  curr_som_dx  (MISSING=0)  .
EXECUTE .

* eating disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (31 thru 32=1)  (81 thru 82=1)  (39=1)  (89=1) (ELSE=0)  INTO  eat1  eat2 eat3 eat4 eat5 eat6 eat7 eat8 eat9 eat10 eat11 eat12.
EXECUTE .
IF (eat1 = 1 OR eat2 = 1 OR eat3 = 1 OR eat4 = 1 OR eat5 = 1 OR eat6 = 1 OR eat7 = 1 OR eat8 = 1 OR eat9 = 1 OR eat10 = 1 OR eat11 = 1 OR eat12 = 1) life_eat_dx = 1 .
EXECUTE .
RECODE
  life_eat_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (31 thru 32=1)  (39=1) (ELSE=0)  INTO  ceat1  ceat2 ceat3 ceat4 ceat5 ceat6 ceat7 ceat8 ceat9 ceat10 ceat11 ceat12.
EXECUTE .
IF (ceat1 = 1 OR ceat2 = 1 OR ceat3 = 1 OR ceat4 = 1 OR ceat5 = 1 OR ceat6 = 1 OR ceat7 = 1 OR ceat8 = 1 OR ceat9 = 1 OR ceat10 = 1 OR ceat11 = 1 OR ceat12 = 1) curr_eat_dx = 1 .
EXECUTE .
RECODE
  curr_eat_dx  (MISSING=0)  .
EXECUTE .

* other disorders .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (33 thru 34=1)  (83 thru 84=1) (ELSE=0)  INTO  oth1  oth2 oth3 oth4 oth5 oth6 oth7 oth8 oth9 oth10 oth11 oth12 .
EXECUTE .
IF (oth1 = 1 OR oth2 = 1 OR oth3 = 1 OR oth4 = 1 OR oth5 = 1 OR oth6 = 1 OR oth7 = 1 OR oth8 = 1 OR oth9 = 1 OR oth10 = 1 OR oth11 = 1 OR oth12 = 1) Life_oth_dx = 1 .
EXECUTE .
RECODE
  life_oth_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (33 thru 34=1)  (ELSE=0)  INTO  coth1  coth2 coth3 coth4 coth5 coth6 coth7 coth8  coth9 coth10 coth11 coth12.
EXECUTE .
IF (coth1 = 1 OR coth2 = 1 OR coth3 = 1 OR coth4 = 1 OR coth5 = 1 OR coth6 = 1 OR coth7 = 1 OR coth8 = 1 OR coth9 = 1 OR coth10 = 1 OR coth11 = 1 OR coth12 = 1) curr_oth_dx = 1 .
EXECUTE .
RECODE
  curr_oth_dx  (MISSING=0)  .
EXECUTE .

* dysthymia .
RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (4=1) (54=1)  (ELSE=0) INTO  dys1  dys2 dys3 dys4 dys5 dys6 dys7 dys8 dys9 dys10 dys11 dys12.
EXECUTE .
IF (dys1 = 1 OR dys2 = 1 OR dys3 = 1 OR dys4 = 1 OR dys5 = 1 OR dys6 = 1 OR dys7 = 1 OR dys8 = 1 OR dys9 = 1 OR dys10 = 1 OR dys11 = 1 OR dys12 = 1) life_dys_dx = 1 .
EXECUTE .
RECODE
  life_dys_dx  (MISSING=0)  .
EXECUTE .

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (4=1)  (ELSE=0)  INTO  cdys1  cdys2 cdys3 cdys4 cdys5 cdys6 cdys7 cdys8 cdys9 cdys10 cdys11 cdys12.
EXECUTE .
IF (cdys1 = 1 OR cdys2 = 1 OR cdys3 = 1 OR cdys4 = 1 OR cdys5 = 1 OR cdys6 = 1 OR cdys7 = 1 OR cdys8 = 1 OR cdys9 = 1 OR cdys10 = 1 OR cdys11 = 1 OR cdys12 = 1) curr_dys_dx = 1 .
EXECUTE .
RECODE curr_dys_dx  (MISSING=0)  .
EXECUTE  . 

*pyschotic 

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (41=1)  (ELSE=0)  INTO  cpychot1  cpychot2 cpychot3 cpychot4 cpychot5 cpychot6 cpychot7 cpychot8 cpychot9 cpychot10 cpychot11 cpychot12.
EXECUTE .
IF (cpychot1 = 1 OR cpychot2 = 1 OR cpychot3 = 1 OR cpychot4 = 1 OR cpychot5 = 1 OR cpychot6 = 1 OR cpychot7 = 1 OR cpychot8 = 1 OR cpychot9 = 1 OR cpychot10 = 1 OR cpychot11 = 1 OR cpychot12 = 1) curr_pychot_dx = 1 .
EXECUTE .
RECODE curr_pychot_dx  (MISSING=0)  .
EXECUTE  . 

RECODE axisI_dx1 axisI_dx2 axisI_dx3 axisI_dx4 axisI_dx5 axisI_dx6 axisI_dx7 axisI_dx8 axisI_dx9 axisI_dx10 axis_dx11 axis_dx12
  (91=1)  (ELSE=0)  INTO  pychot1  pychot2 pychot3 pychot4 pychot5 pychot6 pychot7 pychot8 pychot9 pychot10 pychot11 pychot12.
EXECUTE .
IF (pychot1 = 1 OR pychot2 = 1 OR pychot3 = 1 OR pychot4 = 1 OR pychot5 = 1 OR pychot6 = 1 OR pychot7 = 1 OR pychot8 = 1 OR pychot9 = 1 OR pychot10 = 1 OR pychot11 = 1 OR pychot12 = 1) life_pychot_dx = 1 .
EXECUTE .
RECODE life_pychot_dx  (MISSING=0)  .
EXECUTE  . 

*axis I categories; initiate creation of temp variables for creation of summary variables.
COMPUTE AffCtmp=0 .
EXECUTE .
IF (curr_dep_dx=1 OR curr_mood_dx=1 OR curr_dys_dx=1) AffCtmp=1 .
EXECUTE .
COMPUTE AffLtmp=0 .
EXECUTE .
IF (life_dep_dx=1 OR life_mood_dx=1 OR life_dys_dx=1) AffLtmp=1 .
EXECUTE .
COMPUTE AnxCtmp=0 .
EXECUTE .
IF (curr_anx_dx=1) AnxCtmp=1 .
EXECUTE .
COMPUTE AnxLtmp=0 .
EXECUTE .
IF (life_anx_dx=1) AnxLtmp=1 .
EXECUTE .
COMPUTE SubCtmp=0 .
EXECUTE .
IF (curr_drug_dx=1) SubCtmp=1 .
EXECUTE .
COMPUTE SubLtmp=0 .
EXECUTE .
IF (life_drug_dx=1) SubLtmp=1 .
EXECUTE .

*create hierarchy of categories

COMPUTE AffAnxSubC=0 .
IF (AffCtmp=1 AND AnxCtmp=1 AND SubCtmp=1) AffAnxSubC=1 .
EXECUTE .
COMPUTE AffAnxSubL=0 .
IF (AffLtmp=1 AND AnxLtmp=1 AND SubLtmp=1) AffAnxSubL=1 .
EXECUTE .
COMPUTE AnxSubC=0 .
IF (AnxCtmp=1 AND SubCtmp=1 AND AffCtmp=0) AnxSubC=1 .
EXECUTE .
COMPUTE AnxSubL=0 .
IF (AnxLtmp=1 AND SubLtmp=1 AND AffLtmp=0) AnxSubL=1 .
EXECUTE .
COMPUTE AffSubC=0 .
EXECUTE .
IF (AffCtmp=1 AND SubCtmp=1 AND AnxCtmp=0) AffSubC=1 .
EXECUTE .
COMPUTE AffSubL=0 .
EXECUTE .
IF (AffLtmp=1 AND SubLtmp=1 AND AnxLtmp=0) AffSubL=1 .
EXECUTE .
COMPUTE AffAnxC=0 .
EXECUTE .
IF (AffCtmp=1 AND AnxCtmp=1 AND SubCtmp=0) AffAnxC=1.
EXECUTE .
COMPUTE AffAnxL=0 .
EXECUTE .
IF (AffLtmp=1 AND AnxLtmp=1 AND SubLtmp=0) AffAnxL=1.
EXECUTE .
COMPUTE AnxC=0 .
EXECUTE .
IF (AnxCtmp=1 AND AffAnxC=0 AND AffAnxSubC=0) AnxC=1 .
COMPUTE AnxL=0 .
EXECUTE .
IF (AnxLtmp=1 AND AffAnxL=0 AND AffAnxSubL=0) AnxL=1 .
EXECUTE .
COMPUTE AffC=0 .
EXECUTE .
IF (AffCtmp=1 AND AffAnxC=0 AND AffSubC=0 AND AffAnxSubC=0) AffC=1 .
EXECUTE .
COMPUTE AffL=0 .
EXECUTE .
IF (AffLtmp=1 AND AffAnxL=0 AND AffSubL=0 AND AffAnxSubL=0) AffL=1 .
EXECUTE .
COMPUTE SubC=0 .
EXECUTE .
IF (SubCtmp=1 AND AnxSubC=0 AND AffSubC=0 AND AffAnxSubC=0) SubC=1 .
EXECUTE .
COMPUTE SubL=0 .
EXECUTE .
IF (SubLtmp=1 AND AnxSubL=0 AND AffSubL=0 AND AffAnxSubL=0) SubL=1 .
EXECUTE .
COMPUTE OtherL=0.
EXECUTE .
IF (life_som_dx=1 OR life_eat_dx=1 OR life_oth_dx=1) OtherL=1 .
EXECUTE .
COMPUTE OtherC=0 .
EXECUTE .
IF (curr_som_dx=1 OR curr_eat_dx=1 OR curr_oth_dx=1) OtherC=1 .
EXECUTE .

COMPUTE PychotL=0 .
EXECUTE .
IF (life_pychot_dx=1) PychotL=1 .
EXECUTE .
COMPUTE PychotC=0 .
EXECUTE .
IF (curr_pychot_dx=1) PychotC=1 .
EXECUTE .
COMPUTE NoneC=0 .
EXECUTE .
IF (AnxC=0 AND AffC=0 AND SubC=0 AND AffSubC=0  AND AnxSubC=0 AND AffAnxC=0 AND AffAnxSubC=0 AND OtherC=0 AND PychotC=0)  NoneC=1.
EXECUTE .
COMPUTE NoneL=0 .
IF (AffAnxSubL=0 AND AnxSubL=0 AND AffSubL=0 AND AffAnxL=0 AND AnxL=0 AND AffL=0 AND SubL=0 AND OtherL=0 AND PychotL=0) NoneL=1.
EXECUTE .





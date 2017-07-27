*code written by kimberly nolf on 3.17.14

RECODE
neo1 to neo240 
(1=0) (2=1) (3=2) (4=3) (5=4) .
EXECUTE .

RECODE
 neo1 neo61 neo121 neo181 neo36 neo96 neo156 neo11 neo71 neo46 neo106 neo166 neo21 neo81 neo141 neo231 neo56 neo116
neo176 neo206 neo236 neo32 neo92 neo7 neo67 neo127 neo187 neo42 neo102 neo162 neo222 neo17 neo77 neo137 neo52 neo112
 neo27 neo87 neo147 neo207 neo33 neo93 neo153 neo183 neo213 neo8 neo68 neo128 neo43 neo103 neo163 neo18 neo78 neo138
 neo198 neo228 neo53 neo113 neo173 neo28 neo88 neo148 neo208 neo238 neo4 neo64 neo124 neo39 neo99 neo159 neo189 
neo219 neo14 neo74 neo134 neo49 neo109 neo169 neo199 neo229 neo24 neo84 neo144 neo234 neo59 neo119 neo35 neo95 
neo155 neo10 neo70 neo130 neo190 neo220 neo45 neo105 neo20 neo80 neo140 neo55 neo115 neo175 neo205 neo30 neo90 
neo150 (0=4)  (1=3)  (2=2)  (3=1)  (4=0) INTO neo1r neo61r neo121r neo181r neo36r neo96r neo156r neo11r neo71r neo46r neo106r 
neo166r neo21r neo81r neo141r neo231r neo56r neo116r neo176r neo206r neo236r neo32r neo92r neo7r neo67r neo127r 
neo187r neo42r neo102r neo162r neo222r neo17r neo77r neo137r neo52r neo112r neo27r neo87r neo147r neo207r neo33r neo93r 
neo153r neo183r neo213r neo8r neo68r neo128r neo43r neo103r neo163r neo18r neo78r neo138r neo198r neo228r neo53r neo113r 
neo173r neo28r neo88r neo148r neo208r neo238r neo4r neo64r neo124r neo39r neo99r neo159r neo189r neo219r neo14r neo74r 
neo134r neo49r neo109r neo169r neo199r neo229r neo24r neo84r neo144r neo234r neo59r neo119r neo35r neo95r neo155r neo10r 
neo70r neo130r neo190r neo220r neo45r neo105r neo20r neo80r neo140r neo55r neo115r neo175r neo205r neo30r neo90r neo150r . 
EXECUTE .

Compute ntot = sum.43 (neo1r, neo31, neo61r, neo91, neo121r, neo151, neo181r, neo211, neo6, neo36r, neo66, neo96r, neo126 ,neo156r, 
neo186, neo216, neo11r, neo41, neo71r, neo101, neo131, neo161, neo191, neo221 ,neo16, neo46r, neo76, neo106r, neo136, neo166r,
neo196, neo226, neo21r, neo51, neo81r, neo111, neo141r, neo171, neo201, neo231r, neo26 ,neo56r, neo86, neo116r, neo146 ,neo176r, neo206r, neo236r).
execute.

compute etot =sum.43 (neo2, neo32r, neo62, neo92r, neo122, neo152, neo182, neo212, neo7r, neo37, neo67r, neo97, neo127r, neo157, neo187r, 
neo217, neo12, neo42r, neo72, neo102r, neo132, neo162r, neo192, neo222r, neo17r, neo47, neo77r, neo107, neo137r, neo167, neo197, 
neo227, neo22, neo52r, neo82, neo112r, neo142, neo172, neo202, neo232, neo27r, neo57, neo87r, neo117, neo147r, neo177, neo207r, neo237).
execute.
compute otot = sum.43 (neo3, neo33r, neo63, neo93r, neo123, neo153r, neo183r, neo213r, neo8r, neo38, neo68r, neo98, neo128r, neo158, neo188, neo218
, neo13, neo43r, neo73, neo103r, neo133, neo163r, neo193, neo223, neo18r,  neo48, neo78r, neo108, neo138r, neo168, neo198r, neo228r, 
neo23, neo53r, neo83, neo113r, neo143, neo173r, neo203, neo233, neo28r, neo58, neo88r, neo118, neo148r, neo178, neo208r, neo238r).
execute.
compute atot = sum.43 (neo4r, neo34, neo64r, neo94, neo124r, neo154, neo184, neo214, neo9, neo39r, neo69, neo99r, neo129, neo159r, neo189r, neo219r, 
neo14r, neo44, neo74r, neo104, neo134r, neo164, neo194, neo224, neo19, neo49r, neo79, neo109r, neo139, neo169r, neo199r, neo229r, neo24r, 
neo54, neo84r, neo114, neo144r, neo174, neo204, neo234r, neo29, neo59r, neo89, neo119r, neo149, neo179, neo209, neo239).
execute.
compute ctot = sum.43 (neo5, neo35r, neo65, neo95r, neo125, neo155r, neo185, neo215, neo10r, neo40, neo70r, neo100, neo130r, neo160, 
neo190r, neo220r, neo15, neo45r, neo75, neo105r, neo135, neo165, neo195, neo225, neo20r, neo50, neo80r, neo110, neo140r, neo170, 
neo200, neo230, neo25, neo55r, neo85, neo115r, neo145, neo175r, neo205r, neo235, neo30r, neo60, neo90r, neo120, neo150r, neo180, neo210, 
neo240).
execute.

recode ntot,etot,otot,atot,ctot (sysmis=9999) .
execute .

variable labels ntot 'NEO Neuroticism' etot 'NEO Extraversion' otot 'NEO Openness to Experience' atot 'NEO Agreeableness' ctot 'NEO Conscientiousness' .

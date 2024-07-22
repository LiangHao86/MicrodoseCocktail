;; 1. Based on: RSV_1.1
;; 2. Description: RSV BASEMODEL CMT2

$PROB  BASE STRUCTURE MODEL CMT2
$ABBR DERIV2=NO
$DATA RSV_PK_Full_20240121_GUT_20240411.csv IGNORE=@
$INPUT ID,TIME,DV,AMT,MDV,EVID,CMT,AGE,ALB,CRCL,EGFR,TP,POP,CloXVIII,CloXlVb,RumiBilophila,Clostridium_XlVb,Parabacteroides

$SUBROUTINES  ADVAN13  TOL=5 

$MODEL 
NCOMPARTMENTS=3
COMP=(DEPOT,DEFDOSE)
COMP=(CENTRAL,DEFOBS)
COMP=(PER)

$PK
TVCL     = THETA(1) 
CL       = TVCL* EXP(ETA(1))

TVVC     = THETA(2)
VC       = TVVC * EXP(ETA(2)) 

TVKA     = THETA(3) 
KA       = TVKA* EXP(ETA(3))

TVVP     = THETA(4)
VP       = TVVP * EXP(ETA(4))

TVQ     = THETA(5)
Q       = TVQ * EXP(ETA(5))


;;;;;;;;;;;;;;;;;Initial Value ;;;;;;;;;;;;;;;;;;
;A_0FLG is set by NONMEM for the first record of each subject (i.e. NEWIND.LE.0)
IF (A_0FLG.EQ.1) THEN
A_0(1)   = 0
A_0(2)   = 0
ENDIF


$DES
C1      =  A(2)/VC
C2      =  A(3)/VP
DADT(1) = -KA*A(1)
DADT(2) =  KA*A(1)-C1*(CL+Q)+C2*Q
DADT(3) =  C1*Q-C2*Q


$THETA  
(0,46,500)    ;CL/F, L/h
(0,316,1000)  ;VC/F, L
(0,0.3,3)       ;KA, 1/h
(0,1890,4000)    ;VP/F, L
(0,42,300)     ;Q, L/h


$OMEGA  
0.4              ;ETA_CL
0.2              ;ETA_VC
0 FIX             ;ETA_KA
0 FIX              ;ETA_VP
0 FIX            ;ETA_Q


$ERROR
IF (CMT.EQ.2) IPRED = A(2)/VC 
IRES = DV-IPRED
W=IPRED
DEL = 0

IF(IPRED.EQ.0) DEL = 1
IWRES = IRES/(IPRED+DEL)
Y = IPRED*(1+EPS(1))+EPS(2)


$SIGMA  
0.1 ;pro  
0.00001 FIX ;add

$EST MAXEVAL=9999 METHOD=1 INTER PRINT=5 SIGDIGITS=5 NOABORT
$COV  PRINT=E MATRIX=S
$TABLE ID TIME CMT DV AMT MDV IPRED PRED IWRES RES IRES CWRES WRES ONEHEADER NOPRINT FILE=001.tab
$TABLE ID KA CL VC VP Q ETA1 ETA2 ETA3 ETA4 ETA5 FILE=PAT001.tab ONEHEADER NOPRINT FIRSTONLY

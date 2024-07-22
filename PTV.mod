;; 1. Based on: PTV
;; 2. Description: PTV base model

$PROB  Base model
$ABBR DERIV2=NO
$DATA PTV_PK_Full_20240121_GUT_20240411.csv IGNORE=N
$INPUT  ID,TIME,DV,AMT,MDV,EVID,CMT,CRCL,EGFR,POP,Veil


$SUBROUTINES  ADVAN6  TOL=6 

$MODEL NCOMPARTMENTS=3
COMP=(DEPOT,DEFDOSE)              ;1
COMP=(CENTRAL,DEFOBS)
COMP=PERIPH


$PK
TVCL     = THETA(1) 
CL       = TVCL* EXP(ETA(1))

TVVC     = THETA(2)
VC       = TVVC * EXP(ETA(2)) 

TVKA     = THETA(3) 
KA       = TVKA* EXP(ETA(3))

TVVT     = THETA(4)
VT       = TVVT* EXP(ETA(4))

TVCLD     = THETA(5) 
CLD      = TVCLD* EXP(ETA(5))

;;;;;;;;;;;;;;;;;Initial Value ;;;;;;;;;;;;;;;;;;
;A_0FLG is set by NONMEM for the first record of each subject (i.e. NEWIND.LE.0)
IF (A_0FLG.EQ.1) THEN
A_0(1)   = 0
A_0(2)   = 0
A_0(3)   = 0
ENDIF


$DES
DADT(1) = -KA*A(1)
C1           =  A(2)/VC
C2           =  A(3)/VT

DADT(2) =  KA*A(1)-C1*CL-(C1-C2)*CLD

DADT(3) =   (C1-C2)*CLD



$THETA  

(0,13.8, 40)            ;CL  
(0, 27.4, 100)         ;VC
18.1 FIX          ;KA
(0, 118, 300)        ;VT
(0, 14, 50)            ;CLd



$OMEGA  
0.09               ;ETA_CL
0.09                   ; ETA_VC
0 FIX                   ;ETA_KA
0.09              ; ETA_VT
0 FIX                 ; ETA_CLD




$ERROR  ; CCV ERROR MODEL


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;DU176 ERROR;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF (CMT.EQ.2) IPRED = A(2)/VC 


IRES = DV-IPRED

W    =IPRED

DEL = 0
IF(W.EQ.0) DEL = 1
IWRES = IRES/(W+DEL)
Y = IPRED+EPS(1)*W+EPS(2)


$SIGMA  
0.1058   
0.00001 FIX

$EST MAXEVAL=9999 METHOD=1 INTER PRINT=5 SIGDIGITS=3 NOABORT NOTBT NOOBT NOSBT
$COV MATRIX=R  PRINT=E
$TABLE ID,TIME, DV, CMT, MDV,EVID,IPRED, PRED, WRES CWRES IRES IWRES ONEHEADER NOPRINT FILE=sdtab001.tab
$TABLE ID CL VC KA VT CLD ETA(1) ETA(2) ETA(3) ETA(4) ETA(5) FILE=patab001.tab ONEHEADER NOPRINT FIRSTONLY

#ifdef BFM_NOPOINTERS


#define NOPOINTERS 1

#define O2o(A) D3STATE(A,ppO2o)
#define N1p(A) D3STATE(A,ppN1p)
#define N3n(A) D3STATE(A,ppN3n)
#define N4n(A) D3STATE(A,ppN4n)
#define O4n(A) D3STATE(A,ppO4n)
#define N5s(A) D3STATE(A,ppN5s)
#define N6r(A) D3STATE(A,ppN6r)
#define B1c(A) D3STATE(A,ppB1c)
#define B1n(A) D3STATE(A,ppB1n)
#define B1p(A) D3STATE(A,ppB1p)
#define P1c(A) D3STATE(A,ppP1c)
#define P1n(A) D3STATE(A,ppP1n)
#define P1p(A) D3STATE(A,ppP1p)
#define P1l(A) D3STATE(A,ppP1l)
#define P1s(A) D3STATE(A,ppP1s)
#define P2c(A) D3STATE(A,ppP2c)
#define P2n(A) D3STATE(A,ppP2n)
#define P2p(A) D3STATE(A,ppP2p)
#define P2l(A) D3STATE(A,ppP2l)
#define P3c(A) D3STATE(A,ppP3c)
#define P3n(A) D3STATE(A,ppP3n)
#define P3p(A) D3STATE(A,ppP3p)
#define P3l(A) D3STATE(A,ppP3l)
#define P4c(A) D3STATE(A,ppP4c)
#define P4n(A) D3STATE(A,ppP4n)
#define P4p(A) D3STATE(A,ppP4p)
#define P4l(A) D3STATE(A,ppP4l)
#define Z3c(A) D3STATE(A,ppZ3c)
#define Z3n(A) D3STATE(A,ppZ3n)
#define Z3p(A) D3STATE(A,ppZ3p)
#define Z4c(A) D3STATE(A,ppZ4c)
#define Z4n(A) D3STATE(A,ppZ4n)
#define Z4p(A) D3STATE(A,ppZ4p)
#define Z5c(A) D3STATE(A,ppZ5c)
#define Z5n(A) D3STATE(A,ppZ5n)
#define Z5p(A) D3STATE(A,ppZ5p)
#define Z6c(A) D3STATE(A,ppZ6c)
#define Z6n(A) D3STATE(A,ppZ6n)
#define Z6p(A) D3STATE(A,ppZ6p)
#define R1c(A) D3STATE(A,ppR1c)
#define R1n(A) D3STATE(A,ppR1n)
#define R1p(A) D3STATE(A,ppR1p)
#define R2c(A) D3STATE(A,ppR2c)
#define R3c(A) D3STATE(A,ppR3c)
#define R6c(A) D3STATE(A,ppR6c)
#define R6n(A) D3STATE(A,ppR6n)
#define R6p(A) D3STATE(A,ppR6p)
#define R6s(A) D3STATE(A,ppR6s)
#define O3c(A) D3STATE(A,ppO3c)
#define O3h(A) D3STATE(A,ppO3h)
#define O5c(A) D3STATE(A,ppO5c)

#define PelBacteria(A,B) D3STATE(:,ppPelBacteria(A,B))
#define PhytoPlankton(A,B) D3STATE(:,ppPhytoPlankton(A,B))
#define MesoZooPlankton(A,B) D3STATE(:,ppMesoZooPlankton(A,B))
#define MicroZooPlankton(A,B) D3STATE(:,ppMicroZooPlankton(A,B))
#define PelDetritus(A,B) D3STATE(:,ppPelDetritus(A,B))
#define Inorganic(A,B) D3STATE(:,ppInorganic(A,B))

#define Depth(A) D3DIAGNOS(A,ppDepth)
#define Volume(A) D3DIAGNOS(A,ppVolume)
#define Area(A) D3DIAGNOS(A,ppArea)
#define ETW(A) D3DIAGNOS(A,ppETW)
#define ESW(A) D3DIAGNOS(A,ppESW)
#define ERHO(A) D3DIAGNOS(A,ppERHO)
#define EIR(A) D3DIAGNOS(A,ppEIR)
#define ESS(A) D3DIAGNOS(A,ppESS)
#define EPR(A) D3DIAGNOS(A,ppEPR)
#define DIC(A) D3DIAGNOS(A,ppDIC)
#define CO2(A) D3DIAGNOS(A,ppCO2)
#define pCO2(A) D3DIAGNOS(A,pppCO2)
#define HCO3(A) D3DIAGNOS(A,ppHCO3)
#define CO3(A) D3DIAGNOS(A,ppCO3)
#define ALK(A) D3DIAGNOS(A,ppALK)
#define pH(A) D3DIAGNOS(A,pppH)
#define OCalc(A) D3DIAGNOS(A,ppOCalc)
#define OArag(A) D3DIAGNOS(A,ppOArag)
#define ffCO2(A) D3DIAGNOS(A,ppffCO2)
#define totpelc(A) D3DIAGNOS(A,pptotpelc)
#define totpeln(A) D3DIAGNOS(A,pptotpeln)
#define totpelp(A) D3DIAGNOS(A,pptotpelp)
#define totpels(A) D3DIAGNOS(A,pptotpels)
#define cxoO2(A) D3DIAGNOS(A,ppcxoO2)
#define eO2mO2(A) D3DIAGNOS(A,ppeO2mO2)
#define Chla(A) D3DIAGNOS(A,ppChla)
#define flPTN6r(A) D3DIAGNOS(A,ppflPTN6r)
#define flN3O4n(A) D3DIAGNOS(A,ppflN3O4n)
#define flN4N3n(A) D3DIAGNOS(A,ppflN4N3n)
#define sediR2(A) D3DIAGNOS(A,ppsediR2)
#define sediR3(A) D3DIAGNOS(A,ppsediR3)
#define sediR6(A) D3DIAGNOS(A,ppsediR6)
#define sediO5(A) D3DIAGNOS(A,ppsediO5)
#define xEPS(A) D3DIAGNOS(A,ppxEPS)
#define ABIO_eps(A) D3DIAGNOS(A,ppABIO_eps)
#define BGE(A) D3DIAGNOS(A,ppBGE)
#define BIOALK(A) D3DIAGNOS(A,ppBIOALK)

#define Area2d(A) D2DIAGNOS(A,ppArea2d)
#define ThereIsLight(A) D2DIAGNOS(A,ppThereIsLight)
#define SUNQ(A) D2DIAGNOS(A,ppSUNQ)
#define EWIND(A) D2DIAGNOS(A,ppEWIND)
#define EICE(A) D2DIAGNOS(A,ppEICE)
#define ETAUB(A) D2DIAGNOS(A,ppETAUB)
#define EPCO2air(A) D2DIAGNOS(A,ppEPCO2air)
#define CO2airflux(A) D2DIAGNOS(A,ppCO2airflux)
#define dpco2(A) D2DIAGNOS(A,ppdpco2)
#define dpo2(A) D2DIAGNOS(A,ppdpo2)
#define totsysc(A) D2DIAGNOS(A,pptotsysc)
#define totsysn(A) D2DIAGNOS(A,pptotsysn)
#define totsysp(A) D2DIAGNOS(A,pptotsysp)
#define totsyss(A) D2DIAGNOS(A,pptotsyss)

#define PELSURFACE(A,B) D2DIAGNOS(A,14+B)
#define jsurO2o(B) D2DIAGNOS(B,14+ppO2o)
#define jsurN1p(B) D2DIAGNOS(B,14+ppN1p)
#define jsurN3n(B) D2DIAGNOS(B,14+ppN3n)
#define jsurN4n(B) D2DIAGNOS(B,14+ppN4n)
#define jsurO4n(B) D2DIAGNOS(B,14+ppO4n)
#define jsurN5s(B) D2DIAGNOS(B,14+ppN5s)
#define jsurN6r(B) D2DIAGNOS(B,14+ppN6r)
#define jsurB1c(B) D2DIAGNOS(B,14+ppB1c)
#define jsurB1n(B) D2DIAGNOS(B,14+ppB1n)
#define jsurB1p(B) D2DIAGNOS(B,14+ppB1p)
#define jsurP1c(B) D2DIAGNOS(B,14+ppP1c)
#define jsurP1n(B) D2DIAGNOS(B,14+ppP1n)
#define jsurP1p(B) D2DIAGNOS(B,14+ppP1p)
#define jsurP1l(B) D2DIAGNOS(B,14+ppP1l)
#define jsurP1s(B) D2DIAGNOS(B,14+ppP1s)
#define jsurP2c(B) D2DIAGNOS(B,14+ppP2c)
#define jsurP2n(B) D2DIAGNOS(B,14+ppP2n)
#define jsurP2p(B) D2DIAGNOS(B,14+ppP2p)
#define jsurP2l(B) D2DIAGNOS(B,14+ppP2l)
#define jsurP3c(B) D2DIAGNOS(B,14+ppP3c)
#define jsurP3n(B) D2DIAGNOS(B,14+ppP3n)
#define jsurP3p(B) D2DIAGNOS(B,14+ppP3p)
#define jsurP3l(B) D2DIAGNOS(B,14+ppP3l)
#define jsurP4c(B) D2DIAGNOS(B,14+ppP4c)
#define jsurP4n(B) D2DIAGNOS(B,14+ppP4n)
#define jsurP4p(B) D2DIAGNOS(B,14+ppP4p)
#define jsurP4l(B) D2DIAGNOS(B,14+ppP4l)
#define jsurZ3c(B) D2DIAGNOS(B,14+ppZ3c)
#define jsurZ3n(B) D2DIAGNOS(B,14+ppZ3n)
#define jsurZ3p(B) D2DIAGNOS(B,14+ppZ3p)
#define jsurZ4c(B) D2DIAGNOS(B,14+ppZ4c)
#define jsurZ4n(B) D2DIAGNOS(B,14+ppZ4n)
#define jsurZ4p(B) D2DIAGNOS(B,14+ppZ4p)
#define jsurZ5c(B) D2DIAGNOS(B,14+ppZ5c)
#define jsurZ5n(B) D2DIAGNOS(B,14+ppZ5n)
#define jsurZ5p(B) D2DIAGNOS(B,14+ppZ5p)
#define jsurZ6c(B) D2DIAGNOS(B,14+ppZ6c)
#define jsurZ6n(B) D2DIAGNOS(B,14+ppZ6n)
#define jsurZ6p(B) D2DIAGNOS(B,14+ppZ6p)
#define jsurR1c(B) D2DIAGNOS(B,14+ppR1c)
#define jsurR1n(B) D2DIAGNOS(B,14+ppR1n)
#define jsurR1p(B) D2DIAGNOS(B,14+ppR1p)
#define jsurR2c(B) D2DIAGNOS(B,14+ppR2c)
#define jsurR3c(B) D2DIAGNOS(B,14+ppR3c)
#define jsurR6c(B) D2DIAGNOS(B,14+ppR6c)
#define jsurR6n(B) D2DIAGNOS(B,14+ppR6n)
#define jsurR6p(B) D2DIAGNOS(B,14+ppR6p)
#define jsurR6s(B) D2DIAGNOS(B,14+ppR6s)
#define jsurO3c(B) D2DIAGNOS(B,14+ppO3c)
#define jsurO3h(B) D2DIAGNOS(B,14+ppO3h)
#define jsurO5c(B) D2DIAGNOS(B,14+ppO5c)

#define PELBOTTOM(A,B) D2DIAGNOS(A,65+B)
#define jbotO2o(B) D2DIAGNOS(B,65+ppO2o)
#define jbotN1p(B) D2DIAGNOS(B,65+ppN1p)
#define jbotN3n(B) D2DIAGNOS(B,65+ppN3n)
#define jbotN4n(B) D2DIAGNOS(B,65+ppN4n)
#define jbotO4n(B) D2DIAGNOS(B,65+ppO4n)
#define jbotN5s(B) D2DIAGNOS(B,65+ppN5s)
#define jbotN6r(B) D2DIAGNOS(B,65+ppN6r)
#define jbotB1c(B) D2DIAGNOS(B,65+ppB1c)
#define jbotB1n(B) D2DIAGNOS(B,65+ppB1n)
#define jbotB1p(B) D2DIAGNOS(B,65+ppB1p)
#define jbotP1c(B) D2DIAGNOS(B,65+ppP1c)
#define jbotP1n(B) D2DIAGNOS(B,65+ppP1n)
#define jbotP1p(B) D2DIAGNOS(B,65+ppP1p)
#define jbotP1l(B) D2DIAGNOS(B,65+ppP1l)
#define jbotP1s(B) D2DIAGNOS(B,65+ppP1s)
#define jbotP2c(B) D2DIAGNOS(B,65+ppP2c)
#define jbotP2n(B) D2DIAGNOS(B,65+ppP2n)
#define jbotP2p(B) D2DIAGNOS(B,65+ppP2p)
#define jbotP2l(B) D2DIAGNOS(B,65+ppP2l)
#define jbotP3c(B) D2DIAGNOS(B,65+ppP3c)
#define jbotP3n(B) D2DIAGNOS(B,65+ppP3n)
#define jbotP3p(B) D2DIAGNOS(B,65+ppP3p)
#define jbotP3l(B) D2DIAGNOS(B,65+ppP3l)
#define jbotP4c(B) D2DIAGNOS(B,65+ppP4c)
#define jbotP4n(B) D2DIAGNOS(B,65+ppP4n)
#define jbotP4p(B) D2DIAGNOS(B,65+ppP4p)
#define jbotP4l(B) D2DIAGNOS(B,65+ppP4l)
#define jbotZ3c(B) D2DIAGNOS(B,65+ppZ3c)
#define jbotZ3n(B) D2DIAGNOS(B,65+ppZ3n)
#define jbotZ3p(B) D2DIAGNOS(B,65+ppZ3p)
#define jbotZ4c(B) D2DIAGNOS(B,65+ppZ4c)
#define jbotZ4n(B) D2DIAGNOS(B,65+ppZ4n)
#define jbotZ4p(B) D2DIAGNOS(B,65+ppZ4p)
#define jbotZ5c(B) D2DIAGNOS(B,65+ppZ5c)
#define jbotZ5n(B) D2DIAGNOS(B,65+ppZ5n)
#define jbotZ5p(B) D2DIAGNOS(B,65+ppZ5p)
#define jbotZ6c(B) D2DIAGNOS(B,65+ppZ6c)
#define jbotZ6n(B) D2DIAGNOS(B,65+ppZ6n)
#define jbotZ6p(B) D2DIAGNOS(B,65+ppZ6p)
#define jbotR1c(B) D2DIAGNOS(B,65+ppR1c)
#define jbotR1n(B) D2DIAGNOS(B,65+ppR1n)
#define jbotR1p(B) D2DIAGNOS(B,65+ppR1p)
#define jbotR2c(B) D2DIAGNOS(B,65+ppR2c)
#define jbotR3c(B) D2DIAGNOS(B,65+ppR3c)
#define jbotR6c(B) D2DIAGNOS(B,65+ppR6c)
#define jbotR6n(B) D2DIAGNOS(B,65+ppR6n)
#define jbotR6p(B) D2DIAGNOS(B,65+ppR6p)
#define jbotR6s(B) D2DIAGNOS(B,65+ppR6s)
#define jbotO3c(B) D2DIAGNOS(B,65+ppO3c)
#define jbotO3h(B) D2DIAGNOS(B,65+ppO3h)
#define jbotO5c(B) D2DIAGNOS(B,65+ppO5c)

#define PELRIVER(A,B) D2DIAGNOS(A,0+B)
#define jrivO2o(B) D2DIAGNOS(B,0+ppO2o)
#define jrivN1p(B) D2DIAGNOS(B,0+ppN1p)
#define jrivN3n(B) D2DIAGNOS(B,0+ppN3n)
#define jrivN4n(B) D2DIAGNOS(B,0+ppN4n)
#define jrivO4n(B) D2DIAGNOS(B,0+ppO4n)
#define jrivN5s(B) D2DIAGNOS(B,0+ppN5s)
#define jrivN6r(B) D2DIAGNOS(B,0+ppN6r)
#define jrivB1c(B) D2DIAGNOS(B,0+ppB1c)
#define jrivB1n(B) D2DIAGNOS(B,0+ppB1n)
#define jrivB1p(B) D2DIAGNOS(B,0+ppB1p)
#define jrivP1c(B) D2DIAGNOS(B,0+ppP1c)
#define jrivP1n(B) D2DIAGNOS(B,0+ppP1n)
#define jrivP1p(B) D2DIAGNOS(B,0+ppP1p)
#define jrivP1l(B) D2DIAGNOS(B,0+ppP1l)
#define jrivP1s(B) D2DIAGNOS(B,0+ppP1s)
#define jrivP2c(B) D2DIAGNOS(B,0+ppP2c)
#define jrivP2n(B) D2DIAGNOS(B,0+ppP2n)
#define jrivP2p(B) D2DIAGNOS(B,0+ppP2p)
#define jrivP2l(B) D2DIAGNOS(B,0+ppP2l)
#define jrivP3c(B) D2DIAGNOS(B,0+ppP3c)
#define jrivP3n(B) D2DIAGNOS(B,0+ppP3n)
#define jrivP3p(B) D2DIAGNOS(B,0+ppP3p)
#define jrivP3l(B) D2DIAGNOS(B,0+ppP3l)
#define jrivP4c(B) D2DIAGNOS(B,0+ppP4c)
#define jrivP4n(B) D2DIAGNOS(B,0+ppP4n)
#define jrivP4p(B) D2DIAGNOS(B,0+ppP4p)
#define jrivP4l(B) D2DIAGNOS(B,0+ppP4l)
#define jrivZ3c(B) D2DIAGNOS(B,0+ppZ3c)
#define jrivZ3n(B) D2DIAGNOS(B,0+ppZ3n)
#define jrivZ3p(B) D2DIAGNOS(B,0+ppZ3p)
#define jrivZ4c(B) D2DIAGNOS(B,0+ppZ4c)
#define jrivZ4n(B) D2DIAGNOS(B,0+ppZ4n)
#define jrivZ4p(B) D2DIAGNOS(B,0+ppZ4p)
#define jrivZ5c(B) D2DIAGNOS(B,0+ppZ5c)
#define jrivZ5n(B) D2DIAGNOS(B,0+ppZ5n)
#define jrivZ5p(B) D2DIAGNOS(B,0+ppZ5p)
#define jrivZ6c(B) D2DIAGNOS(B,0+ppZ6c)
#define jrivZ6n(B) D2DIAGNOS(B,0+ppZ6n)
#define jrivZ6p(B) D2DIAGNOS(B,0+ppZ6p)
#define jrivR1c(B) D2DIAGNOS(B,0+ppR1c)
#define jrivR1n(B) D2DIAGNOS(B,0+ppR1n)
#define jrivR1p(B) D2DIAGNOS(B,0+ppR1p)
#define jrivR2c(B) D2DIAGNOS(B,0+ppR2c)
#define jrivR3c(B) D2DIAGNOS(B,0+ppR3c)
#define jrivR6c(B) D2DIAGNOS(B,0+ppR6c)
#define jrivR6n(B) D2DIAGNOS(B,0+ppR6n)
#define jrivR6p(B) D2DIAGNOS(B,0+ppR6p)
#define jrivR6s(B) D2DIAGNOS(B,0+ppR6s)
#define jrivO3c(B) D2DIAGNOS(B,0+ppO3c)
#define jrivO3h(B) D2DIAGNOS(B,0+ppO3h)
#define jrivO5c(B) D2DIAGNOS(B,0+ppO5c)

#define qpcPPY(A,B) D3DIAGNOS(A,ppqpcPPY(B))
#define qncPPY(A,B) D3DIAGNOS(A,ppqncPPY(B))
#define qscPPY(A,B) D3DIAGNOS(A,ppqscPPY(B))
#define qlcPPY(A,B) D3DIAGNOS(A,ppqlcPPY(B))
#define qccPPY(A,B) D3DIAGNOS(A,ppqccPPY(B))
#define qpcMEZ(A,B) D3DIAGNOS(A,ppqpcMEZ(B))
#define qncMEZ(A,B) D3DIAGNOS(A,ppqncMEZ(B))
#define qpcMIZ(A,B) D3DIAGNOS(A,ppqpcMIZ(B))
#define qncMIZ(A,B) D3DIAGNOS(A,ppqncMIZ(B))
#define qpcOMT(A,B) D3DIAGNOS(A,ppqpcOMT(B))
#define qncOMT(A,B) D3DIAGNOS(A,ppqncOMT(B))
#define qscOMT(A,B) D3DIAGNOS(A,ppqscOMT(B))
#define qpcPBA(A,B) D3DIAGNOS(A,ppqpcPBA(B))
#define qncPBA(A,B) D3DIAGNOS(A,ppqncPBA(B))
#define sediPPY(A,B) D3DIAGNOS(A,ppsediPPY(B))
#define sediMIZ(A,B) D3DIAGNOS(A,ppsediMIZ(B))
#define sediMEZ(A,B) D3DIAGNOS(A,ppsediMEZ(B))
#define sunPPY(A,B) D3DIAGNOS(A,ppsunPPY(B))
#define eiPPY(A,B) D3DIAGNOS(A,ppeiPPY(B))
#define ELiPPY(A,B) D3DIAGNOS(A,ppELiPPY(B))


#ifdef INCLUDE_SEAICE






#endif

#define Q1c(A) D2STATE_BEN(A,ppQ1c)
#define Q1n(A) D2STATE_BEN(A,ppQ1n)
#define Q1p(A) D2STATE_BEN(A,ppQ1p)
#define Q11c(A) D2STATE_BEN(A,ppQ11c)
#define Q11n(A) D2STATE_BEN(A,ppQ11n)
#define Q11p(A) D2STATE_BEN(A,ppQ11p)
#define Q6c(A) D2STATE_BEN(A,ppQ6c)
#define Q6n(A) D2STATE_BEN(A,ppQ6n)
#define Q6p(A) D2STATE_BEN(A,ppQ6p)
#define Q6s(A) D2STATE_BEN(A,ppQ6s)
#define Q16c(A) D2STATE_BEN(A,ppQ16c)
#define Q16n(A) D2STATE_BEN(A,ppQ16n)
#define Q16p(A) D2STATE_BEN(A,ppQ16p)
#define Q16s(A) D2STATE_BEN(A,ppQ16s)

#define BenDetritus(A,B) D2STATE_BEN(:,ppBenDetritus(A,B))

#define totbenc(A) D2DIAGNOS_BEN(A,pptotbenc)
#define totbenn(A) D2DIAGNOS_BEN(A,pptotbenn)
#define totbenp(A) D2DIAGNOS_BEN(A,pptotbenp)
#define totbens(A) D2DIAGNOS_BEN(A,pptotbens)

#define qpcBOM(A,B) D2DIAGNOS_BEN(A,ppqpcBOM(B))
#define qncBOM(A,B) D2DIAGNOS_BEN(A,ppqncBOM(B))
#define qscBOM(A,B) D2DIAGNOS_BEN(A,ppqscBOM(B))



#endif



#if defined INCLUDE_BENCO2 && !defined BENTHIC_FULL
#undef INCLUDE_BENCO2
#endif
#if defined INCLUDE_BENPROFILES && !defined BENTHIC_FULL
#undef INCLUDE_BENPROFILES
#endif

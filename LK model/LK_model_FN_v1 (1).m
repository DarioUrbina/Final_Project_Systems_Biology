function dydt = LK_model_FN_v1(t,y,params)

GLY = y(1);
G1P = y(2);
G6P = y(3);
F6P = y(4);
FBP = y(5);
DHAP = y(6);
GAP = y(7);
BPG_13 = y(8);
PG_3 = y(9);
PG_2 = y(10);
PEP = y(11);
PYR = y(12);
LAC = y(13);
P_i = y(14);
ADP = y(15);
ATP = y(16);
AMP = y(17);
PCr = y(18);
Cr = y(19);
NADH = y(20);
NAD = y(21);


%% define parameters
% Glycogen Phosphorylase A
KeqGP1 = params(1,1);     %Keq 
KGLYf = params(2,1);      %KGLYf: this parameter is applied in the eqns in GPB part, not in A part.
KPi1 = params(3,1);        %KPi 
KiGLYf = params(4,1);        %KiGLY : we need to know the KiGLYf and KiGLYb for the 'Glycogen Phosphorylase A' part.
KiGLYb = params(4,1);
KiPi1 = params(5,1);      %KiPi
KGLYb = params(6,1);     %KGLYb
%KG1P1 = params(7,1);      %KG1P : there's no KG1P in the eqns of 'Glycogen Phosphorylase A'
KiG1P1 = params(8,1);     %KiG1P
%params(9:15,1) = 0;

% Glycogen Phosphorylase B
KPi2 = params(1,2);      % KPi    
KiPi2 = params(2,2);      % Kipi
KiGLYf2 = params(3,2);       % KiGLYf
KG1P2 = params(4,2);      % KG1P  
KiG1P2 = params(5,2);      % KiG1P
KiGLYb2 = params(6,2);      % KiGLYb
K_AMP = params(7,2);   % K_AMP
nH = params(8,2);     % nH 
KeqGP2 = params(9,2);    % Keq
%params(10:15,2) = 0;

% Phosphoglucomutase
KG1P3 = params(1,3);     % KG1P
KG6P3 = params(2,3);      % KG6P
KeqPGLM = params(3,3);      % Keq
%params(4:15,3) = 0;

% Phosphoglucoisomerase
KG6P4 = params(1,4);       % KG6P
KF6P4 = params(2,4);         % KF6P
KeqPGI = params(3,4);         % Keq  
%params(4:15,4)=0;

% Phosphofructokinase
KF6P5 = params(1,5);     %KF6P     
K_F6P = params(2,5);       %K_F6P     
KATP5 = params(3,5);     %KATP  
K_ATP = params(4,5);     %K_ATP
KFBP5 = params(5,5);     %KFBP     
K_FBP = params(6,5);     %K_FBP 
KADP5 = params(7,5);      %KADP
K_ADP = params(8,5);      %K_ADP
KiATP5 = params(9,5);     %KiATP 
KaAMP = params(10,5);    %KaAMP
d = params(11,5);    %d
e = params(12,5);    %e
L0 = params(13,5);      %Lo
%KeqPFK = params(14,5);  %Keq
%params(15,5) = 9.5e-5; 

% Aldolase
KFBP6 = params(1,6);      % KFBP
KDHAP6 = params(2,6);         % KDHAP
KGAP6 = params(3,6);         % KGAP
KeqALD = params(4,6);     % Keq
%params(5:15,6) = 0

% Triose Phosphate Isomerase
KGAP7 = params(1,7); %KGAP
KDHAP7 = params(2,7); %KDHAP
KeqTPI = params(3,7); %Keq
%params(4:15,7)= 0;

% Glyceraldehyde-3-Phosphate Dehydrogenase
KGAP8 = params(1,8);     %KGAP
KNAD8 = params(2,8);       %KNAD
KPi8 = params(3,8);       %KPi
K13BPG8 = params(4,8);     %K13BPG
KNADH8 = params(5,8);     %KNADH
KeqGAPDH = params(6,8);      %Keq

% Phosphoglycerate Kinase
K13BPG9 = params(1,9); %K13BPG
KADP9 = params(2,9); %KADP
K3PG9 = params(3,9); %K3PG
KATP9 = params(4,9); %KATP
KeqPGK = params(5,9);%Keq

% Phosphoglycerate Mutase
K3PG10 = params(1,10); %K3PG
K2PG10 = params(2,10); %K2PG
KeqPGM = params(3,10); %Keq

%Enloase
K2PG11 = params(1,11); %K2PG
KPEP11 = params(2,11); %KPEP
KeqENOL = params(3,11); %Keq

%Pyruvate Kinase
KPEP12 = params(1,12); %KPEP
KADP12 = params(2,12); %KADP
KPYR12 = params(3,12); %KPYR
KATP12 = params(4,12); %KATP
KeqPK = params(5,12);%Keq

%Lactate Dehydrogenase
KPYR13 = params(1,13); %KPYR
KNADH13 = params(2,13); %KNADH
KLAC = params(3,13); %KLAC
KNAD13 = params(4,13); %KNAD
KeqLDH = params(5,13);%Keq

%Creatine Kinase
KPCr = params(1,14); %KPCr
KiATP14 = params(2,14); %KiATP
KiADP = params(3,14); %KiADP
KiPCr = params(4,14); %KiPCr
KCr = params(5,14);%KCr
KeqCK = params(6,14);%Keq

%Adenylate Kinase
KAMP = params(1,15); %KAMP
KATP15 = params(2,15); %KATP
KADP15 = params(3,15); %KADP
KeqADK = params(6,14);%Keq

%Vmax
Vmax2 = params(1,16); %GPP B
Vmax1 = params(2,16); %GPP A
Vmax3 = params(3,16); %PGLM
%Vmax4 = params(4,16); %PGI
Vmaxr4 = params(4,16); %PGI
Vmax5 = params(5,16);%PFK
Vmax6 = params(6,16);%ALD
Vmax7 = params(7,16); %TPI
Vmax8 = params(8,16); %GAPDH
%Vmax9 = params(9,16);%PGK
Vmaxr9 = params(9,16);%PGK
Vmax10 = params(10,16);%PGM 
Vmax11 = params(11,16); %EN
Vmax12 = params(12,16); %PK
Vmax13 = params(13,16);%LDH
Vmax15 = params(14,16);%ADK
%Vmax14 = params(15,16);%CK
Vmaxr14 = params(15,16);%CK


%% define reaction rates

%Glycogen Phosphorylase
Vmaxr1 = Vmax1*KGLYb*KiG1P1/(KiGLYf*KPi1*KeqGP1);
Vmaxr2 = Vmax2*KiGLYb2*KG1P2/(KGLYf*KiPi2*KeqGP2);
VGPa = (Vmax1*(GLY*P_i/(KiGLYf*KPi1))-Vmaxr1*(GLY*G1P/(KGLYb*KiG1P1)))/(1+GLY/KiGLYf+P_i/KiPi1+GLY/KiGLYb+G1P/KiG1P1+GLY*P_i/(KGLYf*KiPi1)+GLY*G1P/(KGLYb*KiG1P1));
VGPb = (Vmax2*(GLY*P_i/(KiGLYf2*KPi2))-Vmaxr2*(GLY*G1P/(KiGLYb2*KG1P2)))/(1+GLY/KiGLYf2+P_i/KiPi2+GLY/KiGLYb2+G1P/KiG1P2+GLY*P_i/(KiGLYf2*KPi2)+GLY*G1P/(KiGLYb2*KG1P2))*(((AMP^nH)/K_AMP)/(1+(AMP^nH)/K_AMP));
fraca = 0.4;
fracb = 0.6;
fluxGP = fraca*VGPa+fracb*VGPb;


%Phosphoglucomutase
Vmaxr3 = Vmax3*KG6P3/(KG1P3*KeqPGLM);
VPGLM = ((Vmax3*G1P/KG1P3)-(Vmaxr3*G6P/KG6P3))/(1+G1P/KG1P3+G6P/KG6P3);


%Phosphoglucoisomerase
Vmax4 = Vmaxr4*KG6P4*KeqPGI/KF6P4;
%Vmaxr4 = Vmax4*KF6P4/(KG6P4*KeqPGI);
VPGI = ((Vmax4*G6P/KG6P4)-(Vmaxr4*F6P/KF6P4))/(1+G6P/KG6P4+F6P/KF6P4);


%Phosphofructokinase
delta = (1+F6P/KF6P5)*(1+ATP/KATP5)+ADP/KADP5+FBP/KFBP5*(1+ADP/KADP5);
delta1 = (1+F6P/K_F6P)*(1+ATP/K_ATP)+ADP/K_ADP+FBP/K_FBP*(1+ADP/K_ADP);
alpha = KF6P5*KATP5/(K_F6P*K_ATP);
L = L0*(((1+ATP/KiATP5)/(1+d*ATP/KiATP5))*((1+e*AMP/KaAMP)/(1+AMP/KaAMP)))^4;
Vmaxr5 = Vmax5*KADP5*KFBP5/(KATP5*KF6P5);
VPFK = (Vmax5*ATP*F6P/(KATP5*KF6P5)-Vmaxr5*ADP*FBP/(KADP5*KFBP5))/delta*((1+alpha*L*(delta1/delta)^3)/(1+L*(delta1/delta)^4));


%Aldolase
Vmaxr6 = Vmax6*KDHAP6*KGAP6/(KFBP6*KeqALD);
VALD = (Vmax6*FBP/KFBP6-Vmaxr6*DHAP*GAP/(KDHAP6*KGAP6))/(1+FBP/KFBP6+DHAP/KDHAP6+GAP/KGAP6);


%Triose Phosphate Isomerase
Vmaxr7 = Vmax7*KDHAP7/(KGAP7*KeqTPI);
VTPI = (Vmax7*GAP/KGAP7-Vmaxr7*DHAP/KDHAP7)/(1+GAP/KGAP7+DHAP/KDHAP7);


%Glyceraldehyde-3-Phosphate Dehydrogenase (GAPDH)
DGAPDH = 1+GAP/KGAP8+NAD/KNAD8+P_i/KPi8+GAP*NAD/(KGAP8*KNAD8)+GAP*NAD*P_i/(KGAP8*KNAD8*KPi8)+BPG_13/K13BPG8+NADH/KNADH8+BPG_13*NADH/(K13BPG8*KNADH8);
Vmaxr8 = Vmax8*K13BPG8*KNADH8/(KGAP8*KNAD8*KPi8*KeqGAPDH);
VGAPDH = (Vmax8*GAP*NAD*P_i/(KGAP8*KNAD8*KPi8)-Vmaxr8*BPG_13*NADH/(K13BPG8*KNADH8))/DGAPDH;


%Phosphoglycerate Kinase
Vmax9 = Vmaxr9*K13BPG9*KADP9*KeqPGK/(K3PG9*KATP9);
%Vmaxr9 = Vmax9*K3PG9*KATP9/(K13BPG9*KADP9*KeqPGK);
VPGK = (Vmax9*BPG_13*ADP/(K13BPG9*KADP9)-Vmaxr9*PG_3*ATP/(K3PG9*KATP9))/(1+BPG_13/K13BPG9+ADP/KADP9+BPG_13*ADP/(K13BPG9*KADP9)+PG_3/K3PG9+ATP/KATP9+PG_3*ATP/(K3PG9*KATP9));


%Phosphoglyceromutase
Vmaxr10 = Vmax10*K2PG10/(K3PG10*KeqPGM);
VPGM = (Vmax10*PG_3/K3PG10-Vmaxr10*PG_2/K2PG10)/(1+PG_3/K3PG10+PG_2/K2PG10);


%Enolase
Vmaxr11 = Vmax11*KPEP11/(K2PG11*KeqENOL);
VENOL = (Vmax11*PG_2/K2PG11-Vmaxr11*PEP/KPEP11)/(1+PG_2/K2PG11+PEP/KPEP11);


%Pyruvate Kinase
Vmaxr12 = Vmax12*KATP12*KPYR12/(KPEP12*KADP12*KeqPK);
VPK = (Vmax12*PEP*ADP/(KPEP12*KADP12)-Vmaxr12*PYR*ATP/(KPYR12*KATP12))/(1+PEP/KPEP12+ADP/KADP12+PEP*ADP/(KPEP12*KADP12)+PYR/KPYR12+ATP/KATP12+PYR*ATP/(KPYR12*KATP12));


%Lactate Dehydrogenase (LDH)
Vmaxr13 = Vmax13*KLAC*KNAD13/(KPYR13*KNADH13*KeqLDH);
VLDH = (Vmax13*PYR*NADH/(KPYR13*KNADH13)-Vmaxr13*LAC*NAD/(KLAC*KNAD13))/(1+PYR/KPYR13+NADH/KNADH13+PYR*NADH/(KPYR13*KNADH13)+LAC/KLAC+NAD/KNAD13+LAC*NAD/(KLAC*KNAD13));


%ATP Reactions:Creatine Kinase
Vmax14 = Vmaxr14*KiATP14*KCr*KeqCK/(KiADP*KPCr);
%Vmaxr14 = Vmax14*KiADP*KPCr/(KiATP14*KCr*KeqCK);
VCK = (Vmaxr14*ATP*Cr/(KiATP14*KCr)-Vmax14*ADP*PCr/(KiADP*KPCr))/(1+ADP/KiADP+PCr/KiPCr+ADP*PCr/(KiADP*KPCr)+ATP/KiATP14+ATP*Cr/(KiATP14*KCr));


%ATP Reactions:Adenylate Kinase
Vmaxr15 = Vmax15*KADP15*KADP15/(KATP15*KAMP*KeqADK);

VADK = (Vmax15*ATP*AMP/(KATP15*KAMP)-Vmaxr15*ADP*ADP/(KADP15^2))/(1+ATP/KATP15+AMP/KAMP+ATP*AMP/(KATP15*KAMP)+2*ADP/KADP15+(ADP^2)/(KADP15^2));

%ATP Reactions:ATPase
k = 7.5;  % k has three values, 0.075 0.75 7.5. I'm not sure which one we should use.
VATPase = k*ATP;
output = 0.2*LAC;  %  didnt find exactly the eqn for 'output' in the paper, not 100% sure the output would equal to 0.2*LAC 

%% define ODEs
dydt(1,1)= -fluxGP; %equation for GLY
dydt(2,1)= fluxGP-VPGLM; %equation for G1P
dydt(3,1)= VPGLM-VPGI; %equation for G6P
dydt(4,1)= VPGI-VPFK; %equation for F6P
dydt(5,1)= VPFK-VALD; %equation for FBP
dydt(6,1)= VALD+VTPI; %equation for DHAP
dydt(7,1)= VALD-VTPI-VGAPDH; %equation for GAP
dydt(8,1)= VGAPDH-VPGK; %equation for BPG_13
dydt(9,1)= VPGK-VPGM; %equation for PG_3
dydt(10,1)= VPGM-VENOL; %equation for PG_2
dydt(11,1)= VENOL-VPK; %equation for PEP
dydt(12,1)= VPK-VLDH; %equation for PYR 
dydt(13,1)= VLDH-output; %equation for LAC
dydt(14,1)= -fluxGP-VGAPDH+VATPase; %equation for P_i
dydt(15,1)= VPFK-VPGK-VPK+2*VADK+VCK+VATPase; %equation for ADP
dydt(16,1)= -VPFK+VPGK+VPK-VADK-VCK-VATPase; %equation for ATP
dydt(17,1)= -VADK; %equation for AMP
dydt(18,1)= VCK; %equation for PCr
dydt(19,1)= -VCK; %equation for Cr
dydt(20,1)= VGAPDH-VLDH; %equation for NADH
dydt(21,1)= -(VGAPDH-VLDH); %equation for NAD

return

%% ---- Diver
% Specify parameters
% Glycogen Phosphorylase A
params(1,1) = 0.42;     %Keq 
params(2,1) = 1.7;      %KGLYf 
params(3,1) = 4;        %KPi 
params(4,1) = 2;        %KiGLY 
params(5,1) = 4.7;      %KiPi
params(6,1) = 0.15;     %KGLYf
params(7,1) = 2.7;      %KG1    P 
params(8,1) = 10.1;     %KiG1P
%params(9,1) = 0.42;     %Keq
params(9:15,1) = 0;

% Glycogen Phosphorylase B
params(1,2) = 0.2;      % KPi    
params(2,2) = 4.6;      % Kipi
params(3,2) = 15;       % KiGLYf
params(4,2) = 1.5;      % KG1P  
params(5,2) = 7.4;      % KiG1P
params(6,2) = 4.4;      % KiGLYb
params(7,2) = 1.9e-6;   % K_AMP
params(8,2) = 1.75;     % nH 
params(9,2) = 16.62;    % Keq
params(10:15,2) = 0;

% Phosphoglucomutase

params(1,3) = 0.063;     % KG1P
params(2,3) = 0.03;      % KG6P
params(3,3) = 0.45;      % Keq
params(4:15,3) = 0;
% Phosphoglucoisomerase

params(1,4) = 0.48;       % KG6P
params(2,4) = 0.119;         % KF6P
params(3,4) = 242;         % Keq  
params(4:15,4)=0;


% Phosphofructokinase
params(1,5) = 0.18;     %KF6P     
params(2,5) = 20;       %K_F6P     
params(3,5) = 0.08;     %KATP  
params(4,5) = 0.25;     %K_ATP
params(5,5) = 4.02;     %KFBP     
params(6,5) = 4.02;     %K_FBP 
params(7,5) = 2.7;      %%KADP
params(8,5) = 2.7;      %K_ADP
params(9,5) = 0.87;     %KiATP 
params(10,5) = 0.06;    %KaAMP
params(11,5) = 0.01;    %d
params(12,5) = 0.01;    %e
params(13,5) = 13;      %Lo
params(14,5) = 9.5e-5;  %Keq
params(15,5) = 9.5e-5; 

% Aldolase

params(1,6) = 0.05;      % KFBP
params(2,6) = 2;         % KDHAP
params(3,6) = 1;         % KGAP
params(4,6) = 0.052;     % Keq
params(5:15,6) = 0;

% Triose Phosphate Isomerase
params(1,7)=0.32;   %kGAP
params(2,7)=0.61;   %kDHAP
params(3,7)=0.089;  %keq
params(4:15,7)=0; 

% Glyceraldehyde-3-Phosphate Dehydrogenase
params(1,8)=0.0025;   %kGAP
params(2,8)=0.09;     %kNAD
params(3,8)=0.29;     %kPi
params(4,8)=0.0008;   %k13BPG
params(5,8)=0.0033;   %kNADH
params(6,8)=57109;    %keq
params(7:15,8)=0; 

% Phosphoglycerate Kinase
params(1,9)=0.002;   %k13BPG
params(2,9)=0.008;   %kADP
params(3,9)=1.2;     %k3PG
params(4,9)=0.35;    %kATP
params(5,9)=0.18;    %keq
params(6:15,9)=0; 

% Phosphoglycerate Mutase
params(1,10)=0.2;    %k3PG
params(2,10)=0.014;  %k2PG
params(3,10)=0.49;   %keq
params(4:15,10)=0; 

% Enloase  
params(1,11)=0.1;    %k2PG
params(2,11)=0.37;   %kPEP
params(3,11)=10304;  %keq
params(4:15,11)=0; 

% Pyruvate Kinase 
params(1,12)=0.08;   %kPEP
params(2,12)=0.3;    %kADP
params(3,12)=7.05;   %kPYR
params(4,12)=1.13;   %kATP
params(5,12)=16198;  %keq
params(6:15,12)=0; 

% Lactate Dehydrogenase
params(1,13)=0.335;   %kPYR
params(2,13)=0.002;   %kNADH
params(3,13)=17;      %kLAC
params(4,13)=0.849;   %kNAD
params(5,13)=233;     %keq
params(6:15,13)=0;

% Creatine Kinase
params(1,14)=1.11;   %kPcr
params(2,14)=3.5;    %kiATP
params(3,14)=0.135;  %kiADP
params(4,14)=3.9;    %kiPCR
params(5,14)=3.8;    %kCR
params(6,14)=2.21;   %keq
params(7:15,14)=0; 

% Adenylate Kinase
params(1,15)=0.32;   %kAMP
params(2,15)=0.27;   %kATP
params(3,15)=0.35;   %kADP
params(4:15,15)=0; 

% TABLE 2
params(1,16)=0.03e3;  %GPP B
params(2,16)=0.02e3;  %GPP A
params(3,16)=0.48e3;  %PGLM
params(4,16)=0.88e3;  %PGI
params(5,16)=0.056e3; %PFK
params(6,16)=0.104e3; %ALD
params(7,16)=12e3;    %TPI
params(8,16)=1.265e3; %GAPDH
params(9,16)=1.12e3;  %PGK
params(10,16)=1.12e3; %PGM 
params(11,16)=0.192e3;%EN
params(12,16)=1.44e3; %PK
params(13,16)=1.92e3; %LDH
params(14,16)=0.88e3; %ADK
params(15,16)=0.5e3;  %CK

% specify solver details
step = 0.01;
stop = 3000;
tspan = [0:step:stop];
maxY = 350;
diffeqn = 21;
options = odeset('RelTol',1e-6);

% Specify initial values
initvalue = zeros(diffeqn,1);
initvalue(1,1) = 112; %GLY
initvalue(2,1) = 0.0589; %G1P
initvalue(3,1) = 0.75;  %G6P
initvalue(4,1) = 0.228;  %F6P
initvalue(5,1) = 0.0723;  %FBP
initvalue(6,1) = 0.0764;  %DHAP
initvalue(7,1) = 0.0355;  %GAP
initvalue(8,1) = 0.065;  %13BPG
initvalue(9,1) = 0.052;  %3PG
initvalue(10,1) = 0.005;  %2PG
initvalue(11,1) = 0.0194;  %PEP
initvalue(12,1) = 0.0994;  %PYR
initvalue(13,1) = 1e3;  %LAC
initvalue(16,1) = 1e3;  %ATP
initvalue(15,1) = 1e3;  %ADP
initvalue(17,1) = 2e-5;  %AMP
initvalue(14,1) = 1e3;  %Pi
initvalue(18,1) = 34.67;  %PCr
initvalue(21,1) = 0.5;  %NAD
initvalue(19,1) = 5.33;  %Cr
initvalue(20,1) = 0.5e-3;  %NADH


% Refer to file with ODEs and call solver
[tsim_3a, results_3a] = ode15s(@LK_model_FN_v1,tspan,initvalue,options,params);
%[tsim_3a, results_3a] = ode15s(@LK_model_FN,tspan,initvalue,options,params);
results = results_3a.*10^-3;

% plot results
figure(1); 
plot(tsim_3a,log(results(:,15)),'linewidth',4);
%axis([0 500 -8 1]);
%xlim([0 100])
legend('ADP')
xlabel('time(min)')
ylabel('concentration(M)')
title(['When V1 is ',num2str(params(1,1)),'nM/sec']);

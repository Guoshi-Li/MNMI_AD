
% Analyze intra-regional and inter-regoinal connection strengths for NC,
% MCI and AD subjects
% The connections strengths are estimated by MNMI and stored in the OUTPUT
% directory
% Set "Flag" to 1 to compare between NC and MCI; 2 to compare between NC
% and AD; 3 to compare between MCI and AD
% Written by Guoshi Li (guoshi_li@med.unc.edu)

clc;
close all;
clear all;

Flag = 3; % 1: NC vs. MCI; 2: NC vs AD; 3: MCI vs AD

load AD_SC10;        % structural connectivity
load ADLabelNet2;    % network arrangement of the DK atlas

NS = 48;  % Number of subjects in each group

Threshold = 0.0;

NIC = NP;     % Number of inter-regional connections

NR = 46;
NP = 2*NR+1+NIC;        % Total number of parameters
NN = 4;                 % Number of networks

q = 0.05;


% load NC/MCI/AD estimation output
for i = 1:NS
    
  % load NC output        
    subject = ['NC' int2str(i)];        
    datafile = ['..\OUTPUT\' subject];
    load (datafile);
    FVAL_NC(i) = -fval;
    W0_NC(i,:) = x;

  % load MCI output  
    subject = ['MCI' int2str(i)];        
    datafile = ['..\OUTPUT\' subject];
    load (datafile);
    FVAL_MCI(i) = -fval;
    W0_MCI(i,:) = x;

  % load AD output  
    subject = ['AD' int2str(i)];        
    datafile = ['..\OUTPUT\' subject];
    load (datafile);
    FVAL_AD(i) = -fval;
    W0_AD(i,:) = x;    

end



%% ========================================================================
% Switch the position of PCC (17,18) with cACC (15, 16) so PCC goes to the
% DMN

% For excitatory connections
W1_NC  = W0_NC(:,1:14);
W1_MCI = W0_MCI(:,1:14);
W1_AD  = W0_AD(:,1:14);

CACC1 = W0_NC(:,15:16);
CACC2 = W0_MCI(:,15:16);
CACC3 = W0_AD(:,15:16);

PCC1 = W0_NC(:,17:18);
PCC2 = W0_MCI(:,17:18);
PCC3 = W0_AD(:,17:18);

W1_NC  = [W1_NC  PCC1 CACC1 W0_NC(:,19:end)];
W1_MCI = [W1_MCI PCC2 CACC2 W0_MCI(:,19:end)];
W1_AD  = [W1_AD  PCC3 CACC3 W0_AD(:,19:end)];


% For inhibitory connections
W2_NC  = W1_NC(:,1:(NR+14));
W2_MCI = W1_MCI(:,1:(NR+14));
W2_AD  = W1_AD(:,1:(NR+14));

CACC1 = W0_NC(:,(NR+15):(NR+16));
CACC2 = W0_MCI(:,(NR+15):(NR+16));
CACC3 = W0_AD(:,(NR+15):(NR+16));

PCC1 = W0_NC(:,(NR+17):(NR+18));
PCC2 = W0_MCI(:,(NR+17):(NR+18));
PCC3 = W0_AD(:,(NR+17):(NR+18));

W_NC  = [W2_NC  PCC1 CACC1 W1_NC(:,(NR+19):end)];
W_MCI = [W2_MCI PCC2 CACC2 W1_MCI(:,(NR+19):end)];
W_AD  = [W2_AD  PCC3 CACC3 W1_AD(:,(NR+19):end)];


%% ========================================================================
% Flag = 1: NC vs. MCI; Flag = 2: NC vs. AD; Flag = 3: MCI vs. AD

if (Flag == 1)
   
   W_NM =  W_NC;
   W_MD =  W_MCI;
   FVAL_NM = FVAL_NC;
   FVAL_MD = FVAL_MCI;
   W_ALL = [W_NC; W_MCI];

elseif (Flag == 2)

   W_NM =  W_NC;
   W_MD =  W_AD;
   FVAL_NM = FVAL_NC;
   FVAL_MD = FVAL_AD;
   W_ALL = [W_NC; W_AD];

else

   W_NM =  W_MCI;
   W_MD =  W_AD;
   FVAL_NM = FVAL_MCI;
   FVAL_MD = FVAL_AD;
   W_ALL = [W_MCI; W_AD];       

end

%==========================================================================

Wee_NM = W_NM(:,1:NR);
Wie_NM = W_NM(:,NR+1:2*NR);
Wgc_NM = W_NM(:,2*NR+1:NP-1);
SPI_NM = W_NM(:,NP);

Wee_MD = W_MD(:,1:NR);
Wie_MD = W_MD(:,NR+1:2*NR);
Wgc_MD = W_MD(:,2*NR+1:NP-1);
SPI_MD = W_MD(:,NP);


% Group into different networks

WeeN_NC = zeros(NS, NN);
WieN_NC = zeros(NS, NN);

WeeN_MDD = zeros(NS, NN);
WieN_MDD = zeros(NS, NN);

for i = 1:NN
    
    idx = NID{i};
    
    % For NC/MCI
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_NM(:,idx);
    WTi = Wie_NM(:,idx);
    
    WeeN_NM(:,i) = mean(WTe, 2);
    WieN_NM(:,i) = mean(WTi, 2);
    
   % For MCI/AD 
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_MD(:,idx);
    WTi = Wie_MD(:,idx);
    
    WeeN_MD(:,i) = mean(WTe, 2);
    WieN_MD(:,i) = mean(WTi, 2);
    
        
end

WgcN_NM  = mean(Wgc_NM, 2);
WgcN_MD  = mean(Wgc_MD, 2);

WN_NM = [WeeN_NM  WieN_NM ];
WN_MD = [WeeN_MD  WieN_MD ];

% WN_NM = [WeeN_NM  WieN_NM  SPI_NM];
% WN_MD = [WeeN_MD  WieN_MD  SPI_MD];

% Average across subjects
AWeeN_NM = mean(WeeN_NM);
AWeeN_MD = mean(WeeN_MD);

AWieN_NM = mean(WieN_NM);
AWieN_MD = mean(WieN_MD);

AWeeN_NM_MD = [AWeeN_NM; AWeeN_MD];
AWieN_NM_MD = [AWieN_NM; AWieN_MD];

%==========================================================================
% Maximal correlation 
MAXC_NM = max(FVAL_NM);
MAXC_MD = max(FVAL_MD);

% Average correlation and error
AC_NM = mean(FVAL_NM);
AC_MD = mean(FVAL_MD);


% Average weight
AW_NM = mean(W_NM);
AW_MD = mean(W_MD);

STD_NM = std(W_NM);
STD_MD = std(W_MD);

% Put NC and MDD together
AW_NM_MD  = [AW_NM; AW_MD];
STD_NM_MD = [STD_NM; STD_MD];

% Divide the parameters into four sets
WEE = AW_NM_MD(:, 1:NR);
WIE = AW_NM_MD(:, NR+1:2*NR);
WGC = AW_NM_MD(:, 2*NR+1:NP-1);
WPI = AW_NM_MD(:, NP);

%==========================================================================
m = 1;

GC01 = zeros(NR,NR);
GC02 = zeros(NR,NR);

WGC1 = WGC(1,:);
WGC2 = WGC(2,:);

for j=1:NR
    for i=1:NR
       
       if (MAP(i,j)==1)
           GC01(i,j) = WGC1(m);
           GC02(i,j) = WGC2(m);
           m = m+1;
       end      
        
    end  
end


% =========================================================================
% Exchange cACC (15-16) with PCC (17-18)
% GC1
CACC1 = GC01(15:16,:);
PCC1  = GC01(17:18,:);

GC1 = [GC01(1:14, :); PCC1; CACC1; GC01(19:end,:)];

CACC2 = GC1(:,15:16);
PCC2  = GC1(:,17:18);

GC1 = [GC1(:, 1:14) PCC2 CACC2 GC1(:,19:end)];

% GC2
CACC1 = GC02(15:16,:);
PCC1  = GC02(17:18,:);

GC2 = [GC02(1:14, :); PCC1; CACC1; GC02(19:end,:)];

CACC2 = GC2(:,15:16);
PCC2  = GC2(:,17:18);

GC2 = [GC2(:, 1:14) PCC2 CACC2 GC2(:,19:end)];


%% ========================================================================
%                Statistical Tests        
% =========================================================================

TPC = [];
TPCR = [];
TPCL = [];

TPC_STAT = [];
TPCR_STAT = [];
TPCL_STAT = [];

TPN = [];
TPNR = [];
TPNL = [];

TPN_STAT = [];
TPNR_STAT = [];
TPNL_STAT = [];

% For each parameter
for k = 1:NP

CV = W_NM(:, k);
DV = W_MD(:, k);

[h,p,ci,stats1]  = ttest2(CV, DV);
[h,pr,ci,stats2] = ttest2(CV, DV,'Tail', 'right');
[h,pl,ci,stats3] = ttest2(CV, DV,'Tail', 'left');

TPC = [TPC; p];
TPCR = [TPCR; pr];
TPCL = [TPCL; pl];

TPC_STAT = [TPC_STAT; stats1.tstat];
TPCR_STAT = [TPCR_STAT; stats2.tstat];
TPCL_STAT = [TPCL_STAT; stats3.tstat];


end


%% ========================================================================
% For network comparison (network-averaged recurrent excitation and
% inhibition)
% =========================================================================

for k =1:(2*NN)

CV = WN_NM(:, k);
DV = WN_MD(:, k);

[h,p,ci,stats1]  = ttest2(CV, DV);
[h,pr,ci,stats2] = ttest2(CV, DV,'Tail', 'right');
[h,pl,ci,stats3] = ttest2(CV, DV,'Tail', 'left');

TPN = [TPN; p];
TPNR = [TPNR; pr];
TPNL = [TPNL; pl];

TPN_STAT = [TPN_STAT; stats1.tstat];
TPNR_STAT = [TPNR_STAT; stats2.tstat];
TPNL_STAT = [TPNL_STAT; stats3.tstat];


end

IC = find(TPC<0.05);
IN = find(TPN<0.05);

% [h,pgc,ci,stats1] = ttest2(WgcN_NM, WgcN_MD);

%% ========================================================================

if (Flag == 1)
    disp('Statistical comparison between NC and MCI')
elseif(Flag == 2)
    disp('Statistical comparison between NC and AD')
else
    disp('Statistical comparison between MCI and AD')
end


% Divide the T-statistics and p values into (1) Recurrent excitation; 
% (2) Recurrent inhibition; and (3) Inter-regional EC

TPC1 = TPC(1:NR);
TPC2 = TPC(NR+1:2*NR);
TPC3 = TPC(2*NR+1:NP-1);

TVAL1 = TPC_STAT(1:NR);
TVAL2 = TPC_STAT(NR+1:2*NR);
TVAL3 = TPC_STAT(2*NR+1:NP-1);

% Bonferroni correction
BON1 = TPC1*NR;
BON2 = TPC2*NR;
BON3 = TPC3*length(TPC3);

% FDR correction
FDR1 = mafdr(TPC1, 'BHFDR',1);
FDR2 = mafdr(TPC2, 'BHFDR',1);
FDR3 = mafdr(TPC3, 'BHFDR',1);


%% ========================================================================
% % Bonferroni Correction
% disp('Bonferroni correction:')
% B_thresh1 = 0.05/NR;
% B_thresh2 = 0.05/NR;
% B_thresh3 = 0.05/length(TPC3);
% 
% BT1 = find(TPC1<B_thresh1);
% BT2 = find(TPC2<B_thresh1);
% BT3 = find(TPC3<B_thresh3);
% 
% BT1'
% BT2'
% BT3'


% Find out the corresponding ROIs for significant ECs 

% Intra-reginal excitation & Inhibition
EI = find(IC<=NR);
II = find(IC>NR & IC<=2*NR);

RE1 = IC(EI);     % ROIs within 46 network
RI1 = IC(II)-NR;  % ROIs within 46 network

RE2 = R4(RE1);    % ROIs within 84 network
RI2 = R4(RI1);    % ROIs within 84 network

RE1 = RE1';
RI1 = RI1';

PS1 = TPC1(RE1);
PS2 = TPC2(RI1);

TV1 = TVAL1(RE1);
TV2 = TVAL2(RI1);

BP1 = BON1(RE1);
BP2 = BON2(RI1);

CP1 = FDR1(RE1);
CP2 = FDR2(RI1);

% Recurrent excitation
fprintf('\n')
disp('Comparision Results for Recurrent Excitation:')
fprintf('\n')
disp('The ROIs with significant recurrent excitation are (p<0.05):');
RE1
% RE2

disp('The T-Values for significant intra-regional recurrent excitation:')
TV1'

disp('The uncorrected significant p values for intra-regional recurrent excitation:')
PS1'

disp('The BON-CORRECTED significant p values for intra-regional recurrent excitation:')
BP1'

disp('The FDR-CORRECTED significant p values for intra-regional recurrent excitation:')
CP1'

disp('FDR Correction for intra-regional excitatory connections:');
[pID1,pN1] = gretna_FDR(TPC1,q)


%% ========================================================================
% Recurrent inhibition
fprintf('\n')
disp('Comparision Results for Recurrent Inhibition:')
fprintf('\n')

disp('The ROIs with significant recurrent inhibition are (p<0.05):');
RI1
% RI2

disp('The T-Values for significant intra-regional recurrent inhibition:')
TV2'

disp('The uncorrected significant p values for intra-regional recurrent inhibition:')
PS2'

disp('The BON-CORRECTED significant p values for intra-regional recurrent inhibition:')
BP2'

disp('The FDR-CORRECTED significant p values for intra-regional recurrent inhibition:')
CP2'

disp('FDR Correction for intra-regional inhibitory connections:');
[pID2,pN2] = gretna_FDR(TPC2,q)



%% ========================================================================
% Network-averaged recurrent excitation and inhibition
%==========================================================================
fprintf('\n')
disp('Comparision Results for Network-averaged recurrent excitaiton and inhibition:')
fprintf('\n')

FDRN = mafdr(TPN, 'BHFDR',1);
FDRN';

TVAL_NET = TPN_STAT(IN);

disp('The index of significant parameters are:');
IN'

disp('The T-value for the significant network connections are:')
TVAL_NET'

disp('The uncorrected p-values less than 0.05 are:');
TPN(IN)'

disp('The FDR-corrected p values for network connections are:')
FDRN(IN)'

disp('FDR Correction or network connections:');
[pID4, pN4] = gretna_FDR(TPN,q)





%% ========================================================================
%%  Inter-regional Connections
%% ========================================================================

% Significant inter-reginal connections
IX = find(IC>2*NR & IC<=NP-1);
IR = IC(IX)-2*NR;
IR = IR';

m = 0;

IRX = [];
IRY = [];

for j = 1:NR
    for i = 1:NR
   
         if (MAP(i,j)==1)
             m = m+1;
             
             ch = find(IR==m);
             
             if ~isempty(ch)
                  
                 % ======================================

                 tmpj = j;
                 tmpi = i;

                 % Exchange cACC (15-16) with PCC (17-18)

                 if(tmpj==15)
                    tmpj=17;
                 elseif(tmpj==16)
                    tmpj=18;  
                 elseif(tmpj==17) 
                    tmpj=15;
                 elseif(tmpj==18)
                    tmpj=16;
                 else
                    tmpj = j; 
                 end
                      

                 if(tmpi==15)
                    tmpi=17;
                 elseif(tmpi==16)
                    tmpi=18;  
                 elseif(tmpi==17) 
                    tmpi=15;
                 elseif(tmpi==18)
                    tmpi=16;
                 else
                    tmpi = i; 
                 end

                 % =======================================

                 IRX = [IRX tmpj];
                 IRY = [IRY tmpi];
             
             end
             
         end
        
        
    end
end

IRX2 = R4(IRX);
IRY2 = R4(IRY);

fprintf('\n')
disp('Comparision Results for Inter-regional Connections:')
fprintf('\n')

disp('The indices for significant inter-regional connections are (among 212 parameters):');
IR

disp('The indices for significant inter-regional connections (46 ROIs) are (p<0.05; uncorrected):');
IRX
IRY

% disp('The significant inter-regional connections (84 ROIs) are:');
% IRX2
% IRY2

disp('FDR Correction for inter-regional connections:');
[pID3,pN3] = gretna_FDR(TPC3,q)


% =========================================================================

% Exchange cACC (15-16) with PCC (17-18)
CACC1 = MAP(15:16,:);
PCC1  = MAP(17:18,:);

MAP2 = [MAP(1:14, :); PCC1; CACC1; MAP(19:end,:)];

CACC2 = MAP2(:,15:16);
PCC2  = MAP2(:,17:18);

MAP2 = [MAP2(:,1:14) PCC2 CACC2 MAP2(:,19:end)];


SMAP1 = MAP2;

L = length(IRX);

for k = 1:L
   SMAP1(IRY(k), IRX(k)) = 2;
end

 if (Flag == 1)
    save PARA_NC_MCI.mat W_NM W_MD  FVAL_NM FVAL_MD IR IRX IRY MAP2 SMAP1 TPC1 TPC2 TPC3 RE1 RE2 RI1 RI2 
 elseif (Flag==2)
    save PARA_NC_AD.mat W_NM W_MD   FVAL_NM FVAL_MD IR IRX IRY MAP2 SMAP1 TPC1 TPC2 TPC3 RE1 RE2 RI1 RI2 
 else
    save PARA_MCI_AD.mat W_NM W_MD  FVAL_NM FVAL_MD IR IRX IRY MAP2 SMAP1 TPC1 TPC2 TPC3 RE1 RE2 RI1 RI2  
 end

 
disp('The mean fitness value for NC is:')
mean(FVAL_NM)

disp('The mean fitness value for MCI/AD is:')
mean(FVAL_MD)


%% ========================================================================
    
% Recurrent excitation
figure;
b=bar(WEE', 'FaceColor','flat' );
axis([0, 47, 2.5, 3.5])
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', [])
title('Recurrent Excitation', 'FontSize', 16);
xlabel('ROI');
ylabel('Weight');
box('off');
if(Flag==1)
   legend('NC', 'MCI');
elseif (Flag==2)
   legend('NC', 'AD');       
else
   legend('MCI', 'AD');       
end

% Recurrent Inhibition
figure;
b=bar(WIE', 'FaceColor','flat');
axis([0, 47, 2.5, 3.5])
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', [])
title('Recurrent Excitation', 'FontSize', 16);
xlabel('ROI');
ylabel('Weight');
box('off');
if(Flag==1)
   legend('NC', 'MCI');
elseif (Flag==2)
   legend('NC', 'AD');       
else
   legend('MCI', 'AD');       
end

% Network-averaged recurrent excitation
figure;
bar(AWeeN_NM_MD');
axis([0.5, 4.5, 2.0, 3.5])
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', [])
title('Network-averaged Recurrent Excitation', 'FontSize', 16);
xlabel('ROI');
ylabel('Weight');
box('off');
if(Flag==1)
   legend('NC', 'MCI');
elseif (Flag==2)
   legend('NC', 'AD');       
else
   legend('MCI', 'AD');       
end


figure;
bar(AWieN_NM_MD');
axis([0.5, 4.5, 2.0, 3.5])
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', [])
title('Network-averaged Recurrent Inhibition', 'FontSize', 16);
xlabel('ROI');
ylabel('Weight');
box('off');
if(Flag==1)
   legend('NC', 'MCI');
elseif (Flag==2)
   legend('NC', 'AD');       
else
   legend('MCI', 'AD');       
end


CLIM = [-0.5 1.0];

% Inter-regional EC
figure;
b=bar(WGC', 'FaceColor','flat');
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', [])
title('Inter-regional EC', 'FontSize', 18);
xlabel('Connection');
ylabel('EC');
box('off');
if(Flag==1)
   legend('NC', 'MCI');
elseif (Flag==2)
   legend('NC', 'AD');       
else
   legend('MCI', 'AD');       
end

figure;
imagesc(GC1, CLIM);
colormap jet;   % jet
title('NC', 'FontSize', 14);
set(gca, 'Fontsize', 14);
xlabel('ROI')
ylabel('ROI')
colorbar;
if(Flag==1)
   title('NC', 'FontSize', 18);
elseif (Flag==2)
   title('MCI', 'FontSize', 18);     
else
   title('MCI', 'FontSize', 18);    
end

figure;
imagesc(GC2, CLIM);
colormap jet;  % jet
set(gca, 'Fontsize', 14);
xlabel('ROI')
ylabel('ROI')
colorbar;
if(Flag==1)
   title('MCI', 'FontSize', 18);
elseif (Flag==2)
   title('AD', 'FontSize', 18);     
else
   title('AD', 'FontSize', 18);    
end


% Spontaneous input
figure;
b=bar(WPI', 'FaceColor','flat');
set(gca, 'FontSize', 14);
title('Spontaneous Input', 'FontSize', 16);
box('off');
if(Flag==1)
    set(gca, 'XTickLabel', {'NC', 'MCI'});
elseif (Flag==2)
    set(gca, 'XTickLabel', {'NC', 'AD'});   
else
    set(gca, 'XTickLabel', {'MCI', 'AD'});  
end

% Significant inter-regional EC
figure;
imagesc(SMAP1);
set(gca, 'FontSize', 14);
title('Significant Map', 'Fontsize', 16)
xlabel('ROI')
ylabel('ROI')



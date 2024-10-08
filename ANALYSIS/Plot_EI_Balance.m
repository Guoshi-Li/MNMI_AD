
% Plot intra-regional, inter-regional, and total E-I balance for all 46
% ROIs
% The E-I balance is evaluated as E-I difference (i.e., incoming excitation
% - incoming inhibition)

clc;
close all;
clear all;

Flag = 1; % 1: NC vs. MCI; 2: NC vs AD; 3: MCI vs AD

q = 0.05;

if (Flag == 1)
   load PARA_NC_MCI;
   load GC_NC_MCI.mat
elseif (Flag == 2)
   load PARA_NC_AD;
   load GC_NC_AD.mat;
else
   load PARA_MCI_AD;
   load GC_MCI_AD.mat;
end

load AD_SC10;
load ADLabelNet2;

NR = 46;   % Number of ROIs

CN = NP;   % Number of inter-regional EC
NP = 2*NR+1+CN;  % Number of total estiamted parameters

NM = 4;    % Number of networks

n1 = 48;  % Number of subjects
n2 = 48;  % Number of subjects


MGC_NC = GC_NC;
MGC_MD = GC_MD;


%==========================================================================
%              Within-regional Connections
%==========================================================================

Wee_NC = W_NM(:,1:NR);
Wie_NC = W_NM(:,NR+1:2*NR);
Wgc_NC = W_NM(:,2*NR+1:NP-1);
SPI_NC = W_NM(:,NP);

Wee_MD = W_MD(:,1:NR);
Wie_MD = W_MD(:,NR+1:2*NR);
Wgc_MD = W_MD(:,2*NR+1:NP-1);
SPI_MD = W_MD(:,NP);

for k = 1:n1

   EI_NC(k, :) = Wee_NC(k,:)-Wie_NC(k,:);

end

for k = 1:n2

   EI_MD(k, :) = Wee_MD(k,:)-Wie_MD(k,:);

end


% Calculate mean and STD
Wi_NC = mean(EI_NC);
Wi_MD = mean(EI_MD);

SEi_NC = std(EI_NC)/sqrt(n1); 
SEi_MD = std(EI_MD)/sqrt(n2);

Wi_NC_MD  = [Wi_NC' Wi_MD'];
SEi_NC_MD = [SEi_NC' SEi_MD'];


% T-test

TPC = [];
TPCR = [];
TPCL = [];

TPC_STAT = [];
TPCR_STAT = [];
TPCL_STAT = [];


% For each ROI
for k =1:NR

   CV = EI_NC(:, k);
   DV = EI_MD(:, k);

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

IC = find(TPC<0.05);

% Beferonic Correction
BFT = 0.05/NR;
BI = find(TPC<BFT);

% Corrected Bonferroni value
BON1 = TPC*NR;
% Corrected FDR values
FDR1 = mafdr(TPC,'BHFDR', 1);

TV1 = TPCL_STAT(IC);
BP1 = BON1(IC);
FP1 = FDR1(IC);

disp('Results of Intra-regional E-I balance:')
fprintf('\n')

disp('The ROI index with E/I imbalance are:');
IC'

% disp('The ROI index with E/I imbalance (84 ROIs):');
% R4(IC)

disp('The T-values of significant connections:');
TV1'

disp('The p values less than 0.05 are:');
TPC(IC)'

disp('The Bonferroni-corrected p values:')
BP1'

disp('The FDR-corrected p values:')
FP1'


[pID0,pN0] = gretna_FDR(TPC, q)



%==========================================================================
%              Inter-regional Connections
%==========================================================================
% NC
WGC_NC = zeros(n1, NR);

WGC_POS_NC = zeros(n1, NR);
WGC_NEG_NC = zeros(n1, NR);

WGC_IN_NC = zeros(n1, NR);
WGC_OUT_NC = zeros(n1, NR);

WGC_POS_IN_NC = zeros(n1, NR);
WGC_POS_OUT_NC = zeros(n1, NR);

WGC_NEG_IN_NC = zeros(n1, NR);
WGC_NEG_OUT_NC = zeros(n1, NR);

% MCI/AD
WGC_MD = zeros(n2, NR);

WGC_POS_MD = zeros(n2, NR);
WGC_NEG_MD = zeros(n2, NR);

WGC_IN_MD = zeros(n2, NR);
WGC_OUT_MD = zeros(n2, NR);
  
WGC_POS_IN_MD = zeros(n2, NR);
WGC_POS_OUT_MD = zeros(n2, NR);

WGC_NEG_IN_MD = zeros(n2, NR);
WGC_NEG_OUT_MD = zeros(n2, NR);


% For NC/MCI

for k = 1:n2
    
    PN = [];
    
    IN = [];
    OUT = [];
    
    POS = [];
    NEG = [];
    
    POS_IN = [];
    NEG_IN = [];
    
    POS_OUT = [];
    NEG_OUT = [];
    
    GC = MGC_NC(:,:,k);
    
    for i = 1:NR
                  
      V1 = GC(i,:);
      V2 = GC(:,i);
    
      ind1 = find(V1>0);
      ind2 = find(V1<0);
    
      IN = [IN sum(V1)];
      OUT = [OUT sum(V2)];
      
      POS_IN = [POS_IN sum(V1(ind1))];
      NEG_IN = [NEG_IN sum(V1(ind2))];
    
      ind1 = find(V2>0);
      ind2 = find(V2<0);
    
      POS_OUT = [POS_OUT sum(V2(ind1))];
      NEG_OUT = [NEG_OUT sum(V2(ind2))];
        

      
      V = [V1 V2'];   
      PN = [PN sum((V))];

           
    end
    
    WGC_NC(k,:) = PN;
    
    WGC_IN_NC(k,:)  = IN;
    WGC_OUT_NC(k,:) = OUT;
    
    WGC_POS_NC(k,:) = POS_IN + POS_OUT;
    WGC_NEG_NC(k,:) = NEG_IN + NEG_OUT;    
    
    WGC_POS_IN_NC(k,:) = POS_IN;
    WGC_NEG_IN_NC(k,:) = NEG_IN;
    
    WGC_POS_OUT_NC(k,:) = POS_OUT;
    WGC_NEG_OUT_NC(k,:) = NEG_OUT;
    
end



% For MCI or AD

for k = 1:n2
    
    PN = [];
    
    IN  = [ ];
    OUT = [ ];
    
    POS = [];
    NEG = [];
    
    POS_IN = [];
    NEG_IN = [];
    
    POS_OUT = [];
    NEG_OUT = [];
    
    GC = MGC_MD(:,:,k);
    
    for i = 1:NR
                  
      V1 = GC(i,:);
      V2 = GC(:,i);
         
      ind1 = find(V1>0);
      ind2 = find(V1<0);
      
      IN = [IN sum(V1)];
      OUT = [OUT sum(V2)];
    
      POS_IN = [POS_IN sum(V1(ind1))];
      NEG_IN = [NEG_IN sum(V1(ind2))];
    
      ind1 = find(V2>0);
      ind2 = find(V2<0);
    
      POS_OUT = [POS_OUT sum(V2(ind1))];
      NEG_OUT = [NEG_OUT sum(V2(ind2))];
         
        
      V = [V1 V2'];   
      PN = [PN sum((V))];

           
    end
    
    WGC_MD(k,:) = PN;
    
    WGC_IN_MD(k,:) = IN;
    WGC_OUT_MD(k,:) = OUT;
    
    WGC_POS_MD(k,:) = POS_IN + POS_OUT;
    WGC_NEG_MD(k,:) = NEG_IN + NEG_OUT;      
    
    WGC_POS_IN_MD(k,:) = POS_IN;
    WGC_NEG_IN_MD(k,:) = NEG_IN;
    
    WGC_POS_OUT_MD(k,:) = POS_OUT;
    WGC_NEG_OUT_MD(k,:) = NEG_OUT;
    
end
    

% Calculate mean and STD
W1_NC = mean(WGC_NC);
W2_NC = mean(WGC_POS_NC);
W3_NC = mean(WGC_NEG_NC);
W4_NC = mean(WGC_IN_NC);
W5_NC = mean(WGC_OUT_NC);

W1_MD = mean(WGC_MD);
W2_MD = mean(WGC_POS_MD);
W3_MD = mean(WGC_NEG_MD);
W4_MD = mean(WGC_IN_MD);
W5_MD = mean(WGC_OUT_MD);

SE1_NC = std(WGC_NC)/sqrt(n1);
SE2_NC = std(WGC_POS_NC)/sqrt(n1);
SE3_NC = std(WGC_NEG_NC)/sqrt(n1);
SE4_NC = std(WGC_IN_NC)/sqrt(n1);
SE5_NC = std(WGC_OUT_NC)/sqrt(n1);

SE1_MD = std(WGC_MD)/sqrt(n2);
SE2_MD = std(WGC_POS_MD)/sqrt(n2);
SE3_MD = std(WGC_NEG_MD)/sqrt(n2);
SE4_MD = std(WGC_IN_MD)/sqrt(n2);
SE5_MD = std(WGC_OUT_MD)/sqrt(n2);

% Put NC/MCI and MCI/AD together
W1_NC_MD = [W1_NC' W1_MD'];
W2_NC_MD = [W2_NC' W2_MD'];
W3_NC_MD = [W3_NC' W3_MD'];
W4_NC_MD = [W4_NC' W4_MD'];
W5_NC_MD = [W5_NC' W5_MD'];

SE1_NC_MD = [SE1_NC' SE1_MD'];
SE2_NC_MD = [SE2_NC' SE2_MD'];
SE3_NC_MD = [SE3_NC' SE3_MD'];
SE4_NC_MD = [SE4_NC' SE4_MD'];
SE5_NC_MD = [SE5_NC' SE5_MD'];



%% ========================================================================
% Compute E/I Ratio

% Intra-regional E-I ratio
W_RATIO_INTRA_NC = Wee_NC./Wie_NC;
W_RATIO_INTRA_MD = Wee_MD./Wie_MD;

Wm_RATIO_INTRA_NC = mean(W_RATIO_INTRA_NC);
Wm_RATIO_INTRA_MD = mean(W_RATIO_INTRA_MD);

SE_RATIO_INTRA_NC = std(W_RATIO_INTRA_NC)/sqrt(n1);
SE_RATIO_INTRA_MD = std(W_RATIO_INTRA_MD)/sqrt(n2);

Wm_RATIO_INTRA = [Wm_RATIO_INTRA_NC' Wm_RATIO_INTRA_MD'];
SE_RATIO_INTRA = [SE_RATIO_INTRA_NC' SE_RATIO_INTRA_MD'];

% ===============================================================

% Compute the ratio of total excitation to total inhibition
W_TOTAL_EXC_NC = Wee_NC + WGC_POS_IN_NC;
W_TOTAL_INH_NC = Wie_NC - WGC_NEG_IN_NC;
W_RATIO_NC = W_TOTAL_EXC_NC./W_TOTAL_INH_NC;

W_TOTAL_EXC_MD = Wee_MD + WGC_POS_IN_MD;
W_TOTAL_INH_MD = Wie_MD - WGC_NEG_IN_MD;
W_RATIO_MD = W_TOTAL_EXC_MD./W_TOTAL_INH_MD;

%==============================================================
% TOTAL COMBINED EXCITATION/INHIBITION 
W_TOTAL_NC = EI_NC + WGC_IN_NC;
W_TOTAL_MD = EI_MD + WGC_IN_MD;

Wm_TOTAL_NC = mean(W_TOTAL_NC);
Wm_TOTAL_MD = mean(W_TOTAL_MD);

SE_TOTAL_NC = std(W_TOTAL_NC)/sqrt(n1);
SE_TOTAL_MD = std(W_TOTAL_MD)/sqrt(n2);

Wm_TOTAL_NC_MD = [Wm_TOTAL_NC'  Wm_TOTAL_MD'];
SE_TOTAL_NC_MD = [SE_TOTAL_NC'  SE_TOTAL_MD'];



%===============================================================


Wm_RATIO_NC = mean(W_RATIO_NC);
Wm_RATIO_MD = mean(W_RATIO_MD);

SE_RATIO_NC = std(W_RATIO_NC)/sqrt(n1);
SE_RATIO_MD = std(W_RATIO_MD)/sqrt(n2);

Wm_RATIO = [Wm_RATIO_NC' Wm_RATIO_MD'];
SE_RATIO = [SE_RATIO_NC' SE_RATIO_MD'];


%% ========================================================================
%                  Statistical Tests
% =========================================================================
 
TP1 = [];
TP2 = [];
TP3 = [];
TP4 = [];
TP5 = [];
TP8 = [];
TP9 = [];
TP10 = [];
TP11 = [];

TP4_STAT = [];
TP8_STAT = [];

for k =1:NR

CV1 = WGC_NC(:, k);
DV1 = WGC_MD(:, k);   
    
CV2 = WGC_POS_NC(:, k);
DV2 = WGC_POS_MD(:, k); 

CV3 = WGC_NEG_NC(:, k);
DV3 = WGC_NEG_MD(:, k); 

CV4 = WGC_IN_NC(:, k);
DV4 = WGC_IN_MD(:, k);

CV5 = WGC_OUT_NC(:, k);
DV5 = WGC_OUT_MD(:, k);

CV8 = W_TOTAL_NC(:, k);
DV8 = W_TOTAL_MD(:, k);

CV9 = W_RATIO_NC(:, k);
DV9 = W_RATIO_MD(:, k);

CV10 = W_RATIO_INTRA_NC(:, k);
DV10 = W_RATIO_INTRA_MD(:, k);


[h1,p1,ci,stats1] = ttest2(CV1, DV1);
[h2,p2,ci,stats2] = ttest2(CV2, DV2);
[h3,p3,ci,stats3] = ttest2(CV3, DV3);
[h4,p4,ci,stats4] = ttest2(CV4, DV4);
[h5,p5,ci,stats5] = ttest2(CV5, DV5);

[h8,p8,ci,stats8] = ttest2(CV8, DV8);
[h9,p9,ci,stats9] = ttest2(CV9, DV9);

[h10,p10,ci,stats10] = ttest2(CV10, DV10);


TP1 = [TP1; p1];
TP2 = [TP2; p2];
TP3 = [TP3; p3];
TP4 = [TP4; p4];
TP5 = [TP5; p5];

TP8 = [TP8; p8];
TP9 = [TP9; p9];

TP10 = [TP10; p10];


TP4_STAT = [TP4_STAT stats4.tstat];
TP8_STAT = [TP8_STAT stats8.tstat];


end

% ROI indices in 46 ROIs
I1 = find(TP1<0.05);
I2 = find(TP2<0.05);
I3 = find(TP3<0.05);
I4 = find(TP4<0.05);
I5 = find(TP5<0.05);

I8 = find(TP8<0.05);
I9 = find(TP9<0.05);
I10 = find(TP10<0.05);

% ROI indices in 84 ROIs (DK Atlas)
RI1 = R4(I1);  
RI2 = R4(I2); 
RI3 = R4(I3); 
RI4 = R4(I4); 
RI5 = R4(I5); 
RI8 = R4(I8); 
RI9 = R4(I9); 
RI10 = R4(I10);

% p values
PV1 = TP1(I1);
PV2 = TP2(I2);
PV3 = TP3(I3);
PV4 = TP4(I4);
PV5 = TP5(I5);
PV8 = TP8(I8);
PV9 = TP9(I9);
PV10 = TP10(I10);

TV4 = TP4_STAT(I4);
TV8 = TP8_STAT(I8);

% Transpose
I1 = I1';
I2 = I2';
I3 = I3';
I4 = I4';
I5 = I5';
I8 = I8';
I9 = I9';
I10 = I10';

PV1 = PV1';
PV2 = PV2';
PV3 = PV3';
PV4 = PV4';
PV5 = PV5';
PV8 = PV8';
PV9 = PV9';
PV10 = PV10';


%% ========================================================================
FDR4 = mafdr(TP4,'BHFDR', 1);
FP4  = FDR4(I4);

BON4 = TP4*NR;
BP4  = BON4(I4);

disp('Results of Inter-regional E-I balance:')
fprintf('\n')

disp('The ROI index with inter-regional E/I imbalance are:');
I4

disp('The T-values of significant connections:');
TV4

disp('The p values less than 0.05 are:');
PV4

disp('The Bonferroni-corrected p values:')
BP4'

disp('The FDR-corrected p values:')
FP4'


[pID4,pN4] = gretna_FDR(TP4, q)




%% ========================================================================
% Total E-I balance (E-I difference)

FDR8 = mafdr(TP8,'BHFDR', 1);
FP8  = FDR8(I8);

BON8 = TP8*NR;
BP8  = BON8(I8);

disp('Results of total E-I balance:')
fprintf('\n')

disp('The ROI index with overal E/I imbalance are:');
I8

disp('The T-values of significant connections:');
TV8

disp('The p values less than 0.05 are:');
PV8

disp('The Bonferroni-corrected p values:')
BP8'

disp('The FDR-corrected p values:')
FP8'


[pID8, pN8] = gretna_FDR(TP8, q)





% disp('Total combined E/I Ratio:');
% I9
% RI9
% PV9
% [pID9,pN9] = gretna_FDR(TP9, q)


% disp('Intra-regional E/I Ratio:');
% I10
% RI10
% % RI4
% PV10
% [pID10,pN10] = gretna_FDR(TP10, q)




%% ========================================================================
%                      Ploting
%  ========================================================================

Xn = 1:NR;

% Intra-regional E/I balance
figure;
bar(Xn, Wi_NC_MD);
hold on;
ngroups = size(Wi_NC_MD, 1);
nbars = size(Wi_NC_MD, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, Wi_NC_MD(:,i), SEi_NC_MD(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Intra-regional E-I Balance', 'FontSize', 18);
set(gca, 'XTickLabel', []);
ylabel('E-I Diff.', 'FontSize', 16);
% axis([0.5 4.5 2.0 3.5]);
% set(gca, 'YTick', [2:0.5:3.5]); 
box('off');

if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end



% Inter-regional E-I balance
figure;
bar(Xn, W4_NC_MD);
hold on;
ngroups = size(W4_NC_MD, 1);
nbars = size(W4_NC_MD, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, W4_NC_MD(:,i), SE4_NC_MD(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Inter-regional E-I Balance', 'FontSize', 18);
set(gca, 'XTickLabel', []);
ylabel('E-I Diff.', 'FontSize', 16);
box('off');
if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end



% Total E-I Difference
figure;
bar(Xn, Wm_TOTAL_NC_MD);
hold on;
ngroups = size(Wm_TOTAL_NC_MD, 1);
nbars = size(Wm_TOTAL_NC_MD, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, Wm_TOTAL_NC_MD(:,i), SE_TOTAL_NC_MD(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Total E-I Balance', 'FontSize', 18);
set(gca, 'XTickLabel', []);
ylabel('E-I Diff.', 'FontSize', 16);
% axis([0.5 4.5 2.0 3.5]);
% set(gca, 'YTick', [2:0.5:3.5]); 
box('off');

if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end



% % Total E/I Ratio
% figure;
% bar(Xn, Wm_RATIO);
% hold on;
% ngroups = size(Wm_RATIO, 1);
% nbars = size(Wm_RATIO, 2);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     er = errorbar(x, Wm_RATIO(:,i), SE_RATIO(:,i), '.');
%     er.Color = [0 0 0];
%     er.LineWidth = 1;
%     er.LineStyle = 'none'; 
% end
% hold off
% set(gca, 'FontSize', 14);
% title('Total E-I Balance', 'FontSize', 16);
% set(gca, 'XTickLabel', []);
% ylabel('E/I Ratio', 'FontSize', 16);
% % axis([0.5 4.5 2.0 3.5]);
% % set(gca, 'YTick', [2:0.5:3.5]); 
% box('off');
% xtickangle(45);
% 
% if (Flag==1)
%  legend('NC', 'MCI');
% elseif (Flag==2)
%   legend('NC', 'AD');
% else
%   legend('MCI', 'AD');  
% end
% 
% xtickangle(45);
% 
% 









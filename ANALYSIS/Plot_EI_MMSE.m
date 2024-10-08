
% Find out the regions that show significant (corrected) correlation between E-I
% balance and MMSE score
% Correlations are evaluated for intra-regional E-I balance, inter-regional
% E-I balance and total E-I balance

clc;
close all;
clear all;

load MMSE_SCORE.mat

q = 0.05;

load PARA_NC_MCI;
W_MCI = W_MD;

load PARA_NC_AD;
W_AD = W_MD;

load GC_NC_MCI.mat;
GC_MCI = GC_MD;

load GC_NC_AD.mat
GC_AD = GC_MD;


load AD_SC10;
load ADLabelNet2;

NR = 46;     % Number of ROIs

CN = NP;     % Number of inter-regional EC  
NP = 2*NR+1+CN;  % Number of total estimated parameters

NM = 4;   % Number of networks

n1 = 48;  % Number of subjects
n2 = 48;
n3 = 48;

EC_CON = W_NM;
EC_MCI = W_MCI;
EC_AD  = W_AD;

MGC_NC  = GC_NC;
MGC_MCI = GC_MCI;
MGC_AD  = GC_AD;

%==========================================================================
%             Intra-regional E-I balance
%==========================================================================

Wee_NC = W_NM(:,1:NR);
Wie_NC = W_NM(:,NR+1:2*NR);
Wgc_NC = W_NM(:,2*NR+1:NP-1);
SPI_NC = W_NM(:,NP);

Wee_MCI = W_MCI(:,1:NR);
Wie_MCI = W_MCI(:,NR+1:2*NR);
Wgc_MCI = W_MCI(:,2*NR+1:NP-1);
SPI_MCI = W_MCI(:,NP);

Wee_AD = W_AD(:,1:NR);
Wie_AD = W_AD(:,NR+1:2*NR);
Wgc_AD = W_AD(:,2*NR+1:NP-1);
SPI_AD = W_AD(:,NP);


for k = 1:n1

   EI_NC(k, :) = Wee_NC(k,:)-Wie_NC(k,:);

end

for k = 1:n2

   EI_MCI(k, :) = Wee_MCI(k,:)-Wie_MCI(k,:);

end

for k = 1:n3

   EI_AD(k, :) = Wee_AD(k,:)-Wie_AD(k,:);

end


% Correlation with MMSE

EI_INTRA_NET_ALL = [EI_NC; EI_MCI; EI_AD];
MMSE_ALL  = [MMSE_NC; MMSE_MCI; MMSE_AD];

RV1 = [];
PV1 = [];

for i = 1:NR

    U = EI_INTRA_NET_ALL(:,i);
        
    [r,p] = corr(U, MMSE_ALL);
     
    RV1 = [RV1; r];
    PV1 = [PV1; p];

end


disp('Correlation between Intra-regional E-I difference and MMSE:')

[pID1, pN1] = gretna_FDR(PV1, q)


k1 = find(PV1<=pID1);
RS1 = RV1(k1);
PS1 = PV1(k1);

FDR1 = mafdr(PV1, 'BHFDR', 1);
FP1 = FDR1(k1);

BON1 = PV1*NR;
BP1 = BON1(k1);

disp('The regions that pass FDR correction')
k1'

disp('Corresponding R values:')
RS1'

disp('Uncorrected p-values:')
PS1'

disp('FDR-corrected p-values:')
FP1'

disp('Bonferroni-corrected p-values:')
BP1'


%% ========================================================================
% Compute E/I Ratio

% Intra-regional E-I ratio
W_RATIO_INTRA_NC   = Wee_NC./Wie_NC;
 W_RATIO_INTRA_MCI = Wee_MCI./Wie_MCI;
W_RATIO_INTRA_AD   = Wee_AD./Wie_AD;

EI_INTRA_RATIO_ALL = [W_RATIO_INTRA_NC; W_RATIO_INTRA_MCI; W_RATIO_INTRA_AD];

RV2 = [];
PV2 = [];

for i = 1:NR

    U = EI_INTRA_RATIO_ALL(:,i);
    [r,p] = corr(U, MMSE_ALL);
     
    RV2 = [RV2; r];
    PV2 = [PV2; p];

end


disp('Correlation between Intra-regional E/I ratio and MMSE:')

[pID2, pN2] = gretna_FDR(PV2, q)

if (~isempty(pID2))

   k2  = find(PV2<=pID2)
   RS2 = RV2(k2)
   PS2 = PV2(k2)

   FDR2 = mafdr(PV2, 'BHFDR', 1);
   FP2 = FDR2(k2);

   BON2 = PV2*NR;
   BP2  = BON2(k2);

   disp('The regions that pass FDR correction')
   k2'

   disp('Corresponding R values:')
   RS2'

   disp('Uncorrected p-values:')
   PS2'

   disp('FDR-corrected p-values:')
   FP2'

   disp('Bonferroni-corrected p-values:')
   BP2'

end



%% ========================================================================
%              Inter-regional E-I balance
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


% MCI
WGC_MCI = zeros(n2, NR);

WGC_POS_MCI = zeros(n2, NR);
WGC_NEG_MCI = zeros(n2, NR);

WGC_IN_MCI = zeros(n2, NR);
WGC_OUT_MCI = zeros(n2, NR);

WGC_POS_IN_MCI = zeros(n2, NR);
WGC_POS_OUT_MCI = zeros(n2, NR);

WGC_NEG_IN_MCI  = zeros(n2, NR);
WGC_NEG_OUT_MCI = zeros(n2, NR);


% AD
WGC_AD = zeros(n3, NR);

WGC_POS_AD = zeros(n3, NR);
WGC_NEG_AD = zeros(n3, NR);

WGC_IN_AD = zeros(n3, NR);
WGC_OUT_AD = zeros(n3, NR);
  
WGC_POS_IN_AD = zeros(n3, NR);
WGC_POS_OUT_AD = zeros(n3, NR);

WGC_NEG_IN_AD = zeros(n3, NR);
WGC_NEG_OUT_AD = zeros(n3, NR);


% For NC

for k = 1:n1
    
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


% For MCI

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
    
    GC = MGC_MCI(:,:,k);
    
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
    
    WGC_MCI(k,:) = PN;
    
    WGC_IN_MCI(k,:)  = IN;
    WGC_OUT_MCI(k,:) = OUT;
    
    WGC_POS_MCI(k,:) = POS_IN + POS_OUT;
    WGC_NEG_MCI(k,:) = NEG_IN + NEG_OUT;    
    
    WGC_POS_IN_MCI(k,:) = POS_IN;
    WGC_NEG_IN_MCI(k,:) = NEG_IN;
    
    WGC_POS_OUT_MCI(k,:) = POS_OUT;
    WGC_NEG_OUT_MCI(k,:) = NEG_OUT;
    
end


% For AD

for k = 1:n3
    
    PN = [];
    
    IN  = [ ];
    OUT = [ ];
    
    POS = [];
    NEG = [];
    
    POS_IN = [];
    NEG_IN = [];
    
    POS_OUT = [];
    NEG_OUT = [];
    
    GC = MGC_AD(:,:,k);
    
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
    
    WGC_AD(k,:) = PN;
    
    WGC_IN_AD(k,:) = IN;
    WGC_OUT_AD(k,:) = OUT;
    
    WGC_POS_AD(k,:) = POS_IN + POS_OUT;
    WGC_NEG_AD(k,:) = NEG_IN + NEG_OUT;      
    
    WGC_POS_IN_AD(k,:) = POS_IN;
    WGC_NEG_IN_AD(k,:) = NEG_IN;
    
    WGC_POS_OUT_AD(k,:) = POS_OUT;
    WGC_NEG_OUT_AD(k,:) = NEG_OUT;
    
end
    

WGC_IN_ALL = [WGC_IN_NC; WGC_IN_MCI; WGC_IN_AD];

RV3 = [];
PV3 = [];

for i = 1:NR

    U = WGC_IN_ALL (:,i);
%     V = W_RATIO_INTRA_AD(:,i);
%     UV = [U; V];

    [r,p] = corr(U, MMSE_ALL);
     
    RV3 = [RV3; r];
    PV3 = [PV3; p];

end

fprintf('\n')
disp('Correlation between Inter-regional E-I difference and MMSE:')

[pID3, pN3] = gretna_FDR(PV3, q)

% k3 = find(PV3<=pID3);
% PS3 = PV3(k3)
% RS3 = RV3(k3)



%==========================================================================
%              Combined E-I balance
%==========================================================================

% TOTAL COMBINED EXCITATION/INHIBITION 
W_TOTAL_NC  = EI_NC  + WGC_IN_NC;
W_TOTAL_MCI = EI_MCI + WGC_IN_MCI;
W_TOTAL_AD  = EI_AD  + WGC_IN_AD;


% Compute the ratio of total excitation to total inhibition
W_TOTAL_EXC_NC  = Wee_NC  + WGC_POS_IN_NC;
W_TOTAL_EXC_MCI = Wee_MCI + WGC_POS_IN_MCI;
W_TOTAL_EXC_AD  = Wee_AD  + WGC_POS_IN_AD;

W_TOTAL_INH_NC  = Wie_NC  - WGC_NEG_IN_NC;
W_TOTAL_INH_MCI = Wie_MCI - WGC_NEG_IN_MCI;
W_TOTAL_INH_AD  = Wie_AD  - WGC_NEG_IN_AD;


W_RATIO_NC  = W_TOTAL_EXC_NC./W_TOTAL_INH_NC;
W_RATIO_MCI = W_TOTAL_EXC_MCI./W_TOTAL_INH_MCI;
W_RATIO_AD  = W_TOTAL_EXC_AD./W_TOTAL_INH_AD;

W_TOTAL_ALL = [W_TOTAL_NC; W_TOTAL_MCI; W_TOTAL_AD];
W_RATIO_ALL = [W_RATIO_NC; W_RATIO_MCI; W_RATIO_AD];


%% ======================================================================== 
% Correlation using total E-I difference

RV5 = [];
PV5 = [];

RV6 = [];
PV6 = [];

for i = 1:NR

    U = W_TOTAL_ALL (:,i);
%     V = W_RATIO_INTRA_AD(:,i);
%     UV = [U; V];

    [r,p] = corr(U, MMSE_ALL);
     
    RV5 = [RV5; r];
    PV5 = [PV5; p];

end

fprintf('\n')
disp('Correlation between Total E-I difference and MMSE:')

[pID5, pN5] = gretna_FDR(PV5, q)

k5 = find(PV5<=pID5);
PS5 = PV5(k5);
RS5 = RV5(k5);

FDR5 = mafdr(PV5, 'BHFDR', 1);
FP5 = FDR5(k5);

BON5 = PV5*NR;
BP5 = BON5(k5);

disp('The regions that pass FDR are:')
k5'

disp('The corresponding R-values are:')
RS5'

disp('The uncorrected p-values are:')
PS5'

disp('The FDR-corrected p-values are:')
FP5'

disp('The Bonferroni-corrected p-values are:')
BP5'


EI21_T = W_TOTAL_ALL(:,21);
EI37_T = W_TOTAL_ALL(:,37);
EI38_T = W_TOTAL_ALL(:,38);
EI41_T = W_TOTAL_ALL(:,41);




%% ======================================================================== 
% Correlation using total E/I Ratio

for i = 1:NR

    U = W_RATIO_ALL (:,i);
%     V = W_RATIO_INTRA_AD(:,i);
%     UV = [U; V];

    [r,p] = corr(U, MMSE_ALL);
     
    RV6 = [RV6; r];
    PV6 = [PV6; p];

end

fprintf('\n')
disp('Correlation between Total E/I ratio and MMSE:')

[pID6, pN6] = gretna_FDR(PV6, q)

k6 = find(PV6<=pID6);
PS6 = PV6(k6);
RS6 = RV6(k6);


FDR6 = mafdr(PV6, 'BHFDR', 1);
FP6 = FDR6(k6);

BON6 = PV6*NR;
BP6 = BON6(k6);


disp('The regions that pass FDR are:')
k6'

disp('The corresponding R-values are:')
RS6'

disp('The uncorrected p-values are:')
PS6'

disp('The FDR-corrected p-values are:')
FP6'

disp('The Bonferroni-corrected p-values are:')
BP6'



TOTAL_RATIO_EI41 = W_RATIO_ALL(:,41);


%% ======================================================================
% Correlation between total E-I differnece and MMSE

figure;
scatter(EI21_T, MMSE_ALL,'LineWidth',1);
axis([-8, 8, 18, 31]);
set(gca, 'XTick', [-8:4:8])
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
set(gca, 'FontSize', 16);
title('L.IN', 'FontSize', 18);
xlabel('Total E-I Diff.');
ylabel('MMSE Score');
txt1 = [ 'r = 0.25' ];
txt2 = ['p* = 0.04' ];
text(-7.5, 30, txt1, 'FontSize', 18, 'Color', 'r');
text(-7.5, 29, txt2, 'FontSize', 18, 'Color', 'r');

figure;
scatter(EI37_T, MMSE_ALL,'LineWidth',1);
axis([-10, 10, 18, 31]);
set(gca, 'XTick', [-10:5:10])
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
set(gca, 'FontSize', 16);
title('L.PUT', 'FontSize', 18);
xlabel('Total E-I Diff.');
ylabel('MMSE Score');
txt1 = [ 'r = 0.24' ];
txt2 = ['p* = 0.045' ];
text(6, 27, txt1, 'FontSize', 18, 'Color', 'r');
text(6, 26, txt2, 'FontSize', 18, 'Color', 'r');


figure;
scatter(EI38_T, MMSE_ALL,'LineWidth',1);
axis([-8, 8, 18, 31]);
set(gca, 'XTick', [-8:4:8])
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
set(gca, 'FontSize', 16);
title('R.PUT', 'FontSize', 18);
xlabel('Total E-I Diff.');
ylabel('MMSE Score');
txt1 = [ 'r = 0.25' ];
txt2 = ['p* = 0.04' ];
text(4.8, 27, txt1, 'FontSize', 18, 'Color', 'r');
text(4.8, 26, txt2, 'FontSize', 18, 'Color', 'r');


figure;
scatter(EI41_T, MMSE_ALL,'LineWidth',1);
axis([-6, 6, 18, 31]);
set(gca, 'XTick', [-6:3:6])
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
set(gca, 'FontSize', 16);
title('L.HPC', 'FontSize', 18);
xlabel('Total E-I Diff.');
ylabel('MMSE Score');
txt1 = [ 'r = 0.31' ];
txt2 = ['p** = 0.007' ];
text(3.3, 27, txt1, 'FontSize', 18, 'Color', 'r');
text(3.3, 26, txt2, 'FontSize', 18, 'Color', 'r');




% Correlation between total E/I ratio and MMSE

figure;
scatter(TOTAL_RATIO_EI41, MMSE_ALL,'LineWidth',1);
axis([0, 3, 18, 31]);
h1 = lsline();
h1.LineWidth = 2;
h1.Color = 'r';
set(gca, 'FontSize', 16);
title('L.HPC', 'FontSize', 18);
xlabel('Total E/I Ratio');
ylabel('MMSE Score');
txt1 = [ 'r = 0.29' ];
txt2 = ['p** = 0.017' ];
text(2.2, 27, txt1, 'FontSize', 18, 'Color', 'r');
text(2.2, 26, txt2, 'FontSize', 18, 'Color', 'r');



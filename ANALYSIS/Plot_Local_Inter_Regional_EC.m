
% Plot local recurrent excitation/inhibition and inter-regional EC
% Flag = 1: plot NC and MCI; Flag = 2: plot NC and AD; Flag = 3: Plot MCI and AD 

clc;
close all;
clear all;

Flag = 1; % 1: NC vs. MCI; 2: NC vs AD; 3: MCI vs AD

if (Flag == 1)
   load PARA_NC_MCI;
elseif (Flag == 2)
   load PARA_NC_AD;
else
   load PARA_MCI_AD;
end

load AD_SC10;
load ADLabelNet2;


NR = 46;   % Number of ROIs
CN = NP;   % Number of inter-regional EC

NP = 2*NR+1+CN;  % Number of total parameters
NM = 4;          % Number of networks
Nb = N4;

n1 = 48;  % number of subjects
n2 = 48;  % number of subjects


%==========================================================================
%              Individual Connection Analysis
%==========================================================================

% Average weight
Wm_NC = mean(W_NM);
Wm_MD = mean(W_MD);

STD_NC = std(W_NM);
STD_MD = std(W_MD);

SE_NC = STD_NC/sqrt(n1);
SE_MD = STD_MD/sqrt(n2);

% Put NC and MDD together
Wm_NC_MD  = [Wm_NC;   Wm_MD];
STD_NC_MD = [STD_NC; STD_MD];
SE_NC_MD  = [SE_NC; SE_MD];



% Divide the parameters into four sets
Wm_EE = Wm_NC_MD(:, 1:NR);
Wm_IE = Wm_NC_MD(:, NR+1:2*NR);
Wm_GC = Wm_NC_MD(:, 2*NR+1:NP-1);
Wm_PI = Wm_NC_MD(:, NP);

Wm_GC1 = Wm_GC(:, 1:CN/2);
Wm_GC2 = Wm_GC(:, CN/2+1:end);

% Significnat connections
SGC_NC_MD = Wm_GC'; 
SGC_NC_MD = SGC_NC_MD(IR,:);


% Standard Error
EE_ER = SE_NC_MD(:, 1:NR);
IE_ER = SE_NC_MD(:, NR+1:2*NR);
GC_ER = SE_NC_MD(:, 2*NR+1:NP-1);

GC_ER_NC = GC_ER(1,:);
GC_ER_MD = GC_ER(2,:);

% Significant connections
SSE_NC_MD = GC_ER';
SSE_NC_MD = SSE_NC_MD(IR,:);


% Re-arrange the ROIs according to the fucntional network

Wm_EE_SORT = [];
Wm_IE_SORT = [];

EE_ER_SORT = [];
IE_ER_SORT = [];

for L = 1:NM
    
    idx = NID{L};
    Wm_EE_SORT = [Wm_EE_SORT Wm_EE(:, idx)];
    Wm_IE_SORT = [Wm_IE_SORT Wm_IE(:, idx)];
    
    EE_ER_SORT = [EE_ER_SORT EE_ER(:, idx)];
    IE_ER_SORT = [IE_ER_SORT IE_ER(:, idx)];

end
       

Wm_EE_SORT = Wm_EE_SORT';
Wm_IE_SORT = Wm_IE_SORT';

EE_ER_SORT = EE_ER_SORT';
IE_ER_SORT = IE_ER_SORT';



% Reconstruct the inter-regional NR*NR matrix

m = 1;

GC1 = Wm_GC(1,:);
GC2 = Wm_GC(2,:);

MGC_NC = zeros(NR, NR);
MGC_MD = zeros(NR, NR);

MGC_NC_UNSORTED = zeros(NR, NR);
MGC_MD_UNSORTED = zeros(NR, NR);

TEMP1_GC_NC = [];
TEMP2_GC_NC = [];

TEMP1_GC_MD = [];
TEMP2_GC_MD = [];

%  STE

MSE_NC = zeros(NR, NR);
MSE_MD = zeros(NR, NR);

MSE_NC_UNSORTED = zeros(NR, NR);
MSE_MD_UNSORTED = zeros(NR, NR);

TEMP1_SE_NC = [];
TEMP2_SE_NC = [];

TEMP1_SE_MD = [];
TEMP2_SE_MD = [];


for j=1:NR
    for i=1:NR
       
       if (MAP(i,j)==1)
           MGC_NC_UNSORTED(i,j) = GC1(m);
           MGC_MD_UNSORTED(i,j) = GC2(m);
           
           MSE_NC_UNSORTED(i,j) = GC_ER_NC(m);
           MSE_MD_UNSORTED(i,j) = GC_ER_MD(m);
                      
           m = m+1;
       end      
        
    end  
end


% Weights

for k = 1:NM
    
   idx = NID{k};
   TEMP1_GC_NC = [TEMP1_GC_NC; MGC_NC_UNSORTED(idx,:)];
   TEMP1_GC_MD = [TEMP1_GC_MD; MGC_MD_UNSORTED(idx,:)];
   
   % STE
   TEMP1_SE_NC = [TEMP1_SE_NC; MSE_NC_UNSORTED(idx,:)];
   TEMP1_SE_MD = [TEMP1_SE_MD; MSE_MD_UNSORTED(idx,:)];
     
end

for k = 1:NM
    
   idx = NID{k};
   TEMP2_GC_NC = [TEMP2_GC_NC TEMP1_GC_NC(:, idx)];
   TEMP2_GC_MD = [TEMP2_GC_MD TEMP1_GC_MD(:, idx)];
  
   % STE
   TEMP2_SE_NC = [TEMP2_SE_NC TEMP1_SE_NC(:, idx)];
   TEMP2_SE_MD = [TEMP2_SE_MD TEMP1_SE_MD(:, idx)];
     
end

MGC0_NC = TEMP2_GC_NC;
MGC0_MD = TEMP2_GC_MD;

MSE0_NC = TEMP2_SE_NC;
MSE0_MD = TEMP2_SE_MD;


% =========================================================================
% Exchange cACC (15-16) with PCC (17-18)

MGC_NC = [MGC0_NC(1:14, :); MGC0_NC(17:18, :); MGC0_NC(15:16, :); MGC0_NC(19:end,:)];
MGC_NC = [MGC_NC(:,1:14) MGC_NC(:,17:18) MGC_NC(:,15:16) MGC_NC(:, 19:end)];

MGC_MD = [MGC0_MD(1:14, :); MGC0_MD(17:18, :); MGC0_MD(15:16, :); MGC0_MD(19:end,:)];
MGC_MD = [MGC_MD(:,1:14) MGC_MD(:,17:18) MGC_MD(:,15:16) MGC_MD(:, 19:end)];


MSE_NC = [MSE0_NC(1:14, :); MSE0_NC(17:18, :); MSE0_NC(15:16, :); MSE0_NC(19:end,:)];
MSE_NC = [MSE_NC(:,1:14) MSE_NC(:,17:18) MSE_NC(:,15:16) MSE_NC(:, 19:end)];

MSE_MD = [MSE0_MD(1:14, :); MSE0_MD(17:18, :); MSE0_MD(15:16, :); MSE0_MD(19:end,:)];
MSE_MD = [MSE_MD(:,1:14) MSE_MD(:,17:18) MSE_MD(:,15:16) MSE_MD(:, 19:end)];

% =========================================================================


VGC_NC = MGC_NC(:);
VGC_MD = MGC_MD(:);


% Remove 0 connections
VGC_NC = VGC_NC(VGC_NC~=0);
VGC_MD = VGC_MD(VGC_MD~=0);

VGC_NC_MD = [VGC_NC VGC_MD];

VGC1 = VGC_NC_MD(1:CN/2, :);
VGC2 = VGC_NC_MD(CN/2+1:end, :);

% Standard errors
VSE_NC = MSE_NC(:);
VSE_MD = MSE_MD(:);

VSE_NC = VSE_NC(VSE_NC~=0);
VSE_MD = VSE_MD(VSE_MD~=0);

VSE_NC_MD = [VSE_NC VSE_MD];

VSE1 = VSE_NC_MD(1:CN/2, :);
VSE2 = VSE_NC_MD(CN/2+1:end, :);

%==========================================================================
%                        Network Analysis
%==========================================================================

WNee = zeros(2, NM);
WNie = zeros(2, NM);

% Intra-regional weights
for i = 1:NM
    
    WTe = [];
    WTi = [];
    
    idx = NID{i};
    WTe = Wm_EE(:,idx);
    WTi = Wm_IE(:,idx);
    
    WNee(:,i) = mean(WTe, 2);
    WNie(:,i) = mean(WTi, 2);
   
end


Wee_NC = W_NM(:,1:NR);
Wie_NC = W_NM(:,NR+1:2*NR);
Wgc_NC = W_NM(:,2*NR+1:NP-1);
SPI_NC = W_NM(:,NP);

Wee_MD = W_MD(:,1:NR);
Wie_MD = W_MD(:,NR+1:2*NR);
Wgc_MD = W_MD(:,2*NR+1:NP-1);
SPI_MD = W_MD(:,NP);


WeeN_NC = zeros(n1, NM);
WieN_NC = zeros(n1, NM);

WeeN_MD = zeros(n2, NM);
WieN_MD = zeros(n2, NM);

for i = 1:NM
    
    idx = NID{i};
    
    % For NC
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_NC(:,idx);
    WTi = Wie_NC(:,idx);
    
    WeeN_NC(:,i) = mean(WTe, 2);
    WieN_NC(:,i) = mean(WTi, 2);
    
   % For MDD 
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_MD(:,idx);
    WTi = Wie_MD(:,idx);
    
    WeeN_MD(:,i) = mean(WTe, 2);
    WieN_MD(:,i) = mean(WTi, 2);
    
        
end

WN_NC = [WeeN_NC WieN_NC SPI_NC];
WN_MD = [WeeN_MD WieN_MD SPI_MD];

WeeN_NC_Mean = mean(WeeN_NC);
WeeN_MD_Mean = mean(WeeN_MD);

WieN_NC_Mean = mean(WieN_NC);
WieN_MD_Mean = mean(WieN_MD);

WeeN = [ WeeN_NC_Mean; WeeN_MD_Mean];
WieN = [ WieN_NC_Mean; WieN_MD_Mean];


% Standard errors
STD_NET_NC = std(WN_NC);
STD_NET_MD = std(WN_MD);

SE_NET_NC = STD_NET_NC/sqrt(n1);
SE_NET_MD = STD_NET_MD/sqrt(n2);

SE_NET_NC_MD = [SE_NET_NC; SE_NET_MD];

SE_EE_NET = SE_NET_NC_MD(:,1:NM);
SE_IE_NET = SE_NET_NC_MD(:,NM+1:2*NM);
SE_PI = SE_NET_NC_MD(:,end);


%==========================================================================
%                   Inter-regional weights
%==========================================================================

GCN_NC = zeros(NM, NM);
GCN_EXC_NC = zeros(NM, NM);
GCN_INH_NC = zeros(NM, NM);

GCN_MD = zeros(NM, NM);
GCN_EXC_MD = zeros(NM, NM);
GCN_INH_MD = zeros(NM, NM);

k1 = 1;
l1 = 1;

for i = 1:NM
    
    l1 = 1;
    
    for j = 1:NM
   
      k2 = k1+Nb(i)-1;  
      l2 = l1+Nb(j)-1;
      
      BLOCK = MGC_NC(k1:k2, l1:l2);
      [a, b] = find(BLOCK~=0);
      
      if isempty(a)
          GCN_NC(i, j) = 0;
      else
          GCN_NC(i, j) = sum(sum(BLOCK))/length(a);
      end
            
      PBLOCK = BLOCK(:);
      
      ID1 = find(PBLOCK>0);
      ID2 = find(PBLOCK<0);
      
      if ( length(ID1)~=0 ) 
         TEMP1 = PBLOCK(ID1);
         GCN_EXC_NC(i, j) = sum(TEMP1)/length(ID1);
      end
      
      if ( length(ID2)~=0 )  
         TEMP2 = PBLOCK(ID2);
         GCN_INH_NC(i, j) = sum(TEMP2)/length(ID2);
      end      
        
      l1 = l1+Nb(j);
      
    end
    
    k1 = k1+Nb(i);
end


k1 = 1;
l1 = 1;

for i = 1:NM
    
    k2 = k1+Nb(i)-1;
    l1 = 1;
    
    for j = 1:NM
   
      l2 = l1+Nb(j)-1;
      
      BLOCK = MGC_MD(k1:k2, l1:l2);
      [a, b] = find(BLOCK~=0);
      
      if isempty(a)
          GCN_MD(i, j) = 0;
      else
          GCN_MD(i, j) = sum(sum(BLOCK))/length(a);
      end
            
      PBLOCK = BLOCK(:);
      
      ID1 = find(PBLOCK>0);
      ID2 = find(PBLOCK<0);
      
      if ( length(ID1)~=0 ) 
         TEMP1 = PBLOCK(ID1);
         GCN_EXC_MD(i, j) = sum(TEMP1)/length(ID1);
      end
      
      if ( length(ID2)~=0 )  
         TEMP2 = PBLOCK(ID2);
         GCN_INH_MD(i, j) = sum(TEMP2)/length(ID2);
      end      
        
      l1 = l1+Nb(j);
      
    end
    
    k1 = k1+Nb(i);

end



%==========================================================================

FVm_NC = mean(FVAL_NM);
FVm_MD = mean(FVAL_MD);
FVm    =  [FVm_NC FVm_MD];

SE_FV_NC = std(FVAL_NM)/sqrt(n1);
SE_FV_MD = std(FVAL_MD)/sqrt(n2);
SE_FV = [SE_FV_NC SE_FV_MD];


%==========================================================================

% if (Flag == 1)
%    save MGC_NC_MCI.mat MGC_NC MGC_MD SGC_NC_MD  SSE_NC_MD;
% elseif (Flag == 2)
%    save MGC_NC_AD.mat  MGC_NC MGC_MD SGC_NC_MD  SSE_NC_MD;
% else
%    save MGC_MCI_AD.mat MGC_NC MGC_MD SGC_NC_MD  SSE_NC_MD;
% end


%==========================================================================
%                         Ploting
%==========================================================================

Wm_EE  = Wm_EE';
Wm_IE  = Wm_IE';
Wm_GC  = Wm_GC';
Wm_GC1 = Wm_GC1';
Wm_GC2 = Wm_GC2';

Wm_PI = Wm_PI';

WNee = WNee';
WNie = WNie';

WeeN = WeeN';
WieN = WieN';

FLIP_VGC_NC_MD = flip(VGC_NC_MD);


% STD
EE_ER = EE_ER';
IE_ER = IE_ER';
GC_ER = GC_ER';

GC1_ER = GC_ER(1:CN/2,:);
GC2_ER = GC_ER(CN/2+1:end,:);

SE_EE_NET = SE_EE_NET';
SE_IE_NET = SE_IE_NET';
SE_PI = SE_PI';

FLIP_VSE_NC_MD = flip(VSE_NC_MD);


%% ========================================================================
q = 0.05;

PS1 = TPC1(RE1);
PS2 = TPC2(RI1);

PS1 = PS1';
PS2 = PS2';

% Bonferroni correction
BON1 = TPC1*NR;
BON2 = TPC2*NR;
BON3 = TPC3*length(TPC3);

BP1 = BON1(RE1);
BP2 = BON2(RI1);

% FDR correction
FDR1 = mafdr(TPC1, 'BHFDR',1);
FDR2 = mafdr(TPC2, 'BHFDR',1);
FDR3 = mafdr(TPC3, 'BHFDR',1);

CP1 = FDR1(RE1);
CP2 = FDR2(RI1);



disp('The ROIs with significant recurrent excitation are:');
RE1
% RE2

disp('The significant p values for intra-regional recurrent excitation:')
PS1

disp('The Bonferroni-corrected p-values are:')
BP1'

disp('The FDR-corrected p-values are:')
CP1'

disp('FDR Correction for intra-regional excitatory connections:');
[pID1,pN1] = gretna_FDR(TPC1,q)

fprintf('\n')
disp('The ROIs with significant recurrent inhibition are:');
RI1
% RI2

disp('The significant p values for intra-regional recurrent inhibition:')
PS2

disp('The Bonferroni-corrected p-values are:')
BP2'

disp('The FDR-corrected p-values are:')
CP2'

disp('FDR Correction for intra-regional inhibitory connections:');
[pID2,pN2] = gretna_FDR(TPC2,q)


fprintf('\n')

disp('The significant (uncorrected) inter-regional connections are:');
IR

disp('The corresponding uncorrected p-values are:');
TPC3(IR)'

disp('FDR Correction for inter-regional connections:');
[pID3,pN3] = gretna_FDR(TPC3,q)



%% ========================================================================
% Plot local recurrent excitation and inhibition

ROI_LABEL_NET = {'DMN', 'SAL', 'EXE', 'LIM' };          
                        
Xtick = 1:NR;        

figure;
bar(Xtick, Wm_EE_SORT);
hold on;
ngroups = size(Wm_EE_SORT, 1);
nbars = size(Wm_EE_SORT, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, Wm_EE_SORT(:,i), EE_ER_SORT(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Recurrent Excitation', 'FontSize', 18);
ylabel('Weight', 'FontSize', 16);
xtickangle(45);
axis([0 NR+1 2.5 3.6]);
set(gca, 'YTick', [2.5:0.5:4]); 
set(gca, 'XTickLabel', [ ]); 
box('off');

if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end



figure;
bar(Xtick, Wm_IE_SORT);
hold on;
ngroups = size(Wm_IE_SORT, 1);
nbars = size(Wm_IE_SORT, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, Wm_IE_SORT(:,i), IE_ER_SORT(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Recurrent Inhibition', 'FontSize', 18);
ylabel('Weight', 'FontSize', 16);
xtickangle(45);
axis([0 NR+1 2.5 3.6]);
set(gca, 'YTick', [2:0.5:4]); 
set(gca, 'XTickLabel', [ ]); 
box('off');
if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end



%==========================================================================
% Plot network-averaged local recurrent excitation and inhibition

Xn = 1:NM;

figure;
bar(Xn, WeeN);
hold on;
ngroups = size(WeeN, 1);
nbars = size(WeeN, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, WeeN(:,i), SE_EE_NET(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Recurrent Excitation', 'FontSize', 18);
set(gca, 'XTick', Xn, 'XTickLabel', ROI_LABEL_NET);
ylabel('Weight', 'FontSize', 16);
axis([0.5 4.5 2.0 3.5]);
set(gca, 'YTick', [2:0.5:3.5]); 
box('off');
if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end
xtickangle(45);


figure;
bar(Xn, WieN);
hold on;
ngroups = size(WieN, 1);
nbars = size(WieN, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, WieN(:,i), SE_IE_NET(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Recurrent Inhibition', 'FontSize', 18);
set(gca, 'XTick', Xn, 'XTickLabel', ROI_LABEL_NET);
ylabel('Weight', 'FontSize', 16);
axis([0.5 4.5 2.0 3.5]);
set(gca, 'YTick', [2:0.5:3.5]); 
box('off');
if (Flag==1)
 legend('NC', 'MCI');
elseif (Flag==2)
  legend('NC', 'AD');
else
  legend('MCI', 'AD');  
end
xtickangle(45);



%==========================================================================
% Plot inter-regional EC

CLIM = [-1.0 1.0];

colI = 50;   % 12; color gradations

cmp = [interp1([1 colI],[0.0 0.0 1; 1 1 1],1:colI); interp1([1 colI],[1 1 1; 1 0.0 0.0],2:colI)];

figure;
imagesc(MGC_NC, CLIM);
set(gca, 'Fontsize', 14);
colormap(cmp);  % jet
colorbar;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
if (Flag==1||Flag==2)
     title('Inter-regional EC (NC)', 'FontSize', 18);
else
     title('Inter-regional EC (MCI)', 'FontSize', 18);
end


figure;
imagesc(MGC_MD, CLIM);
set(gca, 'Fontsize', 14);
colormap(cmp); % jet
colorbar;
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
if (Flag==1)
     title('Inter-regional EC (MCI)', 'FontSize', 18);
else
     title('Inter-regional EC (AD)', 'FontSize', 18);
end







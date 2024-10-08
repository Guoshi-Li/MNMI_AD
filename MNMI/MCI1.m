
% MNMI to estiamte the connection parameters (or effective connectivity, EC)
% of a neural mass model (Wilson-Cowan type) based on rs-fMRI

% This code estimates the parameters for MCI 1 (MCI1)
% For MCI subject i, just change "MCI1" to "MCIi" for the "Subject" variable

% This code is supposed to be run in UNC Longleaf cluster
% For local computer or other clusters, the code for parallel computing
% needs to be changed accordingly

% Written by Guosh Li (guoshi_li@med.unc.edu)

tic

clc;
close all;
clear all;

Subject = 'MCI1';
filename = [Subject '.mat'];

load AD_SC10;  % load structural connectivity

%% ========================================================================
% Configuration for parallel computing  

POOLSIZE = 24;

% For UNC Longleaf cluster
filepath = ['/pine/scr/g/u/guoshili/matlab/' Subject];

if ~exist(filepath, 'dir')
    mkdir (filepath);
end

pc = parcluster('local');
pc.JobStorageLocation = filepath;
pc.NumWorkers=24;

pool = parpool(pc, POOLSIZE);


%% ========================================================================

MAX_GEN = 128;     % Maximal generations
FUN_TOL = 1e-3;    % Functional tolerance

rng(66,'twister');   % To ensure results are reproducible

NR = 46;  % Number of regions

% Lower and Upper bound of estimated parameters
MIN_Wee  = 2;   % Recurrent excitation
MAX_Wee  = 4;

MIN_Wie  = 2;   % Recurent inhibition
MAX_Wie  = 4;

MIN_gc =  -2;   % Inter-regional EC 
MAX_gc =   2;    

MIN_Spi = 0.2;  % Spontaenous input
MAX_Spi = 0.4;

MIN_WEE  = MIN_Wee*ones(1, NR);
MAX_WEE  = MAX_Wee*ones(1, NR);

MIN_WIE  = MIN_Wie*ones(1, NR);
MAX_WIE  = MAX_Wie*ones(1, NR);

MIN_GC = MIN_gc*ones(1,NP);  % NP is the nubmer of inter-regional EC
MAX_GC = MAX_gc*ones(1,NP);

MIN_SPI = MIN_Spi;
MAX_SPI = MAX_Spi;

lb = [MIN_WEE MIN_WIE MIN_GC MIN_SPI]; 
ub = [MAX_WEE MAX_WIE MAX_GC MAX_SPI]; 

np = length(lb);  % total number of parameters to be estimated

%% ========================================================================
% Call the Matlab function "ga" to estiamte model parameters

opts = optimoptions('ga', 'MaxGenerations', MAX_GEN, 'FunctionTolerance', FUN_TOL, ...
    'UseParallel',true, 'Display','iter');

% opts = optimoptions('ga', 'MaxGenerations', MAX_GEN, 'FunctionTolerance', FUN_TOL, ...
%     'UseParallel',true, 'Display','iter','PlotFcn', @gaplotbestf);

[x,fval,eflag,outpt, population, scores] = ga(@FCobj, np,...
    [],[],[],[],lb,ub,[],opts);


disp('The optimal solution is:');
x
fval


disp('The maximal correlation is:');
f = -FCobj(x)

% Save the output
save (filename, 'x', 'fval', 'eflag', 'outpt', 'population', 'scores');


% Delete the current parallel pool
poolobj = gcp('nocreate');
delete(poolobj);


disp('Program completed!');

toc

exit;


%% ========================================================================
% Define the objective function to be minimzied by GA

% This code estimates the parameters for MCI 1 (MCI1)
% For MCI subject i, just change "MCI1" to "MCIi" for "Subject" 

function f = FCobj(x)

NR = 46;

rng(66,'twister');

Subject = 'MCI1';

%  Load functional connectivity
ss = ['../DATA/FC/FC_' Subject];     % FC is located at "DATA/FC" directory
load (ss)

EFC  = FC;
VEFC = EFC(:);

% Load Structural connectivity
load AD_SC10;       
SCM = SCN;

Flag_Noise = 1;      % 1: Introduce noise to the neural model
Flag_Mean_BOLD = 0;  % 1: Remove the mean of the neural activity to the Hemodynamic model

DT = 10e-3;   % Integration step (10 ms) 
ST = 200;     % Total simulation time in sec (200 sec)

TR  = 3;      % Scan interval (3 sec)
NTR = TR/DT;  % Number of DT during scan interval

TFC = 20;     % Start time to calculate simulated FC (remove first 20 sec)
NFC = TFC/DT;

W_EI0 = 3.0;  % Fixed E-->I connectoin strength
Wei = W_EI0*ones(NR, 1);

% Assign estimated parameters (x) to the model
Wee = x(1:NR);   
Wie = x(NR+1:2*NR);
Wgc = x(2*NR+1:2*NR+NP);
SPI = x(end);

% Initialize the inter-regional EC matrix
GC = zeros(NR,NR); 

% Construct the inter-regional EC matrix based on estimated parameters and
% the SC mask

m = 1;

for j=1:NR
    for i=1:NR
       
       if (MAP(i,j)==1)
           GC(i,j) = Wgc(m);
           m = m+1;
       end      
        
    end  
end

% Calculate neural activity using the current estimated parameters
[t, X] = Model_NEURAL(NR, DT, ST, SCM, GC, Wee, Wei, Wie, SPI, Flag_Noise);

% Compute the simulated BOLD 
[BOLD] = Model_HEMO_HRF(NR, DT, ST, X, Flag_Mean_BOLD);

% Remove the first 20 sec activity and adjust the resolution
SBOLD = BOLD(NFC:NTR:end, :);

% Compute simulated FC
SFC   = corrcoef(SBOLD);

% Change to a vector
VSFC = SFC(:);

% Compute the correlation between emprical and simulated FC
CO = corrcoef(EFC, SFC);
COF = CO(1,2);

% Minmizie the opposite of the correlation coefficient
f = -COF;



end




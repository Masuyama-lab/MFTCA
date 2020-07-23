% 
% (c) 2020 Naoki Masuyama
% 
% These are the codes of Multilayer Fast Topological CIM-based ART (MFTCA)
% proposed in "N. Amako, N. Masuyama, C. K. Loo, Y. Nojima, Y. Liu, and H. Ishibuchi,
% Multilayer Clustering Based on Adaptive Resonance Theory for Noisy Environments,
% Proc. of 2020 International Joint Conference on Neural Networks (IJCNN 2020), 
% pp. 1-8, Glasgow, UK, July 19-24, 2020."
% 
% Please contact "masuyama@cs.osakafu-u.ac.jp" if you have any problems.
%   

clc;
clear;
close all;
rng('default');

NR = 0.1; % Noise Rate

V = 0.2; % similarity threshold 
Lambda = 50; % topology construction cycle
LEVEL = 5; % Number of Layers

% load Data
load 2D_ClusteringDATASET; c=6; nData=0000; numD=15000;

DATA = [data(1:end,1) data(1:end,2)];

d1 = DATA(1:numD,:);
d2 = DATA(15001:15000+numD,:);
d3 = DATA(30001:30000+numD,:);
d4 = DATA(45001:45000+numD,:);
d5 = DATA(60001:60000+numD,:);
d6 = DATA(75001:75000+numD,:);
D = [d1; d2; d3; d4; d5; d6];
DATA = D;

% Label setting
l1 = 1*ones(size(d1,1),1);
l2 = 2*ones(size(d2,1),1);
l3 = 3*ones(size(d3,1),1);
l4 = 4*ones(size(d4,1),1);
l5 = 5*ones(size(d5,1),1);
l6 = 6*ones(size(d6,1),1);
L = [l1; l2; l3; l4; l5; l6];
LABEL = L;


normDATA = DATA;
for k=1:size(normDATA,2) % Normalization [0-1]
    mmin = min(normDATA(:,k));
    mmax = max(normDATA(:,k));
    if mmin == mmax
        normDATA(:,k) = 1;
    else
        normDATA(:,k) = (normDATA(:,k)-mmin) ./ (mmax-mmin);
    end
end

DATA = normDATA;

% Parameters of MFTCA =================================================
for l=1:LEVEL
    MFTCAnet.Lay{l}.minCIM = V;   % similarity threshold 
    MFTCAnet.Lay{l}.Lambda = Lambda;    % topology construction cycle
    MFTCAnet.Lay{l}.LEVEL = LEVEL;      % Number of Layers
    
    MFTCAnet.Lay{l}.numNodes = 0;       % Number of node
    MFTCAnet.Lay{l}.weight = [];        % Node position
    MFTCAnet.Lay{l}.CountNode = [];     % Counter for each node
    MFTCAnet.Lay{l}.adaptiveSig = [];   % sigma in each node
    MFTCAnet.Lay{l}.edge = [];          % Edge matrix
    MFTCAnet.Lay{l}.NewEdgedNode = [];  % Node which creates new edge
    MFTCAnet.Lay{l}.CountLabel = [];    % Counter for each cluster
    
    MFTCAnet.Lay{l}.InputSample = [];   % Input samples to layer l
    MFTCAnet.Lay{l}.CountSample = 0;    % Number of input samples to layer l
end
% ====================================================================

time_mftca = 0;

% Noise Setting
noiseDATA=[];
noiseLABEL=[];
if NR > 0
    noiseDATA = rand((numD*NR)*c, size(DATA,2));
    noiseLABEL = randi(c,size(noiseDATA,1),1);
end

% Randamize data
ran = randperm(size(DATA,1));
data = DATA(ran,:);
label = LABEL(ran,:);

% Add noise data
data(1:size(noiseDATA,1),:) = noiseDATA;
label(1:size(noiseLABEL,1),:) = noiseLABEL;

% Randamize data
ran = randperm(size(data,1));
data = data(ran,:);
label = label(ran,:);    

% MFTCA  =====================================
tic
MFTCAnet = MFTCA_Lay1(data, label, c, MFTCAnet);
time_mftca = time_mftca + toc;
myPlotMFTCA(MFTCAnet, LEVEL);    % Plot results

disp([' Total Process Time:  MFTCA: ', num2str(time_mftca)]);  




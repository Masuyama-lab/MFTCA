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
function net = MFTCA_LayN(input,label,net,maxLABEL,level)

numNodes = net.Lay{level}.numNodes;       % Number of nodes
weight = net.Lay{level}.weight;           % Node position
CountNode = net.Lay{level}.CountNode;     % Counter for each node
adaptiveSig = net.Lay{level}.adaptiveSig; % sigma in each node
LEVEL  = net.Lay{level}.LEVEL;

% Parameters for Topology
edge = net.Lay{level}.edge;
NewEdgedNode = net.Lay{level}.NewEdgedNode;
Lambda = net.Lay{level}.Lambda;

minCIM = net.Lay{level}.minCIM;

CountLabel = net.Lay{level}.CountLabel;

InputSample = net.Lay{level}.InputSample;
CountSample = net.Lay{level}.CountSample;
CountSample = CountSample + 1;
InputSample(CountSample, :) = input;

% Set a size of CountLabel
if size(weight) == 0
    CountLabel = zeros(1, maxLABEL);
end

estSig = mean(net.Lay{level - 1}.adaptiveSig);

    
if size(weight,1) < 2 % Initialization Process
    % Add Node
    numNodes = numNodes + 1;
    weight(numNodes,:) = input;
    CountNode(numNodes) = 1;
    NewEdgedNode(1, numNodes) = 0;
    edge(numNodes, :) = 0;
    edge(:, numNodes) = 0;
    adaptiveSig(numNodes) = estSig;

    CountLabel(numNodes,label) = 1;

else
        
    % Calculate CIM between data sample and nodes 
    globalCIM = CIM(input, weight, mean(adaptiveSig));

    % Determine the 1st and 2nd winner nodes
    [Lcim_s1, s1] = min(globalCIM);
    globalCIM(s1) = inf;
    [Lcim_s2, s2] = min(globalCIM);

    if Lcim_s1 > minCIM % Case �T
        % Add Node
        numNodes = numNodes + 1;
        weight(numNodes,:) = input;
        CountNode(numNodes) = 1;
        NewEdgedNode(1, numNodes) = 0;
        edge(numNodes,:) = 0;
        edge(:,numNodes) = 0;
        adaptiveSig(numNodes) = mean(adaptiveSig);

        CountLabel(numNodes,label) = 1;

    else
        % Case II
        % Update the 1st winner node
        weight(s1,:) = weight(s1,:) + (1/CountNode(s1)) * (input - weight(s1,:));
        CountNode(s1) = CountNode(s1) + 1;
        CountLabel(s1,label) = CountLabel(s1, label) + 1;
        % go to next level =======================
        if level < LEVEL
            % Save imformation of current network
            net.Lay{level}.adaptiveSig = adaptiveSig;
            net.Lay{level}.numNodes = numNodes;      
            net.Lay{level}.weight = weight;          
            net.Lay{level}.CountNode = CountNode;    
            net.Lay{level}.adaptiveSig = adaptiveSig;

            net.Lay{level}.edge = edge;
            net.Lay{level}.NewEdgedNode = NewEdgedNode;
            net.Lay{level}.Lambda = Lambda;

            net.Lay{level}.CountLabel = CountLabel;

            net.Lay{level}.InputSample = InputSample;
            net.Lay{level}.CountSample = CountSample;

            net = MFTCA_LayN(input,label,net,maxLABEL, level + 1); % Function of next layer
        end

        if Lcim_s2 <= minCIM % Case III
            % Update neighbor nodes
            s1Neighbors = find(edge(s1,:));
            for k = s1Neighbors
                weight(k,:) = weight(k,:) + (1/(10*CountNode(k))) * (weight(s1,:) - weight(k,:));
            end

            % Create an edge between the 1st and 2nd winner nodes
            NewEdgedNode(1,s1) = 1;
            edge(s1,s2) = 1;
            edge(s2,s1) = 1;

        end
    end
end

 % Topology Reconstruction
if mod(CountSample, Lambda) == 0
    % -----------------------------------------------------------------
    % Delete Node based on number of neighbors
    nNeighbor = sum(edge);
    deleteNodeEdge = (nNeighbor == 0);

    % Delete process
    numNodes = numNodes - sum(deleteNodeEdge);
    weight(deleteNodeEdge, :) = [];
    CountNode(deleteNodeEdge) = [];
    NewEdgedNode(:, deleteNodeEdge) = [];
    edge(deleteNodeEdge, :) = [];
    edge(:, deleteNodeEdge) = [];
    adaptiveSig(deleteNodeEdge) = [];
    CountLabel(deleteNodeEdge, :) = [];

    % -----------------------------------------------------------------
    % Delete Intersections of edge
    [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, mean(adaptiveSig));
        
end % mod(CountSample, Lambda) == 0
             

% Cluster Labeling based on edge (Functions are available above R2015b.)
connection = graph(edge ~= 0);
LebelCluster = conncomp(connection);

net.Lay{level}.numNodes = numNodes;    
net.Lay{level}.weight = weight;         
net.Lay{level}.CountNode = CountNode;    
net.Lay{level}.adaptiveSig = adaptiveSig;
if isempty(LebelCluster) == 1
    LebelCluster = 0;
end
net.Lay{level}.LebelCluster = LebelCluster;
net.Lay{level}.edge = edge;
net.Lay{level}.NewEdgedNode = NewEdgedNode;
net.Lay{level}.Lambda = Lambda;

net.Lay{level}.CountLabel = CountLabel;

net.Lay{level}.InputSample = InputSample;
net.Lay{level}.CountSample = CountSample;

end


% Correntropy induced Metric (Gaussian Kernel based)
function cim = CIM(X,Y,sig)
% X : 1 x n
% Y : m x n
[n, att] = size(Y);
g_Kernel = zeros(n, att);

for i = 1:att
    g_Kernel(:,i) = GaussKernel(X(i)-Y(:,i), sig);
end

ret0 = GaussKernel(0, sig);
ret1 = mean(g_Kernel, 2);

cim = sqrt(ret0 - ret1)';
end

function g_kernel = GaussKernel(sub, sig)
g_kernel = 1/(sqrt(2*pi)*sig) * exp(-sub.^2/(2*sig^2));
end



% Delete intersections of edge
function [weight, edge, NewEdgedNode] = DeleteIntersection(weight, edge, NewEdgedNode, sigma)

% for d = 1:size(weight,1); % Search all nodes
for d = find(NewEdgedNode == 1) % Search only new edged nodes
    
    node1 = find(edge(d,:)); % Neighbors of d-th node
    if size(node1,1) >= 1
       posX1 = weight(d,:); % position of d-th node
        for m = 1:size(node1,2) % Search all neighbors of d-th nodes
            posY1 = weight(node1(m),:); % position of m-th neighbor node of d-th node
            for h = 1:size(node1,2)
                target2 = node1(h);
                node2 = find(edge(target2,:)); % Neighbors of m-th node
                posX2 = weight(target2,:); % position of h-th neighbor node of m-th node
                for k = 1:size(node2,2)
                    posY2 = weight(node2(k),:); % position of k-th neighbor node of h-th node
                    isConvex = findIntersection(posX1, posY1, posX2, posY2); % find intersections
                    if isConvex == 1 % If intersection is exist, delete edge which has larger CIM.
                        cim1 = CIM(weight(d,:), weight(node1(m),:), sigma);
                        cim2 = CIM(weight(target2,:), weight(node2(k),:), sigma);
                        if cim2 >= cim1
                            edge(target2, node2(k)) = 0;
                            edge(node2(k), target2) = 0;
                        else
                            edge(d, node1(m)) = 0;
                            edge(node1(m), d) = 0;
                        end
                    end % end isConvex
                end % end k
            end % end h
        end % end m  
    end

end % end d

NewEdgedNode = zeros(size(NewEdgedNode));

end

% Check intersection of edges
function [isConvex] = findIntersection(A, B, C, D)

F1  = B(:,1)-D(:,1);
F2  = B(:,2)-D(:,2);
M11 = B(:,1)-A(:,1);
M21 = B(:,2)-A(:,2);
M12 = C(:,1)-D(:,1);
M22 = C(:,2)-D(:,2);
deter = M11.*M22 - M12.*M21;
lambda = -(F2.*M12-F1.*M22)./deter;
gamma = (F2.*M11-F1.*M21)./deter;

% E = (lambda*[1 1]).*A + ((1-lambda)*[1 1]).*B;
% isConvex = (0 <= lambda & lambda <= 1)  & (0 <= gamma & gamma <= 1);

isConvex = (0 < lambda & lambda < 1)  & (0 < gamma & gamma < 1) ;
isConvex = isConvex';

end





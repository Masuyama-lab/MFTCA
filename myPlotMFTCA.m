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
function myPlotMFTCA(net, LEVEL)

for l = 1:LEVEL
    w = net.Lay{l}.weight;
    edge = net.Lay{l}.edge;

    [N,D] = size(w);
    
    label = net.Lay{l}.LebelCluster;
    
    figure(l);
   
    hold on;
    
    for i=1:N-1
        for j=i:N
            if edge(i,j)==1
                if D==2
                    plot([w(i,1) w(j,1)],[w(i,2) w(j,2)],'k','LineWidth',1.5);
                elseif D==3
                    plot3([w(i,1) w(j,1)],[w(i,2) w(j,2)],[w(i,3) w(j,3)],'k','LineWidth',1.5);
                end
            end
        end
    end
    
    % Change Node color based on LebelCluster
    color = [
        [1 0 0]; % r    
        [0 1 0]; % g
        [0 0 1]; % b   
        [1 0 1]; %pink   
        [0.8500 0.3250 0.0980];% brown  
        [0.9290 0.6940 0.1250];% yellow
        [0.4940 0.1840 0.5560];
        [0.4660 0.6740 0.1880];
        [0.3010 0.7450 0.9330];
        [0.6350 0.0780 0.1840];
        ];
    m = length(color);
    
    for k = 1:N
        if D==2
            plot(w(k,1),w(k,2),'.','Color',color(mod(label(1,k)-1,m)+1,:),'MarkerSize',35);
        elseif D==3
            plot3(w(k,1),w(k,2),w(k,3),'.','Color',color(mod(label(1,k)-1,m)+1,:),'MarkerSize',35);
        end
    end

    
    % back ground color
    set(gca, 'Color', 'w');
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', 'w');
        
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',30,'XTick',...
        [0.0 0.2 0.4 0.6 0.8 1.0],'YTick',[0.0 0.2 0.4 0.6 0.8 1.0]);
    
    axes1.XAxis.TickLabelFormat = '%,.1f';
    axes1.YAxis.TickLabelFormat = '%,.1f';
    
    
    axis equal
    grid on
    set(gca,'GridColor','k')
    set(gca,'layer','bottom');
    hold off
    axis([0 1 0 1]);

end
end
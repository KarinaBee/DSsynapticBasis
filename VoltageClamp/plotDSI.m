% load dTableB2vert.mat
% % % dTable = dTable(strcmp(dTable.GFPLabel, "DRD4"),:);
% % %change here between genotypes and conditions
% % ihexTable = dTable(strcmp(dTable.holdVol,"12") & strcmp(dTable.drugs,"hex") ,:);
% iTable = dTable(strcmp(dTable.holdVol,"12") & strcmp(dTable.drugs, "none"),:);
% % ehexTable = dTable(strcmp(dTable.holdVol,"-60") & strcmp(dTable.drugs,"hex") ,:);
% eTable = dTable(strcmp(dTable.holdVol,"-60") & strcmp(dTable.drugs, "none"),:);

%% IPSCs
%Plot DSI for IPSCs 
    %for ON and OFF
    data = [iTable.DSIon iTable.DSIoff];


   figure('Name','DSI plot')
    tiledlayout(1,2);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('IPSC DSI')

% nexttile
% %Plot DSI for IPSCs in hex
% %     for ON and OFF
%     data = [ihexTable.DSIon ihexTable.DSIoff];
% 
%     plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
%     'yLabel','Direction Selectivity Index')
%     ylim([0 1])
%     title('IPSC DSI Hex')
% 
%     hold on 
    
%% EPSCs
%Plot DSI for EPSCs
    %for ON and OFF
    nexttile
    data = [eTable.DSIon eTable.DSIoff];
    
    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([-.5 1])
    title('ESPC DSI')

%     nexttile
%     hold on 
% %Plot DSI for EPSCs
% %     %for ON and OFF
%     data = [ehexTable.DSIon ehexTable.DSIoff];
% 
%     plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
%     'yLabel','Direction Selectivity Index')
%     ylim([0 1])
%     title('ESPC DSI Hex')



load dTableWTvert.mat
% WTdTable = dTable(strcmp(dTable.stim,"bars"),:);
WTeTable = dTable(strcmp(dTable.holdVol, "-60") & strcmp(dTable.drugs, "none"),:);
WTiTable = dTable(strcmp(dTable.holdVol, "12") & strcmp(dTable.drugs, "none"),:);
% WTiTable = dTable; 

load dTableB2vert.mat
% B2dTable = dTable(strcmp(dTable.stim,"bars"),:);
B2eTable = dTable(strcmp(dTable.holdVol, "-60") & strcmp(dTable.drugs, "none"),:);
B2iTable = dTable(strcmp(dTable.holdVol, "12") & strcmp(dTable.drugs, "none"),:);
% B2iTable = dTable; 

%plotting normalized vector sum for ON and OFF / WT and KO 

 %for ON and OFF WT exc and inh
  data = [WTiTable.normVSon WTiTable.normVSoff];
   figure('Name','Normalized VS')
    tiledlayout(2,2);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Normalized VS')
    ylim([0 2])
    title('WT inh')


    nexttile

   data = [B2iTable.normVSon B2iTable.normVSoff];
   
    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Normalized VS')
    ylim([0 2])
    title('B2 inh')

    nexttile

   data = [WTeTable.normVSon WTeTable.normVSoff];
   
    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Normalized VS')
    ylim([0 2])
    title('WT exc')

    nexttile

   data = [B2eTable.normVSon B2eTable.normVSoff];
   
    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Normalized VS')
    ylim([0 2])
    title('B2 Exc')


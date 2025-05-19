load dTableWThorzPreHex.mat
dTableWThorzPreHex = dTable;

load dTableWThorzPostHex.mat
dTableWThorzPostHex = dTable;

% load dTableWTvertPreHex.mat
% dTableWTvertPreHex = dTable;
% 
% load dTableWTvertPostHex.mat
% dTableWTvertPostHex = dTable;

WTehexTable = dTableWThorzPostHex(strcmp(dTableWThorzPostHex.holdVol,"-60") & strcmp(dTableWThorzPostHex.drugs,"hex") ,:);
WTeTable = dTableWThorzPreHex(strcmp(dTableWThorzPreHex.holdVol,"-60") & strcmp(dTableWThorzPreHex.drugs, "none"),:);
WTihexTable = dTableWThorzPostHex(strcmp(dTableWThorzPostHex.holdVol,"12") & strcmp(dTableWThorzPostHex.drugs,"hex") ,:);
WTiTable = dTableWThorzPreHex(strcmp(dTableWThorzPreHex.holdVol,"12") & strcmp(dTableWThorzPreHex.drugs, "none"),:);
% 
% WTehexTable = dTableWTvertPostHex(strcmp(dTableWTvertPostHex.holdVol,"-60") & strcmp(dTableWTvertPostHex.drugs,"hex") ,:);
% WTeTable = dTableWTvertPreHex(strcmp(dTableWTvertPreHex.holdVol,"-60") & strcmp(dTableWTvertPreHex.drugs, "none"),:);
% WTihexTable = dTableWTvertPostHex(strcmp(dTableWTvertPostHex.holdVol,"12") & strcmp(dTableWTvertPostHex.drugs,"hex") ,:);
% WTiTable = dTableWTvertPreHex(strcmp(dTableWTvertPreHex.holdVol,"12") & strcmp(dTableWTvertPreHex.drugs, "none"),:);

%%  WT EPSCs
% %Plot peak EPSCs for pref and null of cell in HEX 
 %for ON 
    PostPref = abs(WTehexTable.pCurNullON);
    PostNull = abs(WTehexTable.pCurPrefON); 
    data = [PostPref PostNull];
    null =length(PostPref);
    figure('Name','Peak EPSC Currents WT')
    tiledlayout(2,2);
    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1:size(null)
            % Plotting data points
        plot([1 2], [PostPref PostNull], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
         hold on
            xlim([0 3])
            ylim([0 600])
            ylabel('Peak EPSCs (pA)')
            title('Peak ON EPSCs')
             % Adding x-axis tick labels
             set(gca, 'XTick', [1 2], 'XTickLabel', {'Pref', 'Null'})

    % Drawing lines between data points
    plot([1 2], [PostPref(k) PostNull(k)], 'r-');

    hold off
           
        end

    %for OFF
    PostPref = abs(WTehexTable.pCurNullOFF);
    PostNull = abs(WTehexTable.pCurPrefOFF); 
    data = [PostPref PostNull];
    null =length(PostPref);
  
    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1:size(null)
            % Plotting data points
        plot([1 2], [PostPref PostNull], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
         hold on
            xlim([0 3])
            ylim([0 600])
            ylabel('Peak EPSCs (pA)')
            title('Peak OFF EPSCs')
             % Adding x-axis tick labels
             set(gca, 'XTick', [1 2], 'XTickLabel', {'Pref', 'Null'})

    % Drawing lines between data points
    plot([1 2], [PostPref(k) PostNull(k)], 'r-');

    hold off
           
        end


     
%% WT IPSCs
%Plot peak IPSCs for pref and null of cell in HEX 
%for ON 
     PostPref = abs(WTihexTable.pCurNullON);
    PostNull = abs(WTihexTable.pCurPrefON);  
     data = [PostPref PostNull];
    null =length(PostPref);

   
    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(null)
             % Plotting data points
     plot([1 2], [PostPref PostNull], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on

            xlim([0 3])
            ylim([0 1500])
            ylabel('Peak IPSCs (pA)')
            title('Peak ON IPSCs')
             % Adding x-axis tick labels
   set(gca, 'XTick', [1 2], 'XTickLabel', {'Pref', 'Null'})

    % Drawing lines between data points
    plot([1 2], [PostPref(k) PostNull(k)], 'r-');

    hold off
        end

    %for OFF
     PostPref = abs(WTihexTable.pCurNullOFF);
    PostNull = abs(WTihexTable.pCurPrefOFF);  
     data = [PostPref PostNull];
   
    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1 :size(null)
             % Plotting data points
     plot([1 2], [PostPref PostNull], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on

            xlim([0 3])
            ylim([0 1500])
            ylabel('Peak IPSCs (pA)')
            title('Peak OFF IPSCs')
             % Adding x-axis tick labels
   set(gca, 'XTick', [1 2], 'XTickLabel', {'Pref', 'Null'})

    % Drawing lines between data points
    plot([1 2], [PostPref(k) PostNull(k)], 'r-');

    hold off
        end

  

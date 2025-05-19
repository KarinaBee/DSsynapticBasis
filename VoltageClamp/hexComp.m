load dTableWTvertPreHex.mat
dTableWThorzPreHex = dTable;

% dTableWThorzPreHex = dTable(dTable.expDate == 250407, :);

load dTableWTvertPostHex.mat
dTableWThorzPostHex = dTable;
% dTableWThorzPostHex =dTable(dTable.expDate == 250407, :);


WTehexTable = dTableWThorzPostHex(strcmp(dTableWThorzPostHex.holdVol,"-60") & strcmp(dTableWThorzPostHex.drugs,"hex") ,:);
WTeTable = dTableWThorzPreHex(strcmp(dTableWThorzPreHex.holdVol,"-60") & strcmp(dTableWThorzPreHex.drugs, "none"),:);
WTihexTable = dTableWThorzPostHex(strcmp(dTableWThorzPostHex.holdVol,"12") & strcmp(dTableWThorzPostHex.drugs,"hex") ,:);
WTiTable = dTableWThorzPreHex(strcmp(dTableWThorzPreHex.holdVol,"12") & strcmp(dTableWThorzPreHex.drugs, "none"),:);

numColors = size(WTeTable, 1); 
% Generate 'numColors' unique colors dynamically
colors = parula(numColors); % You can change to jet(numColors), turbo(numColors), etc.
%%  WT EPSCs
% %Plot peak EPSCs for pref in pre and post hex for PREF OF CELL 
 %for ON 
    PrePref = abs(WTeTable.pCurNullON);
    PostPref = abs(WTehexTable.pCurNullON); 
    data = [PrePref PostPref];
    null =length(PrePref);
    figure('Name','Peak EPSC Currents WT')
    tiledlayout(2,2);
    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
            % Plotting data points
        plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
         hold on
            xlim([0 3])
            ylim([0 600])
            ylabel('Peak EPSCs (pA)')
            title('Peak ON EPSCs Pref Dir')
             % Adding x-axis tick labels
             set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})
 
             ePrefChangeON = ((PostPref - PrePref) ./ PrePref) * 100;
             ePrefRatioON = (PostPref./PrePref);
             ePrefPreON = PrePref; 
             ePrefPostON = PostPref; 

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
           
        end

    %for OFF
    PrePref = abs(WTeTable.pCurNullOFF);
    PostPref = abs(WTehexTable.pCurNullOFF); 
    data = [PrePref PostPref];

    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
            % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on
            xlim([0 3])
            ylim([0 600])
            ylabel('Peak EPSCs (pA)')
            title('Peak OFF EPSCs Pref Dir')
                % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

    ePrefChangeOFF = ((PostPref - PrePref)./ PrePref) * 100;
    ePrefRatioOFF = (PostPref./PrePref);
    ePrefPreOFF = PrePref; 
    ePrefPostOFF = PostPref; 

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

    %for ON Null dir of cell 
    PrePref = abs(WTeTable.pCurPrefON);
    PostPref = abs(WTehexTable.pCurPrefON); 
    data = [PrePref PostPref];

    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
             % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on
            xlim([0 3])
            ylim([0 600])
            ylabel('Peak EPSCs (pA)')
            title('Peak ON EPSCs Null Dir')
           % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

     eNullChangeON = ((PostPref - PrePref) ./ PrePref) * 100;
     eNullRatioON = (PostPref./PrePref);
     eNullPreON = PrePref; 
     eNullPostON = PostPref; 

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

%   for OFF Null dir of cell 
    PrePref = abs(WTeTable.pCurPrefOFF);
    PostPref = abs(WTehexTable.pCurPrefOFF); 
    data = [PrePref PostPref];

    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
          % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on
            xlim([0 3])
            ylim([0 600])
            ylabel('Peak EPSCs (pA)')
            title('Peak OFF EPSCs Null Dir')
             % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

    eNullChangeOFF = ((PostPref - PrePref) ./ PrePref) * 100;
    eNullRatioOFF = (PostPref./PrePref);
    eNullPreOFF = PrePref; 
     eNullPostOFF = PostPref; 
    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

     
%% WT IPSCs
%Plot peak IPSCs for pref in pre and post hex for PREF OF CELL 
%for ON 
    PrePref = abs(WTiTable.pCurNullON);
    PostPref = abs(WTihexTable.pCurNullON); 
    data = [PrePref PostPref];

    iPrefChangeON = ((PostPref - PrePref) ./ PrePref) * 100;
    iPrefRatioON = (PostPref./PrePref);

    iPrefPreON = PrePref; 
     iPrefPostON = PostPref; 


    figure('Name','Peak IPSC Currents WT')
    tiledlayout(2,2);
    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
             % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on

            xlim([0 3])
            ylim([0 1500])
            ylabel('Peak IPSCs (pA)')
            title('Peak ON IPSCs Pref Dir')
             % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

    %for OFF
    PrePref = abs(WTiTable.pCurNullOFF);
    PostPref = abs(WTihexTable.pCurNullOFF); 
    data = [PrePref PostPref];

   iPrefChangeOFF = ((PostPref - PrePref) ./ PrePref) * 100;
   iPrefRatioOFF = (PostPref./PrePref);
iPrefPreOFF = PrePref; 
     iPrefPostOFF = PostPref;

    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
              % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on

            xlim([0 3])
            ylim([0 1500])
            ylabel('Peak IPSCs (pA)')
            title('Peak OFF IPSCs Pref Dir')
             % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

     %for ON Null dir of cell 
    PrePref = abs(WTiTable.pCurPrefON);
    PostPref = abs(WTihexTable.pCurPrefON); 
    data = [PrePref PostPref];

    iNullChangeON = ((PostPref - PrePref) ./ PrePref) * 100;
    iNullRatioON = (PostPref./PrePref);
    iNullPreON = PrePref; 
    iNullPostON = PostPref;

    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
              % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on

            xlim([0 3])
            ylim([0 1500])
            ylabel('Peak IPSCs (pA)')
            title('Peak ON IPSCs Null Dir')
             % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

    %for OFF Null dir of cell 
    PrePref = abs(WTiTable.pCurPrefOFF);
    PostPref = abs(WTihexTable.pCurPrefOFF); 
    data = [PrePref PostPref];

    iNullChangeOFF = ((PostPref - PrePref) ./ PrePref) * 100;
    iNullRatioOFF = (PostPref./PrePref);

    iNullPreOFF = PrePref; 
    iNullPostOFF = PostPref;

    nexttile 

   plotSpread(data, 'showMM', 5,'xNames', {'Pre' 'Post'})
        for k = 1:size(null)
              % Plotting data points
    plot([1 2], [PrePref PostPref], 'ro-', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    hold on

            xlim([0 3])
            ylim([0 1500])
            ylabel('Peak IPSCs (pA)')
            title('Peak OFF IPSCs Null Dir')
             % Adding x-axis tick labels
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre', 'Post'})

    % Drawing lines between data points
    plot([1 2], [PrePref(k) PostPref(k)], 'r-');

    hold off
        end

       % change EPSCS 
    figure('Name','Percent Change Hex')
    tiledlayout(2,1);
    nexttile 

      data = [ePrefChangeON eNullChangeON ePrefChangeOFF eNullChangeOFF]; 
      plotSpread(data, 'showMM', 5,'xNames', {'PD ON', 'ND ON', 'PD OFF', 'ND OFF'})
      ylim([-100 150])
      title('EPSC PERCENT CHANGE')
      nexttile 
      data = [iPrefChangeON iNullChangeON iPrefChangeOFF iNullChangeOFF]; 
      plotSpread(data, 'showMM', 5,'xNames', {'PD ON', 'ND ON', 'PD OFF', 'ND OFF'})
        ylim([-100 150])
        title('IPSC PERCENT CHANGE')
        %%Ratio 
figure('Name','Drug/Ctrl Ratio Hex')
    tiledlayout(2,1);
    nexttile 

      data = [ePrefRatioON eNullRatioON ePrefRatioOFF eNullRatioOFF]; 
      plotSpread(data, 'showMM', 5,'xNames', {'PD ON', 'ND ON', 'PD OFF', 'ND OFF'})
      % ylim([-100 150])
      title('EPSC RATIO')
      nexttile 
      data = [iPrefRatioON iNullRatioON iPrefRatioOFF iNullRatioOFF]; 
      plotSpread(data, 'showMM', 5,'xNames', {'PD ON', 'ND ON', 'PD OFF', 'ND OFF'})
        % ylim([-100 150])
        title('IPSC RATIO')

      % correlations 
figure('Name','EPSC HEX data')
tiledlayout(2,2);

% EPSCs ON PREF
nexttile
hold on
plot([0 500], [0 500], '--', 'Color', [0.5 0.5 0.5]);  % Unity line
for i = 1:numColors
    plot(ePrefPreON(i), ePrefPostON(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
% mean ± SD
mX = mean(ePrefPreON);   mY = mean(ePrefPostON);
sX = std(ePrefPreON);    sY = std(ePrefPostON);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 500]);   ylim([0 500])
xlabel('EPSC CTRL');   ylabel('EPSC HEX')
title('EPSCs PREF ON')

% EPSCs OFF PREF
nexttile
hold on
plot([0 500], [0 500], '--', 'Color', [0.5 0.5 0.5]);
for i = 1:numColors
    plot(ePrefPreOFF(i), ePrefPostOFF(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
mX = mean(ePrefPreOFF);   mY = mean(ePrefPostOFF);
sX = std(ePrefPreOFF);    sY = std(ePrefPostOFF);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 500]);   ylim([0 500])
xlabel('EPSC CTRL');   ylabel('EPSC HEX')
title('EPSCs PREF OFF')

% EPSCs ON NULL
nexttile
hold on
plot([0 500], [0 500], '--', 'Color', [0.5 0.5 0.5]);
for i = 1:numColors
    plot(eNullPreON(i), eNullPostON(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
mX = mean(eNullPreON);   mY = mean(eNullPostON);
sX = std(eNullPreON);    sY = std(eNullPostON);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 500]);   ylim([0 500])
xlabel('EPSC CTRL');   ylabel('EPSC HEX')
title('EPSCs PREF ON')

% EPSCs OFF PREF
nexttile
hold on
plot([0 500], [0 500], '--', 'Color', [0.5 0.5 0.5]);
for i = 1:numColors
    plot(eNullPreOFF(i), eNullPostOFF(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
mX = mean(eNullPreOFF);   mY = mean(eNullPostOFF);
sX = std(eNullPreOFF);    sY = std(eNullPostOFF);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 500]);   ylim([0 500])
xlabel('EPSC CTRL');   ylabel('EPSC HEX')
title('EPSCs NULL OFF')


figure('Name','IPSC HEX data')
tiledlayout(2,2);

% IPSCs ON PREF
nexttile
hold on
plot([0 1500],[0 1500],'--','Color',[0.5 0.5 0.5]);  % Unity line
for i = 1:numColors
    plot(iPrefPreON(i), iPrefPostON(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
% add mean ± SD
mX = mean(iPrefPreON);  mY = mean(iPrefPostON);
sX = std(iPrefPreON);   sY = std(iPrefPostON);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 1500]); ylim([0 1500])
xlabel('IPSC CTRL'); ylabel('IPSC HEX')
title('IPSCs PREF ON')

% IPSCs OFF PREF
nexttile
hold on
plot([0 1500],[0 1500],'--','Color',[0.5 0.5 0.5]);
for i = 1:numColors
    plot(iPrefPreOFF(i), iPrefPostOFF(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
mX = mean(iPrefPreOFF);  mY = mean(iPrefPostOFF);
sX = std(iPrefPreOFF);   sY = std(iPrefPostOFF);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 1500]); ylim([0 1500])
xlabel('IPSC CTRL'); ylabel('IPSC HEX')
title('IPSCs PREF OFF')

% IPSCs ON NULL
nexttile
hold on
plot([0 1500],[0 1500],'--','Color',[0.5 0.5 0.5]);
for i = 1:numColors
    plot(iNullPreON(i), iNullPostON(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
mX = mean(iNullPreON);  mY = mean(iNullPostON);
sX = std(iNullPreON);   sY = std(iNullPostON);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 1500]); ylim([0 1500])
xlabel('IPSC CTRL'); ylabel('IPSC HEX')
title('IPSCs PREF ON')

% IPSCs OFF NULL
nexttile
hold on
plot([0 1500],[0 1500],'--','Color',[0.5 0.5 0.5]);
for i = 1:numColors
    plot(iNullPreOFF(i), iNullPostOFF(i), 'o', ...
         'MarkerFaceColor', colors(i,:), ...
         'MarkerEdgeColor', colors(i,:));
end
mX = mean(iNullPreOFF);  mY = mean(iNullPostOFF);
sX = std(iNullPreOFF);   sY = std(iNullPostOFF);
errorbar(mX, mY, sY, sY, sX, sX, 's', ...
         'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'CapSize', 8);
hold off
xlim([0 1500]); ylim([0 1500])
xlabel('IPSC CTRL'); ylabel('IPSC HEX')
title('IPSCs NULL OFF')


%% PREF VS NULL COMPARISONS
% figure('Name','PRE HEX PREF VS NULL')
% tiledlayout(2,2);
% 
% % EPSCs PRE HEX ON
% nexttile
% hold on
% for i = 1:numColors
%     plot(ePrefPreON(i), eNullPreON(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs PRE HEX ON')
% 
% % EPSCs PRE HEX OFF
% nexttile
% hold on
% for i = 1:numColors
%     plot(ePrefPreOFF(i), eNullPreOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs PRE HEX OFF')
% 
% % IPSCs PRE HEX ON
% nexttile
% hold on
% for i = 1:numColors
%     plot(iPrefPreON(i), iNullPreON(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs PRE HEX ON')
% 
% % IPSCs PRE HEX OFF
% nexttile
% hold on
% for i = 1:numColors
%     plot(iPrefPreOFF(i), iNullPreOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs PRE HEX OFF')
% 
% %% POST-DRUG PREF VS NULL
% figure('Name','HEX PREF VS NULL')
% tiledlayout(2,2);
% 
% % EPSCs HEX ON
% nexttile
% hold on
% for i = 1:numColors
%     plot(ePrefPostON(i), eNullPostON(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs HEX ON')
% 
% % EPSCs HEX OFF
% nexttile
% hold on
% for i = 1:numColors
%     plot(ePrefPostOFF(i), eNullPostOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs HEX OFF')
% 
% % IPSCs HEX ON
% nexttile
% hold on
% for i = 1:numColors
%     plot(iPrefPostON(i), iNullPostON(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs HEX ON')
% 
% % IPSCs HEX OFF
% nexttile
% hold on
% for i = 1:numColors
%     plot(iPrefPostOFF(i), iNullPostOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs HEX OFF')
% 
% %% PRE-DRUG ON VS OFF
% figure('Name','PRE HEX ON VS OFF')
% tiledlayout(2,2);
% 
% % EPSCs PRE HEX PREF
% nexttile
% hold on
% for i = 1:numColors
%     plot(ePrefPreON(i), ePrefPreOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs PRE HEX PREF')
% 
% % EPSCs PRE HEX NULL
% nexttile
% hold on
% for i = 1:numColors
%     plot(eNullPreON(i), eNullPreOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs PRE HEX NULL')
% 
% % IPSCs PRE HEX PREF
% nexttile
% hold on
% for i = 1:numColors
%     plot(iPrefPreON(i), iPrefPreOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs PRE HEX PREF')
% 
% % IPSCs PRE HEX NULL
% nexttile
% hold on
% for i = 1:numColors
%     plot(iNullPreON(i), iNullPreOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs PRE HEX NULL')
% 
% %% POST-DRUG ON VS OFF
% figure('Name','HEX ON VS OFF')
% tiledlayout(2,2);
% 
% % EPSCs HEX PREF
% nexttile
% hold on
% for i = 1:numColors
%     plot(ePrefPostON(i), ePrefPostOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs HEX PREF')
% 
% % EPSCs HEX NULL
% nexttile
% hold on
% for i = 1:numColors
%     plot(eNullPostON(i), eNullPostOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 500])
% ylim([0 500])
% title('EPSCs HEX NULL')
% 
% % IPSCs HEX PREF
% nexttile
% hold on
% for i = 1:numColors
%     plot(iPrefPostON(i), iPrefPostOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs HEX PREF')
% 
% % IPSCs HEX NULL
% nexttile
% hold on
% for i = 1:numColors
%     plot(iNullPostON(i), iNullPostOFF(i), 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', colors(i,:));
% end
% hold off
% xlim([0 1500])
% ylim([0 1500])
% title('IPSCs HEX NULL')
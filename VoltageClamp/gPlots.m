load gTableWTvert.mat
% wtTable = LOOK;
wtTable =iTable; 
load gTableB2vert.mat
% b2Table = b2Table;
b2Table = iTable;


num_recs = size(wtTable,1);


% for i = 1:num_recs
% 
%     figure 
% tiledlayout(2,1);
%     nexttile 
%     plot(iTable.gInhND{i})
%      ylim([0 1.8e-8])
%     title('ginhND')
%     nexttile 
%     plot(iTable.gInhPD{i})
%     ylim([0 1.8e-8])
%     title('ginhPD')
% 
%     figure 
% tiledlayout(2,1);
%     nexttile 
%     plot(eTable.gExcND{i})
%     title('gexcPD')
%     ylim([0 1.8e-8])
%     nexttile 
%     plot(eTable.gExcPD{i})
%     title('gexcND')
%     ylim([0 1.5e-8])
% 
%     figure 
% tiledlayout(2,1);
% nexttile
% plot(iTable.vmPD{i})
% title('vmPD')
% ylim([-5 40])
% nexttile
% plot(iTable.vmND{i})
% title('vmND')
% ylim([-5 40])
% 
% end 


%WT VM
% 
% for i = 1:num_recs
% 
%     cellNumb = wtTable.cellNumb(i);
%     expDate = wtTable.expDate(i);
%     acqInfo = sprintf(' WT Neuron %g Date %g',cellNumb,expDate);
% 
%     data = [wtTable.vmPD{i}];
% 
%     %pref and null are flipped, top is PD vm response, bottom is ND vm
%     %response
% 
%     figure('Name','vm pref')
% tiledlayout(2,1);
%     nexttile 
%     plot(data)
%    ylim([-70 0]);
% 
%     title(acqInfo);
% 
%     nexttile
%     data = [wtTable.vmND{i}];
% plot(data);
% ylim([-70 0]);
% end
% 
% %B2 VM
% num_recs = size(b2Table,1);
% for i = 1:num_recs
% 
%     cellNumb = b2Table.cellNumb(i);
%     expDate = b2Table.expDate(i);
%     acqInfo = sprintf(' B2 Neuron %g Date %g',cellNumb,expDate);
% 
%     data = [b2Table.vmPD{i}];
% 
%     %pref and null are flipped, top is ND vm response, bottom is PD vm
%     %response
% 
%     figure('Name','vm pref')
% tiledlayout(2,1);
%     nexttile 
%     plot(data)
%    ylim([-70 0]);
%     title(acqInfo);
% 
%     nexttile
%     data = [b2Table.vmND{i}];
% plot(data);
% ylim([-70 0]);
% end


% % plotting input resistance 
%  wtRin = wtTable.mRin; %Ohms 
%  b2Rin = b2Table.mRin; %Ohms 
% %For WT
% data = wtRin;
% 
% figure('Name',' WT Input resistance')
% tiledlayout(1,2)
% nexttile
% plotSpread(data,'showMM',5,'xNames',{'WT'},...
%     'yLabel','Input resistance')
%     ylim([0 1000])
%     title('WT Rin')
% %FOR B2
% nexttile
%     data = b2Rin;
% 
% plotSpread(data,'showMM',5,'xNames',{'B2'},...
%     'yLabel','Input resistance')
%     ylim([0 1000])
%     title('KO Rin')
% 
% 
    % plottingDSI
     % for ON and OFF WT
  data = [wtTable.vmDSIon wtTable.vmDSIoff];
   figure('Name','DSI plot')
    tiledlayout(1,2);
    nexttile 

    plotSpread(data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('WT Vm DSI')

nexttile

    % for ON and OFF
    b2data = [b2Table.vmDSIon b2Table.vmDSIoff];

    plotSpread(b2data,'showMM',5,'xNames',{'ON','OFF'},...
    'yLabel','Direction Selectivity Index')
    ylim([0 1])
    title('B2 Vm DSI')

% %plotting input vs series resistance 
%     b2Rin = b2Table.mRin; 
%     b2Rs = b2Table.mRs; 
%     b2data = [b2Rin b2Rs];
%     b2percent = (b2Rs/b2Rin)*100; 
% 
%     wtRin = wtTable.mRin; 
%     wtRs = wtTable.mRs; 
%     wtdata = [wtRin wtRs];
%     wtpercent = (wtRs/wtRin)*100; 
% 
% %FOR WT 
% % Generate a color map for plotting
% colors = lines(size(wtdata, 1)); % Adjust the color map based on the number of cells
% 
% % Create a new figure
% figure;
% 
% % Plot Rin on the y-axis and Rs on the x-axis
% for i = 1:size(wtdata, 1)
%     plot(wtdata(i, 2), wtdata(i, 1), 'o', 'MarkerFaceColor', colors(i, :), 'MarkerSize', 5); % Adjust marker size as needed
%     hold on;
% end
% 
% % Add labels and title
% xlabel('Series Resistance (Rs)');
% ylabel('Input Resistance (Rin)');
% title(' WT Rin vs Rs');
% 
% %FOR KO
% figure;
% colors = lines(size(b2data, 1)); % Adjust the color map based on the number of cells
% % Plot Rin on the y-axis and Rs on the x-axis
% for i = 1:size(b2data, 1)
%     plot(b2data(i, 2), b2data(i, 1), 'o', 'MarkerFaceColor', colors(i, :), 'MarkerSize', 5); % Adjust marker size as needed
%     hold on;
% end
% 
% % Add labels and title
% xlabel('Series Resistance (Rs)');
% ylabel('Input Resistance (Rin)');
% title(' B2 Rin vs Rs');
% 
% colors = lines(size(wtdata, 1)); % Adjust the color map based on the number of cells
% 
% % Calculate the percentage of Rs over Rin for each cell
% percent_Rs_Rin = (wtdata(:, 2) ./ wtdata(:, 1)) * 100;
% 
% % Create a new figure
% figure;
% 
% % Plot the percentage of Rs over Rin for each cell with the same color
% for i = 1:size(wtdata, 1)
%     plot(i, percent_Rs_Rin(i), 'o', 'MarkerFaceColor', colors(i, :), 'MarkerSize', 5); % Adjust marker size as needed
%     hold on;
% end
% 
% % Add labels and title
% xlabel('Cell');
% ylabel('% Rs / Rin');
% title(' WT Rs/Rin');
% 
% colors = lines(size(b2data, 1)); % Adjust the color map based on the number of cells
% 
% % Calculate the percentage of Rs over Rin for each cell
% percent_Rs_Rin = (b2data(:, 2) ./ b2data(:, 1)) * 100;
% 
% % Create a new figure
% figure;
% 
% % Plot the percentage of Rs over Rin for each cell with the same color
% for i = 1:size(b2data, 1)
%     plot(i, percent_Rs_Rin(i), 'o', 'MarkerFaceColor', colors(i, :), 'MarkerSize', 5); % Adjust marker size as needed
%     hold on;
% end
% 
% % Add labels and title
% xlabel('Cell');
% ylabel('% Rs / Rin');
% title('B2 Rs/Rin');
% 
% % comparing difference in VM peaks 
% 
% wtVmPDon = (wtTable.vmPDonPeak)*1000;
% wtVmPDoff = (wtTable.vmPDoffPeak)*1000;
% wtVmPD = [wtVmPDon wtVmPDoff]; 
% wtVmPDdiff = abs(wtVmPDon-wtVmPDoff);
% colors = lines(size(wtVmPD, 1));
%    figure('Name','WT Peak VM')
%     tiledlayout(1,4);
%     nexttile 
% 
%     plotSpread(wtVmPD, 'showMM', 5,'xNames', {'ON' 'OFF'})
%         for k = 1 :size(wtVmPD)
%             plot([wtVmPDon(k) wtVmPDoff(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor',colors(k, :));
% %             xlim([0 3])
%             ylim([0 60])
%             ylabel('Peak Vm (mV)')
%             title('Preferred Direction')
%             hold on
%         end
% nexttile 
% for k = 1 :size(wtVmPD)
% plotSpread(wtVmPDdiff, 'showMM', 5,'categoryColors', colors(k,:), 'xNames', {'ON' 'OFF'})
%             ylabel('delta ON/OFF')
%             ylim([0 35])
%             title('Preferred Direction')
% end
% 
% wtVmNDon = (wtTable.vmNDonPeak);
% wtVmNDoff = (wtTable.vmNDoffPeak);
% wtVmND = [wtVmNDon wtVmNDoff];
% wtVmNDdiff = abs(wtVmNDon-wtVmNDoff);
% colors = lines(size(wtVmND, 1));
% nexttile 
% 
%  plotSpread(wtVmND, 'showMM', 5,'xNames', {'ON' 'OFF'})
%         for k = 1 :size(wtVmND)
%             plot([wtVmNDon(k) wtVmNDoff(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor',colors(k, :));
% %             xlim([0 3])
%             ylim([0 60])
%             ylabel('Peak Vm (mV)')
%             title('Null Direction')
%             hold on
%         end
% nexttile 
% for k = 1 :size(wtVmND)
% plotSpread(wtVmNDdiff, 'showMM', 5,'categoryColors', colors(k,:), 'xNames', {'ON' 'OFF'})
%             ylabel('delta ON/OFF')
%             ylim([0 35])
%             title('Null Direction')
% end
% 

% %fOR b2 
% b2VmPDon = (b2Table.vmPDonPeak)*1000;
% b2VmPDoff = (b2Table.vmPDoffPeak)*1000;
% b2VmPD = [b2VmPDon b2VmPDoff]; 
% b2VmPDdiff = abs(b2VmPDon-b2VmPDoff);
% 
% 
% colors = lines(size(b2VmPD, 1));
%    figure('Name','B2 Peak VM')
%     tiledlayout(1,4);
%     nexttile 
% 
%     plotSpread(b2VmPD, 'showMM', 5,'xNames', {'ON' 'OFF'})
%         for k = 1 :size(b2VmPD)
%             plot([b2VmPDon(k) b2VmPDoff(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor',colors(k, :));
% %             xlim([0 3])
%             ylim([0 60])
%             ylabel('Peak Vm (mV)')
%             title('Preferred Direction')
%             hold on
%         end
% nexttile 
% for k = 1 :size(b2VmPD)
% plotSpread(b2VmPDdiff, 'showMM', 5,'categoryColors', colors(k,:), 'xNames', {'ON' 'OFF'})
%             ylabel('delta ON/OFF')
%             ylim([0 35])
%             title('Preferred Direction')
% end
% 
% b2VmNDon = (b2Table.vmNDonPeak)*1000;
% b2VmNDoff = (b2Table.vmNDoffPeak)*1000;
% b2VmND = [b2VmNDon b2VmNDoff];
% b2VmNDdiff = abs(b2VmNDon-b2VmNDoff);
% colors = lines(size(b2VmND, 1));
% nexttile 
% 
%     plotSpread(b2VmND, 'showMM', 5,'xNames', {'ON' 'OFF'})
%         for k = 1 :size(b2VmND)
%             plot([b2VmNDon(k) b2VmNDoff(k)],...
%                 'ro-', 'MarkerSize', 5,'MarkerFaceColor',colors(k, :));
% %             xlim([0 3])
%             ylim([0 60])
%             ylabel('Peak Vm (mV)')
%             title('Null Direction')
%             hold on
%         end
% nexttile 
% for k = 1 :size(b2VmND)
% plotSpread(b2VmNDdiff, 'showMM', 5,'categoryColors', colors(k,:), 'xNames', {'ON' 'OFF'})
%             ylabel('delta ON/OFF')
%             ylim([0 35])
%             title('Null Direction')
% end
% 



%comparing VM peaks across geno 

 %for WT ON 
    vmPrefON = (wtTable.vmPDonPeak);
    vmNullON = (wtTable.vmNDonPeak); 
    data = [vmPrefON vmNullON];

    figure('Name','Peak subthreshold Vm')
    tiledlayout(2,2);
    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1:size(vmNullON)
            plot([vmPrefON(k) vmNullON(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');

            ylim([0 40])
            ylabel('Membrane potential (mV)')
            title('WT ON')
            hold on
        end


  %for B2 ON 
    vmPrefON = (b2Table.vmPDonPeak);
    vmNullON = (b2Table.vmNDonPeak); 
    data = [vmPrefON vmNullON];
nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1:size(vmNullON)
            plot([vmPrefON(k) vmNullON(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');

            ylim([0 40])
            ylabel('Membrane potential (mV)')
            title('B2 ON')
            hold on
        end

  %for WT OFF
    vmPrefOFF= (wtTable.vmPDoffPeak);
    vmNullOFF = (wtTable.vmNDoffPeak); 
    data = [vmPrefOFF vmNullOFF];

    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1:size(vmNullOFF)
            plot([vmPrefOFF(k) vmNullOFF(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');

            ylim([0 40])
            ylabel('Membrane potential (mV)')
            title('WT OFF')
            hold on
        end


        %for B2 OFF
    vmPrefOFF= (b2Table.vmPDoffPeak);
    vmNullOFF =(b2Table.vmNDoffPeak); 
    data = [vmPrefOFF vmNullOFF];

    nexttile 

    plotSpread(data, 'showMM', 5,'xNames', {'Pref' 'Null'})
        for k = 1:size(vmNullOFF)
            plot([vmPrefOFF(k) vmNullOFF(k)],...
                'ro-', 'MarkerSize', 5,'MarkerFaceColor', 'r');

            ylim([0 40])
            ylabel('Membrane potential (mV)')
            title('B2 OFF')
            hold on
        end
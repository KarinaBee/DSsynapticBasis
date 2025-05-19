%load d table
load dTableB2vert.mat
dTable = dTable(strcmp(dTable.stim,"bars"),:);
eTable = dTable(strcmp(dTable.holdVol, "-60") & strcmp(dTable.drugs, "none"),:);
iTable = dTable(strcmp(dTable.holdVol, "12") & strcmp(dTable.drugs, "none"),:);

% load dTableWThorzForG.mat
% dTable = dTable(strcmp(dTable.stim,"bars"),:);
% eTable = dTable(strcmp(dTable.holdVol, "-60") & strcmp(dTable.drugs, "none"),:);
% iTable = dTable(strcmp(dTable.holdVol, "12") & strcmp(dTable.drugs, "none"),:);

num_recs = size(eTable.prefDirTrace,1); 

for i = 1:num_recs

   %calculate conductances from the currents for KO 
   %flipped the directions now ND is Null of CELL and PD is pref of CELL 
    iND = iTable.prefDirTrace{i,:}; % in pA
    iPD = iTable.nullDirTrace{i,:}; % in pA
    iND = (iND.*10^-12);% in A
    iPD = (iPD.*10^-12);%in A
    gInhPD = abs(iPD/72)*10^3; 
    gInhND = abs(iND/72)*10^3;
% Assuming gInhPD contains your conductance data with time information
% Assuming gInhPD is a 1x45000 double array

% Set the first 10000 ms to zero
gInhND = gInhND(10000:end);
gInhPD = gInhPD(10000:end);
figure 
tiledlayout(2,1);
    nexttile 
    plot(gInhND)
     ylim([0 2e-8])
    title('ginhND')
    nexttile 
    plot(gInhPD)
    ylim([0 2e-8])
    title('ginhPD')
% % 
% figure 
% tiledlayout(2,1);
%     nexttile 
%     plot(iND)
%     title('iND')
%     ylim([0 10e-10])
%     nexttile 
%     plot(iPD)
%     title('iPD')
%     ylim([0 10e-10])

    eND = eTable.prefDirTrace{i,:};
    ePD = eTable.nullDirTrace{i,:}; 
    eND = (eND.*10^-12); 
    ePD = (ePD.*10^-12);
    gExcPD = abs(ePD/-72)*10^3; 
    gExcND = abs(eND/-72)*10^3;
gExcPD = gExcPD(10000:end);
gExcND = gExcND(10000:end);

figure 
tiledlayout(2,1);
    nexttile 
    plot(gExcND)
    title('gexcPD')
    ylim([0 2e-8])
    nexttile 
    plot(gExcPD)
    title('gexcND')
    ylim([0 2e-8])

%     figure 
% tiledlayout(2,1);
%     nexttile 
%     plot(abs(eND))
%     title('eND')
%     ylim([0 10e-10])
%     nexttile 
%     plot(abs(ePD))
%     title('ePD')
%     ylim([0 10e-10])

    %get gLeak (inverse of Rin) 
    % Rin = iTable.Rin(i,:); %Ohms 
    % mRin = Rin.*10^-6; %mOhms
%     Rs = iTable.Rs(i,:);
%     mRs = Rs.*10^-6;
%     gLeak = (1/Rin);     
gLeak = 3.3e-9;

    %save those conductances 
    iTable.gInhPD{i} = gInhPD;
    iTable.gInhND{i} = gInhND;
    iTable.gExcPD{i} = gExcPD;
    iTable.gExcND{i} = gExcND;

    %Save gLeak 
    iTable.gLeak(i) = gLeak; 
    eTable.gLeak(i) = gLeak; 

    %save Rin 
    % iTable.Rin(i) = Rin; 
    % eTable.Rin(i) = Rin; 
    % iTable.mRin(i) = mRin; 
    % eTable.mRin(i) = mRin; 
%     iTable.mRs(i) = mRs; 
%     eTable.mRs(i) = mRs; 


%use mathews equation 
%pref and null are flipped for VM!!
[Vm] = runModel(gExcPD,gInhPD,gLeak);

VmPD = Vm; 

[Vm] = runModel(gExcND,gInhND,gLeak);

VmND = Vm; 

y = mean(VmND(1:400,:)); %averaging across the first .3 seconds 
VmND= VmND(:,:) - y;
z = mean(VmPD(1:400,:)); %averaging across the first .3 seconds 
VmPD= VmPD(:,:) - z;
% 
% minVMpd = abs(min(VmPD(5000:20000)));
% VmPD = VmPD + minVMpd; 
% minVMnd = abs(min(VmND(5000:20000)));
% VmND = VmND + minVMnd; 
% 
% 
figure 
tiledlayout(2,1);
nexttile
plot(VmPD)
title('vmPD')
ylim([-5 40])
nexttile
plot(VmND)
title('vmND')
ylim([-5 40])


%save Vm
   iTable.vmPD{i} = VmPD; 
   iTable.vmND{i} = VmND; 
   eTable.vmND{i} = VmND; 
   eTable.vmPD{i} = VmPD; 
% 


 x = size(VmND, 1);
% 
if x == 35001
    vmPDon = VmPD(10000:20000);
    vmPDonPeak = max(vmPDon); 
    vmPDoff = VmPD(20000:end);
    vmPDoffPeak = max(vmPDoff); 

    vmNDon = VmND(10000:20000);
    vmNDonPeak = max(vmNDon); 
    vmNDoff = VmND(20000:end);
    vmNDoffPeak = max(vmNDoff); 
else
    vmPDon = VmPD(1:10000);
    vmPDonPeak = max(vmPDon); 
    vmPDoff = VmPD(10000:end);
    vmPDoffPeak = max(vmPDoff); 

    vmNDon = VmND(1:10000);
    vmNDonPeak = max(vmNDon); 
    vmNDoff = VmND(10000:end);
    vmNDoffPeak = max(vmNDoff); 
end
% 
iTable.vmPDonPeak(i) = vmPDonPeak;
iTable.vmPDoffPeak(i) = vmPDoffPeak;
iTable.vmNDonPeak(i) = vmNDonPeak;
iTable.vmNDoffPeak(i) = vmNDoffPeak;
eTable.vmPDonPeak(i) = vmPDonPeak;
eTable.vmPDoffPeak(i) = vmPDoffPeak;
eTable.vmNDonPeak(i) = vmNDonPeak;
eTable.vmNDoffPeak(i) = vmNDoffPeak;
% 
vmDSIon = (vmPDonPeak - vmNDonPeak) / (vmPDonPeak + vmNDonPeak);
vmDSIoff = (vmPDoffPeak - vmNDoffPeak) / (vmPDoffPeak + vmNDoffPeak);

iTable.vmDSIon(i) = vmDSIon; 
iTable.vmDSIoff(i) = vmDSIoff; 
eTable.vmDSIon(i) = vmDSIon; 
eTable.vmDSIoff(i) = vmDSIoff;
% end

%% Functions
end 
function [Vm] = runModel(gExc,gInh,gLeak)
% Function to run a simple conductance model using forward euler method

%Constant properties 
        eExc = 0; % excitatory reversal potential
        eInh = -60; % inhibitory reversal potential
        eLeak = -55; % leak reversal potential
        cap = 8e-11; % 80pF capacitance %formerly 120 pF
        dt = 1e-4; % 100 microsecond sampling interval

% Intialize model
nPts = size(gInh,2);
Vm = NaN(nPts,1);
Vm(1,:) = eLeak;
dV = NaN;

% Integrate via forward euler method
for k = 1:(nPts - 1)%(paramsModel.nPts - 1)
    % Calculate change in V 
    dV = -( ...
        gExc(k) * (Vm(k) - eExc) + ... %exc conductance term
        gInh(k) * (Vm(k) - eInh) + ... %inh conductance term
        gLeak*(Vm(k) - eLeak)) ... % leak conductance term
        / cap; %capacitance

%  Print debugging
%         disp(['Iteration: ' num2str(k)]);
%         disp(['dV: ' num2str(dV)]);

        % Update membrane potential
        Vm(k+1,1) = Vm(k) + dV*dt; % numerically integrate
% 
% %         Print membrane potential
%         disp(['Vm(' num2str(k+1) '): ' num2str(Vm(k+1))]);


end

end
% 

function [newInd] = eyeCalSOS(stimDir,Eye)

       newInd = zeros(24,1); %new empty list for new stim dirs 
       temp = unique(stimDir(:,1));   

       %Right eye calibration
       right0 =180;
       right45 =135;
       right90 =90;
       right135 =45; 
       right180 =0;
       right225 =315;
       right270 =270;
       right315 =225; 

       if Eye == "R"
        replaceDir = [right0; right45; right90; right135; right180; right225; right270; right315];
       end 

       %Left/ventral eye calibration 
       left0 = 0; 
       left45 = 45;
       left90 = 90; 
       left135 = 135;
       left180 = 180;
       left225 = 225;
       left270 = 270;
       left315 = 315; 

       if Eye == "L"
        replaceDir = [left0; left45; left90; left135; left180; left225; left270; left315];
       end 


        %Left/ventral eye calibration 
       left0 = 180; 
       left45 = 225;
       left90 = 270; 
       left135 = 315;
       left180 = 0;
       left225 = 45;
       left270 = 90;
       left315 = 135; 

       if Eye == "F"
        replaceDir = [left0; left45; left90; left135; left180; left225; left270; left315];
       end 

       for ii = 1:numel(temp)
        indices = stimDir(:,1) == temp(ii);
        newInd(indices) = replaceDir(ii);
       end
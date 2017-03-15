function [ssBadFlag, ssS_B, ssNumDiodes, ssModule] = ss(ssCounts, ...
    data_valid, ssDiodeBoresight, ssCent1Vec, ssCent2Vec, ssCountsMin, ssCountsMax, sensorTimeout, useSensors)
%#codegen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name:    neascout_gnc_fsw/ProcessSensors/ssProcessing
% Project:      NEA Scout
% Author(s):    Chris Becker  
%               NASA Marshall Space Flight Center
% Last Mod:     2017-03-03 (Chris Becker)
%
%---------------------------------------------------------------------------------------------------
% INPUTS:
%   ssCounts- Sun sensor diode count values (uint16, 3x4)
%   data_valid- Hardware data validity flag from FSW (uint8)
%   ssDiodeBoresight- Boresight vectors for sun sensor diodes (double, 3x12)
%   ssCent1Vec- Sun vector to use for single diode illuminations (double, 3x12)
%   ssCent2Vec- Sun vector to use for diode pair illuminations (double, 3x12)
%   ssCountsMin- Lower threshold value of sun sensor measurement (uint16)
%   ssCountsMax- Upper threshold value of sun sensor measurement (uint16)
%   sensorTimeout- FSW execution cycles before sensor timeout (uint8)
%   useSensors- Command to ignore sensor data (uint8, 5x1)
%               
% OUTPUTS:
%   ssBadFlag- Sun sensor measurement status flag (int8)
%       0: Good
%       1: Data not valid from fsw
%       2: No diode count values in sun detection range.
%       3: Sensor disabled.
%       5: Initial value set before any measurements processed.
%       6: More than one sensor shows valid (not possible, so all are ignored)
%       7: All modules errored, but different errors
%       8: Calculation error
%       9: Data from ss out of valid count bounds
%   ssS_B- Sun sensor sun vector in body frame (double, 3x1)
%   ssNumDiodes- Number of illuminated diodes of selected sun sensor module (int8)
%   ssModule- Selected sun sensor module (int8)
%
%---------------------------------------------------------------------------------------------------
% Overview: This function processes data from the sun sensors.
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output variables
ssS_B = zeros(3,1);
ssBadFlag = int8(0);
ssS_B_ = zeros(3,3);
ssBadFlag_ = zeros(1,3,'int8');
ssNumDiodes_ = zeros(1,3,'int8');
measSelected = int8(0);
ssNumDiodes = int8(0);
ssModule = int8(0);
%ssMeas_valid = ssMeas_validIn(1);
%ssMeas = ssMeasIn(:,1);

persistent ssCountsPrev ssCountsAvg ssCountsBuf ssCountsBufSum ssCountsAvgPos ssCountsCurPoints timeSinceValid ssBadFlagPrev ssS_BPrev ssNumDiodesPrev ssModulePrev
avgNumPoints = 50;
% fill persistents and set ssBadFlag
if isempty(ssCountsPrev)
    ssCountsPrev = zeros(3,4,'int16');
    ssCountsAvg = zeros(3,4,'int16');
    ssCountsBuf = zeros(3,4,avgNumPoints);
    ssCountsBufSum = zeros(3,4);
    ssCountsAvgPos = ones(3,4);
    ssCountsCurPoints = ones(3,4);
    timeSinceValid = sensorTimeout;
    ssBadFlag_ = [int8(5), int8(5), int8(5)];
    ssBadFlag = int8(5);
    
    % Previous output data
    ssS_BPrev = zeros(3,1);
    ssBadFlagPrev = int8(5);
    ssNumDiodesPrev = int8(0);
    ssModulePrev = int8(0);
end

%% Check validity flag from fsw
if (data_valid == uint8(1))
    timeSinceValid = uint8(0);
else
    timeSinceValid = timeSinceValid + uint8(1);
    
    if (timeSinceValid > sensorTimeout) % 1 second or longer since valid
        ssBadFlag = int8(1);
        return;
    else % return last returned data
        ssBadFlag = ssBadFlagPrev;
        ssS_B = ssS_BPrev;
        ssNumDiodes = ssNumDiodesPrev;
        ssModule = ssModulePrev;
        return;
    end
end

%% sun sensor data checks
for module = 1:3
    % Check valid flag
%     if ssMeas_valid(module) ~= 1
%         ssBadFlag_(module) = int8(2);
%     end
    
    % Data pre-checks
    if (useSensors(2+module) == uint8(0))
        ssBadFlag_(module) = int8(3);
    elseif (checkBounds(ssCounts(module,:),0,ssCountsMax) == int8(0)) % data within bounds
        ssBadFlag_(module) = int8(9);
    end
%     if sum(isfinite(ssCounts(module,:))) < 4
%         ssBadFlag_(module) = int8(9);
%     end
    
    % stale data check, not on a reboot cycle though
%     if ssBadFlag_(module) ~= 5
%         % difference between previous and current ss measurements
%         %if (ssModulePrev == module)
%             deltaMeas = abs(ssCounts(module,:) - ssCountsPrev(module,:));
% 
%             % check to see if ss measurement is stale
%             if max(deltaMeas) == 0
%                 ssTimeStale(module) = ssTimeStale(module) + 1;
%                 if (ssTimeStale(module) > 5)
%                     ssBadFlag_(module) = int8(7);
%                 end
%             else
%                 ssTimeStale(module) = 0;
%             end
%         %end
%     end
    
    % outputs
    if ssBadFlag_(module) == int8(0) % check if measurement is valid
        % Check which diodes are greater than minimum count value
        diodesHigh = zeros(4,1, 'int8');
        for diode = 1:4
            if (ssCounts(module, diode) >= ssCountsMin) % diode value above lower threshold
                diodesHigh(diode) = int8(1);
                
                % Update moving average
                if (ssCountsPrev(module, diode) < ssCountsMin) % value previously below min
                    % Reset moving average
                    ssCountsAvg(module, diode) = ssCounts(module, diode);
                    ssCountsAvgPos(module, diode) = 1;
                    ssCountsCurPoints(module, diode) = 1;
                    ssCountsBuf(module,diode,:) = zeros(1,avgNumPoints);
                    ssCountsBuf(module,diode,1) = ssCounts(module, diode);
                    ssCountsBufSum(module,diode) = ssCounts(module, diode);
                else
                    if (ssCountsAvgPos(module, diode) >= avgNumPoints)
                        ssCountsAvgPos(module, diode) = 1;
                    else
                        ssCountsAvgPos(module, diode) = ssCountsAvgPos(module, diode) + 1;
                    end
                        
                    ssCountsBufSum(module,diode) = ssCountsBufSum(module,diode) - double(ssCountsBuf(module, diode, ssCountsAvgPos(module,diode)));
                    ssCountsBuf(module,diode,ssCountsAvgPos(module, diode)) = ssCounts(module, diode);
                    ssCountsBufSum(module,diode) = ssCountsBufSum(module,diode) + double(ssCounts(module, diode));
                    
                    if (ssCountsCurPoints(module, diode) < avgNumPoints)
                        ssCountsCurPoints(module, diode) = ssCountsCurPoints(module, diode) + 1;
                    end
                    
                    newVal = ssCountsBufSum(module,diode) / ssCountsCurPoints(module, diode);
                    if (newVal > ssCountsMax)
                        ssCountsAvg(module, diode) = ssCountsMax;
                    elseif (newVal < 0)
                        ssCountsAvg(module, diode) = int16(0);
                    else
                        ssCountsAvg(module, diode) = int16(newVal);
                    end
                end
            else % clear diode average
                ssCountsAvg(module, diode) = uint16(0);
                ssCountsAvgPos(module, diode) = uint16(0);
                ssCountsCurPoints(module, diode) = 1;
                ssCountsBuf(module,diode,:) = zeros(1,avgNumPoints);
                ssCountsBuf(module,diode,1) = uint16(0);
                ssCountsBufSum(module,diode) = uint16(0);
            end
        end
        numValidDiodes = int8(sum(diodesHigh));
        ssNumDiodes_(module) = numValidDiodes;
        
        % Determine resulting sun vector
        if (numValidDiodes > 0)             
            if (numValidDiodes == 1) % one valid diode
                moduleStartLoc = 1;
                if (module == 2)
                    moduleStartLoc = 5;
                elseif (module == 3)
                    moduleStartLoc = 9;
                end
                    
                if (diodesHigh(1) == 1)
                    ssS_B_(:,module) = ssCent1Vec(:, moduleStartLoc);
                elseif (diodesHigh(2) == 1)
                    ssS_B_(:,module) = ssCent1Vec(:, moduleStartLoc + 1);
                elseif (diodesHigh(3) == 1)
                    ssS_B_(:,module) = ssCent1Vec(:, moduleStartLoc + 2);
                elseif (diodesHigh(4) == 1)
                    ssS_B_(:,module) = ssCent1Vec(:, moduleStartLoc + 3);
                end
                
                measSelected = int8(module);
                
            elseif (numValidDiodes == 2) % one valid diode
                moduleStartLoc2 = 1;
                if (module == 2)
                    moduleStartLoc2 = 7;
                elseif (module == 3)
                    moduleStartLoc2 = 13;
                end
                
                moduleStartLoc = 1;
                if (module == 2)
                    moduleStartLoc = 5;
                elseif (module == 3)
                    moduleStartLoc = 9;
                end
                
                upperBound = 1.05;
                lowerBound = 0.95;
                vec = [0.0; 0.0; 1.0];
                if (diodesHigh(1) == 1 && diodesHigh(2) == 1)
                    diodeRatio = double(ssCountsAvg(module,1))/double(ssCountsAvg(module,2));
                    
                    if (diodeRatio > upperBound || diodeRatio < lowerBound) % bias towards higher diode
                        vec = diodeRatio*ssCent1Vec(:, moduleStartLoc) + ssCent1Vec(:, moduleStartLoc+1);
%                     elseif (diodeRatio < lowerBound) % bias towards diode 2
%                         vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 1) + ssCent1Vec(:, moduleStartLoc);
                    else
                        vec = ssCent2Vec(:, moduleStartLoc2);
                    end
                    
                elseif (diodesHigh(1) == 1 && diodesHigh(3) == 1)
                    diodeRatio = double(ssCountsAvg(module,1))/double(ssCountsAvg(module,3));
                    if (diodeRatio > upperBound || diodeRatio < lowerBound) % bias towards higher diode
                        vec = diodeRatio*ssCent1Vec(:, moduleStartLoc) + ssCent1Vec(:, moduleStartLoc+2);
%                     elseif (diodeRatio < lowerBound) % bias towards diode 3
%                         vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 2) + ssCent1Vec(:, moduleStartLoc);
                    else
                        vec = ssCent2Vec(:, moduleStartLoc + 1);
                    end

                elseif (diodesHigh(1) == 1 && diodesHigh(4) == 1)
                    diodeRatio = double(ssCountsAvg(module,1))/double(ssCountsAvg(module,4));
                    if (diodeRatio > upperBound || diodeRatio < lowerBound) % bias towards higher diode
                        vec = diodeRatio*ssCent1Vec(:, moduleStartLoc) + ssCent1Vec(:, moduleStartLoc+3);
%                     elseif (diodeRatio < lowerBound) % bias towards diode 4
%                         vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 3) + ssCent1Vec(:, moduleStartLoc);
                    else
                        vec = ssCent2Vec(:, moduleStartLoc + 2);
                    end
                    
                elseif (diodesHigh(2) == 1 && diodesHigh(3) == 1)
                    diodeRatio = double(ssCountsAvg(module,2))/double(ssCountsAvg(module,3));
                    if (diodeRatio > upperBound || diodeRatio < lowerBound) % bias towards higher diode
                        vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 1) + ssCent1Vec(:, moduleStartLoc+2);
%                     elseif (diodeRatio < lowerBound) % bias towards diode 3
%                         vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 2) + ssCent1Vec(:, moduleStartLoc+1);
                    else
                        vec = ssCent2Vec(:, moduleStartLoc + 3);
                    end
                    
                elseif (diodesHigh(2) == 1 && diodesHigh(4) == 1)
                    diodeRatio = double(ssCountsAvg(module,2))/double(ssCountsAvg(module,4));
                    if (diodeRatio > upperBound || diodeRatio < lowerBound) % bias towards higher diode
                        vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 1) + ssCent1Vec(:, moduleStartLoc+3);
%                     elseif (diodeRatio < lowerBound) % bias towards diode 4
%                         vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 3) + ssCent1Vec(:, moduleStartLoc+1);
                    else
                        vec = ssCent2Vec(:, moduleStartLoc + 4);
                    end
                    
                elseif (diodesHigh(3) == 1 && diodesHigh(4) == 1)
                    diodeRatio = double(ssCountsAvg(module,3))/double(ssCountsAvg(module,4));
                    if (diodeRatio > upperBound || diodeRatio < lowerBound) % bias towards higher diode
                        vec = diodeRatio*ssCent1Vec(:, moduleStartLoc+2) + ssCent1Vec(:, moduleStartLoc+3);
%                     elseif (diodeRatio < lowerBound) % bias towards diode 4
%                         vec = diodeRatio*ssCent1Vec(:, moduleStartLoc + 3) + ssCent1Vec(:, moduleStartLoc+2);
                    else
                        vec = ssCent2Vec(:, moduleStartLoc + 5);
                    end
                end    
                ssS_B_(:,module) = vec/norm(vec);
                                    
                        
                
%                 if (diodesHigh(1) == 1 && diodesHigh(2) == 1)
%                     ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc);
%                 elseif (diodesHigh(1) == 1 && diodesHigh(3) == 1)
%                     ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 1);
%                 elseif (diodesHigh(1) == 1 && diodesHigh(4) == 1)
%                     ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 2);
%                 elseif (diodesHigh(2) == 1 && diodesHigh(3) == 1)
%                     ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 3);
%                 elseif (diodesHigh(2) == 1 && diodesHigh(4) == 1)
%                     ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 4);
%                 elseif (diodesHigh(3) == 1 && diodesHigh(4) == 1)
%                     ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 5);
%                 end  
                
                measSelected = int8(module);
                
            else % 3 or 4 valid diodes
                measSelected = int8(module);
                moduleStartLoc = 1;
                if (module == 2)
                    moduleStartLoc = 5;
                elseif (module == 3)
                    moduleStartLoc = 9;
                end
                
                % Solve system of equations for sun vector
                % (http://www.mathworks.com/help/matlab/ref/qr.html)
                A = zeros(4,3);
                b = zeros(4,1);
                if (diodesHigh(1) == 1)
                    A(1,:) = ssDiodeBoresight(:, moduleStartLoc)';
                end
                if (diodesHigh(2) == 1)
                    A(2,:) = ssDiodeBoresight(:, moduleStartLoc+1)';
                end
                if (diodesHigh(3) == 1)
                    A(3,:) = ssDiodeBoresight(:, moduleStartLoc+2)';
                end
                if (diodesHigh(4) == 1)
                    A(4,:) = ssDiodeBoresight(:, moduleStartLoc+3)';
                end
                
                %b = double(ssCounts(module,:)');
                diodeHalfFOV = 33*pi/180.0;
                %b = double(cos((-double(ssCounts(module,:)') + double(ssCountsMax))*diodeHalfFOV/double(ssCountsMax-1600))); % counts assumed to be linear with angle from boresight
                
                % if (issparse(A))
%                     R = qr(A); 
%                 else
%                     R = triu(qr(A)); 
%                 end
%                 x = R\(R'\(A'*b));
%                 r = b - A*x;
%                 err = R\(R'\(A'*r));
%                 x = x + err;
                
                %x = mldivide(A,b)
                
                b = double(ssCountsAvg(module,:)');
                x = mldivide(A,b);
                
                normX = norm(x);
                if (norm(x) > eps)
                    ssS_B_(:,module) = x/normX;
                else
                    ssBadFlag_(module) = int8(8);
                end
            end
        else
            ssBadFlag_(module) = int8(2);
        end
        
        % Store current measurement value
        ssCountsPrev(module,diode) = ssCounts(module,diode);
    
        
    end
    
    % Store last values
    %if (measSelected ~= 0)
    %    ssCountsPrev(module,:) = ssCounts(module,:);
    %end
    
end

% Outputs and persistent state updates
if (ssBadFlag_(1) == int8(0) && (ssBadFlag_(2) == int8(0) || ssBadFlag_(3) == int8(0))) % check for conflicting sensor values from opposite end sensors
   ssBadFlag = int8(6); % ignore all sensor values
elseif (measSelected ~= 0)
    ssS_B = ssS_B_(:,measSelected);
    ssNumDiodes = ssNumDiodes_(measSelected);
    ssModule = measSelected;
else % all sensors bad
    if (ssBadFlag_(1) == ssBadFlag_(2) && ssBadFlag_(2) == ssBadFlag_(3)) % same bad reason
        ssBadFlag = ssBadFlag_(1); 
    else % just return general bad status
        ssBadFlag = int8(7);
    end
end

% Store output
ssBadFlagPrev = ssBadFlag;
ssNumDiodesPrev = ssNumDiodes;
ssModulePrev = ssModule;
ssS_BPrev = ssS_B;

end

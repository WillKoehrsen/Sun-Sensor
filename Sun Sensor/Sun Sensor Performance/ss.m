function [ssBadFlag, ssS_B, ssNumDiodes, ssModule] = ss(ssMeas_valid, ssCounts, ssDiodeBoresight, ssCent1Vec, ssCent2Vec, ssCountsMin, ssCountsMax)
%#codegen

% Validity checking, unit conversion and frame transformation for sun 
% sensors.
%
% Inputs:
%    ssMeas_valid     - valid flag from FSW GNC wrapper.
%    ssCounts         - sun sensor measurement (photodiode count values).
%    ssDiodeBoresight - central boresight vector of each diode on each sensor.
%    ssCent1Vec       - for each diode, sun vector to use when only that diode is illuminated.
%    ssCent2Vec       - sun vector to return for illuminated diode pairs.
%    ssCountsMin      - minimum count value to declare a diode illuminated.
%    ssCountsMax      - maximum count value.
%
% Outputs:
%    ssBadFlag       - Bad Flag (single summary value reported for status of all sensors)
%                       0: Good
%                       2: Valid flag is false
%                       5: Initial value set before any measurements processed.
%                       6: Both sensors show valid (not possible, so both are ignored)
%                       7: Stale measurement
%                       9: Data from ss not finite
%    ssS_B           - sun vector in body frame
%    ssNumDiodes     - number of diodes illuminated for a valid sensor measurement
%
% Devon Sanders, EV42
% Updated for analog sensors on NEAScout - Chris Becker, EV42, 5/01/2016

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
persistent ssCountsPrev ssModulePrev ssTimeStale

% fill persistents and set ssBadFlag
if isempty(ssCountsPrev)
    ssCountsPrev = zeros(3,4,'int16');
    ssModulePrev = int8(0);
    ssTimeStale = zeros(3,1,'uint8');
    ssBadFlag_ = [int8(5), int8(5), int8(5)];
    ssBadFlag = int8(5);
end

%% sun sensor data checks
for module = 1:3
    
    % Check valid flag
    if ssMeas_valid(module) ~= 1
        ssBadFlag_(module) = int8(2);
    end
    
    % isfinite check for incoming data
    if sum(isfinite(ssCounts(module,:))) < 4
        ssBadFlag_(module) = int8(9);
    end
    
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
    if ssBadFlag_(module) == 0 % only need one measurement
        % Check which diodes are greater than minimum count value
        diodesHigh = zeros(4,1, 'int8');
        for diode = 1:4
            if (ssCounts(module, diode) >= ssCountsMin)
                diodesHigh(diode) = int8(1);
            end
        end
        numValidDiodes = sum(diodesHigh);
        ssNumDiodes_(module) = numValidDiodes;
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
                moduleStartLoc = 1;
                if (module == 2)
                    moduleStartLoc = 7;
                elseif (module == 3)
                    moduleStartLoc = 13;
                end
                if (diodesHigh(1) == 1 && diodesHigh(2) == 1)
                    ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc);
                elseif (diodesHigh(1) == 1 && diodesHigh(3) == 1)
                    ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 1);
                elseif (diodesHigh(1) == 1 && diodesHigh(4) == 1)
                    ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 2);
                elseif (diodesHigh(2) == 1 && diodesHigh(3) == 1)
                    ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 3);
                elseif (diodesHigh(2) == 1 && diodesHigh(4) == 1)
                    ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 4);
                elseif (diodesHigh(3) == 1 && diodesHigh(4) == 1)
                    ssS_B_(:,module) = ssCent2Vec(:, moduleStartLoc + 5);
                    
                end
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
                
                b = double(ssCounts(module,:)');
                x = mldivide(A,b);
                
                ssS_B_(:,module) = x/norm(x);
            end
        end   
    end
    
    % Valid measurement selected (if both 2 and 3 are valid, 3 will be the returned measurement selected)
    if (measSelected ~= 0)
        ssCountsPrev(module,:) = ssCounts(module,:);
    end
    
end

% Outputs and persistent state updates
if (ssBadFlag_(1) == 0 && (ssBadFlag_(2) == 0 || ssBadFlag_(3) == 0)) % check for conflicting sensor values from opposite end sensors
   ssBadFlag = int8(6); % ignore all sensor values
elseif (measSelected ~= 0)
    ssModulePrev = measSelected;
    ssS_B = ssS_B_(:,measSelected);
    ssNumDiodes = ssNumDiodes_(measSelected);
    ssModule = measSelected;
else % all sensors bad
    if (ssBadFlag_(1) == ssBadFlag_(2) == ssBadFlag_(3)) % same bad reason
        ssBadFlag = ssBadFlag_(1); 
    else % just return general bad status
        ssBadFlag = int8(1);
    end
end    


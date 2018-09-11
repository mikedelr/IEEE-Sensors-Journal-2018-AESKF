classdef AESKF <handle
%AESKF is the AHRS algorithm for a 9DOF and or 6DOF MIMU/IMU 
% See 'Computationally Efficient Adaptive Error-State Kalman Filter for 
% Attitude Estimation'
%
% https://ieeexplore.ieee.org/document/8434202/
%
% information from each of the sensors has been decoupled and fused
% according to the tuning parameters specified by the user. The
% algorithm assumes that all sensors inputs have been calibrated
%
% Author(s): Michael Del Rosario
%            Phillip Ngo
%            Stephen Redmond 

properties (Access = public)    
    ctr       = 0; % idx corresponding to the number of samples processed
    ctrM      = 0; % idx corresponding to the number of samples processed
    sigma_gyr = NaN;
    
    delta_b_xy = NaN;
    incAngRef = NaN;
    mCosTheta = NaN;
    mSinTheta = NaN;
    
    BufGamM = [];
    N_VAR_GAM_M = NaN;  
    rSumGamM   = 0;
    rSumSqGamM = 0;
    rSumN2GamM = 0;
    AvgConstGamM = NaN;
    nAvgConstGamM = NaN;
    VarConstGamM = NaN;

    BufGamA = [];    
    rSumGamA   = 0;
    rSumSqGamA = 0;
    N_VAR_GAM_A = NaN;    
    AvgConstGamA = NaN;
    nAvgConstGamA = NaN;
    VarConstGamA = NaN;
    
    BufGamA2 = [];    
    rSumGamA2   = 0;
    rSumSqGamA2 = 0;
    N_VAR_GAM_A2 = NaN;    
    AvgConstGamA2 = NaN;
    nAvgConstGamA2 = NaN;
    VarConstGamA2 = NaN;

    N_VAR_GAM_AS = NaN;
    BufGamAS = [];
    rSumGamAS   = 0;
    rSumSqGamAS = 0;
    AvgConstGamAS =  NaN;   
    nAvgConstGamAS = NaN;
    VarConstGamAS = NaN; 
    
    N_VAR_GAM_AL = NaN;
    BufGamAL = [];
    rSumGamAL   = 0;
    rSumSqGamAL = 0;
    AvgConstGamAL = NaN;   
    nAvgConstGamAL = NaN;
    VarConstGamAL = NaN;          
    
    N_VAR_GAM_MS = NaN;
    BufGamMS = [];
    rSumGamMS   = 0;
    rSumSqGamMS = 0;
    AvgConstGamMS =  NaN;   
    nAvgConstGamMS = NaN;
    VarConstGamMS = NaN; 
    
    N_VAR_GAM_ML = NaN;
    BufGamML = [];
    rSumGamML   = 0;
    rSumSqGamML = 0;
    AvgConstGamML = NaN;   
    nAvgConstGamML = NaN;
    VarConstGamML = NaN;        
    
%     R_gamma_th = NaN;
    
    gRef    = 9.796720; % magnitude of acceleration due to gravity (Sydney)
    cA      = 0.1;      % normalised cutoff frequency
    ext_acc = [0 0 0];
    
    mRef    = 57.0732; % magnitude of the local geomagnetic field (Sydney)
    cM      = 0.1;     % normalised cutoff-frequency 
    gMagRef = [0 0 0];
    mag_dis = [0 0 0];
    mag_now = [0 0 0]; 
    mag_pre = [0 0 0]; 
    mag_now_xy = [0 0];
    mag_pre_xy = [0 0];
    gGyr = [0 0 0];
    gAcc = [0 0 0];
    gAcc_now = [0 0 0];
    gAcc_pre = [0 0 0];
    gMag = [0 0 0];
    
    deltaMagBxy = [NaN];
    
    fs          = 100; % additional constants defined for speed
    dt          = 0.01;
    qAngVelDev  = [1 0 0 0];
    qGlobal     = [1 0 0 0];
    qGlobalPrev = [1 0 0 0];        
    muAcc       = 0.005;
    muMag       = 0.005;
    
    % these are for storing the angles gamma_a and gamma_m
    gammaA      = 0;
    gammaM      = 0;    
    %% New variables for Adaptive gain via indirect kalman filter
    % --- Process Noise Covariance Matrix, in this case 1D --- %
    % Q = 2 * dt^{2} * \sigma_{gyr}^{2} i.e., the variance in the gyro when
    % stationary, user should specify the noise level of the gyro for
    % example, sigma = 0.015
    Q       = NaN; 
    R_GamA  = NaN; % variance in the measurement noise  
    R_GamAS = NaN;
    R_GamAL = NaN;
    AvgGamA = NaN;
    P_pri   = NaN; % A Priori error state covariance matrix
    P_pos   = NaN; % A Posteriori error state covariance matrix, can be initialised to Q
    muA     = 0.005; % adaptive kalman gain
    
    Qmag     = NaN; 
    R_GamM   = NaN; % variance in the measurement noise 
    R_GamMS = NaN;
    R_GamML = NaN;
    AvgGamM  = NaN;
    Pmag_pri = NaN; % initialise a posteriori as Q        
    Pmag_pos = NaN; % initialise a priori as Q
    muM      = 0.005;    
end

methods (Access = public)
    function obj = CAHRS_IKFv2(varargin)
        for i = 1:2:nargin
            if strcmp(varargin{i},'fs'), 
                obj.fs = varargin{i+1};
                obj.dt=1/varargin{i+1};
            elseif strcmp(varargin{i},'cA'),
                obj.cA = varargin{i+1};                
            elseif strcmp(varargin{i},'cM'), 
                obj.cM = varargin{i+1};
            elseif strcmp(varargin{i},'gRef'),
                obj.gRef = varargin{i+1};                
            elseif strcmp(varargin{i},'mRef'), 
                obj.mRef = varargin{i+1};                
            elseif strcmp(varargin{i},'q0'),
                obj.qGlobalPrev = varargin{i+1};
                obj.qGlobal = varargin{i+1};
%% --- CAHRS_Indirect Kalman Filter ---   
% add additional fields to support Indirect Kalman Filter
% constructor for CAHRS with indirect kalman filter for dynamic correction
% rate of mu_a and mu_m based on the variance of gamma_a and gamma_m                
            elseif strcmp(varargin{i},'N_VAR_GAM_AS');
                % number of samples to calculate variance in gamma_a, inclination angle error
                obj.N_VAR_GAM_AS = varargin{i+1};   
                % buffer containing the upward errors gamma_a
                obj.BufGamAS = zeros(obj.N_VAR_GAM_AS,1);
                % weighted average for gamma_a
                obj.AvgConstGamAS =  1/obj.N_VAR_GAM_AS;   
                obj.nAvgConstGamAS = -1/obj.N_VAR_GAM_AS;
                % weighted average for variance of gamma_a
                obj.VarConstGamAS = 1/(obj.N_VAR_GAM_AS-1); 
            elseif strcmp(varargin{i},'N_VAR_GAM_AL');
                obj.N_VAR_GAM_AL = varargin{i+1};
                obj.BufGamAL = zeros(obj.N_VAR_GAM_AL,1);
                % weighted average for gamma_a
                obj.AvgConstGamAL =  1/obj.N_VAR_GAM_AL;   
                obj.nAvgConstGamAL = -1/obj.N_VAR_GAM_AL;
                % weighted average for variance of gamma_a
                obj.VarConstGamAL = 1/(obj.N_VAR_GAM_AL-1);             
            elseif strcmp(varargin{i},'N_VAR_GAM_MS');
                % number of samples to calculate variance in gamma_a, inclination angle error
                obj.N_VAR_GAM_MS = varargin{i+1};   
                % buffer containing the upward errors gamma_a
                obj.BufGamMS = zeros(obj.N_VAR_GAM_MS,1);
                % weighted average for gamma_a
                obj.AvgConstGamMS =  1/obj.N_VAR_GAM_MS;   
                obj.nAvgConstGamMS = -1/obj.N_VAR_GAM_MS;
                % weighted average for variance of gamma_a
                obj.VarConstGamMS = 1/(obj.N_VAR_GAM_MS-1); 
            elseif strcmp(varargin{i},'N_VAR_GAM_ML');
                obj.N_VAR_GAM_ML = varargin{i+1};
                obj.BufGamML = zeros(obj.N_VAR_GAM_ML,1);
                % weighted average for gamma_a
                obj.AvgConstGamML =  1/obj.N_VAR_GAM_ML;   
                obj.nAvgConstGamML = -1/obj.N_VAR_GAM_ML;
                % weighted average for variance of gamma_a
                obj.VarConstGamML = 1/(obj.N_VAR_GAM_ML-1);  
            elseif strcmp(varargin{i},'N_VAR_GAM_M');
                % number of samples to calculate variance in gamma_m, yaw angle error                
                obj.N_VAR_GAM_M = varargin{i+1}; 
                % buffer containing the upward errors gamma_m
                obj.BufGamM = zeros(obj.N_VAR_GAM_M,1); 
                % weighted average for gamma_m
                obj.AvgConstGamM = 1/obj.N_VAR_GAM_M; 
                obj.nAvgConstGamM = -1/obj.N_VAR_GAM_M;
                % weighted average for variance of gamma_m
                obj.VarConstGamM = 1/(obj.N_VAR_GAM_M-1);     
            elseif strcmp(varargin{i},'sigma_gyr')
                obj.sigma_gyr = varargin{i+1};
            elseif  strcmp(varargin{i},'delta_b_xy')
                obj.delta_b_xy = varargin{i+1};
            elseif  strcmp(varargin{i},'incAngRef')
                obj.incAngRef = varargin{i+1}; 
            else
                error(['Invalid argument ',varargin{i}]);
            end
        end;
        
        obj.mCosTheta = obj.mRef*cos(obj.incAngRef);
        obj.mSinTheta = obj.mRef*sin(obj.incAngRef);
        
        if ~isnan(obj.sigma_gyr),
            % assume noise in the gyro is equal on the x/y axis
            obj.Q = 2*obj.dt^(2)*obj.sigma_gyr^(2); 
            obj.P_pos = obj.Q; % initialise a priori as Q
            obj.P_pri = obj.Q; % initialise a posteriori as Q

            obj.Qmag = obj.dt^(2)*obj.sigma_gyr^(2);%2*obj.dt^(2)*obj.sigma_gyr^(2); 
            obj.Pmag_pos = obj.Qmag; % initialise a priori as Q
            obj.Pmag_pri = obj.Qmag; % initialise a posteriori as Q            
        end
    end
   
%%
    function obj = Update(obj,Gyr,Acc,Mag)
        qG1 = obj.qGlobal(1); 
        qG2 = obj.qGlobal(2); 
        qG3 = obj.qGlobal(3); 
        qG4 = obj.qGlobal(4); 
        Gyr1 = Gyr(1);    
        Gyr2 = Gyr(2);    
        Gyr3 = Gyr(3);
        % Correct for gyros - get quaternion from gyro in device frame and rotate
        if (~isnan(Gyr1) && ~isnan(Gyr2) && ~isnan(Gyr3) ) 
            dt2 = 0.5*obj.dt; % x1      
            qW1 = 1; 
            qW2 = Gyr1*dt2;   % x1
            qW3 = Gyr2*dt2;   % x1
            qW4 = Gyr3*dt2;   % x1
            qWnorm = 1/sqrt(1 + qW2*qW2 + qW3*qW3 + qW4*qW4); % add 3, x3, invsqrt(1)
            qW1 = qWnorm; 
            qW2 = qW2*qWnorm; % x1
            qW3 = qW3*qWnorm; % x1
            qW4 = qW4*qWnorm; % x1
            obj.qAngVelDev = [qW1 qW2 qW3 qW4];
            % Convert to back to global frame
            qGi1 = qG1*qW1 - qG2*qW2 - qG3*qW3 - qG4*qW4; % sub 3, x4
            qGi2 = qG1*qW2 + qG2*qW1 + qG3*qW4 - qG4*qW3; % add 2, sub 1, x4
            qGi3 = qG1*qW3 - qG2*qW4 + qG3*qW1 + qG4*qW2; % add 2, sub 1, x4
            qGi4 = qG1*qW4 + qG2*qW3 - qG3*qW2 + qG4*qW1; % add 2, sub 1, x4
        else
            qGi1 = qG1; 
            qGi2 = qG2; 
            qGi3 = qG3; 
            qGi4 = qG4; 
        end
% ------------------------------------------------------------------------        
        % Correct for accelerometer - get next rotation required to align 
        % accelerometer/gravity with 'up'.
        % model acceleration due to gravity as a first order low pass filtered 
        % process, i.e., zAcc = a_k - cA * ext_acc;
        % express ext_acc in GFR, apply disturbance model
        % transform back to SFR
        yAcc = Acc;
        zAcc = yAcc - obj.cA * obj.ext_acc; % sub 3, x3 sensor frame
        Acc  = zAcc;
        Acc1 = Acc(1); 
        Acc2 = Acc(2); 
        Acc3 = Acc(3);          
        if ( ~isnan(Acc1) && ~isnan(Acc2) && ~isnan(Acc3) ) 
            % if moving adaptive gain with measurement estimation
            % intermediate calculation in calculating v' = qvq*
            qa1 = -qGi2*Acc1 - qGi3*Acc2 - qGi4*Acc3; % sub 2, x3
            qa2 =  qGi1*Acc1 + qGi3*Acc3 - qGi4*Acc2; % add 1, sub 1, x3
            qa3 =  qGi1*Acc2 - qGi2*Acc3 + qGi4*Acc1; % add 1, sub 1, x3
            qa4 =  qGi1*Acc3 + qGi2*Acc2 - qGi3*Acc1; % add 1, sub 1, x3 
            % Convert acceleration to global frame
            gAcc1 = -qa1*qGi2 + qa2*qGi1 - qa3*qGi4 + qa4*qGi3; % add 2, sub 2, x4
            gAcc2 = -qa1*qGi3 + qa2*qGi4 + qa3*qGi1 - qa4*qGi2; % add 2, sub 2, x4
            gAcc3 = -qa1*qGi4 - qa2*qGi3 + qa3*qGi2 + qa4*qGi1; % add 2, sub 2, x4
            obj.gAcc = [gAcc1,gAcc2,gAcc3];
            % Get fraction of rotation from acc in global to up
            gAcc2Sq = gAcc2*gAcc2; % x1
            gAcc1Sq = gAcc1*gAcc1; % x1
            gAcc12SumSq = gAcc1Sq + gAcc2Sq; % add 1
            axisVecNorm = 1/sqrt(gAcc12SumSq); % invsqrt(1) 
            if (gAcc12SumSq == 0) 
                qUp1 = 1; qUp2 = 0; qUp3 = 0; qUp4 = 0;
                obj.muA = NaN; % purely for visualization
            else
                nGamA      = arctan2_approx(axisVecNorm*gAcc12SumSq,gAcc3); % x1
                nGamASq    = nGamA * nGamA; % x1
                obj.gammaA = nGamA;     
% -------------------------------------------------------------------------                
                % Estimate measurement error, gamma_a using a sliding window
                % Run time calculation of var(gamma_a) used to estimate R_GamA
                obj.ctr = obj.ctr + 1; % add 1
                aCtr    = obj.ctr;
                idxL    = mod(aCtr-1, obj.N_VAR_GAM_AL)+1; % add 1 , sub 1, mod via bitshift because buffer is size(2^X)
                oGamAL  = obj.BufGamAL(idxL);              % the oldest value to subtract from the running sum
                obj.BufGamAL(idxL) = nGamA;                % replace oldest idx in BufGamA with newest value
                obj.rSumGamAL   = obj.rSumGamAL   + nGamA;   % add 1
                obj.rSumSqGamAL = obj.rSumSqGamAL + nGamASq; % add 1 
                if aCtr > obj.N_VAR_GAM_AL
                    obj.rSumGamAL   = obj.rSumGamAL   - oGamAL;          % sub 1
                    obj.rSumSqGamAL = obj.rSumSqGamAL - oGamAL * oGamAL; % sub 1, x1
                end
                nAvgCGamASq = ...
                    obj.nAvgConstGamAL * obj.rSumGamAL * obj.rSumGamAL;  % x2
                obj.R_GamAL = ... 
                    obj.VarConstGamAL * (obj.rSumSqGamAL + nAvgCGamASq); % add 1, x1 
                
                idxS = idxL - (obj.N_VAR_GAM_AS); % sub 1
                % check logic of idxS
                if idxS < 1
                    idxS = obj.N_VAR_GAM_AL + idxS; % add 1
                end
                oGamAS = obj.BufGamAL(idxS);
                obj.rSumGamAS   = obj.rSumGamAS   + nGamA;   % add 1
                obj.rSumSqGamAS = obj.rSumSqGamAS + nGamASq; % add 1 
                if aCtr > obj.N_VAR_GAM_AS
                    obj.rSumGamAS   = obj.rSumGamAS   - oGamAS;          % sub 1
                    obj.rSumSqGamAS = obj.rSumSqGamAS - oGamAS * oGamAS; % sub 1, x1
                end
                nAvgCGamASq_S = ...
                    obj.nAvgConstGamAS * obj.rSumGamAS * obj.rSumGamAS; % x2
                obj.R_GamAS = ... 
                    obj.VarConstGamAS * (obj.rSumSqGamAS + nAvgCGamASq_S); % add 1, x1 
                % set measurement noise based on sliding windows
                % only pick the smaller window variance if the external
                % acceleration is very small
                if ( obj.ext_acc(1) * obj.ext_acc(1) + ...
                        obj.ext_acc(2)*obj.ext_acc(2) + ...
                        obj.ext_acc(3)*obj.ext_acc(3) ) < 1 % add 2, x3
                    obj.R_GamA = obj.R_GamAS;    
                else
                    obj.R_GamA = obj.R_GamAL;    
                end
                obj.P_pri    = obj.P_pos + obj.Q;                  % add 1 
                P_pri_R_GamA = obj.P_pri + obj.R_GamA;             % add 1
                S_a_inv      = 1/sqrt( P_pri_R_GamA*P_pri_R_GamA ); % x1, invsqrt(1)
                mu_a         = obj.P_pri * S_a_inv;                 % x1
                % use a small fixed gain whilst the buffer is not full
                if aCtr < obj.N_VAR_GAM_AL, mu_a = 0.001; end 
                if aCtr < obj.N_VAR_GAM_AS, mu_a = 0.001; end
% -------------------------------------------------------------------------                
                obj.muA   = mu_a;
                obj.P_pos = (1 - mu_a) * obj.P_pri; % sub 1, x1
                muAgamA   = (mu_a*nGamA)*0.5;       % x2 
                qUp1      =  1;
                qUp2      =  gAcc2 * axisVecNorm * muAgamA;    % x2
                qUp3      = -gAcc1 * axisVecNorm * muAgamA;    % sub 1, x2
                qUpNorm   = 1/sqrt(1 + qUp2*qUp2 + qUp3*qUp3); % +2, x2, invsqrt(1)
                % normalise to unit quaternions
                qUp1 = qUpNorm; 
                qUp2 = qUp2*qUpNorm; % x1
                qUp3 = qUp3*qUpNorm; % x1
                qUp4 = 0; 
                % update accumulating error
            end
            % Rotate global frame towards 'up'
            qGii1 = qUp1*qGi1 - qUp2*qGi2 - qUp3*qGi3 ; % sub 2, x3 
            qGii2 = qUp1*qGi2 + qUp2*qGi1 + qUp3*qGi4 ; % add 2, x3
            qGii3 = qUp1*qGi3 - qUp2*qGi4 + qUp3*qGi1 ; % add 1, sub 1, x3
            qGii4 = qUp1*qGi4 + qUp2*qGi3 - qUp3*qGi2 ; % add 1, sub 1, x3
% -----------------------------------------------------------------------            
            % estimate external acceleration in sensor frame
            z1 = 2*(qGii2*qGii4 - qGii1*qGii3);     % sub 1, x3 
            z2 = 2*(qGii3*qGii4 + qGii1*qGii2);     % add 1, x3
            z3 = 2*(qGii1*qGii1 + qGii4*qGii4) - 1; % add 1, sub 1, x3
            upSFR    = [z1 z2 z3];
            obj.ext_acc = yAcc - obj.gRef * upSFR; % sub 3, x3   
        else
            qGii1 = qGi1;  qGii2 = qGi2;  qGii3 = qGi3;  qGii4 = qGi4;
        end
% -----------------------------------------------------------------------                    
        % Correct for magnetometer - get rotation around vertical to align measured global xy 
        % model magnetic field strength as a first order low pass filtered
        % process, i.e., zMag = m_k - cM * mag_dis;
        yMag = Mag;
        obj.mag_now = yMag;
        zMag = yMag - obj.cM(1) * obj.mag_dis;  % sub 3, x3
        Mag  = zMag;
        Mag1 = Mag(1);    
        Mag2 = Mag(2);    
        Mag3 = Mag(3);
        % Transform magnetometer reading into global frame
        if ( ~isnan(Mag1) && ~isnan(Mag2) && ~isnan(Mag3) ) 
            qm1  = -qGii2*Mag1 - qGii3*Mag2 - qGii4*Mag3; % sub 3, x3
            qm2  =  qGii1*Mag1 + qGii3*Mag3 - qGii4*Mag2; % add 1, sub 1, x3
            qm3  =  qGii1*Mag2 - qGii2*Mag3 + qGii4*Mag1; % add 1, sub 1, x3
            qm4  =  qGii1*Mag3 + qGii2*Mag2 - qGii3*Mag1; % add 1, sub 1, x3
            gMag1 = -qm1*qGii2 + qm2*qGii1 - qm3*qGii4 + qm4*qGii3; % add 2, sub 2, x4 
            gMag2 = -qm1*qGii3 + qm2*qGii4 + qm3*qGii1 - qm4*qGii2; % add 2, sub 2, x4
            obj.mag_now_xy = [gMag1 gMag2];
            obj.gMag       = [gMag1 gMag2 0 ];
            
% ---------------------------------------------------------------------------                
% ------------ MAGNETOMETER MEASUREMENT VALIDATION ---------------------------
% ---  change in magnitude of the xy component of the estimated global magnetic field   
            delta_mag_b_XY = abs(... % add 2, sub 1, mult 4
                 (obj.mag_pre_xy(1)*obj.mag_pre_xy(1) + ...
                  obj.mag_pre_xy(2)*obj.mag_pre_xy(2)) - ...
                 (obj.mag_now_xy(1)*obj.mag_now_xy(1) + ...
                  obj.mag_now_xy(2)*obj.mag_now_xy(2)) ) * obj.fs; % do not count * obj.fs, instead obj.delta_b_xy /= obj.fs prior to runtime

            obj.deltaMagBxy = delta_mag_b_XY;              
            obj.mag_pre_xy = obj.mag_now_xy; % set as previous for next time
            obj.mag_pre = yMag;
            if delta_mag_b_XY < obj.delta_b_xy
                % Get fraction of rotation from mag in global to north (xy components only)
                qn2 = 0; qn3 = 0;
                if gMag2 == 0 % y-component is 0 therefore pointing north
                    qn1 = 1;  qn4 = 0;
                else
                    nGamM      = arctan2_approx(abs(gMag2),gMag1);
                    nGamMSq    = nGamM * nGamM; % x1
                    obj.gammaM = nGamM;
% -------------------------------------------------------------------------                
                    % Estimate measurement error, gamma_m using a sliding window
                    % Run time calculation of var(gamma_m) used to estimate R_GamM
                    obj.ctrM = obj.ctrM + 1;                 % add 1
                    mCtr   = obj.ctrM;
                    mIdx   = mod(mCtr-1, obj.N_VAR_GAM_M)+1; % add 1, sub 1
                    oGamM = obj.BufGamM(mIdx);              % keep track of oldest value in the BufGamA
                    obj.BufGamM(mIdx) = nGamM;   % replace oldest idx in BufGamM with newest value
%--------------------------------------------                    
                    obj.rSumGamM   = obj.rSumGamM   + nGamM;    % add 1
                    obj.rSumSqGamM = obj.rSumSqGamM + nGamMSq;  % add 1
                    if mCtr > obj.N_VAR_GAM_M
                        obj.rSumGamM   = obj.rSumGamM   - oGamM;         % sub 1
                        obj.rSumSqGamM = obj.rSumSqGamM - oGamM * oGamM; % sub 1, x1
                    end
                    nAvgCGamMSq = ...
                        obj.nAvgConstGamM * obj.rSumGamM * obj.rSumGamM ;   % x2
                    obj.R_GamM = ...
                        obj.VarConstGamM * (obj.rSumSqGamM + nAvgCGamMSq ); % add 1, x1  
% -----------------------------------------------------                                
                    obj.Pmag_pri = obj.Pmag_pos + obj.Qmag;   % add 1
                    Pmag_pri_R_GamM = obj.Pmag_pri + obj.R_GamM; % add 1
                    S_m_inv = 1/sqrt( Pmag_pri_R_GamM * Pmag_pri_R_GamM ); % x1, invsqrt(1) 
                    mu_m    = obj.Pmag_pri * S_m_inv; % x1
                    if mCtr < obj.N_VAR_GAM_ML, mu_m = 0.001; end 
                    obj.muM      = mu_m;                
                    obj.Pmag_pos = (1 - mu_m) * obj.Pmag_pri; % sub 1, x1
                    muMgamM = mu_m * nGamM * 0.5; % x2
                    qn1 = 1;
                    qn4 = java.lang.Math.signum(-gMag2)*muMgamM; % sub 1, x1
                    qNorm = 1/sqrt(1 + qn4 * qn4);    % add 1, x1, invsqrt(1)
                    qn1 = qNorm; 
                    qn4 = qn4*qNorm;                             % x1
                end
                % Rotate global frame towards 'north'        
                qGiii1 = qn1*qGii1 - qn4*qGii4; % sub 1, x2
                qGiii2 = qn1*qGii2 - qn4*qGii3; % sub 1, x2
                qGiii3 = qn1*qGii3 + qn4*qGii2; % add 1, x2
                qGiii4 = qn1*qGii4 + qn4*qGii1; % add 1, x2
            else
                obj.muM = NaN; % purely for visualization
                qGiii1 = qGii1; 
                qGiii2 = qGii2; 
                qGiii3 = qGii3; 
                qGiii4 = qGii4;
            end
            % -----------------------------------------------------------------------            
            % estimate magnetometer disturbance in sensor frame
            x1 = 2*(qGiii1*qGiii1 + qGiii2*qGiii2) -1; % add 1, sub 1, x3            
            x2 = 2*(qGiii2*qGiii3 - qGiii1*qGiii4);    % sub 1, x3
            x3 = 2*(qGiii2*qGiii4 + qGiii1*qGiii3);    % add 1, x3
            northSFR = [x1 x2 x3];
            obj.mag_dis = yMag - obj.mCosTheta*northSFR ...
                - obj.mSinTheta*upSFR; % sub 6, x6
        else
            qGiii1 = qGii1; 
            qGiii2 = qGii2; 
            qGiii3 = qGii3; 
            qGiii4 = qGii4;
        end
        obj.qGlobal = [qGiii1 qGiii2 qGiii3 qGiii4];        
    end
    
end 

end


function [atan2x]=arctan2_approx(y,x)
    % only works when y is positive i.e. the numerator in atan2(y,x)
    if x >= 0
        if x>=y
            invsqrtxsq = java.lang.Math.signum(x)*(1/sqrt(x*x)); % can be replaced by the fast inverse square root
            %atan2x = arctan_approx(y/x);
            atan2x = arctan_approx(y*invsqrtxsq);
        else % x < y
            invsqrtysq = java.lang.Math.signum(y)*(1/sqrt(y*y)); % can be replaced by the fast inverse square root
            %atan2x = pi/2 minus arctan_approx(x/y);
            atan2x = pi/2 - arctan_approx(x*invsqrtysq);
        end
    else % x<=0
        if y>abs(x)
            %atan2x = pi/2 plus arctan_approx(abs(x)/y);
            invsqrtysq = java.lang.Math.signum(y)*(1/sqrt(y*y)); % can be replaced by the fast inverse square root
            atan2x = pi/2 + arctan_approx(abs(x)*invsqrtysq);
        else 
            %atan2x = pi minus arctan_approx(y/abs(x));                        
            invsqrtxsq = 1/sqrt(x*x); % can be replaced by the fast inverse square root
            atan2x = pi - arctan_approx(y*invsqrtxsq);
        end
    end
end

% add 1, sub 2, x4
function atanx = arctan_approx(x)
    qtr_pi = pi/4; % predefined constant
    atanx = (qtr_pi*x) - x*(abs(x)-1) * (0.2447 + 0.0663 * abs(x));
end


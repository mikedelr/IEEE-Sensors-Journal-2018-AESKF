function [ varargout ] = wrapper_MIMU_AESKF( varargin )
%WRAPPER_MIMU_AESKF Wrapper function designed to process a dataset of
%accelerometer, gyroscope and magnetometer through AESKF
%   [ quatsMIMU ] = ...
%   wrapper_MIMU_AESKF( ... )
%
%   Inputs:: (must be 'string', value pairs, e.g., 'fs',100 .... 
%             e.g., sampling rate of 100 Hz)
%       'q0'          - initial orientation quaternion
%       'fs'          - sampling rate
%       'Acc'         - accelerometer signal size = [N x 3]
%       'Gyr'         - gyroscope     signal size = [N x 3]
%       'Mag'         - magnetometer  signal size = [N x 3]
%       'N_VAR_GAM_AS' - number of samples over which to calculate the 
%                       short-term variance of the angle between the up vector and the 
%                       accelerometer measurement vector expressed in the
%                       estimated GFR. Only used when device appears
%                       stationary
%       'N_VAR_GAM_AL' - number of samples over which to calculate the 
%                       variance of the angle between the up vector and the 
%                       accelerometer measurement vector expressed in the
%                       estimated GFR. Used when device is appears in
%                       motion
%       'N_VAR_GAM_M' - number of samples over which to calculate the 
%                       variance of the angle between the north vector and
%                       the magnetometer measurement vector projected to
%                       the horizontal plane of the estimated GFR
%       'sigma_gyr'   - standard deviation in the gyro when stationary
%   Outputs::
%       quatsMIMU - quaternions from a 9DOF MIMU
%       R_gamma_a - adaptive measurement noise based on the variance of the
%                   angle between the up vector and the accelerometer 
%                   measurement vector expressed in the estimated GFR over
%                   a sliding window of size 'N_VAR_GAM_A'
%       gammaA - angle between the up vector and the accelerometer 
%                measurement vector expressed in the estimated GFR
%       Pk_pri - a priori error
%       Pk_pos - a posteriori error
%       Kk     - Kalman gain


q0          = [];
N_VAR_GAM_M = NaN;
sigma_gyr   = NaN;
% sigma_acc   = NaN;
for i = 1:2:nargin
    if strcmp(varargin{i}, 'q0'),
        q0 = varargin{i+1};      
    elseif strcmp(varargin{i},'fs'), 
        fs = varargin{i+1};   
    elseif strcmp(varargin{i},'Acc')
        ACC = varargin{i+1};
        [NUM_ACC,COLS_A]=size(ACC);
        if COLS_A ~= 3
            error('Accelerometer data must be size [N x 3]');
        end        
    elseif strcmp(varargin{i},'Gyr')
        GYRO = varargin{i+1};
        [NUM_GYR,COLS_G]=size(GYRO);
        if COLS_G ~= 3
            error('Gyroscope data must be size [N x 3]');
        end
    elseif strcmp(varargin{i},'Mag')
        MAGNO = varargin{i+1};        
        [NUM_MAG,COLS_M]=size(MAGNO);
        if COLS_M ~= 3
            error('Magnetometer data must be size [N x 3]');
        end
    elseif strcmp(varargin{i},'N_VAR_GAM_AS');
        N_VAR_GAM_AS = varargin{i+1};
    elseif strcmp(varargin{i},'N_VAR_GAM_AL');
        N_VAR_GAM_AL = varargin{i+1};              
    elseif strcmp(varargin{i},'N_VAR_GAM_M');
        N_VAR_GAM_M = varargin{i+1};        
    elseif strcmp(varargin{i},'sigma_gyr');
        sigma_gyr = varargin{i+1};           
    elseif strcmp(varargin{i},'cA');
        cA = varargin{i+1};
    elseif strcmp(varargin{i},'cM');
        cM = varargin{i+1};       
    elseif strcmp(varargin{i},'gRef');
        gRef = varargin{i+1};
    elseif strcmp(varargin{i},'mRef');
        mRef = varargin{i+1};      
    elseif  strcmp(varargin{i},'delta_b_xy')
        delta_b_xy = varargin{i+1};
    elseif  strcmp(varargin{i},'incAngRef')
        incAngRef = varargin{i+1};        
    else error('Invalid argument');
    end    
end

if (NUM_ACC ~= NUM_GYR) && (NUM_ACC ~= NUM_MAG) && (NUM_GYR ~= NUM_MAG)
    error('Number of samples in IMU data are unequal');
end

ACC_VAR_THRESH = 9; % WIN = 1*fs
WIN_SECS = 0.25;
VAR_WIN  = floor(fs*WIN_SECS);

% define storage before use for speed
FRAMES    = NUM_ACC;
quatsMIMU = zeros(FRAMES,4);

gam_a     = nan(FRAMES,1);
avg_gam_a = nan(FRAMES,1);
Pa_pri    = zeros(FRAMES,1);
Pa_pos    = zeros(FRAMES,1);
K_a       = zeros(FRAMES,1);
R_gam_a   = zeros(FRAMES,1);
R_gam_as  = zeros(FRAMES,1);
R_gam_al  = zeros(FRAMES,1);
ext_acc   = zeros(FRAMES,3);
a_mag_var = zeros(FRAMES,1);


gam_m     = nan(FRAMES,1);
avg_gam_m = nan(FRAMES,1);
Pm_pri    = zeros(FRAMES,1);
Pm_pos    = zeros(FRAMES,1);
K_m       = zeros(FRAMES,1);
R_gam_m   = zeros(FRAMES,1);
R_gam_ms  = zeros(FRAMES,1);
R_gam_ml  = zeros(FRAMES,1);
mag_dis   = zeros(FRAMES,3);

gAcc   = zeros(FRAMES,3);
gMag   = zeros(FRAMES,3);
deltaMagXY = zeros(FRAMES,1);

if ~exist('incAngRef','var')
    incAngRef = NaN;
end

aeskf_Obj = AESKF('fs',fs,...
                  'q0',q0,...
                  'sigma_gyr',sigma_gyr,...
                  'cA',cA,...
                  'cM',cM,...
                  'gRef',gRef,...
                  'mRef',mRef ,...
                  'delta_b_xy',delta_b_xy,...
                  'N_VAR_GAM_AS',N_VAR_GAM_AS,...
                  'N_VAR_GAM_AL',N_VAR_GAM_AL,...
                  'N_VAR_GAM_M',N_VAR_GAM_M,...
                  'incAngRef',incAngRef);              
              
aeskf_Obj.gMagRef = quatrotate(q0, MAGNO(1,:));
              
for frame = 1:FRAMES 
    aeskf_Obj.Update(GYRO(frame,:), ACC(frame,:),MAGNO(frame,:));
    quatsMIMU(frame, :) = aeskf_Obj.qGlobal;
    % upward angle components
    gam_a(frame, :)     = aeskf_Obj.gammaA;
    avg_gam_a(frame, :) = aeskf_Obj.AvgGamA;
    Pa_pri(frame, :)    = aeskf_Obj.P_pri;
    Pa_pos(frame, :)    = aeskf_Obj.P_pos;
    K_a(frame, :)       = aeskf_Obj.muA;
    R_gam_a(frame,:)    = aeskf_Obj.R_GamA;
    R_gam_as(frame,:)    = aeskf_Obj.R_GamAS;
    R_gam_al(frame,:)    = aeskf_Obj.R_GamAL;
    ext_acc(frame,:)    = aeskf_Obj.ext_acc;
    % northward angle components
    gam_m(frame, :)     = aeskf_Obj.gammaM;
    avg_gam_m(frame, :) = aeskf_Obj.AvgGamM;
    Pm_pri(frame, :)    = aeskf_Obj.Pmag_pri;
    Pm_pos(frame, :)    = aeskf_Obj.Pmag_pos;
    K_m(frame, :)       = aeskf_Obj.muM;
    R_gam_m(frame,:)    = aeskf_Obj.R_GamM; 
    R_gam_ms(frame,:)   = aeskf_Obj.R_GamMS;
    R_gam_ml(frame,:)   = aeskf_Obj.R_GamML;    
    mag_dis(frame,:)    = aeskf_Obj.mag_dis;
    gAcc(frame,:)       = aeskf_Obj.gAcc; 
    gMag(frame,:)       = aeskf_Obj.gMag;
    deltaMagXY(frame,:) = aeskf_Obj.deltaMagBxy;
end

% variable output assignment
if nargout >= 1, varargout{1} = quatsMIMU; end
if nargout >= 2, varargout{2} = R_gam_a;   end    
if nargout >= 3, varargout{3} = avg_gam_a; end    
if nargout >= 4, varargout{4} = gam_a;     end
if nargout >= 5, varargout{5} = Pa_pri;    end    
if nargout >= 6, varargout{6} = Pa_pos;    end
if nargout >= 7, varargout{7} = K_a;       end

if nargout >= 8,  varargout{8} = R_gam_m;   end    
if nargout >= 9,  varargout{9} = avg_gam_m; end    
if nargout >= 10, varargout{10} = gam_m;     end
if nargout >= 11, varargout{11} = Pm_pri;    end    
if nargout >= 12, varargout{12} = Pm_pos;    end
if nargout >= 13, varargout{13} = K_m;       end

if nargout >= 14, varargout{14} = ext_acc;   end
if nargout >= 15, varargout{15} = mag_dis;   end
if nargout >= 16, varargout{16} = gAcc;   end
if nargout >= 17, varargout{17} = gMag;   end
if nargout >= 18, varargout{18} = deltaMagXY;   end
if nargout >= 19, varargout{19} = R_gam_as; end
if nargout >= 19, varargout{20} = R_gam_al; end

end





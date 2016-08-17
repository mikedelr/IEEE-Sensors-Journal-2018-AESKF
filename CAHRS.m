classdef CAHRS <handle
%CAHRS is the AHRS algorithm for a 9DOF and or 6DOF MIMU/IMU 
% See 'Quaternion based Complementary Filter for Attitude Determination of a
% Smartphone'
% information from each of the sensors has been decoupled and fused
% according to the tuning parameters specified by the user. The
% algorithm assumes that all sensors inputs have been calibrated
%
% Author(s): Stephen Redmond & Michael Del Rosario (m.delrosario@unsw.edu.au)

properties (Access = public)
    fs          = 100;               % additional constants defined for speed
    dt          = 0.01;
    qAngVelDev  = [1 0 0 0];
    qGlobal     = [1 0 0 0];
    qGlobalPrev = [1 0 0 0];        
    muAcc = 0.005;
    muMag = 0.005;
    qGlobalGyr = NaN;
    AccGlobal = NaN;
    qRotTowardsUp = NaN;
    MagGlobal = NaN;
    qRotTowardsNorth = NaN;
    qGlobalPrevReal  = NaN;
    qGlobalRotatedUp = NaN;
    
    % these are for storing the angles gamma_a and gamma_m
    gammaA         = NaN;
    gammaM         = NaN;    
    
    % store global acceleration vector
    gAcc = [];
end

methods (Access = public)
    %% default constructor
    function obj = CAHRS(varargin)
        for i = 1:2:nargin
            if strcmp(varargin{i},'fs'), 
                obj.fs = varargin{i+1};
                obj.dt=1/varargin{i+1};
            elseif strcmp(varargin{i},'muAcc'),
                obj.muAcc = varargin{i+1};
            elseif strcmp(varargin{i},'muMag'), 
                obj.muMag = varargin{i+1};
            elseif strcmp(varargin{i},'qInit')
                obj.qGlobalPrev = varargin{i+1};
            else error('Invalid argument');
            end
        end;
        obj.muAcc  = obj.muAcc*obj.dt; 
        obj.muMag  = obj.muMag*obj.dt;
    end        

    %% UpdateNoCallsArcTan
    % 128 multiplications, 43 additions, 43 subtractions, 6 calls to
    % invSqrt (4 multiplications, 2 subtractions)
    % 2 calls to arctan_approx (multiply (x4), addition (x1), subtraction (x2))
    function obj=UpdateNoCallsArcTan(obj,Gyr,Acc,Mag)
        Acc1=Acc(1);    Acc2=Acc(2);    Acc3=Acc(3);
        Gyr1=Gyr(1);    Gyr2=Gyr(2);    Gyr3=Gyr(3);
        Mag1=Mag(1);    Mag2=Mag(2);    Mag3=Mag(3);
        
        muA2 = obj.muAcc/2; 
        muM2 = obj.muMag/2; 
        qG1 = obj.qGlobal(1); 
        qG2 = obj.qGlobal(2); 
        qG3 = obj.qGlobal(3); 
        qG4 = obj.qGlobal(4);            
        %% Correct for gyros - get quaternion from gyro in device frame and rotate
        if ~isnan(Gyr)
            dt2 = 0.5*obj.dt;
            qW1 = 1; qW2 = Gyr1*dt2; qW3 = Gyr2*dt2; qW4 = Gyr3*dt2;
            qWnorm = 1/java.lang.Math.sqrt(qW1*qW1 + qW2*qW2 + qW3*qW3 + qW4*qW4); % can be replaced by the fast inverse square root
            qW1 = qW1*qWnorm; qW2 = qW2*qWnorm; qW3 = qW3*qWnorm; qW4 = qW4*qWnorm;  
            obj.qAngVelDev = [qW1 qW2 qW3 qW4];
            % Convert to back to global frame
            qGi1 = qG1*qW1 - qG2*qW2 - qG3*qW3 - qG4*qW4;
            qGi2 = qG1*qW2 + qG2*qW1 + qG3*qW4 - qG4*qW3;
            qGi3 = qG1*qW3 - qG2*qW4 + qG3*qW1 + qG4*qW2;
            qGi4 = qG1*qW4 + qG2*qW3 - qG3*qW2 + qG4*qW1; 
        else
            qGi1 = qG1; qGi2 = qG2; qGi3 = qG3; qGi4 = qG4; 
        end
        %% Correct for accelerometer - get next rotation required to align 
        % accelerometer/gravity with 'up'.
        if ~isnan(Acc)
            % intermediate calculation in calculating v' = qvq*
            qa1  = -qGi2*Acc1 - qGi3*Acc2 - qGi4*Acc3;
            qa2  =  qGi1*Acc1 + qGi3*Acc3 - qGi4*Acc2;
            qa3  =  qGi1*Acc2 - qGi2*Acc3 + qGi4*Acc1;
            qa4  =  qGi1*Acc3 + qGi2*Acc2 - qGi3*Acc1;  
            % Convert acceleration to global frame
            gAcc1 = -qa1*qGi2 + qa2*qGi1 - qa3*qGi4 + qa4*qGi3;
            gAcc2 = -qa1*qGi3 + qa2*qGi4 + qa3*qGi1 - qa4*qGi2;
            gAcc3 = -qa1*qGi4 - qa2*qGi3 + qa3*qGi2 + qa4*qGi1;
            obj.gAcc = [gAcc1,gAcc2,gAcc3];
            % Get fraction of rotation from acc in global to up
            if muA2<0,    
                muA2 = 0;    
                warning('Capping mu at 0');
            elseif muA2>0.5,
                muA2 = 0.5;    
                warning('Capping mu at 0.5');
            end
            gAcc2Sq = gAcc2*gAcc2;
            gAcc1Sq = gAcc1*gAcc1;
            gAcc12SumSq = gAcc1Sq+gAcc2Sq;
            axisVecNorm = 1/java.lang.Math.sqrt((gAcc12SumSq)); % can be replaced by the fast inverse square root
            if axisVecNorm == 0
                qUp1 = 1; qUp2 = 0; qUp3 = 0; qUp4 = 0;
            else
                gamma_a        = arctan2_approx(axisVecNorm*gAcc12SumSq,gAcc3);
                obj.gammaA     = gamma_a;
                mu_a_gamma_a   = muA2*gamma_a;
                
                cos_mu_a_gamma_a = 1;            % small angle approximation
                sin_mu_a_gamma_a = mu_a_gamma_a; % small angle approximation
                obj.cos_muA_gammaA = cos_mu_a_gamma_a;
                obj.sin_muA_gammaA = sin_mu_a_gamma_a;
                
                qUp1    = cos_mu_a_gamma_a;
                qUp2    =  gAcc2 * axisVecNorm * sin_mu_a_gamma_a;
                qUp3    = -gAcc1 * axisVecNorm * sin_mu_a_gamma_a; 
                qUpNorm = 1/java.lang.Math.sqrt((qUp1*qUp1 + qUp2*qUp2 + qUp3*qUp3)); % can be replaced by the fast inverse square root
                % normalise to unit quaternions
                qUp1 = qUp1*qUpNorm; qUp2 = qUp2*qUpNorm; qUp3 = qUp3*qUpNorm; qUp4 = 0;
            end
            % Rotate global frame towards 'up'
            qGii1 = qUp1*qGi1 - qUp2*qGi2 - qUp3*qGi3 - qUp4*qGi4;
            qGii2 = qUp1*qGi2 + qUp2*qGi1 + qUp3*qGi4 - qUp4*qGi3;
            qGii3 = qUp1*qGi3 - qUp2*qGi4 + qUp3*qGi1 + qUp4*qGi2;
            qGii4 = qUp1*qGi4 + qUp2*qGi3 - qUp3*qGi2 + qUp4*qGi1;
        else
            qGii1 = qGi1;  qGii2 = qGi2;  qGii3 = qGi3;  qGii4 = qGi4;
        end
        %% Correct for magnetometer - get next rotation around vertical to align measured global xy 
        % Transform magnetometer reading into global frame
        if ~isnan(Mag)
            qm1  = -qGii2*Mag1 - qGii3*Mag2 - qGii4*Mag3;
            qm2  =  qGii1*Mag1 + qGii3*Mag3 - qGii4*Mag2;
            qm3  =  qGii1*Mag2 - qGii2*Mag3 + qGii4*Mag1;
            qm4  =  qGii1*Mag3 + qGii2*Mag2 - qGii3*Mag1;  
            gMag1 = -qm1*qGii2 + qm2*qGii1 - qm3*qGii4 + qm4*qGii3;
            gMag2 = -qm1*qGii3 + qm2*qGii4 + qm3*qGii1 - qm4*qGii2;
            % Get fraction of rotation from mag in global to north (xy components only)
            if muM2<0,    
                muM2 = 0;    warning('Capping mu at 0');
            elseif muM2>0.5,
                muM2 = 0.5;  warning('Capping mu at 0.5');
            end
            qn2 = 0; qn3 = 0;
            if gMag2 == 0 % y-component is 0 therefore pointing north
                qn1 = 1;  qn4 = 0;
            else
                gamma_m      = arctan2_approx(abs(gMag2),gMag1); % absolute value unsure if always correct
                obj.gammaM   = gamma_m;
                mu_m_gamma_m = muM2*gamma_m;

                cos_mu_m_gamma_m = 1;
                sin_mu_m_gamma_m = mu_m_gamma_m;
                
                obj.cos_muM_gammaM = cos_mu_m_gamma_m;
                obj.sin_muM_gammaM = sin_mu_m_gamma_m;
                
                qn1 = cos_mu_m_gamma_m;
                qn4 = java.lang.Math.signum(-gMag2)*sin_mu_m_gamma_m;
                qNorm = 1/java.lang.Math.sqrt(qn1*qn1+qn4*qn4);  % can be replaced by the fast inverse square root
                qn1 = qn1*qNorm;
                qn4 = qn4*qNorm; % normalise to unit quaternions
            end
            % Rotate global frame towards 'up'        
            qGiii1 = qn1*qGii1 - qn2*qGii2 - qn3*qGii3 - qn4*qGii4;
            qGiii2 = qn1*qGii2 + qn2*qGii1 + qn3*qGii4 - qn4*qGii3;
            qGiii3 = qn1*qGii3 - qn2*qGii4 + qn3*qGii1 + qn4*qGii2;
            qGiii4 = qn1*qGii4 + qn2*qGii3 - qn3*qGii2 + qn4*qGii1;
        else
            qGiii1 = qGii1; qGiii2 = qGii2; qGiii3 = qGii3; qGiii4 = qGii4;
        end
        obj.qGlobal = [qGiii1 qGiii2 qGiii3 qGiii4];
    end % END UpdateNoCallsArcTan
 
end

end

function [atan2x]=arctan2_approx(y,x)
    % only works when y is positive i.e. the numerator in atan2(y,x)
    if x >= 0
        if x>=y
            invsqrtxsq = java.lang.Math.signum(x)*(1/java.lang.Math.sqrt(x*x)); % can be replaced by the fast inverse square root
            %atan2x = arctan_approx(y/x);
            atan2x = arctan_approx(y*invsqrtxsq);
        else % x < y
            invsqrtysq = java.lang.Math.signum(y)*(1/java.lang.Math.sqrt(y*y)); % can be replaced by the fast inverse square root
            %atan2x = pi/2 - arctan_approx(x/y);
            atan2x = pi/2 - arctan_approx(x*invsqrtysq);
        end
    else % x<=0
        if y>abs(x)
            %atan2x = pi/2 + arctan_approx(abs(x)/y);
            invsqrtysq = java.lang.Math.signum(y)*(1/java.lang.Math.sqrt(y*y)); % can be replaced by the fast inverse square root
            atan2x = pi/2 + arctan_approx(abs(x)*invsqrtysq);
        else 
            %atan2x = pi - arctan_approx(y/abs(x));                        
            invsqrtxsq = 1/java.lang.Math.sqrt(x*x); % can be replaced by the fast inverse square root
            atan2x = pi - arctan_approx(y*invsqrtxsq);
        end
    end
end

function [atanx]=arctan_approx(x)
qtr_pi = pi/4;
% multiply (x4)
% addition (x1)
% subtraction (x2)
atanx = (qtr_pi.*x) - x.*(abs(x)-1) .* (0.2447 + 0.0663 .* abs(x));
end


% 1) synchronise the imu and robot data streams,
% 2) express robot orientation in GFR using the IMU data during static periods
% to derive a quaternion
% 3) simulate the accelerometer data use the real gyroscope data
% 4) validate that the euler angles derived from the quaternions are
% equivalent
function [] = main()
clear all;clc;
addpath('.\CAHRS'); % path to CAHRS functions and class files
addpath('.\AESKF_functions'); % path to AESKF functions and class files
addpath('.\quaternion_functions');
addpath('.\Data');
rng(1);
sigma_gyr = 0.01;
gValue = 9.796720;

cA = 0.1;
cM = 0.99;

N_VAR_GAM_AS = 8;
N_VAR_GAM_AL = 50;

N_VAR_GAM_M = 8;
fsMIMU      = 100;
delta_mag_b_xy = (2.5*10^-3) * fsMIMU;

%%
% Good expDays, 1,3, 5,6
% EXP_DAYS =[1 3 5 ];
% EXP_DAYS = 1;

% load('Up-Down.mat')
% load('Left-Right.mat')
load('ZigZag.mat')

close all;

%%
ACC_VAR_THRESH = 12; % WIN = 1*fs
muAcc = 0.003;
muMag = 0.001;
acc_cal = table2array(imu_real_sim_data(:,{'ImuAccX','ImuAccY','ImuAccZ'}));
gyr_cal = table2array(imu_real_sim_data(:,{'ImuGyrX','ImuGyrY','ImuGyrZ'}));
mag_cal = NaN(size(acc_cal));
WIN_SECS = 0.25;
VAR_WIN  = floor(fsMIMU*WIN_SECS); % NUM_SAMPLES

MOV_VAR_ACC_2 = movingvar(sum(acc_cal.^2,2),VAR_WIN);
bIsDynamic = MOV_VAR_ACC_2 > ACC_VAR_THRESH;
bIsDynamic(1:VAR_WIN-1) = true; % assume moving during start-up transient of variance calculation

figure;
subplot(2,1,1);plot([acc_cal,bIsDynamic]);
subplot(2,1,2);plot(MOV_VAR_ACC_2);

q_mea = table2array(imu_real_sim_data(:,{'qw','qx','qy','qz'}));

%% --- Real measurements from the smartphone
% CAHRS, fixed gain, 
std_gyr = mean(std(acc_cal));
std_acc = mean(std(gyr_cal));
[quats_cahrs] = ...
opMimuDynamicCAHRS_ArcTan('qInit',q_mea(1,:),...
'bDynamicMu',bIsDynamic,...
'dynamic_mu_acc',muAcc,'dynamic_mu_mag',muMag,...
'static_mu_acc',0.75,'static_mu_mag',0,...
'fs',fsMIMU,'Acc',acc_cal,'Mag',mag_cal,'Gyr',gyr_cal);


N_TIME = length(imu_real_sim_data.timestamps); 

i_static = (~bIsDynamic) & ([1:N_TIME]' < 1000);
mRef = 1;
gRef = mean(vecnormalize(acc_cal(i_static,:)));

[quats_aeskf, R_gamma_A_IKF, avg_gamma_A_IKF, ...
 gamma_A_IKF, Pk_pri_IKF, Pk_pos_IKF, K_gain_IKF,...
 R_gamma_M_IKF, avg_gamma_M_IKF, ...
 gamma_M_IKF, Pm_pri, Pm_pos, Km, ext_acc,mag_dis, gAcc,gMag] = ...
    wrapper_MIMU_AESKF(... 
        'q0',q_mea(1,:),'fs',fsMIMU,'N_VAR_GAM_M',N_VAR_GAM_M,...
        'N_VAR_GAM_AS',N_VAR_GAM_AS,'N_VAR_GAM_AL',N_VAR_GAM_AL,...
        'sigma_gyr',sigma_gyr,'cA',cA,'gRef',gRef,'cM',cM,'mRef',mRef,...
        'Acc',acc_cal,'Mag',mag_cal,'Gyr',gyr_cal,...
        'delta_b_xy',delta_mag_b_xy*mRef*mRef); 

est_eul_cahrs = rad2deg(quaternion2nautical(quats_cahrs));  % ahrs estimate 
est_eul_aeskf = rad2deg(quaternion2nautical(quats_aeskf)); % ahrs estimate2
tru_eul_mea   = rad2deg(quaternion2nautical(q_mea)); % truth measurement
        
tru_eul_mea = tru_eul_mea;
real_eul_fix_sync = est_eul_cahrs;
real_eul_IKF_sync = est_eul_aeskf;

for e = 1:3
    eulerAngErr = tru_eul_mea(:,e)-real_eul_fix_sync(:,e); % order not important since error is abs
    greater180 = eulerAngErr > 180;  
    real_eul_fix_sync(greater180,e) = real_eul_fix_sync(greater180,e) + 360;
    % if difference is less than -180, add 360
    lesserN180 = eulerAngErr < -180;
    real_eul_fix_sync(lesserN180,e) = real_eul_fix_sync(lesserN180,e) - 360;
    eulerAngErr2 = tru_eul_mea(:,e)-real_eul_IKF_sync(:,e); % order not important since error is abs
    greater180 = eulerAngErr2 > 180;  
    real_eul_IKF_sync(greater180,e) = real_eul_IKF_sync(greater180,e) + 360;
    % if difference is less than -180, add 360
    lesserN180 = eulerAngErr2 < -180;
    real_eul_IKF_sync(lesserN180,e) = real_eul_IKF_sync(lesserN180,e) - 360;
end
%%
updateFigureContents(['Real Phone Data Acc/Gyr CAHRS vs CAHRS_IKF',nexusFile]);hsub=[];
hsub(1)=subplot(4,1,1);hold on;lgstr={};h=[];hz=zoom;set(hz,'Enable','on')
title('Roll');
ptime = 1:length(phonedata.time(estIdx));%phonedata.time;
h(end+1)=plot(ptime,tru_eul_mea(:,1),'k','LineWidth',2);lgstr{end+1}='REF';
h(end+1)=plot(ptime,real_eul_fix_sync(:,1),'oc','LineStyle','none');lgstr{end+1}='CAHRS';
h(end+1)=plot(ptime,real_eul_IKF_sync(:,1),'xb','LineStyle','none');lgstr{end+1}='CAHRS-IKF';
legend(h,lgstr,'Orientation','Horizontal');

sync_segs = reshape(phonedata.syncSegments',[],1);
stat_segs = reshape(phonedata.statSegments',[],1);
% adjust stationary period for convergence
sync_segs(2:end-1) = sync_segs(2:end-1);
% adjust for algorithm convergence
sync_segs(1) = sync_segs(1);
% adjust for alignment between reference systems
sync_segs(end) = sync_segs(end)+ t21; 
stat_segs(end) = stat_segs(end)+ t21;
sub_sync_segs = [sync_segs(1:6:end);sync_segs(6:6:end)];
sub_stat_segs = [stat_segs(1:6:end);stat_segs(6:6:end)];
mid_segs = floor(mean([sync_segs(1:6:end),sync_segs(6:6:end)] +1,2));
line(ptime(repmat(sub_sync_segs +1,1,2))',repmat(ylim,18,1)','Color','k');

hsub(2)=subplot(4,1,2);hold on;lgstr={};h=[];
title('Pitch');
h(end+1)=plot(ptime,tru_eul_mea(:,2),'k','LineWidth',2);lgstr{end+1}='REF';
h(end+1)=plot(ptime,real_eul_fix_sync(:,2),'om','LineStyle','none');lgstr{end+1}='CAHRS';
h(end+1)=plot(ptime,real_eul_IKF_sync(:,2),'xr','LineStyle','none');lgstr{end+1}='CAHRS-IKF';
legend(h,lgstr,'Orientation','Horizontal');

hsub(3)=subplot(4,1,3);hold on;lgstr={};h=[];
title('Error');
h(end+1)=plot(ptime,tru_eul_mea(:,1)-real_eul_fix_sync(:,1),'--c');lgstr{end+1}='\phi_{CAHRS}';
h(end+1)=plot(ptime,tru_eul_mea(:,1)-real_eul_IKF_sync(:,1),'--b');lgstr{end+1}='\phi_{CAHRS-IKF}';
h(end+1)=plot(ptime,tru_eul_mea(:,2)-real_eul_fix_sync(:,2),'--m');lgstr{end+1}='\theta_{CAHRS}';
h(end+1)=plot(ptime,tru_eul_mea(:,2)-real_eul_IKF_sync(:,2),'--r');lgstr{end+1}='\theta_{CAHRS-IKF}';        
legend(h,lgstr,'Orientation','Horizontal');

hsub(4)=subplot(4,1,4);hold on;lgstr={};h=[];
title('Error State Kalman Filter');
h(end+1)=plot(ptime,(R_gamma_A_IKF(estIdx)),'k');lgstr{end+1}='\gamma_{a}^{2} -real';
h(end+1)=plot(ptime,((gamma_A_IKF(estIdx))),'c');lgstr{end+1}='\gamma_{a} -real';
h(end+1)=plot(ptime, K_gain_IKF(estIdx),'g','LineWidth',2);lgstr{end+1}='k-real';
h(end+1)=plot(ptime,Pk_pos_IKF(estIdx),'r','LineWidth',2);lgstr{end+1}='Pk-real';
legend(h,lgstr,'Orientation','Horizontal');        
linkaxes(hsub,'x');    
       
%% Calculate the RMSE of the segments
colorscheme = 'cool';
xlabels={'A','B','C','D','E','F','G','H','I'};
ylabels={'CAHRS','AESKF'};
ylabels={'',''};
blogmap = true;
NumTicks = 9;

rmse_segs = [sync_segs(1:6:end),sync_segs(6:6:end)];
alphaSegs = cell(length(estIdx),1);
rmse_segs(1,1) = 1000;
[N_RUNS,~]=size(rmse_segs);

rse_real_cahrs = sqrt( (tru_eul_mea-real_eul_fix_sync).^2 );
rse_real_erkf = sqrt( (tru_eul_mea-real_eul_IKF_sync).^2 );

rmse_real_cahrs = nan(N_RUNS,2);
rmse_real_erkf  = nan(N_RUNS,2);

for r = 1:N_RUNS
    seg = [rmse_segs(r,1):rmse_segs(r,2)]';
    alphaSegs(seg) = xlabels(r);
    rmse_real_cahrs(r,:) = sum(rse_real_cahrs(seg,1:2))./length(seg);
    rmse_real_erkf(r,:) = sum(rse_real_erkf(seg,1:2))./length(seg);
end

rollHeatMapReal = [rmse_real_cahrs(:,1),...
                   rmse_real_erkf(:,1)]';

pitchHeatMapReal = [rmse_real_cahrs(:,2),...
                    rmse_real_erkf(:,2)]';

mincolor =  0;
maxcolor = 20;

if expDays == 1
    titleStr = sprintf('Up-Down');
    filesuffix = 'UpDown.csv';
elseif expDays == 2
    titleStr = sprintf('Up-Down with 90%s Rotation',char(176));
    filesuffix = 'UpDown_90_deg_rot.csv';
elseif expDays == 3
    titleStr = sprintf('Left-Right');
    filesuffix = 'LeftRight.csv';    
elseif expDays == 4
    titleStr = sprintf('Left-Right with 90%s Rotation',char(176));
    filesuffix = 'LeftRight_90_deg_rot.csv';

elseif expDays == 5
    titleStr = sprintf('Zig-Zag');
    filesuffix = 'ZigZag.csv';    
elseif expDays == 6
    titleStr = sprintf('Zig-Zag + 45%s ',char(176));
    filesuffix = 'ZigZag_45_deg_rot.csv';        
else
    
end           
%           

[p_imu_roll] = signrank(rmse_real_cahrs(:,1),rmse_real_erkf(:,1),'alpha',0.05,'tail','both');
[p_imu_pitch]= signrank(rmse_real_cahrs(:,2),rmse_real_erkf(:,2),'alpha',0.05,'tail','both');
fprintf('%s\n', titleStr);
fprintf('imu: %s: $p$ = %0.3f, %s: $p$= %0.3f \n','$\phi$',p_imu_roll,'$\theta$', p_imu_pitch);


h_fig = updateFigureContents(['Heat Map RMSE - Larger Font - ',titleStr]);
set(h_fig,'Units', 'Centimeters','Position', [0.5, 0, 9, 6], ...
    'PaperUnits', 'Centimeters', 'PaperSize', [9, 6]);
% subplot = @(m,n,p) ...
%     subtightplot (m, n, p, [0.075 0.01], [0.1 0.1], [0.06 0.1]);

subplot(2,1,1)
himg = heatmap(rollHeatMapReal,xlabels,[],'%1.2f',...
    'Colormap',colorscheme,'UseLogColormap', blogmap, ...
    'MinColorValue', mincolor, 'MaxColorValue', maxcolor,...
    'GridLines', ':','FontSize',10);

[RotX,RotY,Axis,XTicks,XTickLabels,YTicks,YTickLabels] =...
    XYrotalabel(0,90,gca,2:2:8,xlabels(2:2:8),1:2,ylabels,[],[]);
set(gca,'XTickLabel',[]);


subplot(2,1,2)
himg = heatmap(pitchHeatMapReal,xlabels,[],'%1.2f',...
    'Colormap',colorscheme,'UseLogColormap', blogmap, ...
    'MinColorValue', mincolor, 'MaxColorValue', maxcolor,...
    'GridLines', ':','FontSize',10);

[RotX,RotY,Axis,XTicks,XTickLabels,YTicks,YTickLabels] =...
    XYrotalabel(0,90,gca,2:2:8,xlabels(2:2:8),1:2,ylabels,[],[]);
set(gca,'XTickLabel',[]);

set(gcf,'NextPlot','add');
axes;


h = title(titleStr,'FontSize',12);
set(gca,'Visible','off');
set(h,'Visible','on');

%% data for export
% break
% %%
% imu_real_sim_data = table(phonedata.time(estIdx),...
% ...    sim_acc(truIdx,1),sim_acc(truIdx,2),sim_acc(truIdx,3),...
% ...    sim_gyr(truIdx,1),sim_gyr(truIdx,2),sim_gyr(truIdx,3),...
%     q_mea(:,1),q_mea(:,2),q_mea(:,3),q_mea(:,4),...
%     alphaSegs,...
%     phonedata.acc(estIdx,1),phonedata.acc(estIdx,2),phonedata.acc(estIdx,3),...
%     phonedata.gyr(estIdx,1),phonedata.gyr(estIdx,2),phonedata.gyr(estIdx,3),...
%     'VariableNames',...
%     {'timestamps',...
% ...    'SimAccX','SimAccY','SimAccZ',...
% ...    'SimGyrX','SimGyrY','SimGyrZ',...
%     'qw','qx','qy','qz','segment',...
%     'ImuAccX','ImuAccY','ImuAccZ',...
%     'ImuGyrX','ImuGyrY','ImuGyrZ'});
%               
% writetable(imu_real_sim_data,...
%     fullfile(savedir,[getCurrentDateTime,'_',filesuffix]),...
%     'FileType','text','Delimiter',';');
% 
% save(fullfile(savedir,...
%     [getCurrentDateTime,'_',strrep(filesuffix,'.csv','.mat')]),...
%     'imu_real_sim_data');
end



clear all;clc;close all;

opal_data = importdata('opal_test_data.mat' );


bIsDynamic = true(length(opal_data.acc(:,1)),1);
bIsDynamic(1:800) = false;


quats = opMimuDynamicCAHRS_ArcTan(...
    'bDynamicMu',bIsDynamic,'fs',opal_data.fsMimu,...
    'dynamic_mu_acc',0.003,'dynamic_mu_mag',0.001,...
    'static_mu_acc',0.5,'static_mu_mag',0.5,...
    'Acc',opal_data.acc,'Mag',opal_data.mag,'Gyr',opal_data.gyr);

figure;plot(quats,'DisplayName','quats')
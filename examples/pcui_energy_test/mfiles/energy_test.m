%%
%2D energy test results
clear all; close all; clc;

%dEbdt vs. phi_d
figure;
load /home/barthur/Desktop/energy_128.mat; 
tdt = t(2:end-1);
subplot(2,1,1)
plot(tdt,dEbdt,'b',tdt,dEpdt,'r--')
subplot(2,1,2)
plot(tn,1e-6*phi_d,'k',tn,phi_i,'k--');

%dEbdt vs. phi_d
figure;
load /home/barthur/Desktop/energy_256.mat; 
tdt = t(2:end-1);
subplot(2,1,1)
plot(tdt,dEbdt,'b',tdt,dEpdt,'r--')
subplot(2,1,2)
plot(tn,1e-6*phi_d,'k',tn,phi_i,'k--');

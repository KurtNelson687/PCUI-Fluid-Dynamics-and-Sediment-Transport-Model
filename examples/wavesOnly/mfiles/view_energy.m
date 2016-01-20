%view energy data for pcui internal seiche case
clc; clear all; close all;

%Low pass filter of energy
dt = 0.01;
Wn = 1/dt;
[B,A] = butter(2,1/(Wn/2));

figure;
hold on;
load energy_512x512_k-6.mat
% dEbdt = filter(B,A,dEbdt);
plot(tdt,dEbdt,'k',tn,phi_i,'m',t,phi_d,'m--');
load energy_256x256_k-6.mat
% dEbdt = filter(B,A,dEbdt);
plot(tdt,dEbdt,'b');
load energy_64x64_k-6.mat
% dEbdt = filter(B,A,dEbdt);
plot(tdt,dEbdt,'r');
load energy_32x32_k-6.mat
% dEbdt = filter(B,A,dEbdt);
plot(tdt,dEbdt,'g');
box on;
hold off;

figure;
hold on;
load energy_512x512_k-6.mat
plot(t,phi_d,'m',tn,phi_i,'k--');
load energy_256x256_k-6.mat
plot(t,phi_d,'k');
load energy_128x128_k-6.mat
plot(t,phi_d,'b');
load energy_64x64_k-6.mat
plot(t,phi_d,'r');
load energy_32x32_k-6.mat
plot(t,phi_d,'g');
box on;
hold off;

figure;
hold on;
load energy_512x512_k0.mat
plot(tdt,dEbdt,'k');
load energy_512x512_k-6.mat
plot(tdt,dEbdt,'b',t,phi_d,'b--',tn,phi_i,'b:');
box on;
hold off;

figure;
hold on;
load energy_256x256_k0.mat
plot(tdt,dEbdt,'kx-');
load energy_256x256_k0_harold.mat
plot(tdt,dEbdt,'b');
box on;
hold off;

figure;
hold on;
load energy_512x512_k-6.mat
plot(tdt,dEbdt,'k',t,phi_d,'k--',tn,phi_i,'k:');
load energy_256x256_k-6.mat
plot(tdt,dEbdt,'b',t,phi_d,'b--',tn,phi_i,'b:');
load energy_128x128_k-6.mat
plot(tdt,dEbdt,'r',t,phi_d,'r--',tn,phi_i,'r:');
box on;
hold off;



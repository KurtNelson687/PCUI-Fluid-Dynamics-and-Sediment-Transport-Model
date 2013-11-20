clc; clear all; close all;
load energy.mat

% -------------------------------------------------------------------------
% Plot energy quantities
% ------------------------------------------------------------------------- 
Et = Ek + Ep(25:25:14000);
dEtdt = (Et(3:end)-Et(1:end-2))/2/0.005/25;
dEkdt = (Ek(3:end)-Ek(1:end-2))/2/0.005/25;
phi_d_eff = phi_d(25:25:14000) - phi_i;

figure(1);
plot(tn,Et,'k',t,Ep,'k--',t,Eb,'b',t,Ea,'r',tn,Ek,'k:');

figure(2);
plot(tn(2:end-1),dEtdt,'k',tdt,dEpdt,'k--',tdt,dEbdt,'b',tdt,dEadt,'r',tn(2:end-1),dEkdt,'k:',tn,phi_z,'r:',[0 tn(end)],[0 0],'k')

figure(3);
plot(t,phi_d,'b',tn,eps_tot,'k',tn,phi_i,'b:',tn,phi_d_eff,'r',tdt,smooth(smooth(dEbdt)),'b');

figure(4);
eta = phi_d_eff./(phi_d_eff+eps_tot);
eta_cum = cumtrapz(tn,phi_d_eff)./(cumtrapz(tn,phi_d_eff)+cumtrapz(tn,eps_tot));
plot(tn,eta,'b',tn,eta_cum,'k');

figure(5);
% plot(tn,eps_top,'m',tn,eps_int,'b',tn,eps_low,'g',tn,eps_bot,'r',tn,eps_tot,'k')
bar(tn,[eps_top,eps_low,eps_int,eps_bot],'stacked');
hold on;
plot(tn,eps_tot,'k','linewidth',2)
hold off;
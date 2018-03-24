function [zo] = getZo(Cd,H,guess)
% Script: getZo
%
% Author: Kurt Nelson
%
% Purpose: This function computes the bed roughness, zo, assuming a log law
% velocity profile.
%
% Input:
% 1) Cd - drag coefficient
% 2) H - water depth
% 3) guess - guess for zo
kappa = 0.41;
f = @(zo)Cd- (1/kappa*(log(H/zo)+zo/H-1))^-2;
[zo,~] = fzero(f,guess);
end


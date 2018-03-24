function [dataWC] = computeCurrentSedFlux(dataWC,c1,c2,runsCup,aveStart,num_steps_period)
% Function: computeCurrentSedFlux
%
% Author: Kurt Nelson
%
% Purpose: This function computes diffusion, depostion, and turbulent
% sedimet fluxes associated with currents by period averaging.
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) c1,c2 -  hold coefficients for 1st and 2nd order derivative
%    respectively
% 3) runsCup - simulaitons to compute fluxes for
% 4) aveStart - start index for period-averaging
% 5) num_steps_period - number of steps per period
%
% Outputs:
% 1) dataWC.CsedSteady - sediment concentration associated with currents
% 2) dataWC.Fdc - sediment diffusion flux
% 3) dataWC.vCsedSteady - vertical turbulent sediment flux associated with
%    currents
% 4) dataWC.steadySedTurb - divergence of vertical turbulent sediment
%    flux associated with currents
% 5) dataWC.steadyDeposition - divergence of sediment depostion flux
%    associated with currents
% 6) dataWC.steadyDiffusion - divergence of sediment diffusion flux
%    associated with currents
%%
for plotNum = runsCup
    [m,n] = size(dataWC(plotNum).Csed);
    endAve = floor(n/num_steps_period)*num_steps_period; % index to end averaging
    aveIndex = aveStart:endAve; % vector containing period average indicies
    dataWC(plotNum).CsedSteady = mean(dataWC(plotNum).Csed(:,aveIndex),2)'; % sediment concentration
    dataWC(plotNum).vCsedSteady = mean(dataWC(plotNum).vCsed(:,aveIndex),2)'; % vertical turbulent sediment flux
    
    dataWC(plotNum).Ftc = mean(dataWC(plotNum).vCsed(:,aveIndex),2)'; % turbulent flux
    dataWC(plotNum).Fwc = mean(dataWC(plotNum).Csed(:,aveIndex),2)'*dataWC(plotNum).ws; % depositional flux
    
    % Intialize arrays
    dataWC(plotNum).steadyDeposition = zeros(m,1);
    dataWC(plotNum).steadyDiffusion = zeros(m,1);
    dataWC(plotNum).steadySedTurb = zeros(m,1);
    dataWC(plotNum).Fdc = zeros(m,1);
    
    % Compute sediment flux divergences
    for j=1:m
        if j ==1 % at the bed
            dataWC(plotNum).steadySedTurb(j) = c1(j,:)*dataWC(plotNum).vCsedSteady(j:j+4)';
            dataWC(plotNum).steadyDeposition(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j:j+4)'*dataWC(plotNum).ws;
            dataWC(plotNum).steadyDiffusion(j) = c2(j,:)*dataWC(plotNum).CsedSteady(j:j+4)'*dataWC(plotNum).ak;
            dataWC(plotNum).Fdc(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j:j+4)'*dataWC(plotNum).ak;
        elseif j ==2 % one cell up from the bed
            dataWC(plotNum).steadySedTurb(j) = c1(j,:)*dataWC(plotNum).vCsedSteady(j-1:j+3)';
            dataWC(plotNum).steadyDeposition(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-1:j+3)'*dataWC(plotNum).ws;
            dataWC(plotNum).steadyDiffusion(j) = c2(j,:)*dataWC(plotNum).CsedSteady(j-1:j+3)'*dataWC(plotNum).ak;
            dataWC(plotNum).Fdc(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-1:j+3)'*dataWC(plotNum).ak;
        elseif j == m % at the free-surface
            dataWC(plotNum).steadySedTurb(j) = c1(j,:)*dataWC(plotNum).vCsedSteady(j-4:j)';
            dataWC(plotNum).steadyDeposition(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-4:j)'*dataWC(plotNum).ws;
            dataWC(plotNum).steadyDiffusion(j) = c2(j,:)*dataWC(plotNum).CsedSteady(j-4:j)'*dataWC(plotNum).ak;
            dataWC(plotNum).Fdc(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-4:j)'*dataWC(plotNum).ak;
        elseif j ==m-1 % one cell down from the free-surface
            dataWC(plotNum).steadySedTurb(j) = c1(j,:)*dataWC(plotNum).vCsedSteady(j-3:j+1)';
            dataWC(plotNum).steadyDeposition(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-3:j+1)'*dataWC(plotNum).ws;
            dataWC(plotNum).steadyDiffusion(j) = c2(j,:)*dataWC(plotNum).CsedSteady(j-3:j+1)'*dataWC(plotNum).ak;
            dataWC(plotNum).Fdc(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-3:j+1)'*dataWC(plotNum).ak;
        else % remaining internal cells
            dataWC(plotNum).steadySedTurb(j) = c1(j,:)*dataWC(plotNum).vCsedSteady(j-2:j+2)';
            dataWC(plotNum).steadyDeposition(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-2:j+2)'*dataWC(plotNum).ws;
            dataWC(plotNum).steadyDiffusion(j) = c2(j,:)*dataWC(plotNum).CsedSteady(j-2:j+2)'*dataWC(plotNum).ak;
            dataWC(plotNum).Fdc(j) = c1(j,:)*dataWC(plotNum).CsedSteady(j-2:j+2)'*dataWC(plotNum).ak;
        end
    end
    
    dataWC(plotNum).CsedSteady = dataWC(plotNum).CsedSteady';
    dataWC(plotNum).vCsedSteady = dataWC(plotNum).vCsedSteady';
end

end


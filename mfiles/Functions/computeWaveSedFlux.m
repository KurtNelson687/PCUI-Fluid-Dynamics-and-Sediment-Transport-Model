function [dataWC] = computeWaveSedFlux(dataWC,c1,c2,runsCup,aveStart,num_steps_period,phases)
% Function: computeWaveSedFlux
%
% Author: Kurt Nelson
%
% Purpose: This function computes diffusion, depostion, and turbulent
% sedimet fluxes associated with waves by phase-averaging.
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) c1,c2 -  hold coefficients for 1st and 2nd order derivative
%    respectively
% 3) runsCup - simulaitons to compute fluxes for
% 4) aveStart - start index for phase-averaging
% 5) num_steps_period - number of steps per period
% 6) phases - index containing phases to average over
%
% Outputs:
% 1) dataWC.CsedWave - sediment concentration associated with waves
% 2) dataWC.waveDepositionFlux - sediment diffusion flux associated with
%    waves
% 3) dataWC.vCsedWave - vertical turbulent sediment flux associated with
%    waves
% 4) dataWC.waveSedTurb - divergence of vertical turbulent sediment
%    flux associated with waves
% 5) dataWC.waveDeposition - divergence of sediment depostion flux
%    associated with waves
% 6) dataWC.waveDiffusion - divergence of sediment diffusion flux
%    associated with waves
%%
%%%%% Computing the wave components

for plotNum = runsCup
    count = 1;
    [m,n] = size(dataWC(plotNum).Csed);
    
    % initialize arrays
    dataWC(plotNum).CsedWave= zeros(m,length(phases));
    dataWC(plotNum).waveDeposition = zeros(m,length(phases));
    dataWC(plotNum).waveDiffusion= zeros(m,length(phases));
    dataWC(plotNum).waveSedTurb= zeros(m,length(phases));
    
    for phase = phases
        aveIndex = aveStart-1+phase:num_steps_period:n; % indicies for phase-averaging
        
        % sediment concentration from waves - phase average then subtract
        % out period-average
        dataWC(plotNum).CsedWave(:,count) = mean(dataWC(plotNum).Csed(:,aveIndex),2);
        dataWC(plotNum).CsedWave(:,count) = dataWC(plotNum).CsedWave(:,count)-dataWC(plotNum).CsedSteady; 
        
        dataWC(plotNum).vCsedPhase(:,count) = mean(dataWC(plotNum).vCsed(:,aveIndex),2); % wave turbulent sediment flux
        
        %divergence of fluxes
        for j=1:m
            if j ==1 % at the bed
                dataWC(plotNum).waveDeposition(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j:j+4,count)*dataWC(plotNum).ws;
                dataWC(plotNum).waveDiffusion(j,count) = c2(j,:)*dataWC(plotNum).CsedWave(j:j+4,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveDiffusionFlux(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j:j+4,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveSedTurb(j,count) = c1(j,:)*dataWC(plotNum).vCsedPhase(j:j+4,count);
            elseif j ==2 % one cell up from the bed
                dataWC(plotNum).waveDeposition(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-1:j+3,count)*dataWC(plotNum).ws;
                dataWC(plotNum).waveDiffusion(j,count) = c2(j,:)*dataWC(plotNum).CsedWave(j-1:j+3,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveDiffusionFlux(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-1:j+3,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveSedTurb(j,count) = c1(j,:)*dataWC(plotNum).vCsedPhase(j-1:j+3,count);
            elseif j ==m % at the free-surface
                dataWC(plotNum).waveDeposition(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-4:j,count)*dataWC(plotNum).ws;
                dataWC(plotNum).waveDiffusion(j,count) = c2(j,:)*dataWC(plotNum).CsedWave(j-4:j,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveDiffusionFlux(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-4:j,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveSedTurb(j,count) = c1(j,:)*dataWC(plotNum).vCsedPhase(j-4:j,count);
            elseif j ==m-1 % one cell down from the free-surface
                dataWC(plotNum).waveDeposition(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-3:j+1,count)*dataWC(plotNum).ws;
                dataWC(plotNum).waveDiffusion(j,count) = c2(j,:)*dataWC(plotNum).CsedWave(j-3:j+1,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveDiffusionFlux(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-3:j+1,count)*dataWC(plotNum).ak;
                dataWC(plotNum).waveSedTurb(j,count) = c1(j,:)*dataWC(plotNum).vCsedPhase(j-3:j+1,count);
            else % remaining internal cells
                dataWC(plotNum).waveDeposition(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-2:j+2,count)*dataWC(plotNum).ws;
                dataWC(plotNum).waveDiffusion(j,count) = c2(j,:)*dataWC(plotNum).CsedWave(j-2:j+2,count)*dataWC(plotNum).ak;       
                dataWC(plotNum).waveDiffusionFlux(j,count) = c1(j,:)*dataWC(plotNum).CsedWave(j-2:j+2,count)*dataWC(plotNum).ak;    
                dataWC(plotNum).waveSedTurb(j,count) = c1(j,:)*dataWC(plotNum).vCsedPhase(j-2:j+2,count);
            end
        end
        dataWC(plotNum).waveSedTurb(:,count) = dataWC(plotNum).waveSedTurb(:,count)-dataWC(plotNum).steadySedTurb;
        dataWC(plotNum).vCsedWave(:,count) = dataWC(plotNum).vCsedPhase(:,count)-dataWC(plotNum).vCsedSteady;
        count = count+1;
    end
end


end


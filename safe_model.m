function [pns_total,pns_axis,sim_struct] = safe_model(gwf, dt, hw)
% [pns_total,pns_frac,sim_struct] = safe_model(gwf, dt, hw)
% function adapted from safe_gwf_to_pns() for IDEA simulation
%INPUTS
% gwf (nx3) in T/m
% dt  (1x1) in s
% hw  (struct) is structure that describes the hardware configuration and PNS
% response. Example: hw = safe_example_hw().
% OUTPUTS:
% pns_total: (nx1) in percentage of total stimulution limit
% pns_axis: (nx3) in percentage of stimulation limit for individual axis
% Sim_struct: siemens simulation predictions of the following
% 'SFX', 'SFY', 'SFZ' - Stimulation caused by invidual axis( in T/s)
% SLL: Total stimulation limit(deny level)(fraction)
% SLT: Total stimulation threshold(first level)(fraction)
%
% This PNS model is based on the SAFE-abstract;
% SAFE-Model - A New Method for Predicting Peripheral Nerve Stimulations in MRI
% by Franz X. Herbank and Matthias Gebhardt. Abstract No 2007. 
% Proc. Intl. Soc. Mag. Res. Med. 8, 2000, Denver, Colorado, USA
% https://cds.ismrm.org/ismrm-2000/PDF7/2007.PDF
% 
% The main SAFE-model was coded by Thomas Witzel @ Martinos Center,
% MGH, HMS, Boston, MA, USA.
% 
% The code was adapted/expanded by Filip Szczepankiewicz @ LMI
% BWH, HMS, Boston, MA, USA.
%
%praveenivp: shifting impulse response by 0.5 gradient raster brings us 
%even closer to IDEA simulation results.

dgdt = padarray(diff(gwf,1),[1 0],0,'pre') / dt;
pns = zeros(size(dgdt));

for i = 1:size(dgdt,2)
    
    switch i
        case 1
            p = hw.x;
        case 2
            p = hw.y;
        case 3
            p = hw.z;
    end
    
    pns(:,i) = safe_pns_model(dgdt(:,i), dt, p); %T/s
    
end
%IDEA simulation: SLL(deny level)
stim_limit=[hw.x.stim_limit hw.y.stim_limit hw.z.stim_limit]; %T/s
pns_frac=bsxfun(@rdivide,pns,stim_limit); %fraction
SLL=sqrt(sum(pns_frac.*pns_frac,2)); %eucledian norm of all axis(>1 exceeds safety limit)


pns_axis=pns_frac*100; %percentage
pns_total=SLL(:)*100; %percentage

%IDEA simulation: SLT(first level)
stim_thres=[hw.x.stim_thresh hw.y.stim_thresh hw.z.stim_thresh]; %T/s
pns_frac=bsxfun(@rdivide,pns,stim_thres); %fraction
SLT=sqrt(sum(pns_frac.*pns_frac,2)); %eucledian norm of all axis(>1 exceeds safety limit)
sim_struct=struct('SFX',pns(:,1),'SFY',pns(:,2),'SFZ',pns(:,3),'SLL',SLL(:));%,'SLT',SLT(:));

end

function stim = safe_pns_model(dgdt, dt, hw)
% function stim = safe_pns_model(dgdt, dt, hw)
%
% dgdt (nx3) is in T/m/s
% dt   (1x1) is in s
% All time coefficients (a1 and tau1 etc.) are in ms.
% 
% OUTPUT:
% stim: (nx3) is  in T/s for individual axis 
%
% from:https://github.com/filip-szczepankiewicz/safe_pns_prediction.git
stim1 = hw.a1 * abs( RCfilter_timedomain(dgdt     , hw.tau1, dt * 1000) );
stim2 = hw.a2 *      RCfilter_timedomain(abs(dgdt), hw.tau2, dt * 1000)  ;
stim3 = hw.a3 * abs( RCfilter_timedomain(dgdt     , hw.tau3, dt * 1000) );
stim = ((stim1 + stim2 + stim3)*hw.g_scale );

end


function fw = RCfilter_timedomain(dgdt, tau, dt)
% function fw = RCfilter_timedomain(dgdt, tau, dt)
mode='IIR'; 
switch mode
    case 'FIR' % bit slow only for testing
        t=0:dt:(10*tau); %5*tau -10*tau
        psf=exp(-1*(t+dt*0.5)/tau);
        psf=psf./(tau/dt);
        fw=conv(dgdt,psf);
        fw((end-length(psf)+2):end)=[];
    case 'IIR'
%         from:https://github.com/filip-szczepankiewicz/safe_pns_prediction.git
        alpha = dt / (tau + 0.5*dt);
        fw = zeros(size(dgdt));
        fw(1) = dgdt(1);       
        for s = 2:length(dgdt)
            fw(s) = alpha * dgdt(s) + (1-alpha) * fw(s-1);
        end 
end

end







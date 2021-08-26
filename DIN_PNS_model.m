function [pns_total,pns_axis,sim_struct] = DIN_PNS_model(gwf, dt, stim_limit,tc)
% [pns_total,pns_axis,sim_struct] = DIN_PNS_model(gwf, dt, stim_limit,tc)
% PNS-model based on DIN document
%INPUTS
% gwf (nx3) in T/m
% dt  (1x1) in s
% stim_limit - stimulation limit(rb) in T/s (default: 20 T/s)
% tc - Time constant in s(deafult: 0.36ms)
% OUTPUTS:
% pns_total: (nx1) in percentage of total stimulution limit
% pns_axis: (nx3) in percentage of stimulation limit for individual axis
% Sim_struct: siemens simulation predictions of the following
% 'SFX', 'SFY', 'SFZ' - Stimulation caused by invidual axis( in T/s)
% SLL: Total stimulation limit(deny level)(fraction)
% SLT: Total stimulation threshold(first level)(fraction)
%
% 

switch nargin
    case 1
        dt =10e-6; %s
    case 2
        stim_limit=20; %T/s
        tc=0.36e-3; %s
    case 3
        tc=0.36e-3; %s
end

%slew rate
dgdt = padarray(diff(gwf,1),[1 0],0,'pre') / dt;

% calculate R(t): see page 110 of DIN-VDE document
t=0:dt:dt*(length(dgdt)-1);
R=zeros(size(dgdt));
for j=1:size(R,2)
    for i=1:size(R,1)
        R(i,j)=sum((dgdt(1:i,j)*tc)./col((tc+max(t(1:i))-t(1:i)).^2),1)*dt;
    end
end

%IDEA simulation: SLL(deny level)
R_scaled=abs(R)/pi; % pi scaling gives agreement with IDEA SIM
SLL=sqrt(sum(R_scaled.*R_scaled,2))/stim_limit;
sim_struct=struct('SFX',R_scaled(:,1),'SFY',R_scaled(:,2),'SFZ',R_scaled(:,3),'SLL',SLL);

%same outputs as 
pns_frac=bsxfun(@rdivide,R_scaled,stim_limit); %fraction
pns_axis=pns_frac*100; %percentage
pns_total=SLL(:)*100; %percentage

% %IDEA simulation: SLT(first level)
% stim_thres=0.8*stim_limit;
% pns_frac=bsxfun(@rdivide,R_scaled,stim_thres); %fraction
% SLT=sqrt(sum(pns_frac.*pns_frac,2)); %eucledian norm of all axis(>1 exceeds safety limit)
% sim_struct.SLT=SLT;
end


 function x = col(x)
x = x(:);
 end







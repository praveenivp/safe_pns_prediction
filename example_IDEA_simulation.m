%% read DSV files and simple plot
fn=uigetdir(userpath);%'folder path to all dsv files';
[traces_grad,header_grad]=readDSVFolder(fn,{'GRX','GRY','GRZ'});
[traces_slew,header_slew]=readDSVFolder(fn,{'D0X','D0Y','D0Z'});
[traces_stim,header_stim]=readDSVFolder(fn,{'SFX', 'SFY', 'SFZ','SLL'},true);
% figure,plot_linked(traces_grad)
% figure,plot_linked(traces_slew)
% figure,plot_linked(traces_stim)
%% load Gradient Specifications
[fn2,fp]=uigetfile(fullfile(fn,'\*.asc')); %select asc file
hw = safe_hw_from_asc(fullfile(fp,fn2), true);

%% predict stimulation limit with SAFE model
gwf=cat(1,traces_grad.GRX ,traces_grad.GRY,traces_grad.GRZ)'*1e-3;% T/m
dt=10e-6; %Raster time in s
[pns_total,pns_frac,SAFE_struct]  = safe_model(gwf, dt, hw);
figure,plot_linked(traces_stim),
plot_linked(SAFE_struct),
legend('IDEA Simulation','SAFE model')

%% DIN model vs safe model
tc=0.36e-3; %s
stim_limit_DIN=20; %rb in T/s when effective stimulation time is <1 ms
[pns_total,pns_axis,DIN_struct] = DIN_PNS_model(gwf, dt, stim_limit_DIN,tc);
figure,plot_linked(traces_stim),
plot_linked(DIN_struct)
legend('SAFE model','DIN-PNS model')

%% (old) predict stimulation limit with SAFE model ( https://github.com/praveenivp/safe_pns_prediction.git )
doPadding=false;
rf=[];
[pns, res] = safe_gwf_to_pns(gwf, rf, dt, hw, doPadding);

%pns percentage to SFX,SFY,SYZ,SLL
pns=padarray(pns,[1 0],0,'pre');
stim_limit_siemens=[hw.x.stim_limit hw.y.stim_limit hw.z.stim_limit]; %T/s
SF_XYZ=bsxfun(@times,pns,stim_limit_siemens/100);
figure,plot_linked(traces_stim),
plot_linked(struct('SFX',SF_XYZ(:,1),'SFY',SF_XYZ(:,2),'SFZ',SF_XYZ(:,3),'SLL',sqrt(sum(pns.*pns,2))./100))
legend('IDEA Simulation','SAFE model')
%%
function plot_linked(traces)
names=fieldnames(traces);
% figure,
for i=1:length(names)
ax{i}=subplot(length(names),1,i);
hold on,plot(getfield(traces,names{i}));
title(names{i})

if(any(strcmpi(names{i},{'SFX','SFY','SFZ'}))) 
    ylabel('Stimulation (T/s)'),xlabel('Time(x10us)')
elseif(any(strcmpi(names{i},{'SLL','SLT'}))) 
    ylabel('Fraction'),xlabel('Time(x10us)')
elseif(any(strcmpi(names{i},{'GRX','GRY','GRZ'})))
    ylabel('Amp (mT/m)'),xlabel('Time(x10us)')
elseif(any(strcmpi(names{i},{'D0X','D0Y','D0Z'})))
    ylabel('Slewrate (mT/m/ms)'),xlabel('Time(x10us)')
end
end
 linkaxes([ax{:}],'x')
end



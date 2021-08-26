function [traces,header]=readDSVFolder(FolderPath,paraList,plotON)
% [traces,header]=readDSVFolder(FolderName,paraList,plotON)
% function to read IDEA simulation Folder in MATLAB
%
%INPUTS:
%Foldepath: full folder path or []
%paraList - List of list files to read(default: {'ADC','GRX','GRY','GRZ'})
%'RF2', 'RFD' -RF amplitude and phase
% 'GRX', 'GRY', 'GRZ' -Gradients 
% 'D0X', 'D0Y', 'D0Z' -slew rate
% 'M0X', 'M0Y', 'M0Z' - 0th order gradient moment
% 'M1X', 'M1Y', 'M1Z' -1st order gradient moment
% 'MXX', 'MXY', 'MXZ' - Maxwell gradient moment
% 'SFX', 'SFY', 'SFZ' - Stimulation caused by invidual axis
% 'PXA', 'PYA', 'PZA', - Gradient duty Cycle(time constast 1ms)
% 'PXB', 'PYB', 'PZB', - Gradient duty Cycle(time constast 10 ms)
% 'PXC', 'PYC', 'PZC' -  Gradient duty Cycle(time constast 100 ms)
% 'SAP'- Level of stimualtion 
%'SLL' - Total stimulation limit(deny level)
%'SLT' -Total stimualtion threshold (warn level)
% 'ADC', 'ERR', 'GAM', 'GSR', 'INF', 'NC1', 
%
%plotON: Booleean to enable plotting(deafult: false)
%
%OUTPUT:
%traces- Struct with simulation traces
%header- struct with the header information of about simulation traces
%
%praveen.ivp@gmail.com


if(~exist('FolderPath','var')|| ~isfolder(FolderPath))
        FolderPath=uigetdir(userpath,'Select the DSV folder to read');
end
if(~exist('paraList','var'))
        paraList={'ADC','GRX','GRY','GRZ'};
end
if(~exist('plotON','var'))
        plotON=false;
end


fileList=dir(fullfile(FolderPath,'*.dsv'));
    for i=1:length(fileList)
        for currPara=1:length(paraList)
        if(contains(fileList(i).name,paraList{currPara}))
            [waveform,hd]=dsv_read(fullfile(FolderPath,fileList(i).name));
            eval(sprintf('traces.%s=waveform;',paraList{currPara}));
            eval(sprintf('header.%s=hd;',paraList{currPara}));
        end
        end
    end
    
    if(plotON)
        plot_linked(traces)
    end
end


function [data, header]= dsv_read(filename)
%[data, header]=function dsv_read(filename)
% reads DSV_V100 files(type of run lengthencoding file)

%Author: praveen.ivp@gmail.com
%Error checking

if nargin==0
    warning('No input file'); 
    [filename,path]=uigetfile('*.dsv','Select a dsv file');
    filename=fullfile(path,filename);
end
if(~exist(filename,'file'))
    error('Input file not found');
end

 %% DSV scanning     
 fid=fopen(filename);
frewind(fid);
idx=1;

%read the header file
while(~feof(fid))
    line=fgetl(fid);

if(contains(line,'='))
    temp=regexp(line,'=','split');
    if(~isnan(str2double(temp{2})))
    evalc(strcat('header.',temp{1},'=',temp{2},';'));  
    else
    evalc(strcat('header.',temp{1},'=''',temp{2},''';'));   
    end
end
if(strcmp(line,'[VALUES]'))
    break;
end
    idx=idx+1;
end    

if(~(strcmp(header.FORMAT,'DSV_V0100')))
    error('Invalid file format')
end
%read encoded data
rle_en=cell2mat(textscan(fid,'%d'));
fclose(fid);

%Converting: DSV_0V100 is a kind of runlength encoding scheme 
%DSV_v100: 0,0,0,0,1,5,5,5,6,7 -> 0,0,2,1,5,5,1,6,7
data=zeros(1,header.SAMPLES);
c=zeros(1,length(rle_en));
idx=1;  %index for encoded array
idx1=1; %index for decoded array
while(idx<length(rle_en))
    a=double(rle_en(idx));
    b=rle_en(idx+1);
    if(a==b)
         data(1,idx1:(idx1+rle_en(idx+2)+1))=a.*ones(1,rle_en(idx+2)+2);
        idx1=(idx1+rle_en(idx+2)+2);
        c(idx)=rle_en(idx+2)+2;
        idx=idx+3;
    else
         data(1,idx1)=a;
        idx=idx+1;
       idx1=idx1+1;
    end
end

%data post processing
%the first data is the absolute value all others are signal change
%so that it can compress ramps/linear signal change
%all data points are scaled by VERTFACTOR
data=cumsum(data./header.VERTFACTOR);
%Error checking
if(~(length(data)==header.SAMPLES))
    warning('decoded data samples didn''t match the header specification(#samples)'); 
end
if((max(data)>(header.MAXLIMIT))||(round(min(data))<(header.MINLIMIT)))
    warning('decoded data samples didn''t match the header specification(MAXLIMIT||MINLIMIT)'); 
end
end

function plot_linked(traces)
names=fieldnames(traces);
 figure,
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
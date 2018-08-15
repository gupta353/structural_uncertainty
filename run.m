% save_dir=directory which contains rainfall, streamflow and other relevant
% data
% wt_matrix=mat-file of watershed topology matrix

clear
close all
clc
fclose('all');

global GEOMORPH GLOBAL_DATA

addpath(['D:\Research\Thesis_work\Structural_'...
    'uncertainty\MatLab_codes\20180222']);
%% input data
save_dir=['D:\Research\Thesis_work\Structural_uncertainty'...
    '\MatLab_codes\20180222\huc_0712000408'];
wt_matrix=['D:\Research\Thesis_work\Structural_uncertainty\MatLab'...
     '_codes\20180222\huc_0712000408\wt_0712000408'];

arealen_filename=fullfile(save_dir,'area_length.txt');
fid=fopen(arealen_filename,'r');
area_length=textscan(fid,'%s%s%s');
fclose(fid);
ddrain_area=area_length{2}(2:end);
length_link=area_length{3}(2:end);
ddrain_area=cellfun(@str2num,ddrain_area);
length_link=cellfun(@str2num,length_link);

load(wt_matrix);
wt=wt_0712000408;
nodeobs=1;
support=0:1015200;           % support for instantaneous unit-hydrograph computation
R_time_scale=60;             % time-scale (in mins) of rainfall data (use streamflow data at hourly scale)
S_time_scale=60;             % time-scale (in mins) of streamflow data   
param=[3,2,1,1,1,1];
%% Run the mainScript
GEOMORPH=mainScript(wt,ddrain_area,length_link,save_dir);
%% read rainfall data
p_filename=fullfile(save_dir,'rainfall_excess.txt');
fid=fopen(p_filename,'r');
R=textscan(fid,'%s%s%s');
fclose(fid);
R=R{3}(2:end);
R=cellfun(@str2num,R);     % daily rainfall data

%% read the text-file containing the list of links draining into nodes
fid=fopen(['D:\Research\Thesis_work\Structural_uncertainty'...
    '\MatLab_codes\20180222\huc_0712000408\links_draining'...
    '_into_each_node.txt'],'r');
node_link=textscan(fid,'%s%s');
links=node_link{2}{nodeobs+1};
links=strsplit(links,',');
links=cellfun(@str2num,links);

%% read area-length data
fid=fopen(['D:\Research\Thesis_work\Structural_uncertainty'...
    '\MatLab_codes\20180222\huc_0712000408\area_length.txt'],'r');
area_length=textscan(fid,'%s%s%s');
fclose(fid);
ddarea=area_length{2}(2:end);
ddarea=cellfun(@str2num,ddarea);
darea=sum(ddarea(links));
%% read baseflow separated surface flow data (in m^3/s)
strm_filename=fullfile(save_dir,'streamflow_data',...
    strcat('surfaceflow_event_',num2str(nodeobs),'.txt'));
fid=fopen(strm_filename,'r');
strmobs=textscan(fid,'%s%s%s');
time_steps_date=strmobs{1}(2:end);
time_steps_time=strmobs{2}(2:end);
wrapper=@(x1,x2)strcat(x1,32,x2);
time_steps=cellfun(wrapper,time_steps_date,time_steps_time,...
    'UniformOutput',false);
wrapper_1=@(x)datenum(x,'dd-mmm-yyyy HH:MM');
time_steps=cellfun(wrapper_1,time_steps);
strmobs=strmobs{3}(2:end);
strmobs=0.0283168*cellfun(@str2num,strmobs); % 0.0283168 for conversion of units from ft^3/s to m^3/s
% sampling of streamflow data at hourly-time scale (average streamflow of each is measured)
scale=24*60*(time_steps(2)-time_steps(1));  % in mins
skip_tsteps=60/scale;
count=0;
for i=1:round(skip_tsteps):length(strmobs)
    count=count+1;
    if length(strmobs)>i+3
        strm(count)=mean(strmobs(i:i+3));
        tsteps(count)=time_steps(i);
    end
end
strmobs=strm';
time_steps=tsteps';
%% data to pass
GLOBAL_DATA.save_dir=save_dir;
GLOBAL_DATA.R=R;
GLOBAL_DATA.darea=darea;
GLOBAL_DATA.nodeobs=nodeobs;
GLOBAL_DATA.strmobs=strmobs;
GLOBAL_DATA.support=support;
GLOBAL_DATA.R_time_scale=R_time_scale;
GLOBAL_DATA.S_time_scale=S_time_scale;
GLOBAL_DATA.time_steps=time_steps;
break
%% objective function computation
%
% L= objectiveComputation(param,GLOBAL_DATA,GEOMORPH);
%}
%% optimization
%
loss=@(param)objectiveComputation(param,...     % loss function
    GLOBAL_DATA,GEOMORPH);           
parent=[0.35,0.1,1,1,1,1];                     % initial guess
lb=[0.00001,0.00001,0,0,0.5,0.00001];           % lower bound
ub=[1,1,1,1,1,10];                              % upper bound
% [minimum,fval] = anneal(loss,parent);
[bestx,bestf] = sceua(parent,lb,ub,100000,5,0.01,[],20,2000,0);
%}
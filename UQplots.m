% plot of the observed, simulated streamflow with 95% credible-region


global GLOBAL_DATA GEOMORPH

%% input data
load(['D:\Research\Thesis_work\Structural_uncertainty\MatLab'...
    '_codes\20180222\results\optimization_06_29_2018']);

save_dir=GLOBAL_DATA.save_dir;
wt=GEOMORPH.wt;
nodeobs=GLOBAL_DATA.nodeobs;
support=0:2000000;           % support for instantaneous unit-hydrograph computation
time_scale='hourly';        % time-scale of rainfall data
param=bestx;
%% read rainfall data
R=GLOBAL_DATA.R;

%% read the text-file containing the list of links draining into nodes
fid=fopen(['D:\Research\Thesis_work\Structural_uncertainty'...
    '\MatLab_codes\20180222\huc_0512011115\links_draining'...
    '_into_each_node.txt'],'r');
node_link=textscan(fid,'%s%s');
links=node_link{2}{nodeobs+1};
links=strsplit(links,',');
links=cellfun(@str2num,links);
%% read area-length data
fid=fopen(['D:\Research\Thesis_work\Structural_uncertainty'...
    '\MatLab_codes\20180222\huc_0512011115\area_length.txt'],'r');
area_length=textscan(fid,'%s%s%s');
fclose(fid);
ddrain_area=area_length{2}(2:end);
ddrain_area=cellfun(@str2num,ddrain_area);
darea=sum(ddrain_area(links));
length_link=area_length{3}(2:end);
length_link=cellfun(@str2num,length_link);
%% read baseflow separated surface-flow data (in m^3/s)
strmobs=GLOBAL_DATA.strmobs;
time_steps=GLOBAL_DATA.time_steps;
%% Run the mainScript
% mainScript;
%% plots for structural uncertainty
theta_opt=bestx(1:2);       % optimal model-parameter
Cp_opt=diag(bestx(3:4));    % optimal parameter-perturbation variances
beta_opt=bestx(5);          % optimal shape-parameter
sigmaobs2_opt=bestx(6);     % optimal observation-noise variance
mu=hydrograph(theta_opt,... % mean vector
    GLOBAL_DATA,GEOMORPH)';
n=10000;                      % number of samples to be drawn

Ct_opt=CovarianceMatrixEstimation(Cp_opt,theta_opt,sigmaobs2_opt,...
    GLOBAL_DATA,GEOMORPH);

% surface plot of correlation matrix
%
figure;
h=surface(corrcov(Ct_opt));
colormap gray
colorbar
xlim([1 25])
ylim([1 25])
xlabel('time-step','fontname','arial','fontsize',14)
ylabel('time-step','fontname','arial','fontsize',14)
set(gca,'fontname','arial','fontsize',14)
set(h,'EdgeColor','none')
box.linewidth=2;
set(gca,box)
save_fig_name=fullfile(save_dir,...
    strcat('correlatin_matrix_',num2str(nodeobs)));
print(save_fig_name,'-r300','-dtiff')
close all
    %}
% structural uncertainty plot
%{
r=mulgennormrnd(n,mu,Ct_opt,beta_opt);

r975=prctile(r,97.5);
% r025=prctile(r,2.5);
r025=prctile(r,2.5);
r025(r025<0)=0;

figure;
f2=area(time_steps,r975,'FaceColor',[0.5,0.8,0.5]);
hold on
f3=area(time_steps,r025,'FaceColor',[1,1,1]);
f1=plot(time_steps,mu,'color',[0,0,0],'linewidth',3);
f4=plot(time_steps,strmobs,'o','MarkerFaceColor',...
    [0.80,0,0],'color',[0.80,0,0],'MarkerSize',9);
datetick('x','mmmm dd')
xlim([time_steps(1) time_steps(end)])
xlabel('time (date)','fontname','arial','fontsize',14)
ylabel('streamflow (m^{3} s^{-1})','fontname','arial',...
    'fontsize',14,'interpreter','Tex')
legend([f1,f2,f4],{'mean simulated',...
    '95% credible-region','observed'},...
    'fontname','arial','fontsize',14);
legend('boxoff')
box.linewidth=2;
set(gca,'fontsize',14,'fontname','arial',...
    'plotBoxAspectRatio',[1.2,1,1],box)
print(fullfile(save_dir,strcat('node',...
    num2str(nodeobs))),'-dtiff','-r300')
%}
%%% Estimation of streamflow by convolution of rainfall and instantanous unit hydrograph
% inputs: uh=instantanous unit hydrograph at point of observation/simulation (values at every second)
%         R=rainfall time-series (in mm)
%         darea=draiange area of the point of observation/simulation (in hectares)
%         time_scale=a string with value either 'hourly' (for hourly
%         rainfall data) or 'daily' (for daily rainfall data)
%         support=time-length at which instantaneous unit-hydrograph would
%         be computed (at seconds-scale)
% output: streamflow=daily average streamflow in m^3/s

function [streamflow,uh]=hydrograph(theta,GLOBAL_DATA,GEOMORPH)

R=GLOBAL_DATA.R;
darea=GLOBAL_DATA.darea;
node=GLOBAL_DATA.nodeobs;
support=GLOBAL_DATA.support;
R_time_scale=GLOBAL_DATA.R_time_scale;
S_time_scale=GLOBAL_DATA.S_time_scale;

vel_stream=theta(1);
vel_hillslope=theta(2);
uh=unitHydrograph(node,vel_stream,vel_hillslope,support,GEOMORPH);

% if strcmpi(time_scale,'hourly')
%     Tr=3600; % for hourly rainfall data 
% elseif strcmpi(time_scale,'daily')
%     Tr=24*3600; % for daily rainfall data
% else
%     error('time_scale should be either ''hourly'' or ''daily''');
% end
Tr=60*min(R_time_scale,S_time_scale);
% convsersion of uh-scale to minimum of rainfall time-scale and stteamflow scale
l_uh=length(uh);
addit=Tr-Tr*(l_uh/Tr-floor(l_uh/Tr));
uh=[uh;zeros(ceil(addit),1)];
uh=reshape(uh,Tr,length(uh)/Tr);
uh=mean(uh);
uh=[0;uh'];

% re-arranging rainfall data
% for i=1:length(R)
%     Rs(1+Tr*(i-1):i*Tr)=R(i)/Tr;
% end

streamflow_s=10*darea*conv(uh,R); % 1-mm x 1-hectare=10-m^3
streamflow=streamflow_s';
%convsersion of streamflow from hourly to daily scale
% count=0;
% Ts=24;
% for ii=1:Ts:length(streamflow_s)
%     count=count+1;
%     ind=ii:ii+Ts-1;
%     if ind(end)<length(streamflow_s)
%         streamflow(count)=mean(streamflow_s(ind));
%     end
% end
% streamflow=[0,streamflow];
end
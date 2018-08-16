% WATERSHED GEOMORPHOLOGY
% User inputs: wt=watershed-topology matrix
%              ddrain_area=direct-drainage area of each link (in Hectare)
%              length_link=length of each of the links (in meters)
%              save_dir=directory in which the output text-files have to be
%              saved
% wt: User has to define a watershed topology matrix (wt), Each row of wt
% corresponds to a unique path in watershed and each column in wt
% corresponds to a link in the watershed.Total number of columns in wt 
% should be equal to total number of links in the watershed. wt is a 
% zero-one matrix with ones defining the links in the watershed. 
% For example, suppose path-i consists of link-numbers 1, 2 and 4, then 
% ith row of wt would contain ones in 1st, 2nd and 4th columns and rest 
% of the entries of the wt would be zero.
% Outputs: GEOMORPH=a structure array of that contains
%                   wt=user input watershed-topology matrix
%                   RTD=residence time-distribution in each reach which is
%                   assumed to 'expoenntial' (other RTDs not allowed)
%                   net=river network as a tree
%                   node_net=netwrok of nodes as a tree
%                   link_net=network of links as a tree
%                   strahler_order=strahler order of each link
%                   ddrain_area=direct-drainage area of each link (in Hectares)
%                   pathPoints=
%                   length_link=length of each of the links (in meters)
%                   length_num_mat:
%                   length_num_mat(:,1)=strahler-order 
%                   length_num_mat(:,2)=number of streams of that order;
%                   length_num_mat(:,3)=total length of streams of that order;
%                   length_num_mat(:,4)=average length;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Except above mentioned outputs this routine
% (1) plots the tree-representation of the link-network of the watershed
% and asigns a number to each of it's nodes
% (2) saves a text file conatining the output length_num_mat into directory 'save_dir'
% (4) Computes Horton's bifurcation ratio (Rb) and save it in a text-file into directory 'save_dir'
% (5) saves a text-file containing links draining into each node into directory 'save_dir'
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes and caution: 
% (1) tempsum in section 'tree structure' should be either 'numpaths' or
%     'numpaths-1'; automate the selection
% (2) temp_net in section 'link_network' should either be 'net(:,2:end)' or
%     'net(:,1:end)'; automate the selection
% (3) 'wt' matrix contains 'links' in the path, NOT 'nodes'

function [GEOMORPH]=mainScript(wt,ddrain_area,length_link,save_dir)

addpath(['D:\Research\Thesis_work\Structural'...
    '_uncertainty\MatLab_codes\20180222\uniquecell'])

%% user-inputs
%%%%%%%%%% Data required %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix of watershed-topology (rows correspond to # of paths, columns correspond to # links)
%{
% hypothetical wt matrix
 wt=[1,1,0,0,0,0,0;
     0,1,1,0,1,1,0;
     0,1,1,1,0,1,0;
     0,1,1,0,0,0,1];

ddrain_area=10*[1:7];                    % direct-draiange area of link

length_link=5:5:35;
vel_stream=1;
vel_hillslope=0.5;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
numPaths=size(wt,1);                % compute # of paths from watershed-topology matrix
numLinks=sum(wt,2);                 % compute # of links in each path from watershed-topology matrix
RTD='exponential';                  % residence-time distribution (RTD)
numNodes=numLinks;                  % # of nodes in each path
totnumNodes=sum(numLinks)+1;        % total # of nodes
for pathind=1:numPaths
    pathPoints{pathind}=...         % link indices in each path
        find(wt(pathind,:)==1);
end
pathPoints=pathPoints';

for linkind=1:size(wt,2)
    linkPaths{linkind}=...         % paths common to each link
        find(wt(:,linkind)==1);
end
linkPaths=linkPaths';
%% error messages
wt_err=wt;
wt_err(wt_err==1)=0;
sum_err=sum(wt(:,1:end-1),1);
sum_err(sum_err<numPaths)=0;
num_unique_paths=size(unique(wt,'rows'),1);
if ~isempty(find(wt_err))
    error('atleast one of the entries in wt matrix is neither zero nor one');
elseif ~isempty(find(sum_err))
    error(['atleast one of the links except the most downstreamis common' ...
    'to all the paths which is physically impossible']);
elseif num_unique_paths<numPaths
    error('atleast one of the path has been specified more than once in wt matrix');
end
%% tree structure (output variable=net)
s=sum(wt,1);                         % number of paths shared by each link
tempsum=numPaths;
net=repmat({1:numPaths},numPaths,1);
col_ind=ones(numPaths,1);

while tempsum>0
    ind=find(s==tempsum);
    if isempty(ind)
        tempsum=tempsum-1;
    else
        for i=1:length(ind)
            tempcell=linkPaths{ind(i)};
            count=0;
            for row=1:numPaths                
                col=col_ind(row);
                if isempty(MY_intersect(tempcell,net{row,col}))
                else
                    count=count+1;
                    if count<=tempsum
                        col_ind(row)=col_ind(row)+1;
                        net{row,col+1}=tempcell;
                    else
                        break
                    end
                end
            end
        end
        tempsum=tempsum-1;
    end
end
%% node network (output variable=node_net) (node_net assigns upstream node-number of a link)
%
treedepth=size(net,2);
node=0;
[Au,idx,idx2]=uniquecell(net);
idx2_res=reshape(idx2,size(net,1),size(net,2));
for i=1:treedepth
    un=unique(idx2_res(:,i));
    for j=1:length(un)
        if un(j)~=1 % if the corresponding element of net(:,i) doesn't contain empty matrix
            node=node+1;
            ind=find(idx2_res(:,i)==un(j));
            node_net(ind,i)=node;
        end
    end
end
parents=zeros(numPaths,1);
for i=2:treedepth
    for j=1:numPaths
        parents(j,i)=node_net(j,i-1);
    end
end
max_node=max(node_net(:));
for i=1:max_node
    [r,c]=find(node_net==i);
    nodes(i)=parents(r(1),c(1));
end
treeplot(nodes)
set(gca,'Ydir','reverse')
[x,y]=treelayout(nodes);
for i=1:max_node
    text(x(i),y(i),num2str(i),'VerticalAlignment',...
        'top','fontSize',16,'fontName','arial');
end
%}
%% link network (output_variable=link_net)
%
temp_net=net;
temp_wt=wt;
temp_treedepth=size(temp_net,2);
link_net=size(wt,2)*ones(numPaths,1);
for i=2:temp_treedepth
    temp_link=[];
    for j=1:numPaths
        paths=temp_net{j,i};
        
        if isempty(paths)
            link_net(j,i)=0;
        else
            temp_s=sum(temp_wt(paths,:),1);
            link=find(temp_s==length(paths));
            temp_link=[temp_link,link];
            if isempty(link)
                link_net(j,i)=0;
            else
                link_net(j,i)=link;
            end
        end
       
    end
    
     temp_wt(:,temp_link)=0;
end
%}
%% define strahler order of links (output_variable=strahler_order)
link_net_depth=size(link_net,2);
ind=find(s==1);
strahler_order(ind)=1;
for i=link_net_depth-1:-1:2         % last number should be 2 in one tempsum=numPaths (L76) and 1 in case tempsum=numPaths-1 (L76)
    uni_links=unique(link_net(:,i));
    
    if isempty(uni_links==0)       
    else
        uni_links(uni_links==0)=[];
    end
    
    for j=1:length(uni_links)
        ind=find(link_net(:,i)==uni_links(j));
        ddrain_links=unique(link_net(ind,i+1));
        
        if ddrain_links~=0
            so1=strahler_order(ddrain_links(1));
            so2=strahler_order(ddrain_links(2));
            if so1==so2
                strahler_order(uni_links(j))=so1+1;
            else
                strahler_order(uni_links(j))=max(so1,so2);
            end
        end
        
    end
    
end
% strahler_order=strahler_order';
% filename=fullfile(save_dir,'strahler_order.txt');
% fid=fopen(filename,'wt');
% fprintf(fid,'%s %s\n','subbasin','strahler_order');
% data_write=[(1:size(wt,2))',(strahler_order)];
% dlmwrite(filename,data_write,'-append','delimiter','\t');
% fclose(fid);
%% mean length of the streams of each strahler-order in entire network (output_variable=length_num_mat)
% length_num_mat(:,1)=strahler-order; length_num_mat(:,2)=number of streams of that order;
% length_num_mat(:,3)=total length of streams of that order;
% length_num_mat(:,4)=average length;

max_order=max(strahler_order);
length_num_mat=[(1:max_order)',...
    zeros(max_order,1),zeros(max_order,1)];
links_used=[];
for order=1:max_order
    links=find(strahler_order==order);                        % links of order 'order'
    
    for path=1:numPaths
        link_path_int{path}=intersect(links,pathPoints{path});
        numLinks_int(path)=length(link_path_int{path});
    end
    
    while ~isempty(links)      
        [max_int,ind]=max(numLinks_int);
        links_curr_iter=link_path_int{ind};
        length_num_mat(order,3)=length_num_mat(order,3)...
            +sum(length_link(links_curr_iter));
        length_num_mat(order,2)=length_num_mat(order,2)+1;
        links=MY_setdiff(links,links_curr_iter); 
        for path=1:numPaths
            if ~isempty(MY_intersect(links_curr_iter,link_path_int{path}))
                link_path_int{path}=[];
                numLinks_int(path)=0;
            end
        end
    end
   
end
length_num_mat(:,4)=length_num_mat(:,3)./length_num_mat(:,2);
filename=fullfile(save_dir,'mean_length.txt');
fid=fopen(filename,'wt');
fprintf(fid,'%s %s %s %s\n','strahler_order',...
    'number_of_streams','total_length(meters)','mean_length(meters)');
dlmwrite(filename,length_num_mat,'-append','delimiter','\t');
fclose(fid);
%}
%% determination of Horton's bifurcation ratio
%
% plot(length_num_mat(:,1),log(length_num_mat(:,2)))
mdl=fitlm(length_num_mat(:,1),log(length_num_mat(:,2)));
Rb=exp(-mdl.Coefficients(2,1).Estimate);
filename=fullfile(save_dir,'hortons_bifurcation.txt');
fid=fopen(filename,'wt');
fprintf(fid,'%s %s\n','Horton''s bifurcation ratio (Rb) =',num2str(Rb));
fclose(fid);
%}
%% determination of Horton's Length ratio
%
% plot(length_num_mat(:,1),log(length_num_mat(:,2)))
mdl=fitlm(length_num_mat(:,1),log(length_num_mat(:,4)));
Rl=exp(mdl.Coefficients(2,1).Estimate);
filename=fullfile(save_dir,'hortons_bifurcation.txt');
fid=fopen(filename,'a');
fprintf(fid,'%s %s\n','Horton''s length ratio (Rl) =',num2str(Rl));
fclose(fid);
%}
%% total drainage area of each node
max_node=max(node_net(:));
filename=fullfile(save_dir,'links_draining_into_each_node.txt');
fid=fopen(filename,'wt');
fprintf(fid,'%s %s\n','node','links_draining');
for node=1:max_node
    [node_r,node_c]=find(node_net==node);
    links=link_net(node_r,node_c+1:treedepth);
    links=unique(links);
    links(links==0)=[];
    
    if isempty(links) % if the node is a source node
        link=link_net(node_r,node_c);
        drain_area(node)=ddrain_area(link);
        links_node{node}=link_net(node_r(1),node_c(1));
        link_node_str{node}=sprintf([repmat('%d,',1,length(links_node{node})-1)...
            ,'%d'],links_node{node}(1:end));
        fprintf(fid,'%s %s\n',num2str(node),link_node_str{node});
    else % if the node is not a source node
        drain_area(node)=sum(ddrain_area(links));
        links_node{node}=links;
        link_node_str{node}=sprintf([repmat('%d,',1,length(links_node{node})-1)...
            ,'%d'],links_node{node}(1:end));
        fprintf(fid,'%s %s\n',num2str(node),link_node_str{node});
    end
    
end
fclose(fid);

% a structure to be passed
GEOMORPH.wt=wt;
GEOMORPH.RTD=RTD;
GEOMORPH.net=net;
GEOMORPH.node_net=node_net;
GEOMORPH.link_net=link_net;
GEOMORPH.strahler_order=strahler_order;
GEOMORPH.ddrain_area=ddrain_area;
GEOMORPH.pathPoints=pathPoints;
GEOMORPH.length_link=length_link;
GEOMORPH.length_num_mat=length_num_mat;
end
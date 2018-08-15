function nodeGeo = geomorphology(wt,length_link)

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


    numPaths=size(wt,1);                % compute # of paths from watershed-topology matrix
    
    pathPoints=cell(numPaths,1);
    for pathind=1:numPaths
        pathPoints{pathind}=...         % link indices in each path
            find(wt(pathind,:)==1);
    end
    pathPoints=pathPoints';

    linkPaths=cell(size(wt,2),1);
    for linkind=1:size(wt,2)
        linkPaths{linkind}=...         % paths common to each link
            find(wt(:,linkind)==1);
    end
    linkPaths=linkPaths';

    % tree structure
    s=sum(wt);                         % number of paths shared by each link
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

    % node network
    %
    treedepth=size(net,2);
    node=0;
    [~,~,idx2]=uniquecell(net);
    idx2_res=reshape(idx2,size(net,1),size(net,2));
    
    for i=1:treedepth
        un=unique(idx2_res(:,i));
        for j=1:length(un)
            if un(j)~=1
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
    %{
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
    %% link network
    %
    temp_net=net(:,2:end);
    temp_wt=wt;
    temp_treedepth=size(temp_net,2);
    
    link_net=zeros(numPaths,temp_treedepth);
    for i=1:temp_treedepth
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
    %% define strahler order of streams
    %
    link_net_depth=size(link_net,2);
    ind=find(s==1);
    strahler_order(ind)=1;
    for i=link_net_depth-1:-1:1
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
    strahler_order=strahler_order';
    %}
    %% mean length of the streams of each strahler-order in entire network
    %
    max_order=max(strahler_order);
    length_num_mat=[(1:max_order)',...
        zeros(max_order,1),zeros(max_order,1)];

    for order=1:max_order
        links=find(strahler_order==order);                        % links of order 'order'
        
        for path=1:numPaths
            link_path_int{path}=MY_intersect(links,pathPoints{path});
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
    %}
    
    nodeGeo{1}=net;
    nodeGeo{2}=link_net;
    nodeGeo{3}=strahler_order;
    nodeGeo{4}=length_num_mat;
end


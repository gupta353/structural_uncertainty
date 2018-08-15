% Computation of Jacobian of a fucntion given that the function is defined
% only in positive quadrant
% think about precision  problem and optimal step-size
% inputs: fun=funcional handle of which jacobian is to be computed
%         a=variable value at which the jacobian is to be computed
%         h=step-size
% output: J=Jacobian of 'fun' at 'a' 

function J=jacobiancomp(fun,a,h)

h0=h;   
nx=length(a);
hfrac=0.5;
prec_thresh=10^(-16);       % floating-point precision threshold
count_thresh=100;           % threshold for number of time step-size is adjusted

today=date;
fid=fopen(fullfile(['D:\Research\Thesis_work\Structural_'...
    'uncertainty\MatLab_codes\20180222\'...
    'jacobiancomp_func'],strcat('report_',today,'.txt')),'a');

% determination of an appropriate step-size
for variter=1:nx
    xc=[a(variter)-2*h,a(variter)-h,a(variter)+h,a(variter)+2*h];
    count=0;
    while sum(xc>0)<length(xc) && count<count_thresh
        h=h*hfrac;
        xc=[a(variter)-2*h,a(variter)-h,a(variter)+h,a(variter)+2*h];
        count=count+1;
    end
    par{variter}=xc;
end

% choice of differentiation method: two-point or four-point
for variter=1:nx
    
    xc=par{variter};
    funp=@(x)fun(a-[zeros(1,variter-1)...
    ,a(variter)-x,zeros(1,length(a)-variter)]);

    if sum(xc>0)==length(xc)&& h>prec_thresh
        J(variter,:)=numderfourp(funp,par{variter});
        fprintf(fid,'%s%s %s%s %s\n','variable=',...
            num2str(variter),'h=',num2str(h),'four-point method');
    else
        h=h0; % replace with optimal value
        xc=[a(variter),a(variter)+h];
        J(variter,:)=numdertwop(funp,xc);
        fprintf(fid,'%s%s %s%s %s\n','variable=',...
            num2str(variter),'h=',num2str(h),'two-point method');
    end
end
fclose(fid);
end

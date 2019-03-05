function [y_1,y_2,y_3,y_or]=MWHA(yi, xi, x, ylu, nf, dm, HiLo, thr, fet)
%MWHA
%
% INPUTS
% xi  : rowvec, 已有节点坐标
% x   : rowvec, 待插值（滤波）的节点坐标
% nf  : 谐波个数（默认值为1）
% dm  : 支持域半径（默认值为100day/2, for SPOT 10-day NDVI is 5, for MODIS 16-day NDVI is 3）
% HiLo: indicating rejection of -1: low, 1: high outliers. 
% thr : 判断低值点的阈值(NDVI设为0.1)
% fet : 拟合误差(0.01)
% high: 数据值域最高值(NDVI=1）
% low : 数据值域最低值(NDVI=-1)
% yi  : colvec, 单像素时间序列（列向量）
%
% OUTPUTS
% y_1:  iteration 1
% y_2:  iteration 2
% y_3:  iteration 3
% y_or: final adjusted smoothed value
%
% Author
% Yang, G. (2015).
%
% References
% Yang, G., Shen, H., Zhang, L., He, Z., Li, X., 2015. 
%   A Moving Weighted Harmonic Analysis Method for Reconstructing High-Quality 
%   SPOT VEGETATION NDVI Time-Series Data. IEEE Trans. Geosci. Remote Sens. 53. 
%   https://doi.org/10.1109/TGRS.2015.2431315
low  = ylu(1);
high = ylu(2);

nnodes=length(xi);
num_0=find(yi,1);

HiLo = -HiLo; % keep consistent with HANTS

if isempty(num_0)
    y_or=zeros(nnodes,1);
else
    p_yi=ones(1,nnodes);
    diff_yi=diff(yi);
    dm2=10;
    
    if HiLo == 1
    % if (strcmp(HiLo,'Hi'))
        lo_py1=find(diff_yi<=-thr);
        lo_py2= diff_yi>=thr;
        p_yi(lo_py1+1)=0;
        p_yi(lo_py2)=0;
    elseif HiLo == -1 %(strcmp(HiLo,'Lo'))
        lo_py1=find(diff_yi>=thr);
        lo_py2= diff_yi<=-thr;
        p_yi(lo_py1+1)=0;
        p_yi(lo_py2)=0;
    end
    
    lo=find(p_yi==0);
    x1=x(lo);
    
    %%
    npoints=length(x1);
    
    V=zeros(2*nf+1,npoints);
    w=zeros(npoints,nnodes);
    p=ones(2*nf+1,nnodes);
    px=ones(2*nf+1,npoints);
    
    noutmax=nnodes-2*nf-2;

    ang=(2.*pi*(1:nf)/nnodes)';
    cs1=cos(ang*xi);
    sn1=sin(ang*xi);
    
    cs2=cos(ang*x1);
    sn2=sin(ang*x1);
    
    for m=1:nf
        p(2*m,:)=cs1(m,:);
        p(2*m+1,:)=sn1(m,:);
        px(2*m,:)=cs2(m,:);
        px(2*m+1,:)=sn2(m,:);
    end
    
    d1=ones(nnodes,npoints)*diag(x1);
    d2=ones(npoints,nnodes)*diag(xi);
    di=transpose(d1)-d2;
    
    p_num=ones(1,nnodes);
    p_num(yi<low | yi>high)=0;
    p_num(yi==0)=0;
    
    for j = 1 : npoints
        dmi=dm2;
        nout=nnodes;
        while  (nout>noutmax&&dmi<5*dm2)
            r=abs(di(j,:))/dmi;
            lo1= r>1.0;
            lo2=find(r<=0.5);
            lo3=find(r>0.5&r<=1);
            
            wj(lo1)=0;
            wj(lo2)=2/3 - 4*r(lo2).^2 + 4*r(lo2).^3;
            wj(lo3)= 4/3 - 4*r(lo3) + 4*r(lo3).^2 - 4*r(lo3).^3/3;
            
            w(j,:)=wj.*p_num.*p_yi;
            nout=sum(w(j,:)==0);
            dmi=dmi+1;
        end
        
        if (nout>noutmax)
            % disp('Not enough data points')
            return
        end
        
        B=p*diag(w(j,:));
        
        A=B*p';
        v=A\B*yi;
        V(:,j)=v;
    end
    YY=V'*px;
    Y_mw=diag(YY);
    yi(lo)=Y_mw;
    
    % %%
    npoints=length(x);
    %  Y=zeros(npoints,1);
    V=zeros(2*nf+1,npoints);
    %  Y_mw=zeros(npoints,1);
    
    w=zeros(npoints,nnodes);
    p=ones(2*nf+1,nnodes);
    px=ones(2*nf+1,npoints);
    
    noutmax=nnodes-2*nf-2;
        
    ang=(2.*pi*(1:nf)/nnodes)';
    cs1=cos(ang*xi);
    sn1=sin(ang*xi);
    
    cs2=cos(ang*x);
    sn2=sin(ang*x);
    
    for m=1:nf
        p(2*m,:)=cs1(m,:);
        p(2*m+1,:)=sn1(m,:);
        px(2*m,:)=cs2(m,:);
        px(2*m+1,:)=sn2(m,:);
    end
    
    d1=ones(nnodes,npoints)*diag(x);
    d2=ones(npoints,nnodes)*diag(xi);
    di=transpose(d1)-d2;
    
    p_num=ones(1,nnodes);
    p_num(yi<low | yi>high)=0;
    p_num(yi==0)=0;
    
    y_or=yi;
    max_y_err=fet+1;
    
    % ready=~strcmp(HiLo,'none');
    i_num=0;
    while max_y_err>fet&&i_num<1000
        i_num=i_num+1;
        for j = 1 : npoints
            dmi=dm;
            nout=nnodes;
            while  (nout>noutmax&&dmi<5*dm)
                r=abs(di(j,:))/dmi;
                lo1= r>1.0;
                lo2=find(r<=0.5);
                lo3=find(r>0.5&r<=1);
                
                wi(lo1)=0;
                wi(lo2)=2/3 - 4*r(lo2).^2 + 4*r(lo2).^3;
                wi(lo3)= 4/3 - 4*r(lo3) + 4*r(lo3).^2 - 4*r(lo3).^3/3;
                
                w(j,:)=wi.*p_num;
                nout=sum(w(j,:)==0);
                dmi=dmi+1;
            end
            
            if (nout>noutmax)
                % disp('Not enough data points')
                return
            end
            
            B=p*diag(w(j,:));
            A=B*p';
            v=A\B*y_or;
            V(:,j)=v;
        end
        YY=V'*px;
        Y_mw=diag(YY);
        % p_yi=ones(1,nnodes);
        p_num=ones(1,nnodes);
        
        %%
        if i_num==1
            y_1=Y_mw;
        end
        %%
                
        % if (strcmp(HiLo,'Hi'))   
        if HiLo == 1
            y_err=y_or-Y_mw;
            diffVec=p_num'.*y_err;
            lo_diff1=find(diffVec<0);
            y_or(lo_diff1)=Y_mw(lo_diff1);
            
            [~, rankVec]=sort(diffVec,'ascend');
            max_y_err=y_err(rankVec(nnodes));
            
            %    elseif (strcmp(HiLo,'Lo'))
            %        y_err=Y_mw-yi;
            %        p_y_err=y_err<0.05;
            %        diffVec=p_yi'.*p_num'.*p_y_err.*y_err;
            %         lo_diff=find(diffVec>0);
            %         y_or=Y_mw;
            %         y_or(lo_diff)=yi(lo_diff);
            %
            %
            %        [~, rankVec]=sort(diffVec,'ascend');
            %        max_y_err=y_err(rankVec(nnodes));
        else
            max_y_err=fet-1;
            y_or=Y_mw;
        end
        % p_yi=ones(1,nnodes);
        %%
        if i_num==1||i_num==2
            y_2=y_or;
        end
        
        p_num=ones(1,nnodes);
    end
    
    y_3=y_or;
    % % ;
    
    %% 计算三个均值
    yi1=y_1;
    % yi2=(y_2+yi)/2;
    yi2=yi;
    mat_mean1=ones(nnodes,1)*mean(yi);
    mat_mean2=ones(nnodes,1)*mean(yi(yi>mean(yi)));
    mat_mean3=ones(nnodes,1)*mean(yi(yi<mean(yi)));
    
    %% 点同时位于均值和上均值之间的情况
    
    diff11=y_or-mat_mean1;
    diff12=yi-mat_mean1;
    diff13=mat_mean2-y_or;
    diff14=mat_mean2-yi;
    
    p_y_or2=diff11>=0;
    p_yi2=diff12>=0;
    p_y_or3=diff13>=0;
    p_yi3=diff14>=0;
    
    p2=p_y_or2.*p_yi2.*p_y_or3.*p_yi3;
    lo_p2=find(p2==1);
    w_new=(diff11-diff12)./diff11;
    w_or=diff12./diff11;
    y_or(lo_p2)=w_new(lo_p2).*y_or(lo_p2)+w_or(lo_p2).*yi2(lo_p2);
    
    %% 点同时位于上均值以上的情况
    diff21=y_or-mat_mean2;
    diff22=yi-mat_mean2;
    
    p_y_or1=diff21>=0;
    p_yi1=diff22>=0;
    p1=p_y_or1.*p_yi1;
    lo_p1=find(p1==1);
    
    w_new=(diff21-diff22)./diff21;
    w_or=diff22./diff21;
    
    y_or(lo_p1)=w_new(lo_p1).*y_or(lo_p1)+w_or(lo_p1).*yi2(lo_p1);
    
    %% 点同时位于均值和下均值之间的情况
    diff31=y_or-mat_mean3;
    diff32=yi-mat_mean3;
    
    p_y_or4=diff11<=0;
    p_yi4=diff12<=0;
    p_y_or5=diff31>=0;
    p_yi5=diff32>=0;
    p_yi6=diff32<=0;
    
    p3=p_y_or4.*p_yi4.*p_y_or5.*p_yi5;
    lo_p3=find(p3==1);
    w_new=(diff31-diff32)./diff31;
    w_or=diff32./diff31;
    y_or(lo_p3)=w_new(lo_p3).*y_or(lo_p3)+w_or(lo_p3).*yi2(lo_p3);
    
    %% 点跨越上均值两边的情况
    p4=p_y_or1.*p_yi3.*p_yi2;
    lo_p4=find(p4==1);
    w_new=diff21./(diff21+diff14);
    w_or=diff14./(diff21+diff14);
    temp=w_new;
    w_new=max(w_new,w_or);
    w_or=min(temp,w_or);
    y_or(lo_p4)=w_new(lo_p4).*y_or(lo_p4)+w_or(lo_p4).*yi1(lo_p4);
    
    %% 点跨越均值两边的情况
    p5=p_y_or2.*p_yi4.*p_yi5.*p_y_or4;
    lo_p5=find(p5==1);
    w_new=diff11./(diff11+abs(diff12));
    w_or=abs(diff12)./(diff11+abs(diff12));
    temp=w_new;
    w_new=max(w_new,w_or);
    w_or=min(temp,w_or);
    y_or(lo_p5)=w_new(lo_p5).*y_or(lo_p5)+w_or(lo_p5).*yi1(lo_p5);
    
    %% 点跨越下均值两边
    p6=p_y_or5.*p_yi6.*p_y_or4;
    lo_p6=find(p6==1);
    w_new=diff31./(diff31+abs(diff32));
    w_or=abs(diff32)./(diff31+abs(diff32));
    temp=w_new;
    w_new=max(w_new,w_or);
    w_or=min(temp,w_or);
    y_or(lo_p6)=w_new(lo_p6).*y_or(lo_p6)+w_or(lo_p6).*yi1(lo_p6);
end

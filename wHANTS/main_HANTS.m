clear, clc
close all
%%

file_y = 'TSM_CANS6_y.txt';
file_w = 'TSM_CANS6_w.txt';

y = importdata(file_y);
w = importdata(file_w);

n = length(y);

x    = (1:n);
xi   = x;

nf   = 3;
dm   = 7;
HiLo = 1;
thr  = 0.1;
fet  = 0.02;
high = 1.0;
low  = -0.1;
yi   = y';

[y_1, y_2, y_3, y_or]=MWHA(x,xi,nf,dm,HiLo,thr,fet,high,low,yi);

plot([yi, y_1, y_2, y_3, y_or])
legend({'yi', 'y-1', 'y-2', 'y-3', 'y-or'})

ys = table(yi, y_1, y_2, y_3);

%% HANTS
nptperyear = 23;
ylu = [low, high];
noutmax = 10/nptperyear*n;
delta = 0.1;

[yr, amp, phi] = HANTS(yi, x', HiLo, nf, ylu, nptperyear, fet, noutmax, delta);

close all
plot([yi, yr]);

% function [y_1,y_2,y_3,y_or]=MWHA(x,xi,nf,dm,HiLo,thr,fet,high,low,yi)

% xi  : 已有节点坐标
% x   : 待插值（滤波）的节点坐标
% nf  : 谐波个数（默认值为1）
% dm  : 支持域半径（默认值为7）
% HiLo: 'Hi','Lo','none'.
% thr : 判断低值点的阈值（NDVI设为0.1）
% fet : 拟合误差（0.01）
% high: 数据值域最高值（NDVI=1）
% low : 数据值域最低值（NDVI=-1）
% yi  : 单像素时间序列（列向量）

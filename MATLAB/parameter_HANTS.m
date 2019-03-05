% G:/Github/GEE/gee_whittaker/MATLAB/parameter_HANTS

%% PARAMETERS
ylu        = [-0.1, 1.0];
HiLo = -1;

% HANTS
nf_HANTS   = 4; % The number of frequency
nptperyear = 23;

noutmax    = 9/nptperyear*nptperyear; % DOD=5
delta      = 0.5;  % damping factor, 0.5
fet        = 0.02; % fit error tolerance, 0.05 (Zhou 2015)

% WMHA
nf_WMHA    = 1;
dm         = 3;   % D = (2m+1)*16 = 100,
thr        = 0.1; % MWHA, threshold for low value

t = (1:nptperyear);

%% 

file_y = 'TSM_CANS6_y.txt';
file_w = 'TSM_CANS6_w.txt';

y = importdata(file_y)';
w = importdata(file_w)';

n    = length(y);
x    = (1:n);

[y_1, y_2, y_3, y_or]=MWHA(y, x, x, ylu, nf_WMHA, dm, HiLo, thr, fet);

figure
subplot(211)
plot([y, y_1, y_2, y_3, y_or])
legend({'yi', 'y-1', 'y-2', 'y-3', 'y-or'})
ys = table(y, y_1, y_2, y_3);

subplot(212)
[yr, amp, phi] = HANTS(y, x', ylu, nf_HANTS, HiLo, nptperyear, fet, noutmax, delta);
plot([y, yr]);
legend({'y', 'yr'})
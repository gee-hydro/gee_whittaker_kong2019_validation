function res = HANTS_multiYears(file_y)
%% This function only works for TSF input file. 
% USED to validate gee_Whittaker performance
%
% EXAMAPLES:
% HANTS_multiYears TSM_rep11_SG_y.txt;

%% PARAMETERS
ylu        = [-0.1, 1.0];
HiLo       = -1;

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
read_input = @(file) importdata(file, '\t', 1);

ys = read_input(file_y);

info = textscan(ys.textdata{1}, '%d %d %d'); 
info = cell2mat(info);

ys = ys.data;
nyear = info(1);
% nptperyear = info(2);
ngrid = info(3);

%% 
t = (1:nptperyear);

res = nan(nptperyear*nyear, ngrid);

for i = 1:ngrid
    y = ys(i, :)';
    ysim = y*nan;
    for j = 1:nyear
        I = ((j - 1)*nptperyear+1):(j*nptperyear);
        yj = y(I);
        [yr, amp, phi] = HANTS(yj, t', ylu, nf_HANTS, HiLo, nptperyear, fet, noutmax, delta);
        ysim(I) = yr';
    end
    res(:, i) = ysim;
end

outfile = sprintf('HANTS_%s', file_y);
csvwrite(outfile, res)

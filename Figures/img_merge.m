
files = dirR('.\', ['.*sites.', 'jpg']);
files = files([2, 1]);

FUN  = @(fs) cellfun(@(x) x(85:end-50, 180:end-180, :), fs, 'UniformOutput', false);
% FUN  = [];
nrow = 2; ncol = 1; byrow = true;
outfile = ['Fig1_sites_dist', '.tif'];
layout_image(files, outfile, nrow, ncol, FUN, byrow)
% system(outfile)


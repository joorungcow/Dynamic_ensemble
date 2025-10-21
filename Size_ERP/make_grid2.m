function[gridx, gridy] = make_grid2(rect) % input = dp.resolution

xsize = rect(1);
ysize = rect(2);

slotsize = 145; % adjust the value depending on the display resolution

numx = floor(xsize / slotsize); 
numy = floor(ysize / slotsize);

xx = (0 : slotsize : (slotsize * (numx - 1))) + (slotsize/2); % x-centres
xx = xx(4:end-3); % set a margin 

yy = (0 : slotsize : (slotsize * (numy - 1))) + (slotsize/2); % y-centers
yy = yy(2:end-1);

gridx = repmat(xx, [length(yy), 1]);
gridy = repmat(yy', [1, length(xx)]);

fixrow = (size(gridx, 1)+1)/2; % empty slot for the fixation cross
fixcol = (size(gridx, 2)+1)/2;

gridx(fixrow, fixcol) = NaN;
gridy(fixrow, fixcol) = NaN;

end
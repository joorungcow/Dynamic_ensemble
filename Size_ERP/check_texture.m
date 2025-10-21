function[Simg, Timg] = check_texture(SindivR, TindivR, ppd)


%%%%%%%%%%%%%%%%%%%% grain size varies
% %% standard set
% Simg = struct();
% for i = 1:length(SindivR)
%     imsize = SindivR(i)*2; % diameter
%     [x, y] = meshgrid(linspace(-1, 1, imsize), linspace(-1, 1, imsize)); % canvas
%     [~, r] = cart2pol(x, y);
%     check = sign(sin(2*pi*sf*x) .* sin(2*pi*sf*y));
%     check(r > 1) = 0; 
% 
%     Simg(i).img1 = round((check + 1) * 0.5 * 255); % * 255 for PTB
%     Simg(i).img2 = round((-check + 1) * 0.5 * 255); % phase reversed
% end
% 
% %% test set
% Timg = struct();
% for j = 1:length(TindivR)
%     imsize = TindivR(j)*2;
%     [x, y] = meshgrid(linspace(-1, 1, imsize), linspace(-1, 1, imsize)); % canvas
%     [~, r] = cart2pol(x, y);
%     check = sign(sin(2*pi*sf*x) .* sin(2*pi*sf*y));
%     check(r > 1) = 0; 
% 
%     Timg(j).img1 = round((check + 1) * 0.5 * 255);
%     Timg(j).img2 = round((-check + 1) * 0.5 * 255);
% end

%%%%%%%%%%%%%%%%%%%% constant grain size (cpd)
imsize = max([SindivR, TindivR]) * 2;

nCycleperDeg = 3; % cycle per deg
ppc = ppd / nCycleperDeg; % pixel per cycle in 1-deg canvas
thisncycle = imsize / ppc; 
sf = floor(thisncycle)/2; % results sf*2 B/W reps 

% standard set
Simg = struct();
for i = 1:length(SindivR)
    [x, y] = meshgrid(linspace(-1, 1, imsize), linspace(-1, 1, imsize));
    [~, r] = cart2pol(x, y);
    if mod(size(r, 1), 2) == 1 % odd
        cc = floor(size(r, 1)/2) + 1;
    else % even
        cc = floor(size(r, 1)/2);
    end
    thisr = r(cc, cc + floor(SindivR(i))); 
    check = sign(sin(2*pi*sf*x) .* sin(2*pi*sf*y));
    check(r > thisr) = 0;

    Simg(i).img1 = round((check + 1) * 0.5 * 255); % * 255 for PTB
    Simg(i).img2 = round((-check + 1) * 0.5 * 255); % phase reversed
end

% test set
Timg = struct();
for j = 1:length(TindivR)
    [x, y] = meshgrid(linspace(-1, 1, imsize), linspace(-1, 1, imsize));
    [~, r] = cart2pol(x, y);
    if mod(size(r, 1), 2) == 1 % odd
        cc = floor(size(r, 1)/2) + 1;
    else % even
        cc = floor(size(r, 1)/2);
    end
    thisr = r(cc, cc + floor(TindivR(j))); 
    check = sign(sin(2*pi*sf*x) .* sin(2*pi*sf*y));
    check(r > thisr) = 0;

    Timg(j).img1 = round((check + 1) * 0.5 * 255);
    Timg(j).img2 = round((-check + 1) * 0.5 * 255);
end


end
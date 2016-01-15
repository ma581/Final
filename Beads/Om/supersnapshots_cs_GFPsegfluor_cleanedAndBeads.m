
%% Manoj edits
% 1. -g replaced with -r

function supersnapshots_cs_GFPsegfluor_cleanedAndBeads(varargin);

% function analsnapshots(varargin);
% NOTE:
%     the best way to call this is to run it with the following options:
%
%     >> supersnapshots('noflat','naturalback','sfile');
%
% options: 
%     noflat = turn off flatfield corrections.
%     naturalback = take background values from images
%     sfile = save all the S's at the end to sfile.mat
%     edge -- uses an edge filter to clean up the phase image.
%     beads -- goes into "bead mode" -- for pictures of 0.5 um fluorescent beads...

defaultctimestr = '1';
defaultcgainstr = 'med'
defaultcexptime = 1;

beadmode = 0;
useedge = 0;
useflat = 1;
naturalback = 0;
for i = 1:length(varargin),
    theparam = lower(varargin{i});
    switch(strtok(theparam)),
    case 'edge',
        disp('edge turned on');
        useedge = 1;
    case 'noflat',
        disp('flatfield turned off');
        useflat = 0;
    case 'naturalback',
        naturalback = 1;
        disp('natural backgrounds');
    case 'beads'
        beadmode = 1;
        disp ('bead mode on');                    
    end;
end;

colorcode = 'cygrx';
Ncolors = length(colorcode);

shift = zeros(5,2);

iptsetpref('imshowborder','tight');
erodesize = 5;

% the qualfunc is a filter function that filters out which cells
% count and which get thrown away.
if ~exist('qual'),    
    qualfunc = @qual;
else
    qualfunc = [];
end;

%set these params before using:
suffix = input(['choose your suffix [default = -1]'],'s');
if isempty(suffix),
    suffix = '-1';
end;

prefixes = findprefixes;

for i = 1:length(prefixes),
    x = input(['do you want to do this one: [Y/n] ',char(prefixes(i))],'s');
    if isempty(x),
        x = 'Y';
    end;
    if upper(x) == 'N',
        doit(i) = 0;
    else
        doit(i) = 1;
    end;
end;
prefixes = prefixes(logical(doit));
prefixes

ucolors = [];
colorused = zeros(Ncolors,1);
for i = 1:Ncolors,
    D = dir(['*-',colorcode(i),'.tif']);
    
    if length(D)>0, 
        ucolors = cat(2,ucolors, colorcode(i));
        colorused(i) = 1;
    end;    
end;
ucolors,

if useflat,
    flat = [];
    for i = 1:length(ucolors),
        flatname = ['flat-',ucolors(i),'.tif'];
        if exist(flatname),
            flat = imread(flatname);
        else,
            disp(['can''t find ',flatname,'...quitting']);
        end;
        flatbackname = ['flatback-',ucolors(i),'.tif'];
        if exist(flatbackname),
            flatbac = imread(flatbackname);
        else,
            disp(['can''t find ',flatbackname,'...quitting']);
        end;
        flat = imsubtract(flat,flatbac);
        myflat{i} = double(medfilt2(flat,[5 5]));
    end;    
end;

border = 10;

%%
%Giant for loop
L=[];
for p = 1:length(prefixes),
    
    myname = char(prefixes(p));
    mynamesimple = myname;
    mynamesimple(find(mynamesimple=='-'))=[];
    mynamesimple(find(mynamesimple=='+'))=[];
    mynamesimple(find(mynamesimple=='='))=[];  
    
    L=[];
    D = dir([myname,'-*-r.tif']);
    for i = 1:length(D),
        L(i) = length(char(D(i).name));
    end;
    L0 = length(myname)+9;
    D = D(find(L==L0));
    
    uucolors = [];
    for i = 1:length(ucolors),
        D1 = D(1).name;
        
        f = findstr(D1,'-r.tif');
        
        D1(f(1)+1) = ucolors(i);
        if exist(D1),
            uucolors = cat(2,uucolors,ucolors(i));
            disp(['using color ',ucolors(i)]);
        end;
    end;
    exptime = [];
    for i = 1:length(uucolors),
        
        [Timestr,Gainstr,Exptime,Cube] = imsettings(D(1).name,uucolors(i));
        if isempty(Timestr),
            disp('using default values for c...')
            Timestr = defaultctimestr;
            Gainstr = defaultcgainstr;
            Exptime = defaultcexptime;
        end;
        
        %         timestr(i) = cellstr(Timestr);
        %         gainstr(i) = cellstr(Gainstr);
        %         exptime(i) = exptime; 
        %         cubecolor(i) = colorcode(Cube+1);
        
        if naturalback==0,
            backimg = ['back-',Timestr,'-',Gainstr,'-',uucolors(i),'.tif']
            %         defaultback = cellstr(['back-default-',uucolors(i),'.tif']);                    
            defaultback = ['back-default-',uucolors(i),'.tif'];
            if ~exist(backimg) & exist(defaultback);
                backimg = defaultback;
                disp(['Warning: we''re using the default back image:', backimg]);
            end;
            back = imread(backimg);
            myback{i} = mean2(double(back));	
        end;
        exptime(i) = Exptime;
    end;
    
    
    D(1).name;
    
    disp(['Starting ',myname])
    
    
    %	Nimgs = length(dir([myname,'-*-r.tif']))	
    Nimgs = length(D);
      
    S = [];
    
    for i = 1:Nimgs,
        disp(['Working on ' myname '-' num2str(i)]);
       %padding with zeros
        xx = num2str(i);
        if length(xx) < 2,
            xx = ['0',xx];
        end;
        xxx = num2str(i,'%03d');
        
        
        pname = ['',myname,'-',xx,'-r.tif'];
                
        if exist(pname),
            imp = imread(pname);
            d = dir(['',myname,'-',xx,'-r.tif']);
            %[dstr,tstr,dn] = imdatetime(pname);
            imf = [];
            for j = 1:length(uucolors),
                f = find(colorcode==uucolors(j));
                
                im = double(imshift(imread(['',myname,'-',xx,'-',uucolors(j),'.tif']),shift(f,:)));
                im=imresize(im,[1040,1396],'bilinear'); %changed by CPS 27.10.14 
                    
                if naturalback,
                    % 				ims = sort(im(:));
                    % 				imb = im-ims(floor(0.01*length(ims)));
                    %                 disp('sliding...');
                    %                 im2 = colfilt(im, [80 80], 'sliding', @localback);
                    imb = im;
                else
                    imb = (im-myback{j});
                end;
                
                if useflat,
                    imf{j} = imb./myflat{j};
                else,
                    imf{j} = imb;
                end;
            end; %not used
            
%==========================================================================          
            [L,maskThr] = segfluorToSnapScript(imp,beadmode);             
%========================================================================== 
             
            [s,Lq,Lo] = collectstats(imf,L,beadmode,qualfunc);

            if ~isempty(s),
                eval(['Lq',mynamesimple,'',xx,'=Lq;']);
                if naturalback,
                    for i = 1:length(s),
                        lowbox = floor(max([ones(1,2); s(i).Centroid - 64]));
                        highbox = floor(min([size(imp'); s(i).Centroid + 64]));
                        for kk = 1:length(imf),
                            box = imf{kk}(lowbox(2):highbox(2),lowbox(1):highbox(1));
                            sbox = sort(box(:));
                            bval = sbox(floor(length(sbox)*.005));
                            s(i).bval(kk) = bval;
                        end;
                    end;
                end;
                
                for i = 1:length(s),
                    %s(i).date = dn;
                    s(i).fname = pname;
                    s(i).exptime = exptime;
                end;
                
                if ~isempty(s),
                    if isempty(S),
                        S = s;
                        %			eval(['S = s',mynamesimple,'',xx,';']);
                    else
                        %			eval(['S = cat(2,S,s',mynamesimple,'',xx,');']);
                        S = cat(2,S,s);
                    end;
                end;
                imshowlabel(Lq);%,maskThr);
            end;
        else
            disp(['skipping ',pname]);
        end;
    end;
    
    eval(['S',mynamesimple,'=S;']);
    if ~isempty(S),
        eval(['save S',mynamesimple,suffix,' exptime S',mynamesimple,' Lq',mynamesimple,'*']);
    else
        eval(['save S',mynamesimple,suffix,' exptime ']);
    end;
    eval(['clear S',mynamesimple,' Lq',mynamesimple,'*']);
end;

% finally, create the S file, if specified:

clear S*;

D = dir(['S*',suffix,'.mat']);
for i = 1:length(D), load(D(i).name,'S*'); end;
save('sfile.mat','S*');
%%
%everythign else

function [mystats,Lq,Lo] = collectstats(fimgs,L,beadmode,qualfunc);

if nargin >=3,
    if ~isempty(qualfunc),
        ourqualfunc = qualfunc;
    else
        if beadmode
            ourqualfunc = @qualbeads;
        else
            ourqualfunc = @qual;
        end
    end
else
	ourqualfunc = @qual;
end;

sz = size(L);

Lo = L;
n = max(max(L)) % num objects
if n>0,
    %stats = imfeature(L,'MinorAxisLength','MajorAxisLength','Eccentricity','Area','Centroid','ConvexArea','EulerNumber','Solidity',4);
    %     stats = imfeature(L,'Area','Centroid',4);
    stats = regionprops(L,'MinorAxisLength','MajorAxisLength','Eccentricity','Area','Centroid','ConvexArea','EulerNumber','Solidity'); %changed by CPS 27.10.14 
    qualifies = feval(ourqualfunc,stats);
    Lq = zeros(size(L));
else
    disp(['no cells found on this img']);
    qualifies = [];
    Lq = [];
    mystats = [];
    return;
end;

k = 1;
for i = 1:n,
    if qualifies(i)
        %		p = cat(1,stats(i).PixelList(:,1),stats(i).PixelList(:,2));
        %		if (~ismember(1,p)) & (~ismember(sz(1),p))
        
        Lq(find(Lo==i))=i;
        
        ff = find(L==i);
        
        for cnum = 1:length(fimgs),
            fimg = fimgs(cnum);
            fimg = fimg{1};
            newstruc.mean(cnum) = mean(fimg(ff));
            newstruc.max(cnum) = max(fimg(ff));
            newstruc.min(cnum) = min(fimg(ff));
            newstruc.median(cnum) = median(double(fimg(ff)));
            newstruc.total(cnum) = sum(fimg(ff));
            newstruc.std(cnum) = std(fimg(ff));
            s = sort(fimg(ff));
            %             newstruc.melowmeas(cnum) = mean(s(end-2:-1:end-9));
            
        end;
        
        %newstruc(k).stats = stats(i);
        % copy stats:
        statfields = fieldnames(stats(i));
        for g = 1:length(statfields),
            fn = char(statfields(g));
            fv = getfield(stats(i),fn);
            newstruc = setfield(newstruc,fn,fv);
        end;
        
        mystats(k) = newstruc;
        
        k = k + 1;
    end;
end;
if all(qualifies==0),
    disp(['no cells found on this img']);
    qualifies = [];
    Lq = [];
    mystats = [];
    return;
end;
k

function b = qualbeads(stats);
if length(stats) == 0,
    b = [];
else,
    b = ones(length(stats),1);
%     b = b & ([stats.MinorAxisLength]' <25);
%     b = b & ([stats.Area]' <= 25);
%     b = b & ([stats.MajorAxisLength]' <25);
%     b = b & ([stats.EulerNumber]' == 1); % no holes.
%     b = b & ([stats.Eccentricity]' < 0.3);
end;

function b = qual(stats);
if length(stats) == 0,
    b = [];
else,
    b = ones(length(stats),1);
    b = b & ([stats.MinorAxisLength]' <25);
    b = b & ([stats.Area]' >= 75);
    b = b & ([stats.MajorAxisLength]' >18);
    b = b & ([stats.EulerNumber]' == 1); % no holes.
    b = b & ([stats.Eccentricity]' > 0.3);
end;
% for i = 1:length(b),
% 	if ~b(i),
% 		stats(i);
% 	end;
% end;

function [exptimestr, gainstr,exptime,cube] = imsettings(pname, color);

pos = findstr('-r', pname);
pname(pos+1) = color;
pname,

if exist(pname),

    iminfo = imfinfo(pname);

    if isfield(iminfo,'ImageDescription'),

        descrip = iminfo.ImageDescription
    else
        descrip = [];
    end;
else
    exptimestr = '';
    gainstr = '';
    exptime = 0;
    cube=-1;
    return;
end;

if isempty(descrip),
    exptimestr = '';
    gainstr = '';
    exptime = 0;
    cube=-1;
    return;
end;


exptimepos = findstr('Exptime=',descrip) + length('Exptime=');
exptime = sscanf(descrip(exptimepos:end),'%f');
exptimestr = num2str(exptime);
cubestrpos = findstr('Cube=',descrip) + length('Cube=');
cube = sscanf(descrip(cubestrpos:end),'%d');

% if exptime < 1,
%     exptimestr = ['0',exptimestr];
% end;
exptimestr(exptimestr=='.')=[];

gainpos = findstr('Gain=',descrip) + length('Gain=');
gain = sscanf(descrip(gainpos:end), '%f');
switch(gain),
    case 1,
        gainstr = 'low';
    case 2,
        gainstr = 'med';
    case 3,
        gainstr = 'high';
    otherwise,
        disp('can''t find gain setting -- using high');
        gainstr = 'high';
end;
%--------------------------------------------------------------------------
% ... and ends.
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% New function starts ...
% It returns the characteristics of the image 'name'
%--------------------------------------------------------------------------
function [datestr,timestr,daten] = imdatetime(name);

if ~exist(name),
    datestr = '';
    timestr = '';
    daten = 0;
    disp('Error: can''t find file -- exiting');
    return;
end;

iminfo = imfinfo(name);
if ~isfield(iminfo,'ImageDescription'),
    datestr = '';
    timestr = '';
    daten = 0;
    disp('Error: No description information for date');
    return;
end;
descrip = iminfo.ImageDescription;
keystring = 'Acquired at';
f = findstr(keystring,descrip);
if isempty(f),
    datestr = '';
    timestr = '';
    daten = 0;
    disp(['Error: Couldn''t find keystring=',keystring]);
    return;
end;

d = descrip(f(1):end);
A = sscanf(d,[keystring,'%d-%d-%d %d:%d:%d']);
datestr = [num2str(A(1)),'-',num2str(A(2)),'-',num2str(A(3))];
timestr = [num2str(A(4)),':',num2str(A(5)),':',num2str(A(6))];
daten = datenum(A(1),A(2),A(3),A(4),A(5),A(6));

%--------------------------------------------------------------------------
% ... and ends
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% New function...
%--------------------------------------------------------------------------
function prefixes = findprefixes;
D = dir('*-r.tif');
prefixes = [];
for i = 1:length(D),
    name = D(i).name;
    L(i) = length(name);
    p = findstr('-r.tif',name);
    %       disp(name(1:p-4));
    prefixes{i} = name(1:p-4);
end;

delme = [];
prefixes = unique(prefixes');
for i = 2:length(prefixes),
    f = strmatch(upper(prefixes(i)),upper(prefixes(1:i-1)),'exact');
    if ~isempty(f),
        delme = [delme i];
    end;
end;
if ~isempty(delme),
    disp(['found case duplicates -- removing: ']);
    prefixes(delme),
    prefixes(delme) = [];
end;


j = 0;

for i = 1:length(prefixes),
    D = dir([char(prefixes(i)),'-*-r.tif']);
    L(i) = length(D);
    if L(i) > 1,
        j = j + 1;
        prefixes2(j) = prefixes(i);
    end;
end;

%    prefixes = prefixes2;
%--------------------------------------------------------------------------
% The 'findprefixes' function has ended.
%--------------------------------------------------------------------------

function y = minmin(x);
y = min(min(x));

function y = max2(x);
y = max(max(double(x)));

function y = maxmax(x);
y = max(max(x));


function dL = carefuldilate(L,blob,niter,refimage);
% function dL = carefuldilate(L, blob, niter, refimage)
%
% dilates L by strel blob niter times
% if regfimage is included, points are only added that exist in refimage

if nargin < 4,
    refimage = ones(size(L));
end;


for i = 1:niter,

    LD = imdilate(L,blob);
    BW0 = (L>0);
    BW1 = imdilate(BW0,blob);
    NewPix = BW1 & ~BW0 & refimage;

    L(NewPix) = LD(NewPix);

end;

dL = L;


function outim = imshowlabel(L,bwscreen);

if islogical(L),
    L = double(L);
end;

L = mod(L,256);
M = maxmax(L)+1;

mymap = melowhsv(M);
mymap2 = mymap;
[s,I] = sort(rand(M-1,1));
mymap2(2:end,:) = mymap2(I+1,:);

if nargin>=2,
    rgb = 0.5 * ind2rgb(L,mymap2);
    bwscreen = double(bwscreen);
    bwscreen = 0.5 * bwscreen / maxmax(bwscreen);
    rgb(:,:,1) = rgb(:,:,1) + bwscreen;
    rgb(:,:,2) = rgb(:,:,2) + bwscreen;
    rgb(:,:,3) = rgb(:,:,3) + bwscreen;
    figure(1); %changed by CPS 27.10.14 
    imshow(rgb);
    outim = rgb;
else
    imshow(L+1, mymap2);
    if nargout == 1,
        outim = L;
    end;
end;

function map = melowhsv(m)
%HSV    Hue-saturation-value color map.
%   HSV(M) returns an M-by-3 matrix containing an HSV colormap.
%   HSV, by itself, is the same length as the current colormap.
%
%   An HSV colormap varies the hue component of the hue-saturation-value
%   color model.  The colors begin with red, pass through yellow, green,
%   cyan, blue, magenta, and return to red.  The map is particularly
%   useful for displaying periodic functions.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(hsv)
%
%   See also GRAY, HOT, COOL, BONE, COPPER, PINK, FLAG, PRISM, JET,
%   COLORMAP, RGBPLOT, HSV2RGB, RGB2HSV.

%   See Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78.
%   C. B. Moler, 8-17-86, 5-10-91, 8-19-92, 2-19-93.
%   Copyright 1984-2000 The MathWorks, Inc.
%   $Revision: 5.7 $  $Date: 2000/06/01 02:53:43 $

if nargin < 1, m = size(get(gcf,'colormap'),1); end
h = (0:m-1)'/max(m,1);
if isempty(h)
    map = [];
else
    map = hsv2rgb([h ones(m,2)]);
end

map(1,:) = 0;

function y = imshift(x,shift);

y = x;
y(:,:,:,:) = minmin(y);

% first shift x:
shx = shift(1);
if shx>0,
    y(1+shx:end,:,:,:) = x(1:end-shx,:,:,:);
elseif shx<0,
    shx = abs(shx);
    y(1:end-shx,:,:,:) = x(1+shx:end,:,:,:);
end;

% then shift y:
shy = shift(2);
if shx~=0,
    x = y;
end;
if shy>0,
    y(:,1+shy:end,:,:) = x(:,1:end-shy,:,:);
elseif shy<0,
    shy = abs(shy);
    y(:,1:end-shy,:,:) = x(:,1+shy:end,:,:);
end;

if max(abs(shift))==0,
    y = x;
end;



function redseg(doit)


%Default values for parameters.
%Useflat specifies whether you have background corrected files.
%Naturalback 
useflat = 0;
naturalback = 0;

%Assumes first image saved in series starts on 1 (eg 'mytif-01-p.tif) 
firstimnum = 1;

%In case images are taken not through Michael's program.
defaultctimestr = '1';
defaultcgainstr = 'med'
defaultcexptime = 1;

%implement quickreg or fixed shift as an option here
%Each row represents registration for separate channel
shift = zeros(5,2); %No registration

global doreg %Currently we have no clue what this does

%Going through user arguments
%for i = 1:length(varargin),
%    theparam = lower(varargin{i});
%   switch(strtok(theparam)),
%       case 'flat',
%            disp('flatfield turned on');
%            useflat = 1;
%        case 'naturalback',
%            naturalback = 1;
%            disp('natural backgrounds');
%    end;
%end

naturalback=1;


colorcode = 'cygrx'; %All possible colors, notably 't' from Olympus scopes is not included - YET
Ncolors = length(colorcode); %Ncolors is used later
shift = zeros(5,2); %Used for image registration later (see fx imshift)

% This is a setting which changes how images are displayed.
iptsetpref('imshowborder','tight'); %What this changes is unlcear....
erodesize = 5;

%Why?
if exist('qual'),    
    qualfunc = @qual; %Used in fx collectstats
else
    qualfunc = [];
end;

%suffix is used during saving
%suffix = input(['choose your suffix [default = -1]'],'s');
%if isempty(suffix),
%    suffix = '-1';
%end;
suffix='-1';

prefixes = findprefixes; %prefixes is an array of filename strings e.g. EthoNacl001

%for i = 1:length(prefixes),
%    doit(i) = 1;
%end

%User specifes which prefixes to segment here
%for i = 1:length(prefixes),
 %   x = input(['do you want to do this one: [Y/n] ',char(prefixes(i))],'s');
 %   if isempty(x),
 %       x = 'Y';
 %   end;
 %   if upper(x) == 'N',
 %       doit(i) = 0;
 %   else
 %       doit(i) = 1;
 %   end;
%end;
prefixes = prefixes(logical(doit)); %turns doit into an logical array and then only keepts entries with 1
prefixes

ucolors = [];
colorused = zeros(Ncolors,1);

%Specifies which channels are present in snapshots
for i = 1:Ncolors,
        D = dir(['*-',colorcode(i),'.tif']);
    if length(D)>0, 
        ucolors = cat(2,ucolors, colorcode(i)); % returns a string with used colors e.g. yr
        colorused(i) = 1; % anaotes positions in the matrix which where there is a hit with 1
    end;    
end;
ucolors,

%----------------------------------------------------------------------------------
%flat and flatback file need for this. So far I have not seen these files
%--------------------------------------------------------------------------

%Allows for user flatfield and background substraction
%'flat' is a user created fit (eg parabola) of image flatness in each channel
%eg flat-y.tif
%'flatback' is a empty image in each channel (no cells)
if useflat,
    flat = [];
    for i = 1:length(ucolors),
        flatname = ['flat-',ucolors(i),'.tif'];
        if exist(flatname),
            flat = imread(flatname);
        else
            disp(['can''t find ',flatname,'...quitting']);
        end;
        flatbackname = ['flatback-',ucolors(i),'.tif'];
        if exist(flatbackname),
            flatbac = imread(flatbackname);
        else
            disp(['can''t find ',flatbackname,'...quitting']);
        end;
        flat = imsubtract(flat,flatback);  %subtraction of image
        myflat{i} = double(medfilt2(flat,[5 5])); % calc of median value to reduce noise
    end;    
end;



%Giant for loop
for p = 1:length(prefixes), % over the number of different set pads
    L=[];
    
%removing +,' ',+,=,. from file name
    myname = char(prefixes(p));
    mynamesimple = myname;
    mynamesimple(find(mynamesimple=='-'))=[];
    mynamesimple(find(mynamesimple==' '))=[];
    mynamesimple(find(mynamesimple=='+'))=[];
    mynamesimple(find(mynamesimple=='='))=[];
    mynamesimple(find(mynamesimple=='.'))=[];

    D = dir([myname,'-*-r.tif']);
    
%Get the length of the Files which comply with dir 
    for i = 1:length(D),
        L(i) = length(char(D(i).name));
    end;
% Checks if the naming is correct. Length of prefixes + 9 hast to be equal
%to the length of of the files found in dir 
    L0 = length(myname)+9;
     D = D(find(L==L0));

    uucolors = [];
    %Checking file integrity of saved files, again.
    for i = 1:length(ucolors),
        D1 = D(1).name;
        f = findstr(D1,'-r.tif'); %position at which '-r.tif' starts
        D1(f(1)+1) = ucolors(i);  %replaces position f with ucolors eg. r with y -> -y.tif

        %Generating a string with used or existing filenames if the filename
        %exists.
        if exist(D1),
            uucolors = [uucolors,ucolors(i)];
            disp(['using color ',ucolors(i)]);
        end;
    end;
% Asigning value to Timestr, Gainstr, Exptime if they are empty
    exptime = [];
    for i = 1:length(uucolors),
        [Timestr,Gainstr,Exptime,Cube] = imsettings(D(1).name,uucolors(i)); %should get imagesettings but does nothing really
        if isempty(Timestr),
            disp('using default values for c...')
            Timestr = defaultctimestr;
            Gainstr = defaultcgainstr;
            Exptime = defaultcexptime;
        end;
    end

    disp(['Starting ',myname]) %Displays which prefix is currently running
    Nimgs = length(D);
    S = [];

    %Adding leading zeros to image numbers.
    for i = 1:Nimgs,
        disp(['Working on ' myname '-' num2str(i)]);
        xx = num2str(i+firstimnum-1);
        if length(xx) < 2,
            xx = ['0',xx];
        end;


        pname = ['',myname,'-',xx,'-r.tif'];
        pname2= ['',myname,'-',xx]
        if exist(pname),
            imp = imread(pname);
            psize = size(imp);
            d = dir(['',myname,'-',xx,'-r.tif']);
%             [dstr,tstr,dn] = imdatetime(pname); %CHECK:  Do we really need imdatetime? [Jon]
            imf = [];
            ph3 = imp;

            %Rescales image values to the 25th biggest value, maybe not
            %necesary
            for i= 1:size(ph3,3),
                ph3 = ph3;
                ph3(:,:,i)= medfilt2(ph3(:,:,i),[3 3]); %operates on image
                x= double(ph3(:,:,i));
                s= sort(x(:)); %reduces matrix x to a vector
                small= s(25);  %returns the 25th smalles value
                %ALT CODE:  Subtract the background prior to rescaling
                %small1000 = s(1000);
                %ph3(ph3<small1000)=0;
                big= s(end-25); %returns the 25th biggest value
                rescaled=(x - small)/(big - small);%Thresholding
                rescaled(rescaled<0)= 0; %setting values smaller than 0 to 0
                ph3(:,:,i)= uint16(10000*rescaled);%rescaling pixels
            end;
%------------------------------------------------------------------------------
%bis hier her 25.11.2013
%--------------------------------------------------------------------------------

            for j = 1:length(uucolors),
                f = find(colorcode==uucolors(j));

                %Registers each channel image before reading it in as 'im'
                %manipulates image with out any effect
                im = double(imshift(imread(['',myname,'-',xx,'-',uucolors(j),'.tif']),shift(f,:))); %shifts something but has now function as shift is 0
                im = imresize(im,psize,'bilinear'); %why??? resizes images of different colors to same size
                if useflat,
                    imf{j} = im./myflat{j};
                else
                    imf{j} = im;
                end;
            end;

            T=progthreshfluor(ph3); %Main segmentation function
            T = imclearborder(T,4); %Removes cells near border

            %Smooth cell edges by eroding and dilating
            t1= imerode(T, strel('disk',4));
            t2= bwlabel(t1);
            L= imdilate(t2, strel('disk', 4));

            %Goes through and collects cell data and stores it into structure 's'.
            [s,Lq,Lo] = collectstats(imf,L,qualfunc,pname2);% imf:flat image, L=segmented image, for now qualfunc=[]

            %Makes colormap image and displays it
            imagesc(makergb(T,Lq,imp));
            drawnow;

            %Draws a box around each cell in image and takes a low value to
            %stores in background value (b.val) in the 's' structure.
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

                if isempty(S),
                    S = s;
                    %			eval(['S = s',mynamesimple,'',xx,';']);
                else
                    %			eval(['S = cat(2,S,s',mynamesimple,'',xx,');']);
                    S = cat(2,S,s);
                end;
            end;
        else
            disp(['skipping ',pname]);
        end;
    end;

    
    %Saving single 'S' files.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUBFUNCTIONS BELOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% This is the key function here.
%--------------------------------------------------------------------------
function progfluor = progthreshfluor(imt)
p = []; %Initialize p

%Parameters below influence segmentation
p.edge_lapofgauss_sigma = 3;
p.minCellArea = 50;
p.minCellLengthConservative= 12; % JCR observed 2005-02-16
p.minCellLength = 20;            % JCR observed 2005-02-16
p.maxCellWidth = 20;             % JCR observed 2005-02-16
p.minNumEdgePixels = 215;        % ME set to 215 in mail to JCR 2005-08-18
p.maxThreshCut = 0.3;
p.maxThreshCut2 = 0.2;
p.maxThresh = 0.1;
p.minThresh = 0.1;
p.imNumber1 = 2;
p.imNumber2 = 1;
p.radius = 5;
p.angThresh = 2.7;

imt = imt - 1200; %Subtract the bottom 1.2% of image
e = edge(imt,'log',0); %transforms the image so that the edges become visible
f = imfill(e,'holes'); %filles holes
L = bwlabel(f); %lables connected cells

%Find and kill small cells
r = regionprops(L,'Area');
flittle = find([r.Area]<50);
disp(['Removing small (crud) cells(', num2str(length(flittle)),').']);
goodones = setdiff(1:max(max(L)), flittle);
bw2 = ismember(L, goodones);
L2 = bwlabel(bw2);
% L2= renumberimage(L2); %Remove for speed?

%Find cells with good Solidity
r2 = regionprops(L2,'Solidity');
fsolid = find([r2.Solidity]>0.2);
bw3 = ismember(L2,fsolid);
L3 = bwlabel(bw3);

%Remove short cells
r= regionprops(L3,'majoraxislength');
fshort= find([r.MajorAxisLength] < p.minCellLengthConservative);
disp(['Removing short cells(', num2str(length(fshort)),').']);
longones = setdiff(1:max(max(L3)), fshort);
bw4 = ismember(L3, longones);
L4 = bwlabel(bw4);

%Remove unfilled cells on edge
r= regionprops(L4, 'area');
fpts= find([r.Area]<200);
fpts = find([r.Area]<100); % JHL Changed 2011-06-26 for B. subtilis in SMS 20% glucose
edgy = setdiff(1:max(max(L4)), fpts);
bw5 = ismember(L4, edgy);
L5= bwlabel(bw5);

L6= L5;
L7=L6;
label = max2(L7);

% BREAKING UP BIG CELLS
% Goes along the THINned cell and looks for places where there is a change 
% in the fluor value - i.e. where there could be a space between cells.
r= regionprops(L6, 'majoraxislength');
fbiggies= find(([r.MajorAxisLength]>50));
disp(['Breaking up big cells(', num2str(length(fbiggies)),').']);
for i = 1:length(fbiggies),
    Lcell= +(L6 == fbiggies(i)); % + converts logical to double
    Lcell(Lcell == 1)= fbiggies(i);
    [cutcell, num]= breakcellfluorsnaps(Lcell, imt(:,:,1), ...
        p.maxThresh, p.minThresh, p.minCellLength);
    L7(L6 == fbiggies(i))= 0;
    % place cutcell
    label= label+length(fbiggies); %max2(L7);
    for j = 1:num,
        L7(cutcell==j)= label+j;
    end;
end;
L7= renumberimage(L7);

% DEKINKING CELLS with low solidity
% Tries again to look for phase extrema along the thin, using another phase
% slice (if available). Then uses "dekinker" to look for sharp angles along
% the bacterial spine (thin) - if those are found, the cell is cut at the
% corner.
% uses PHSUB(:,:,p.imNumber2) (default p.imNumber2=1)
L8= L7;
r= regionprops(L7, 'solidity');
fkinkies= find(([r.Solidity] < 0.85));
disp(['Breaking kinky cells(',num2str(length(fkinkies)),').']);
label = max2(L8);
for i = 1:length(fkinkies)
    %disp([' ',num2str(i),': dekinking cell number ',num2str(fkinkies(i))]);
    Lcell= +(L7 == fkinkies(i)); % + converts logical to double
    Lcell(Lcell == 1)= fkinkies(i);

    [cutcell, num] = breakcellfluorsnaps(Lcell, imt(:,:,1), ...
        p.maxThresh, p.minThresh, p.minCellLength);
    if num == 1
        cutcell= dekinker(Lcell, p.radius, p.minCellLength, p.angThresh);
    end;

    L8(L7 == fkinkies(i))= 0;
    % place cutcell
    for j = 1:num,
        L8(cutcell==j)= label+length(fkinkies);
    end;
end;

L8= renumberimage(L8);

% FINAL IMAGE
disp('Almost done.');

%Final check for small cells:  probably useless
r= regionprops(L8, 'area');
flittles= find([r.Area] < p.minCellArea);
for i= 1:length(flittles),
    % delete cell
    L8(L8 == flittles(i))= 0;
end;
%L8= renumberimage(L8);

%Plump the cells for final mask
L9= carefuldilate(+L8, strel('diamond',1), 1);
progfluor = L9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End Progthreshfluor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mystats,Lq,Lo] = collectstats(fimgs,L,qualfunc,myname); %fimgs: flat colored images normally (rfp and yfp)

if nargin >=4,
    if ~isempty(qualfunc),
        ourqualfunc = qualfunc;
    else,
        ourqualfunc = @qual;
    end;
else
	ourqualfunc = @qual;
end;

sz = size(L); %size of segmented image
Lo = L; %segmented image
n = max(max(L)); % num objects (Largest value in segmented image -> cell number in Image)

if n>0,
    stats = regionprops(L,'MinorAxisLength','MajorAxisLength','Eccentricity','Area','Centroid','ConvexArea','EulerNumber','Solidity'); %stats
    qualifies = feval(ourqualfunc,stats); %creats vector of the size of the number cells in the image with 1s.
    Lq = zeros(size(L)); %creats zero matrix of size L
else
    disp(['no cells found on this img']);
    qualifies = [];
    Lq = [];
    mystats = [];
    return;
end;

%register the images with cells
if(sum(qualifies))
    for i=1:length(fimgs)
        %resize the fluors if binned
        fimgs{i} = imresize(fimgs{i},sz,'bilinear'); %making sure rfp and yfp have the same image size
        [subreg,shift,back]=regsnaps(L(100:424,100:424),fimgs{i},[100 100 424 424],3);
        fimgs{i}=imshift(fimgs{i},-shift);
        disp(['shift' num2str(i) ' = ' num2str(shift(1)) ' ' num2str(shift(2))]);
    end
end

k = 1;
for i = 1:n,
    if qualifies(i) 
        Lq(find(Lo==i))=i;%sets all cells except for the ith cell 0

        ff = find(L==i);%finds cells which are none 0 -> 1 and gives their position a one dim index
        %--------------------------------------------------------------------------------------
        %Extrating information on the fluorscence in the cells eg. rfp and
        %yfp
        %----------------------------------------------------------------------------------------
        for cnum = 1:length(fimgs),
            fimg = fimgs(cnum); % gets fimgs of specific color
            fimg = fimg{1}; %sets fimg to the imge
            newstruc.mean(cnum) = mean(fimg(ff)); % fimg(ff) retruns the pixes vectors for one cell as vector and then mean pixel val is calc
            newstruc.max(cnum) = max(fimg(ff));
            newstruc.min(cnum) = min(fimg(ff));
            newstruc.median(cnum) = median(double(fimg(ff)));
            newstruc.total(cnum) = sum(fimg(ff));
            newstruc.std(cnum) = std(fimg(ff));
            newstruc.cellno=i;
            newstruc.imname=myname;
            
            s = sort(fimg(ff));
            %             newstruc.melowmeas(cnum) = mean(s(end-2:-1:end-9));

        end;

        %newstruc(k).stats = stats(i);
        % copy stats:
        
        %------------------------------------------------------------------------
        %Merging the two structures from stats and the ones calc from
        %fluorscence 
        %---------------------------------------------------------------
        statfields = fieldnames(stats(i));
        for g = 1:length(statfields),
            fn = char(statfields(g));
            fv = getfield(stats(i),fn);
            newstruc = setfield(newstruc,fn,fv);
        end;

        mystats(k) = newstruc;  %putting structure into my stats

        k = k + 1;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsubreg, xshift, xback] = regsnaps(L, imx, rect, deltamax)
% function [xsubreg, xshift, xback] = regsnaps(L, imx, rect, deltamax)
%  calculates translation between phase and fluorescent images
% load images if necessary
if min(size(imx))==1,
    if isempty(findstr('.tif', imx)),
        imx = [imx,'.tif'];
    end;
    if exist(imx)==2
        imx= imread(imx);
    end
end
if min(size(imx))>1
% get subimages
imx1= double(imx(rect(1):rect(3), rect(2):rect(4)));
LL= +(L(1+deltamax:end-deltamax, 1+deltamax:end-deltamax) > 0);
% try all possible translations
for di= -deltamax:deltamax,
    for dj= -deltamax:deltamax,
        imx2= imx1(deltamax+di+1:end-deltamax+di, deltamax+dj+1:end-deltamax+dj);
        xsum(di+deltamax+1,dj+deltamax+1)= sum(imx2(LL > 0));
    end;
end;
% best translation is the one with largest csum (the most white pixels in LL)
[xbesti, xbestj]= find(xsum == max2(xsum));
xbesti= xbesti - deltamax - 1;
xbestj= xbestj - deltamax - 1;
% record translation
xshift= [xbesti(1) xbestj(1)];
% record translated images
xsubreg= imx((rect(1):rect(3))+xbesti(1), (rect(2):rect(4))+xbestj(1));
% calculate background fluorescence
imxb = double(imx);
imxb(rect(1):rect(3), rect(2):rect(4))=0;
imxbvect=imxb(imxb>0);
xback=median(imxbvect);
else
    xsubreg=[];
    xshift=[];
    xback=[];
end;

%%%%%%%%%%%%%%%%%% End Regsnaps

%--------------------------------------------------------------------------
% This function is related to displaying the images as the segmentation
% code runs.
%--------------------------------------------------------------------------
function rgb = makergb(r,g,b);

if nargin<2,
	disp('MelowError: I need at least 2 images here!');
	return;
end;

if size(r,1) == 1,
	% these must be filenames;
	r = imread(r);
	g = imread(g);
	if nargin>=3,
		b = imread(b);
	end;
end;

if nargin==2,
	b = r;
	b(:) = 0;
end;

r = double(r(:,:,1,1));
g = double(g(:,:,1,1));
b = double(b(:,:,1,1));

rgb(:,:,1) = (r-minmin(r)) / (maxmax(r)-minmin(r));
rgb(:,:,2) = (g-minmin(g)) / (maxmax(g)-minmin(g));
rgb(:,:,3) = (b-minmin(b)) / (maxmax(b)-minmin(b));

rgb(isnan(rgb)) = 0;
rgb = uint8(rgb * 255);
%--------------------------------------------------------------------------
%makergb ends here.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% The function imshift is not used as shift =0
%-------------------------------------------------------------------------
function y = imshift(x,shift);
y = x;
y(:,:,:,:) = minmin(y); % why more dimensions than the image has

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

%Qualifies each cell in an image for data gathering
function b = qual(stats);
if length(stats) == 0,
    b = [];
else,
    b = ones(length(stats),1);
end;

%--------------------------------------------------------------------------
% The function 'imsettings' starts...
%--------------------------------------------------------------------------
function [exptimestr, gainstr,exptime,cube] = imsettings(pname, color);

pos = findstr('-p', pname);
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

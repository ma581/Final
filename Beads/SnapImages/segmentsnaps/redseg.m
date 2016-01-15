function redsegbeads(varargin)

%% Beginning stuff
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
for i = 1:length(varargin),
    theparam = lower(varargin{i});
    switch(strtok(theparam)),
        case 'flat',
            disp('flatfield turned on');
            useflat = 1;
        case 'naturalback',
            naturalback = 1;
            disp('natural backgrounds');
    end;
end

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
suffix = input(['choose your suffix [default = -1]'],'s');
if isempty(suffix),
    suffix = '-1';
end;

prefixes = findprefixes; %prefixes is an array of filename strings e.g. EthoNacl001

%User specifes which prefixes to segment here
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

%--------------------------------------------------------------------------
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


%------------------------------------------------------------------------------
%% Giant for loop
%------------------------------------------------------------------------------
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
                im = imresize(im,psize,'bilinear'); %why???
                if useflat,
                    imf{j} = im./myflat{j};
                else
                    imf{j} = im;
                end;
            end;
%---------------------------------------------------------------------------       
        %Main segmentation function  
%---------------------------------------------------------------------------
            [T] = progthreshfluor(ph3); %Main segmentation function
%            T = imclearborder(T,4); %Removes cells near border

            %Smooth cell edges by eroding and dilating
            t1= imerode(T, strel('disk',4));
            t2= bwlabel(t1);
            L= imdilate(t2, strel('disk', 4));

            %Goes through and collects cell data and stores it into structure 's'.
            [s,Lq,Lo] = collectstats(imf,L,qualfunc);% imf:flat image, L=segmented image, for now qualfunc=[]
            
%---------------------------------------------------------------------------       
% Getting rid of low intensity beads- Manoj
%---------------------------------------------------------------------------
%             mean_threshold = 1400;
%             for j = 1:length(s)
%                 mean(j) = s(j).mean(:,2);
%             end
%              mean_usefull = find(mean< mean_threshold);
%             
%             goodbeads = setdiff(1:max(max(L)), mean_usefull);
%             %     C = SETDIFF(A,B) for vectors A and B, returns the values in A that 
%             %     are not in B with no repetitions. C will be sorted.
% 
%             goodbeads_bool = ismember(L, goodbeads);
%             %   LIA = ISMEMBER(A,B) for arrays A and B returns an array of the same
%             %   size as A containing true where the elements of A are in B and false
%             %   otherwise.    
% 
%             L = bwlabel(goodbeads_bool);
%             [s,Lq,Lo] = collectstats(imf,L,qualfunc);
%---------------------------------------------------------------------------

            
               
            
              %Makes colormap image and displays it
%             imagesc(makergb(T,Lq,imp));
%             drawnow;

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
function [varargout] = progthreshfluor(imt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Manoj edit 01/07/15 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
% %  Thresholding, Edge detection and fill
%----------------------------------------------------------------------

imt = imt - 1200; %Subtract the bottom 1.2% of image
e = edge(imt,'log',0); %transforms the image so that the edges become visible
f = imfill(e,'holes'); %filles holes
L = bwlabel(f); %lables connected cells
L4 = L;

%----------------------------------------------------------------------
%Remove unfilled beads on edge 
%----------------------------------------------------------------------
r= regionprops(L4, 'area');
fpts = find([r.Area]<50); 
edgy = setdiff(1:max(max(L4)), fpts);
bw5 = ismember(L4, edgy);
L5= bwlabel(bw5);

%----------------------------------------------------------------------
% Kill beads that are too close 
%----------------------------------------------------------------------
a = imdilate(L5,strel('disk',15));%dilates so close beads merge
f = imfill(a,'holes'); %filles holes
a = bwlabel(f); %labels connected cells

r = regionprops(a, 'area');
flarge = find([r.Area]>1200);%Identifies merged beads
goodbeads = setdiff(1:max(max(a)), flarge);
    %     C = SETDIFF(A,B) for vectors A and B, returns the values in A that 
    %     are not in B with no repetitions. C will be sorted.

goodbeads_bool = ismember(a, goodbeads);
    %   LIA = ISMEMBER(A,B) for arrays A and B returns an array of the same
    %   size as A containing true where the elements of A are in B and false
    %   otherwise.    

L6 = bwlabel(goodbeads_bool);

%----------------------------------------------------------------------
% Kill beads that are small/close to edge 
%----------------------------------------------------------------------
fsmall = find([r.Area]<1090);%Identifies small beads
goodbeads = setdiff(1:max(max(L6)), fsmall);
goodbeads_bool = ismember(L6, goodbeads);
L6 =  bwlabel(goodbeads_bool);

L7 = imerode(L6,strel('disk',15)); %erode back to original bead size
goodbeads_bool = logical(L7);    %convert to logical
goodbeads_bool = double(goodbeads_bool); % convert so we can multiply to image
% L8 = goodbeads_bool.*L7; %get original intensities only where we want

%----------------------------------------------------------------------
%Output
%----------------------------------------------------------------
varargout{1} = double(L7); %Final segmented beads
varargout{2} = goodbeads_bool; % Boolean of segmented beads
varargout{3} = L; % all segemented beads

%-------------------------------------------End Edit %-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End Progthreshfluor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mystats,Lq,Lo] = collectstats(fimgs,L,qualfunc);

if nargin >=3,
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
n = max(max(L)) % num objects (Largest value in segmented image -> cell number in Image)

if n>0,
    stats = regionprops(L,'MinorAxisLength','MajorAxisLength','Eccentricity','Area','Centroid','ConvexArea','EulerNumber','Solidity'); %stats
    qualifies = feval(ourqualfunc,stats); %No clue what this is doing
    Lq = zeros(size(L));
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
        fimgs{i} = imresize(fimgs{i},sz,'bilinear');
        [subreg,shift,back]=regsnaps(L(100:424,100:424),fimgs{i},[100 100 424 424],3);
        fimgs{i}=imshift(fimgs{i},-shift);
        disp(['shift' num2str(i) ' = ' num2str(shift(1)) ' ' num2str(shift(2))]);
    end
end

k = 1;
for i = 1:n,
    if qualifies(i) 
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

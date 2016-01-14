function p = segmoviefluorH (p,varargin)


% Modified version of SEGMOVIEFLUORC to incorporate the HYBRID method for segmentation. 
% Example :
%
% p = segfluormovieH(p);
%
% It is modified to read RFP images into the variable X and phase images into the variable Y. 
% It is modified to call segfluorH(X,Y,p) 
% - Manoj - 






% SEGMOVIEPHASE   Segment movie from phase images.
%
%   SEGMOVIEPHASE(P) segments the movie from phase images according to the
%   parameters and controls described in the schnitzcells parameter structure P,
%   typically generated using INITSCHNITZ.
%
%   SEGMOVIEPHASE(P,'Field1',Value1,'Field2',Value2,...) segments the movie
%   from phase images according to the parameters and controls described in the
%   schnitzcells parameter structure P, but only after P has been updated by
%   setting P.Field1 = Value1, P.Field2 = Value2, etc.  Thus any schnitzcells
%   parameters can be updated (overridden) in the function call via these
%   optional field/value pairs.  (This is in the style of setting MATLAB
%   properties using optional property/value pairs.)  For a complete list of
%   schnitzcells parameter structure fields, see INITSCHNITZ.
%
%   SEGMOVIEPHASE returns a struct array, the updated schnitzcells parameter
%   structure p that reflects fields updated via any optional field/value
%   arguments, as well as fields updated to reflect the status of the
%   segmentation at the time the segmentation exited.
%
%   For example, P = SEGMOVIEPHASE(P,'segRange',[5:10]) would perform
%   segmentation on frames 005 thru 010 of the movie described in P, and
%   return the updated schnitzcells parameter structure P.
%
%-------------------------------------------------------------------------------
% Schnitzcells Fields / Parameters that you can adjust to control segmoviephase
%
%   segRange                range of frame numbers to segment; by default
%                           all image frames in imageDir will be segmented
%
%   prettyPhaseSlice        the phase image number with best focus; with more
%                           than one phase slice this defaults to 2, otherwise
%                           defaults to 1
%
%   segmentationPhaseSlice  the phase image number used when performing cell
%                           boundary segmentation; defaults to 1
%
%-------------------------------------------------------------------------------
%

%-------------------------------------------------------------------------------
% Parse the input arguments, input error checking
%-------------------------------------------------------------------------------

numRequiredArgs = 1;
if (nargin < 1) | ...
        (mod(nargin,2) == 0) | ...
        (~isSchnitzParamStruct(p))
    errorMessage = sprintf ('%s\n%s\n%s\n',...
        'Error using ==> segmoviephase:',...
        '    Invalid input arguments.',...
        '    Try "help segmoviephase".');
    error(errorMessage);
end

%-------------------------------------------------------------------------------
% Override any schnitzcells parameters/defaults given optional fields/values
%-------------------------------------------------------------------------------

% Loop over pairs of optional input arguments and save the given fields/values
% to the schnitzcells parameter structure
numExtraArgs = nargin - numRequiredArgs;
if numExtraArgs > 0
    for i=1:2:(numExtraArgs-1)
        if (~isstr(varargin{i}))
            errorMessage = sprintf ('%s\n%s%s%s\n%s\n',...
                'Error using ==> segmoviephase:',...
                '    Invalid property ', num2str(varargin{i}), ...
                ' is not (needs to be) a string.',...
                '    Try "help segmoviephase".');
            error(errorMessage);
        end
        fieldName = schnitzfield(varargin{i});
        p.(fieldName) = varargin{i+1};
    end
end

% Define segmentation parameters that vary with Movie Kind here:
switch lower(p.movieKind)
    
    case 'e.coli'
        disp('Movie kind: e.coli');
        if ~isfield(p,'edge_lapofgauss_sigma')
            p.edge_lapofgauss_sigma = 2;
        end
        if ~isfield(p,'minCellArea')
            p.minCellArea = 30;
        end
        if ~isfield(p,'minCellLengthConservative')
            p.minCellLengthConservative = 20;
        end
        if ~isfield(p,'minCellLength')
            p.minCellLength = 30;
        end
        if ~isfield(p,'maxCellWidth')
            p.maxCellWidth = 7;
        end
        if ~isfield(p,'minNumEdgePixels')
            %     p.minNumEdgePixels = 300;  % JCR: was dropped to 250 from original 300
            p.minNumEdgePixels = 250;  % JCR: value used for nature methods paper
        end
        if ~isfield(p,'maxThreshCut')
            p.maxThreshCut = 0.3;
        end
        if ~isfield(p,'maxThreshCut2')
            p.maxThreshCut2 = 0.2;
        end
        if ~isfield(p,'maxThresh')
            p.maxThresh = 0.25;
        end
        if ~isfield(p,'minThresh')
            p.minThresh = 0.25;
        end
        if ~isfield(p,'imNumber1')
            p.imNumber1 = 2;
        end
        if ~isfield(p,'imNumber2')
            p.imNumber2 = 1;
        end
        if ~isfield(p,'radius')
            p.radius = 5;
        end
        if ~isfield(p,'angThresh')
            p.angThresh = 2.7;
        end
        
    case 'bacillus'
        disp('Movie kind: bacillus');
        if ~isfield(p,'edge')
            p.edge_lapofgauss_sigma = 3; %3
        end
        if ~isfield(p,'minCellArea')
            p.minCellArea = 100; %60
        end
        if ~isfield(p,'minCellLengthConservative')
            p.minCellLengthConservative= 12;   %12  % JCR observed 2005-02-16
        end
        if ~isfield(p,'minCellLength')
            p.minCellLength = 20;  %20          % JCR observed 2005-02-16
        end
        if ~isfield(p,'maxCellWidth')
            p.maxCellWidth = 20;      %20       % JCR observed 2005-02-16
        end
        if ~isfield(p,'minNumEdgePixels')
            p.minNumEdgePixels = 215;        % ME set to 215 in mail to JCR 2005-08-18
        end
        if ~isfield(p,'maxThreshCut')
            p.maxThreshCut = 0.3; %0.3
        end
        if ~isfield(p,'maxThreshCut2')
            p.maxThreshCut2 = 0.2;  %0.2
        end
        if ~isfield(p,'maxThresh')
            p.maxThresh = 0.050; % 0.050  % changed from 0.05 10/12/08
        end
        if ~isfield(p,'minThresh')
            p.minThresh = 0.050 ; %0.050   % changed from 0.05 5/2/08
        end
        if ~isfield(p,'imNumber1')
            p.imNumber1 = 2; %2
        end
        if ~isfield(p,'imNumber2')
            p.imNumber2 = 1; %1
        end
        if ~isfield(p,'radius')
            p.radius = 5; %1
        end
        if ~isfield(p,'angThresh')
            p.angThresh = 2.7; %2.7
        end
        
    otherwise
        errorMessage = sprintf ('%s\n%s\n',...
            'Error using ==> segmoviephase:',...
            '    Invalid movie kind setting internal movie kind parameters.');
        error(errorMessage);
end

disp ('using schnitzcells parameter structure:');
disp (p);

% figure out (automatically) if we have just 1 phase or multiple phase images
% This affects/selects which phase images to use... at least initially
% for determining number and range of frames in the movie.
mnamePhase2   = [p.movieName,'-y-2-*.tif'];
mnamePhaseAll = [p.movieName,'-y-*.tif'];
%mnamePhaseAll = [p.movieName,'-t-*.tif']; % if just RFP channel

Dphase2   = dir([p.imageDir, mnamePhase2]);
DphaseAll = dir([p.imageDir, mnamePhaseAll]);

if ~isfield(p,'numphaseslices')
    % Derive the number of phase slices by examining the image directory files.
    % We understand 3 cases:
    %  1) if one   phase slice,  Dphase2 is empty
    %  2) if two   phase slices, DphseAll is 2x size of Dphase2
    %  3) if three phase slices, DphseAll is 3x size of Dphase2

    if isempty(DphaseAll)
        error(['Can''t find any images in directory ' p.imageDir]);
    end

    if isempty(Dphase2)
        p.numphaseslices = 1;
        disp('You appear to have 1 phase image per frame: setting p.numphaseslices = 1;');
    elseif length(DphaseAll) == length(Dphase2)*2
        p.numphaseslices = 2;
        disp('You appear to have 2 phase images per frame: setting p.numphaseslices = 2;');
    elseif length(DphaseAll) == length(Dphase2)*3
        p.numphaseslices = 3;
        disp('You appear to have 3 phase images per frame: setting p.numphaseslices = 3;');
    else
        error(['Can''t figure out how many phase slices are in your movie!']);
    end
end

% This code is not robust: problem if there are 2 or 3 phases and you set
% numphaseslices = 1, D is too big.  Smarter to figure out which D is
% smaller, and use the smaller one  (BTW, same code's in expand_segimage)
if p.numphaseslices == 1
    D = DphaseAll;
    if length(D) == 0
        errorMessage = sprintf('%s%s\n%s%s\n',...
            '    No images found in directory ', p.imageDir, ...
            '    matching pattern ', mnamePhaseAll);
        error(errorMessage);
    end
else
    D = Dphase2;
    if length(D) == 0
        errorMessage = sprintf('%s%s\n%s%s\n',...
            '    No images found in directory ', p.imageDir, ...
            '    matching pattern ', mnamePhase2);
        error(errorMessage);
    end
end

% Look at phase image file names to figure out how many frames we have
[s,I] = sort({D.name}');
D = D(I);
numpos= findstr(D(1).name, '.tif')-3;

% if user didn't specify a range of frames to segment, figure it out
% given the image names
if ~isfield(p,'segRange')
    imageNameStrings = char(s);
    p.segRange = str2num(imageNameStrings(:,numpos:numpos+2))';
end

%For test run, do not save segmentation
if isfield(p,'testrun')
    SAVESEG = 0;
else
    SAVESEG = 1;
end

% set some defaults if they don't exist
if exist('outprefix')~=1
    outprefix = [p.movieName 'seg'];
end
if ~isfield(p,'prettyPhaseSlice');
    if p.numphaseslices>1
        p.prettyPhaseSlice = 2;
    else
        p.prettyPhaseSlice = 1;
    end
end
if ~isfield(p,'segmentationPhaseSlice');
    p.segmentationPhaseSlice = 1;
end
regsize= 0;% max translation (pixels) between phase and fl. images - obsolete.

for i= p.segRange
    
    mynum = str3(i);
    % Manoj edit
    Dframe_RFP = dir([p.imageDir p.movieName '*-t*-' str3(i) '.tif']);
    Dframe_phase = dir([p.imageDir p.movieName '*-p*-' str3(i) '.tif']);
    
    
    disp(['reading ',p.imageDir,Dframe_RFP(1).name]);
    pname_RFP = Dframe_RFP(1).name;
    pname_phase = Dframe_phase(1).name;

    if p.numphaseslices==1
        X(:,:,1) = imread([p.imageDir,pname_RFP]);
        Y(:,:,1) = imread([p.imageDir,pname_phase]);
    else
        fstr = findstr(pname_RFP,'-p-');
        for islice = 1:p.numphaseslices
            pname_RFP(fstr+3) = num2str(islice);
            X(:,:,islice) = imread([p.imageDir,pname_RFP]);
            Y(:,:,islice) = imread([p.imageDir,pname_phase]);

        end
    end
    
    [phsub,LNsub,rect]= segfluorH(X,Y,p);
    phaseFullSize = size(LNsub);
    if size(phsub,3)>1, % nitzan June25: this used to be: phsub = phsub(:,:,2);
        phsub = phsub(:,:,p.prettyPhaseSlice);
    end;
    savelist=['''phsub'',''LNsub'',''rect'',''timestamp'',''phaseFullSize'''];
    % nitzan 2005June25 added the following to save both phases if different:
    if p.segmentationPhaseSlice~=p.prettyPhaseSlice
        phseg = phsub(:,:,p.segmentationPhaseSlice);
        savelist = [savelist,',''phseg'''];
    end
    % end nitzan 2005June25
    %     [TSTa,TSTb,TSTc,TSTd,timestamp] = imsettings([p.imageDir,Dframe(1).name],'p');
    %     clear TST*
    
    Lname= [outprefix,mynum];
    cname= [p.imageDir,p.movieName,'-c-',mynum,'.tif'];
    yname= [p.imageDir,p.movieName,'-y-',mynum,'.tif'];
    gname= [p.imageDir,p.movieName,'-g-',mynum,'.tif'];
    rname= [p.imageDir,p.movieName,'-t-',mynum,'.tif'];
    if exist(cname)==2
        disp('found CFP image');
        [creg, cshift, cback, cbinning]= quicknoreg(LNsub,cname,rect,regsize,phaseFullSize);
        %       [exptcstr, gainc, exptc]= imsettings(cname);
        %       savelist=[savelist,',''creg'',''cshift'',''exptc'',''gainc'',''cback'',''cbinning'''];
        savelist=[savelist,',''creg'',''cshift'',''cback'',''cbinning'''];
    end
    if exist(yname)==2
        disp('found YFP image');
        [yreg, yshift, yback, ybinning]= quicknoreg(LNsub,yname,rect,regsize,phaseFullSize);
        % [exptystr, gainy, expty]= imsettings(yname);
        % savelist=[savelist,',''yreg'',''yshift'',''expty'',''gainy'',''yback'',''ybinning'''];
        savelist=[savelist,',''yreg'',''yshift'',''yback'',''ybinning'''];
    end
    
    if exist(gname)==2
        disp('found GFP image');
        [greg, gshift, gback, gbinning]= quicknoreg(LNsub,gname,rect,regsize,phaseFullSize);
        %       [exptgstr, gaing, exptg]= imsettings(gname);
        %       savelist=[savelist,',''greg'',''gshift'',''exptg'',''gaing'',''gback'',''gbinning'''];
        savelist=[savelist,',''greg'',''gshift'',''gback'',''gbinning'''];
    end
    
    if exist(rname)==2
        disp('found RFP image');
        [rreg, rshift, rback, rbinning]= quicknoreg(LNsub,rname,rect,regsize,phaseFullSize);
        %   [exptrstr, gainr, exptr]= imsettings(rname);
        % savelist=[savelist,',''rreg'',''rshift'',''exptr'',''gainr'',''rback'',''rbinning'''];
        savelist=[savelist,',''rreg'',''rshift'',''rback'',''rbinning'''];
    end
    timestamp = [];
    if SAVESEG
        eval(['save(''',p.segmentationDir,Lname,''',',savelist,');']);
        disp(['saved file ',p.segmentationDir,Lname]);
    end
    
    clear L*;
end
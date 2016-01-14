function L7 = breakingcells(bwimage,phase, maxthresh, minthresh, minCellLength)

% Breaks cells 
% Output = breakingcells(I,PHASE,maxthresh, minthresh, minCellLength)

% maxthresh & minthresh:    threshold values (smaller they are the more cuts are included)
% mincelllength:    cuts which create cells smaller than this will be ignored


% BREAKING CELLS
% break up big cells
% Goes along the THINned cell and looks for places where there is a change
% in the phase value - i.e. where there could be a space between cells.
% uses PHSUB(:,:,p.imNumber1) (default p.imNumber1=2)

L6 = bwlabel(bwimage);
L7= bwlabel(bwimage);

% p.maxThresh = 0.050;
% p.minThresh = 0.050 ; %0.050   % changed from 0.05 5/2/08
% p.minCellLength = 30;

r= regionprops(L6, 'majoraxislength');
fbiggies= find(([r.MajorAxisLength]>50));
disp(['Breaking up big cells(', num2str(length(fbiggies)),').']);
imshowlabel(L7);
for z = 1:length(fbiggies),
    %disp([' ',num2str(i),': checking cell number ',num2str(fbiggies(i))]);
    Lcell= +(L6 == fbiggies(z)); % + converts logical to double
    Lcell(Lcell == 1)= fbiggies(z);
    %         cutcell= breakcellfluor(Lcell, imcomplement(combined), ...
    %             0.2, 0.005, p.minCellLength);
%     cutcell= breakcellfluor(Lcell, imcomplement(phase), ...
%         0.02, 0.004, p.minCellLength);
    
    cutcell= breakcellfluor(Lcell, imcomplement(phase), ...
        maxthresh, minthresh, minCellLength);

    L7(L6 == fbiggies(z))= 0;
    % place cutcell
    label= max2(L7);
    for j = 1:max2(cutcell),
        L7(cutcell==j)= label+j;
    end;
end;

L7 = renumberimage(L7);

end
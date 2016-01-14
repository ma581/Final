function L6  = cuttingcells(L4, maxthresh, maxcellwidth, mincelllength)

%% CUTTING INDIVIDUAL CELLS at the narrow waist: makes sure septation events
% are properly identified and that cells are identified seperately. Cuts
% cells at points where both sides of the cell are sufficiently concave.
%   maxthresh:      smaller maxtresh the more points are cut
%   maxcellwidth:   two points must be closer than this for cut to be accepted
%   mincelllength:  cut ignored if it creates `cell' smaller than this
%COmmented out 8/2/08
L6 = L4;
r= regionprops(L6, 'solidity');
fkinks1= find(([r.Solidity] > 0.75 ));
fkinks2= find(([r.Solidity] < 0.90 ));
fkinks=intersect(fkinks1,fkinks2);
disp(['Cutting individual cells(',num2str(length(fkinks)),').']);
for z= 1:length(fkinks)
    %disp([' ',num2str(i),':cutting cell number ',num2str(fkinks(i))]);
    Lcell= (L6 == fkinks(z));
%     cutcell= cutcurv(Lcell, 0.4, 30, 10);
       cutcell= cutcurv(Lcell, maxthresh, maxcellwidth, mincelllength);
    cellnos= unique(cutcell);
    L6(Lcell)= 0;
    label= max2(L6);
    for j= 2:length(cellnos)
        L6(find(cutcell == cellnos(j)))= label+j;
    end;
end;

L6= renumberimage(L6);

end
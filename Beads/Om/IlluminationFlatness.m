function IlluminationFlatness(filename,rad);

% sfile = dir([filename '\Sbeads*']);
sfile = dir([filename '\SJLB021_3_001-1*']);
sname = sfile.name;
smatLoad = load([filename '\' sname]);
p = findstr('-1.mat',sname);
snameshort = sname(1:p-1);

eval(['beadsdata = smatLoad.' snameshort ';']);

beadLoc = cell2mat(arrayfun(@(x)(x.Centroid),beadsdata,'uniformoutput',false));
beadLoc = reshape(beadLoc,[2,length(beadLoc)/2])';

beadVal = cell2mat(arrayfun(@(x)x.mean,beadsdata,'uniformoutput',false));

%%%%%%%%%%%%%%%%%%Manoj Edit %%%%%%%%%%%%
beadVal = beadVal(1:2:length(beadVal)); %Looks only at RFP means

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
beadsCalib = zeros(1040,1392);
% beadsCalib = zeros(104,139);
tint = 1;
for i = 1:tint:1040
    for j = 1:tint:1392
        beadLocDist = ((beadLoc(:,2)-i).^2 + (beadLoc(:,1)-j).^2).^(1/2);
        beadLocClose = beadLocDist <= rad;
        beadValClose = beadVal(beadLocClose);
        beadsCalib(i,j) = median(beadValClose);

    end
end

figure;
imshow(beadsCalib,[]);
eval(['save ', filename ,'\beadsCalib' num2str(rad) ' beadsCalib']); 
toc

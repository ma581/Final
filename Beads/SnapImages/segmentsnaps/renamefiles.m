
name = 'mod_hi_od_50_bioW_y_sucC_c';

for j = 1 : 15
    im_p = imread(['p' num2str(j) '.tif']);
    im_c = imread(['cfp' num2str(j) '.tif']);
    im_r = imread(['rfp' num2str(j) '.tif']);
    im_y = imread(['yfp' num2str(j) '.tif']);
    
    imwrite(im_p, [name '-' str2(j) '-p.tif'],'tif');
    imwrite(im_c, [name '-' str2(j) '-c.tif'],'tif');
    imwrite(im_r, [name '-' str2(j) '-r.tif'],'tif');
    imwrite(im_y, [name '-' str2(j) '-y.tif'],'tif');
end
function crossApEn = calcCrossApEn_noFilt(brainMask,img_4D,m,r,ROI_voxel)
ROI_ts = squeeze(img_4D(ROI_voxel(1),ROI_voxel(2),ROI_voxel(3),:));
std_ROI_ts = std(double(ROI_ts));
r_val = r*std_ROI_ts;
brainVox = find(brainMask==max(brainMask(:)));
img_size = size(img_4D);
% cross Approximate-entropy computation
crossApEn = zeros(img_size(1),img_size(2),img_size(3));
nFail = 0;
msg = 'calculating CrossApEn.';
h = waitbar(0,msg);
for vox = 1:length(brainVox)
    [row,col,sl] = ind2sub([img_size(1) img_size(2) img_size(3)],brainVox(vox));
    tmp_TS1 = squeeze(img_4D(row,col,sl,:));
    std_tmp_TS1 = std(double(tmp_TS1));
    multFactor = std_ROI_ts/std_tmp_TS1;
    TS1 = tmp_TS1.*multFactor;
    TS2 = ROI_ts;
    tmp = cross_approx_entropy(m,r_val,TS1,TS2); %max(xcorr(double(TS1),double(TS2')));
    crossApEn(row,col,sl) = tmp(1);
    nFail = nFail+tmp(2); %0; 
    waitbar(vox/length(brainVox));
end
close(h);
return

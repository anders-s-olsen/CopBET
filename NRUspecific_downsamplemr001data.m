function image_4D = NRUspecific_downsamplemr001data(image_4D)
disp('Resampling to TR=2s')
datatmp = double(reshape(image_4D,[91*109*91,size(image_4D,4)])');
datatmpresamp = resample(datatmp,2,5)';
image_4D = reshape(datatmpresamp,[91,109,91,size(datatmpresamp,2)]);
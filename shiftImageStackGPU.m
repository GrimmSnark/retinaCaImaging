function imgStack = shiftImageStackGPU(imgStack,XYOffsets)

for i = 1:size(imgStack,3)
% % shift image discretely by whole pixel values
%     row_shift = floor(XYOffsets(i,1));
%     col_shift = floor(XYOffsets(i,2));
%     registeredImage = circshift(imgStack(:,:,i),[row_shift col_shift]);
%     
%     % make rows and cols zero where we shifted
%     if(row_shift>0), registeredImage(1:row_shift,:,:) = 0; elseif(row_shift<0), registeredImage((end+row_shift+1):end,:,:) = 0; end
%     if(col_shift>0), registeredImage(:,1:col_shift,:) = 0; elseif(col_shift<0), registeredImage(:,(end+col_shift+1):end,:) = 0; end   
% 
%     imgStack(:,:,i) = registeredImage;

tForm = transltform2d(XYOffsets(i,:));
  imgStack(:,:,d) = imwarp(imgStack(:,:,i),tForm);

end
end
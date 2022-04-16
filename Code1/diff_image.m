function [u] = diff_image(image)

[image_h,image_w]=size(image);
% y1=image(2:image_h,:)-image(1:image_h-1,:);
y1 = diff(image); % the matrix of row differences
% y2=image(:,2:image_w)-image(:,1:image_w-1);
% image_t = image';
% y2 = diff(image_t);
% y2 = y2';
y2 = diff(image,1,2);

y1=[y1;zeros(1,image_w)];
y2=[y2,zeros(image_h,1)];
u=[y1,y2];
end
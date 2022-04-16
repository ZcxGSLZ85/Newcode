function [ u ] = div_image(image)

[image_h,image_w]=size(image);
y1=image(:,1:image_w/2);
y2=image(:,image_w/2+1:image_w);

y1 = y1(1:image_h-1,:);
% p1=y1(2:image_h-1,:)-y1(1:image_h-2,:);
p1 = diff(y1);
p1 = [y1(1,:);p1;-y1(image_h-1,:)];

y2 = y2(:,1:image_w/2-1);
% p2=y2(:,2:image_w/2-1)-y2(:,1:image_w/2-2);
p2 = diff(y2,1,2);
p2 = [y2(:,1),p2,-y2(:,image_w/2-1)];

u=p1+p2;
end


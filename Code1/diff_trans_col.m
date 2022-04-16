function [ yt_vec ] = diff_trans_col( y,m,n )


% the transpose of difference matrix D*x, that is
% D^{T}*x

y_map = reshape(y,m,n);
diff_y = -diff(y_map(:,1:n-1),1,2);
yt_map = [-y_map(:,1),diff_y,y_map(:,end-1)];
yt_vec = yt_map(:);


end


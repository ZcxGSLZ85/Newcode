function [ x_vec ] = diff_two_dimensional( x,m,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% difference matrix D*x

X_map = reshape(x,m,n);
diff_row = [diff(X_map);zeros(1,n)];
diff_col = [diff(X_map,1,2),zeros(m,1)];
x_vec = [diff_col(:);diff_row(:)];




end


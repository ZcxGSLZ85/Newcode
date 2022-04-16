
% original image
 f = double(rgb2gray(imread('goldhill512','png')));  

% blur operator
% mean filter
 psf = fspecial('average', 9);        % Define a mean filter of size 9*9;
% Gaussian filter
% psf = fspecial('gaussian',7,10);    % Define a Gaussian filter of size 7*7;

% noisy image
% random noise
% rng('default');  
% fn = imfilter(f, psf) + 1.5 * randn( size(f) );    % PSF filtering and adding 1.5 random normal distribution noise£»
% fn = imfilter(f, psf) + 3 * randn( size(f) );      % PSF filtering and adding 3 random normal distribution noise£»
 load goldhill_noise3_date.mat       % Mean filtering of size 9*9 and Gaussian noise of standard deviation 1.5£»
% load goldhill_noise4_date.mat      % Mean filtering of size 9*9 and Gaussian noise of standard deviation 3£»
% load goldhill_noise1_date.mat      % Gaussian filtering of size 7*7 with standard deviation 10 and Gaussian noise with standard deviation of 1.5£»
% load goldhill_noise2_date.mat      % Gaussian filtering of size 7*7 with standard deviation 10 and Gaussian noise with standard deviation of 3£»


[m,n] = size(f);   

D = make_derivatives(m,n);   % difference
epsilon = norm(diff_image(f),1);

dy = D*f(:);
ani_tv = sum(abs(dy));  % Anisotropic TV

L = sum(abs(psf(:)));


  tol = 5e-4;
 iter = 1000;

% initial value 
x0 = zeros(size(f)); % The initial iterative guess 
v0 = zeros(size(f,1),2*size(f,2));



%% Compared experiments
% Choosing parameters
mu = 0.1;
alpha=0.8;
gamma = 0.1/sqrt(8);
pro = 1.5;

[u1,snr1,i1,time1] =Forward_backward_half_forward_deblurring(f,x0,v0,fn,psf,gamma,mu,tol,iter);
[u2,snr2,i2,time2] =Relaxed_Inertial_forward_backward_half_forward_deblurring(f,x0,v0,fn,psf,alpha,gamma,mu,pro,tol,iter);



set(0,'defaultfigurecolor','w') 
figure(1)
subplot(2,2,1); imshow(f/255,[]); disp('(a)');  %'Original image'
xlabel('(a)')
subplot(2,2,2); imshow(fn/255,[]); disp('(b)');  %'Noisy image'
xlabel('(b)')
subplot(2,2,3); imshow(u1/255,[]); disp('(c)');  %'FBHF'
xlabel('(c)')
subplot(2,2,4); imshow(u2/255,[]); disp('(d)');  %'RIFBHF'
xlabel('(d)')



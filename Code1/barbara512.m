% original image
 f = double(imread('barbara512.png')); 
 
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
 load barbara_noise3_date.mat     % Mean filtering of size 9*9 and Gaussian noise of standard deviation 1.5£»
% load barbara_noise4_date.mat    % Mean filtering of size 9*9 and Gaussian noise of standard deviation 3£»
% load barbara_noise1_date.mat    % Gaussian filtering of size 7*7 with standard deviation 10 and Gaussian noise with standard deviation of 1.5£»
% load barbara_noise2_date.mat    % Gaussian filtering of size 7*7 with standard deviation 10 and Gaussian noise with standard deviation of 3£»

[m,n] = size(f);   

D = make_derivatives(m,n);     % difference 
epsilon = norm(diff_image(f),1); 

dy = D*f(:); 
ani_tv = sum(abs(dy));   % Anisotropic TV

L = sum(abs(psf(:)));   

 tol=5e-4;
 iter = 500;

% initial value 
x0 = zeros(size(f)); % The initial iterative guess 
v0 = zeros(size(f,1),2*size(f,2));


%% The results of SNR and ||f-f*|| of multiple parameters gamma or lambda

% Choose parameters
mu = 0.1;
alpha=0.1;
gamma = 0.2/sqrt(8);
lambda = [0.1 0.4 0.7 1 1.4 1.7];

% Cellular array
SNR2 = cell(1,length(lambda));
u = cell(1,length(lambda));
RES2 = cell(1,length(lambda));
t2 = cell(1,length(lambda));
k2 = zeros(1,length(lambda));


tic
for i = 1:length(lambda)
[u{i},SNR2{i},k2(i),RES2{i},t2{i}]= Relaxed_Inertial_forward_backward_half_forward_deblurring(f,x0,v0,fn,psf,alpha,gamma,mu,lambda(i),tol,iter);
end
time = toc;

% \alpha=0.1,\gamma = 0.2/sqrt(8); 
 figure(1);
semilogy(t2{1},RES2{1},'b',t2{2},RES2{2},'g',t2{3},RES2{3},'r',t2{4},RES2{4},'c',t2{5},RES2{5},'m',t2{6},RES2{6},'k');
 title('$\alpha=0.1, \gamma=0.2/\sqrt{8}$','Interpreter','latex');
xlabel('Time t'),ylabel('||f-f*||');
 legend('\lambda=0.1','\lambda=0.4','\lambda=0.7','\lambda=1','\lambda=1.4','\lambda=1.7');

 figure(2);
plot(t2{1},SNR2{1},'b',t2{2},SNR2{2},'g',t2{3},SNR2{3},'r',t2{4},SNR2{4},'c',t2{5},SNR2{5},'m',t2{6},SNR2{6},'k');
 title('$\alpha=0.1, \gamma=0.2/\sqrt{8}$','Interpreter','latex');
xlabel('Time t'),ylabel('SNR');
 legend('\lambda=0.1','\lambda=0.4','\lambda=0.7','\lambda=1','\lambda=1.4','\lambda=1.7');


%% 
%% Parametric three-dimensional diagram

% mu = 0.1;
% gamma = 0.1/sqrt(8);
% k1=[];
% figure
% jj=1;
% % % Selection parameters
% alpha_aa = [0 0.05 0.1 0.2 0.3 0.5 0.7 0.9 1];
% lambda_aa  = [0.1 0.2 0.3 0.5 0.7 1 1.3 1.5 1.7];
% 
% 
%     for j=1:9
%         alpha=alpha_aa(j);
%         for k=1:9
%         lambda=lambda_aa(k);   
%        [X1,SNR1,k1,RES,t1] = Relaxed_Inertial_forward_backward_half_forward_deblurring(f,x0,v0,fn,psf,alpha,gamma,mu,lambda,tol,iter);
%          bb(j,k)=k1;  
%         end       
%     end        
%     %%
%     subplot(1,2,2)
%     pcolor(alpha_aa,lambda_aa,bb);
%     suptitle('a ');
%     shading interp
%     xlabel(' $\alpha$ ','Interpreter','latex','Fontsize',12);
%     ylabel(' $\lambda$','Interpreter','latex','Fontsize',12);
%     colorbar
%     %%
%     subplot(1,2,1)
%     surf(alpha_aa,lambda_aa,bb);
%     shading interp
%     xlabel(' $\alpha$ ','Interpreter','latex','Fontsize',12);
%     ylabel(' $\lambda$','Interpreter','latex','Fontsize',12);
%     zlabel('Iter','Interpreter','latex','Fontsize',12);
%      hold on

 
%% Compared experiments
%Choosing parameters
% mu = 0.1;
% alpha=0.8;
% gamma = 0.1/sqrt(8);
% pro = 1.5;
% 
% [u1,snr1,i1,time1] =Forward_backward_half_forward_deblurring(f,x0,v0,fn,psf,gamma,mu,tol,iter);
% [u2,snr2,i2,time2] =Relaxed_Inertial_forward_backward_half_forward_deblurring(f,x0,v0,fn,psf,alpha,gamma,mu,pro,tol,iter);
%  
% 
% set(0,'defaultfigurecolor','w') 
% figure(1)
% subplot(2,2,1); imshow(f/255,[]); disp('(a)');  %'Original image'
% xlabel('(a)')
% subplot(2,2,2); imshow(fn/255,[]); disp('(b)');  %'Noisy image'
% xlabel('(b)')
% subplot(2,2,3); imshow(u1/255,[]); disp('(c)');  %'FBHF'
% xlabel('(c)')
% subplot(2,2,4); imshow(u2/255,[]); disp('(d)');  %'RIFBHF'
% xlabel('(d)')
% 
% 

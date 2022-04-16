 
% clc,clear
% 
% path(path,'./rof/');
% path(path,'./image bechmark/');
% path(path,'./toolbox_optim/toolbox/')
% path(path,'./toolbox_optim/tests/')
% path(path,'./toolbox_optim/')
% path(path,'./noise_data/')
% path(path,'./data_wu/')

% f = double(imread('lena256','png'));
% f = double(imread('house','png'));
% f = double(imread('cameraman512','tif'));
% f = double(imread('vendredi','tif'));
% f = double(imread('fingerprint512','png'));
% f = double(imread('peppers256','png'));
% f = double(imread('house256','png'));
% f = double(imread('house','png'));
 f = double(imread('barbara512','png'));
% f = double(imread('baboon512[1]','jpg'));
%  f = double(imread('untitled','png'));
%  f = double(imread('Pirate','tif'));
% f = double(imread('Lena512','png'));

% f = double(rgb2gray(imread('building_org','png')));
% f = f/255;

% f = double(imread('boat512','png'));
% f = double(imread('goldhill','png'));


% f = double(imread('text','png'));

[m,n] = size(f);
D = make_derivatives(m,n);


% noise level
 sigma = 0.1;
% 
% g = f + sqrt(sigma)*randn(m,n);
% load data_barbara512_30.mat;
% load boat512_30.mat;
% load data_lena512.mat;
% load goldhill512_30.mat;

% load text_001.mat;
% load text_005.mat;
% load text_01.mat;

% load building_15.mat;
load barbara_noise1_date.mat

%% image denoise
 epsilon = 1e-2;
 gamma = 0.11;
 a1 = 1;
 iter = 100;

% lambda =  [5.9:0.1:6.1]; % 0.22 for building 0.1 noise level
 
lambda = [2.7:0.1:3];

% lambda2 = [0.01:0.01:0.3]; % 0.9 for building 0.1 noise level
% 
y0 = zeros(m,2*n);
% 
% 

%%
% [x_update,k,SNR,SSIM,PSNR,t1,t2]= tv_denoise_dual_constrained(f,g,y0,lambda,gamma,iter,epsilon,a1);


%%
SNR = cell(1,length(lambda));
SSIM = zeros(1,length(lambda));
x_update = cell(1,length(lambda));
PSNR = cell(1,length(lambda));
t = cell(1,length(lambda));
k = zeros(1,length(lambda));

tic
for i = 1:length(lambda)
[x_update{i},k(i),SNR{i},SSIM(i),PSNR{i},t{i}]= tv_denoise_dual_constrained(f,fn,y0,lambda(i),gamma,iter,epsilon,a1);
end
time = toc;
% 
% 
% [ms,is] = max(SNR);
% 
% lambda(is)

%%
% % 
%  plot(k2(1),SNR{2},'r-*');
% %     title('Iteration steps of GMAOR with diferent parameters k (m=20)');
%     axis ([0.1 40 0 17]);
%     xlabel('The parameter alpha_msor ','Interpreter','latex','Fontsize',12);
%     ylabel('The iteration steps','Interpreter','latex','Fontsize',12);
%     legend('arctan(x)')
%      hold on;
% 

%  figure; colormap gray;
% subplot(221); imagesc(f); axis image; axis off; title('Original');
% subplot(222); imagesc(fn); axis image; axis off; title('Noisy');
% subplot(223); imagesc(x_update{is});axis off; axis image; 
% subplot(224); imagesc(x_update{length(lambda)});axis off; axis image; 
% title('AIsotropic TV denoising');

figure(1);
plot(1:k(1),SNR{1},'k:',1:k(2),SNR{2},'y-.',1:k(3),SNR{3},'y-.',1:k(4),SNR{4},'g-.')

%subplot(2,2,1);
% plot(1:k(1),SNR{1}),title('SNR vs iteration numbers'),xlabel('iteration numbers'),ylabel('SNR');
% %subplot(2,2,2);
% hold on
% plot(1:k(2),SNR{2}),title('SNR vs iteration numbers'),xlabel('iteration numbers'),ylabel('SNR');
% %subplot(2,2,3);
% hold on
% plot(1:k(3),SNR{3}),title('SNR vs iteration numbers'),xlabel('iteration numbers'),ylabel('SNR');
% %subplot(2,2,4);
% hold on
% plot(1:k(4),SNR{4}),title('SNR vs iteration numbers'),xlabel('iteration numbers'),ylabel('SNR');
% %axis ([0 40 0 20]);
% legend('\lambda=1','\lambda=2','\lambda=3','\lambda=4')
% hold on
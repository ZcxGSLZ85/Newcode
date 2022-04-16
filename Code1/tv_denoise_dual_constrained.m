function [x_update,k,SNR,SSIM,PSNR,t2]= tv_denoise_dual_constrained(f_true,f,y0,lambda,gamma,iter,epsilon,a1)

 [m,n]=size(f);
% y =zeros(m,2*n);

y = y0;
k =1;

SNR = [];
primal_residual = [];


tic;
t1 = clock;

while k <= iter
    
% tic;
% t2 = clock;
    

%     x_update = proj_bound(f + gamma*div_image(y),0,255);
     x_update = (f + gamma*div_image(y));
    
    xd = diff_image(x_update);
    y1 = y + xd;
    
     yp = max(abs(y1)-lambda/gamma,0).*sign(y1);
    
   y_update = (1-a1)*y + a1 *( y1 - yp); % anisotropic total variation
      
      % isotropic total variation
%       y1_sum = sqrt(y1(1:m*n).^2 + y1(m*n+1:end).^2);
%       y11= max(y1_sum - lambda/gamma,0).*(y1(1:m*n)./y1_sum);
%       y12 = max(y1_sum - lambda/gamma,0).*(y1(m*n+1:end)./y1_sum);
%       y_update = (1-a1)*y + a1 * (y1-[y11;y12]);


      % isotropic total variation
%       y1_sum = sqrt(y1(:,1:n).^2 + y1(:,n+1:end).^2);
%       y11= max(y1_sum - lambda/gamma,0).*(y1(:,1:n)./y1_sum);
%       y12 = max(y1_sum - lambda/gamma,0).*(y1(:,n+1:end)./y1_sum);
%       y_update = (1-a1)*y + a1 * (y1-[y11,y12]);
      
      
      
   %  y_update = (1-a1)*y1 + a1 *( y1 - max(abs(y1)-lambda/gamma,0).*sign(y1));   %
   
        SNR(k)   = 20*log10(norm(f_true(:))/norm(f_true(:)-x_update(:)));
        SSIM = ssim(f_true,x_update);
         pr = xd - yp;
        primal_residual(k) = norm(pr(:)); 
      PSNR(k)   = 20*log10(sqrt(m*n)*255/norm(f_true(:)-x_update(:))) ;
   %     fval(k) = 0.5*norm(x_update-f)^2 + lambda*norm(D*x_update,1);
    
   
  %  time1(k) = etime(clock,t2);
   t2(k) = etime(clock,t1);
   
    if  norm(y_update(:)-y(:))/norm(y(:)) <= epsilon
        break;
    else
       y = y_update;   
       k = k+1;
    end
    
 
% u=f+lamda*div(y1_update,y2_update,1);

% mse(k) = norm(f_true - u,'fro')^2/(m*n);

end
end


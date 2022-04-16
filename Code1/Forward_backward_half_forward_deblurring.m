


function [z1_update,snr,i,time] = Forward_backward_half_forward_deblurring(x_true,x0,v0,b,psf,gamma,mu,tol,iter)
%  [m,n] = size(x_true);
%  xaverage = mean(x_true);
%  f = @(x) 0.5*norm(A*x-b)^2 + mu*norm(diff_image(x),1);

% relerror = zeros(iter,1);
% abserror = zeros(iter,1);
% eerror = zeros(iter,1);
%  derror = zeros(iter,1);
% snr = zeros(iter,1);
% psnr = zeros(iter,1);
% mssim = zeros(iter,1);



% data fidelity
A_dir  = @(x) imfilter(x, psf);
A_adj  = @(x) imfilter(x, rot90(psf,2));  % WARNING: 'psf' must be a (2n+1)-by-(2n+1) matrix
g = @(x) A_adj(A_dir(x) - b); %A'(Ax-b)

x = x0;
v = v0;



 i = 1;

 tic;
 while i <= iter
    
     
     nabla_f = g(x);
     u1 = x - gamma*(-div_image(v) + nabla_f);
     x_update = max(u1,0);
     u2 = v + gamma*diff_image(x);
     v_update = u2 - gamma*max(abs(u2/gamma)-mu/gamma,0).*sign(u2/gamma);
    
    
     z1_update = x_update + gamma*(-div_image(v) + div_image(v_update));
     z2_update = v_update + gamma*(-diff_image(x) + diff_image(x_update));
    
 

%z^{1,k}=J_{gamma_{k}f}(x-gamma_{k}(L^{*}y+Cx));
%z^{2,k}=J_{gamma_{k}g^{*}}(y+gamma_{k}Lx);
%x^{k+1}=z^{1,k}+gamma_{k}(L^{*}y-L^{*}z^{2,k});
%y^{k+1}=z^{2,k}+gamma_{k}(-Lx+Lz^{1,k}).


          %derror(i) = norm(z1_update(:)-x_true(:))/norm(x_true(:)); 
        
          snr(i) = 20*log10(norm(x_true(:))/norm(x_true(:)-z1_update(:)));      
        
    if  norm(z1_update(:)-x(:))/norm(x(:)) <= tol
         break;
      else 
        x = z1_update;
        v = z2_update;
        i = i+1;
    end
 end


  time = toc;   


        
        
        
end

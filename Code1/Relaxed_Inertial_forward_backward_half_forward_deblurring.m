
function [z1_update,SNR,i,RES,t2] = Relaxed_Inertial_forward_backward_half_forward_deblurring(x_true,x0,y0,b,psf,alpha,gamma,mu,pro,tol,iter)
%function [z1_update,SNR,i,t2] = Relaxed_Inertial_forward_backward_half_forward_deblurring(x_true,x0,y0,b,psf,alpha,gamma,mu,pro,tol,iter)
 
%[m,n] = size(x_true);
 %xaverage = mean(x_true);
% f1 = 0.5*norm(A*x_true-b)^2 + mu*norm(diff_image(x_true),1);

% data fidelity
A_dir  = @(x) imfilter(x, psf);
A_adj  = @(x) imfilter(x, rot90(psf,2));  % WARNING: 'psf' must be a (2n+1)-by-(2n+1) matrix
g = @(x) A_adj(A_dir(x) - b); %A'(Ax-b)
 
yp = 0.5*norm(A_dir(x_true)-b)^2 + mu*norm(diff_image(x_true),1);
 
%Initial variable
x0 = x0;
x1 = x0;
y0 = y0;
y1 = y0;


SNR = [];
RES = [];

i = 1;

 tic;
t1 = clock;

 while i <= iter
      
     w1 = x1 + alpha*(x1-x0);
     w2 = y1 + alpha*(y1-y0);
     
     nabla_f = g(w1);
     u1 = w1 - gamma*(-div_image(w2) + nabla_f);
     x_update = max(u1,0);
     u2 = w2 + gamma*diff_image(w1);
     y_update = u2 - gamma*max(abs(u2/gamma)-mu/gamma,0).*sign(u2/gamma);
    
    
     z1_update = (1-pro)*w1+pro*(x_update + gamma*(-div_image(w2) + div_image(y_update)));
     z2_update = (1-pro)*w2+pro*(y_update + gamma*(-diff_image(w1) + diff_image(x_update)));
    

     
%w^{1,k}=x^{k}+alpha_{k}(x^{k}-x^{k-1});
%w^{2,k}=y^{k}+alpha_{k}(y^{k}-y^{k-1});
%z^{1,k}=J_{gamma_{k}A}(w^{1,k}-gamma_{k}(L^{*}w^{2,k}+Cw^{1,k}));
%z^{2,k}=J_{gamma_{k}B^{-1}}(w^{2,k}+gamma_{k}Lw^{1,k};
%x^{k+1}=z^{1,k}+gamma_{k}(L^{*}w^{2,k}-L^{*}z^{2,k});
%y^{k+1}=z^{2,k}+gamma_{k}(-Lw^{1,k}+Lz^{1,k}).

%Estimate of SNR
     SNR(i) = 20*log10(norm(x_true(:))/norm(x_true(:)-z1_update(:)));

     t2(i)=etime(clock,t1);   
 

    
     xd = 0.5*norm(A_dir(z1_update)-b)^2 + mu*norm(diff_image(z1_update),1);
     pr = xd - yp;
     RES(i)=norm(pr(:));%%²Ð²î
    
  % Stopping criteria
    if  norm(z1_update(:)-x1(:))/norm(x1(:)) <= tol
         break
    else
        
        x0 = x1;
        y0 = y1;
        x1 = z1_update;
        y1 = z2_update;
        i = i+1;
    end
 end
 %t2 = toc;
 
%end


        
        
        


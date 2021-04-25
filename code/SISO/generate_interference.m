n_chn_realizations = 40;
Nsub = 64;  
interference = zeros(Nsub,Nsub,n_chn_realizations);
whitening = zeros(Nsub,Nsub);
for a = 60:Nsub
    for b = 60:Nsub
        whitening(a,b) = 1;
    end
end
for i = 1:n_chn_realizations
    x1 = 8*(randperm(8)-ones(1,8));
    x2 = 8*(randperm(8)-ones(1,8));
    H_int = LTE_Eva_Model;    
    H_relevant = H_int(1:Nsub,1:8);
    H_int2 = repmat(H_relevant, 1, Nsub/8);
    x = 2*randi(1,Nsub)-1;
    interference(:,:,i) = (sqrt(0.1)*norm(H_org(:,:,i))/norm(H_int2))*H_int2.*x;
    whitening = rand(Nsub,1)>.90; %Change before Jupyter
    interference(:,:,i) = interference(:,:,i).*whitening;
%     for k = 1:2
%         x1_k = x1(k);
%         x2_k = x2(k);
%         whitening(x1_k+1:x1_k+8,x2_k+1:x2_k+8) = ones(8,8);
%     end
%     interference(:,:,i) = interference(:,:,i).*repmat(whitening,1,64); %Change before Jupyter
end
save('interference.mat','interference');
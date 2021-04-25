H_org_noisy = zeros(Nsub,Nsub,Nr,length(SNR_vec),n_chn_realizations);
mse_ls_interference_ci = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse_interference_ci = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse_true_interference_ci = zeros(n_chn_realizations,length(SNR_vec));

for iter = 1:n_chn_realizations
    chn_norm_factor = norm(reshape(H_org(:,1,:,iter),[Nsub,Nr]))/sqrt(Nsub*Nr); 
    noise_matrix = randn(Nsub,Nsub,Nr)+1j*randn(Nsub,Nsub,Nr);

    for k = 1:length(SNR_vec)
        snr = SNR_vec(k);
        std_dev = sqrt(1/power(10,snr/10));
        
        noise = Nr*chn_norm_factor * std_dev * noise_matrix; 
        H_org_noisy(:,:,:,k,iter) = H_org(:,:,:,iter) + noise + interference(:,:,:,iter);

        % Pilots are considered as 1 without loss of any generality
        % LS estimation
        H_ls = H_org_noisy(:,1,:,k,iter);
        diff_ls = H_org(1:Nsub,1,:,iter) - H_ls(1:Nsub,1,:);
        mse_ls_interference_ci(iter,k) = norm(reshape(diff_ls,[Nsub,Nr]))^2/norm(reshape(H_org(1:Nsub,1,:,iter),[Nsub,Nr]))^2;
        
        H_LMMSE = Rhh*((Rhh+chn_norm_factor/power(10,snr/10)*eye(Nsub*Nr))\reshape(H_ls,[Nsub*Nr,1]));
        diff_mmse = reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]) - H_LMMSE;
        mse_mmse_interference_ci(iter,k) = norm(diff_mmse)^2/norm(reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]))^2;
        
        H_2D = reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]);
        interference_2D = reshape(interference(:,1,:,iter),[Nsub*Nr,1]);
        Rhh_true = H_2D*H_2D';
        H_LMMSE_true = Rhh_true*((H_2D*H_2D' + interference_2D*interference_2D'+chn_norm_factor/power(10,snr/10)*eye(Nsub*Nr))\reshape(H_ls,[Nsub*Nr,1]));
        diff_mmse_true = reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]) - H_LMMSE_true;
        mse_mmse_true_interference_ci(iter,k) = norm(diff_mmse_true)^2/norm(reshape(H_org(1:Nsub,1,:,iter),[Nsub*Nr,1]))^2;
    end

noiseless = H_org;
noisy_0dB = H_org_noisy(:,:,:,1,:);
noisy_5dB = H_org_noisy(:,:,:,2,:);
noisy_10dB = H_org_noisy(:,:,:,3,:);  
noisy_15dB = H_org_noisy(:,:,:,4,:);
noisy_20dB = H_org_noisy(:,:,:,5,:);
end

figure
plot(SNR_vec, 10*log10(mean(mse_ls_interference_ci)), '-+','LineWidth',2)
grid
hold
plot(SNR_vec, 10*log10(mean(mse_mmse_interference_ci)), '-+','LineWidth',2)
plot(SNR_vec, 10*log10(mean(mse_mmse_true_interference_ci)), '-+','LineWidth',2)
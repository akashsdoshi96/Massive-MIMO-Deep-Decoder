%clc; clear;
rng('default')
Nsub = 64;
SNR_vec = 0:5:40;
n_chn_realizations = 40;
Rhh = zeros(Nsub);
H_org = zeros(Nsub,Nsub,n_chn_realizations);
H_org_noisy = zeros(Nsub,Nsub,length(SNR_vec),n_chn_realizations);
mse_ls = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse_true = zeros(n_chn_realizations,length(SNR_vec));

for ii = 1:n_chn_realizations
    H_org(:,:,ii) = EPA_Uplink_TDD;
    for i = 1:64
        H_org(i,:,ii) = H_org(i,:,ii)./norm(H_org(i,:,ii));
    end
    %Rhh = Rhh + (1/n_chn_realizations)*H_org(:,1,ii)*H_org(:,1,ii)';
end

for iter = 1:n_chn_realizations 
    chn_norm_factor = norm(H_org(:,1,iter))/sqrt(Nsub); 
    noise_matrix = randn(Nsub,Nsub)+1j*randn(Nsub,Nsub);
    for k = 1:length(SNR_vec)
        snr = SNR_vec(k);
        std_dev = sqrt(1/power(10,snr/10));

        noise = std_dev*(abs(H_org(:,:,iter)).*noise_matrix);
        H_org_noisy(:,:,k,iter) = H_org(:,:,iter) + noise;

        % Pilots are considered as 1 without loss of any generality
        % LS estimation
        
        H_ls = H_org_noisy(:,1,k,iter);
        diff_ls = H_org(:,1,iter) - H_ls(:,1);
        mse_ls(iter,k) = norm(diff_ls)^2/norm(H_org(:,1,iter))^2;

        % MMSE estimation
        Rhh_true = H_org(:,1,iter)*H_org(:,1,iter)';
        H_LMMSE = Rhh*((Rhh+(power(10,-snr/10))*(diag(abs(H_org(:,1,iter))).*eye(Nsub)))\H_ls);
        H_LMMSE_true = Rhh_true*((Rhh_true+(power(10,-snr/10))*(diag(abs(H_org(:,1,iter))).*eye(Nsub)))\H_ls);
        diff_mmse = H_org(1:Nsub,1,iter) - H_LMMSE(1:Nsub,1);
        mse_mmse(iter,k) = norm(diff_mmse)^2/norm(H_org(1:Nsub,1,iter))^2;
        diff_mmse_true = H_org(1:Nsub,1,iter) - H_LMMSE_true(1:Nsub,1);
        mse_mmse_true(iter,k) = norm(diff_mmse_true)^2/norm(H_org(1:Nsub,1,iter))^2;
    end

noiseless = H_org;
noisy_0dB = H_org_noisy(:,:,1,:);
noisy_5dB = H_org_noisy(:,:,2,:);
noisy_10dB = H_org_noisy(:,:,3,:);  
noisy_15dB = H_org_noisy(:,:,4,:);
noisy_20dB = H_org_noisy(:,:,5,:);
noisy_25dB = H_org_noisy(:,:,6,:);
noisy_30dB = H_org_noisy(:,:,7,:);
noisy_35dB = H_org_noisy(:,:,8,:);  
noisy_40dB = H_org_noisy(:,:,9,:);
end

save('H_EVA_64x64_allSNR.mat','noiseless','noisy_0dB','noisy_5dB','noisy_10dB','noisy_15dB','noisy_20dB','noisy_25dB','noisy_30dB','noisy_35dB','noisy_40dB');

figure
semilogy(SNR_vec, mean(mse_ls), '-+','LineWidth',2)
grid
hold
semilogy(SNR_vec, mean(mse_mmse), '-o','LineWidth',2)
semilogy(SNR_vec, mean(mse_mmse_true), '-o','LineWidth',2)
clc;
rng('default')
Nsub = 64;
SNR_vec = 0:5:40;
n_chn_realizations = 40;
Rhh = zeros(Nsub);
H_org = zeros(Nsub,Nsub,n_chn_realizations);
H_org_noisy = zeros(Nsub,Nsub,length(SNR_vec),n_chn_realizations);
Htotal = zeros(1200,12,n_chn_realizations);
mse_ls = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse = zeros(n_chn_realizations,length(SNR_vec));
mse_mmse_true = zeros(n_chn_realizations,length(SNR_vec));

for ii = 1:n_chn_realizations
    Htotal(:,:,ii) = LTE_Eva_Model;
    H(:,:) = Htotal(:,:,ii);
    H_relevant = H(1:Nsub,1:8);
    H_org(:,:,ii) = repmat(H_relevant, 1, Nsub/8);
    Rhh = Rhh + (1/n_chn_realizations)*H_org(1:Nsub,1,ii)*H_org(1:Nsub,1,ii)';
end

% for ii = 1:n_chn_realizations
%     H_org(:,:,ii) = relevant_channel_models;%(1/sqrt(2))*randn(Nsub,Nsub) + sqrt(-1)*(1/sqrt(2))*randn(Nsub,Nsub);
%     Rhh = Rhh + (1/n_chn_realizations)*H_org(1:Nsub,1,ii)*H_org(1:Nsub,1,ii)';
% end

% for ii = 1:n_chn_realizations
%     tap_3_cg = randn(3,1) + 1j*randn(3,1);
%     dft_64 = fft(tap_3_cg,64);
%     H_org(:,:,ii) = repmat(reshape(dft_64,[64,1]),1,64);
%     Rhh = Rhh + (1/n_chn_realizations)*H_org(1:Nsub,1,ii)*H_org(1:Nsub,1,ii)';
% end

for iter = 1:n_chn_realizations 
    chn_norm_factor = norm(H_org(:,1,iter))/sqrt(Nsub); 
    noise_matrix = randn(Nsub,Nsub)+1j*randn(Nsub,Nsub);
    for k = 1:length(SNR_vec)
        snr = SNR_vec(k);
        std_dev = sqrt(1/power(10,snr/10));

        noise = std_dev*chn_norm_factor*noise_matrix;
        H_org_noisy(:,:,k,iter) = H_org(:,:,iter) + noise;

%         % Pilots are considered as 1 without loss of any generality
%         % LS estimation
%         H_ls = H_org_noisy(:,1,k,iter);
%         diff_ls = H_org(1:Nsub,1,iter) - H_ls(1:Nsub,1);
%         mse_ls(iter,k) = norm(diff_ls)^2/norm(H_org(1:Nsub,1,iter))^2;
% 
%         % MMSE estimation
%         Rhh_true = H_org(1:Nsub,1,iter)*H_org(1:Nsub,1,iter)';
%         H_LMMSE = Rhh*((Rhh+(chn_norm_factor/power(10,snr/10))*eye(Nsub))\H_ls);
%         H_LMMSE_true = Rhh_true*((Rhh_true+(chn_norm_factor/power(10,snr/10))*eye(Nsub))\H_ls);
%         diff_mmse = H_org(1:Nsub,1,iter) - H_LMMSE(1:Nsub,1);
%         mse_mmse(iter,k) = norm(diff_mmse)^2/norm(H_org(1:Nsub,1,iter))^2;
%         diff_mmse_true = H_org(1:Nsub,1,iter) - H_LMMSE_true(1:Nsub,1);
%         mse_mmse_true(iter,k) = norm(diff_mmse_true)^2/norm(H_org(1:Nsub,1,iter))^2;
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

save('H_64x64_allSNR.mat','noiseless','noisy_0dB','noisy_5dB','noisy_10dB','noisy_15dB','noisy_20dB','noisy_25dB','noisy_30dB','noisy_35dB','noisy_40dB');

% figure
% semilogy(SNR_vec, mean(mse_ls), '-+','LineWidth',2)
% grid
% hold
% semilogy(SNR_vec, mean(mse_mmse), '-o','LineWidth',2)
% semilogy(SNR_vec, mean(mse_mmse_true), '-o','LineWidth',2)
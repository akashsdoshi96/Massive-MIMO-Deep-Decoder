n_chn_realizations = 20;
Nsub = 64;
Nr = 64;
interference = zeros(Nsub,Nsub,Nr,n_chn_realizations);
whitening1 = zeros(Nsub,Nsub);
whitening2 = zeros(Nsub,Nsub);
for i = 1:n_chn_realizations
    x1_1 = [0 Nsub-8]; x1_2 = [0 Nsub-8];
    x2_1 = [0 Nsub-8]; x2_2 = [0 Nsub-8];
    x1_1_k = randi([1,2],1); x1_2_k = randi([1,2],1);
    x2_1_k = randi([1,2],1);
    equal = 0;
    while equal == 0
        x2_2_k = randi([1,2],1);
        if x1_1_k ~= x2_1_k
            equal = 1;
        elseif x1_2_k ~= x2_2_k
            equal = 1;
        end
    end
    H_int_1 = EPA_Uplink_TDD;
    H_int_2 = EPA_Uplink_TDD;
    x_antenna1 = 2*randi([0,1],Nsub,Nsub)-1;
    x_antenna2 = 2*randi([0,1],Nsub,Nsub)-1;
    whitening1(x1_1(x1_1_k)+1:x1_1(x1_1_k)+8,x1_2(x1_2_k)+1:x1_2(x1_2_k)+8) = ones(8,8);
    whitening2(x2_1(x2_1_k)+1:x2_1(x2_1_k)+8,x2_2(x2_2_k)+1:x2_2(x2_2_k)+8) = ones(8,8);
    x_antenna1 = x_antenna1.*whitening1;
    x_antenna2 = x_antenna2.*whitening2;
    x1 = repmat(reshape(x_antenna1,[Nsub,Nsub,1]),[1,1,Nr]);
    x2 = repmat(reshape(x_antenna2,[Nsub,Nsub,1]),[1,1,Nr]);
    H_norm_org = reshape(H_org(:,:,:,i),[Nsub*Nr,Nsub]);
    H_norm_int1 = reshape(H_int_1,[Nsub*Nr,Nsub]);
    H_norm_int2 = reshape(H_int_2,[Nsub*Nr,Nsub]);
    interference1 = (0.5*norm(H_norm_org)/norm(H_norm_int1))*H_int_1.*x1;
    interference2 = (0.5*norm(H_norm_org)/norm(H_norm_int2))*H_int_2.*x2;
    interference(:,:,:,i) = interference1 + interference2;
end
save('interference.mat','interference');  
function [rho_dynamic, rho_batch, rho_clean, rho_mix, clean_dynamic,clean_batch,clean_mix]...
         =check_cov(check_corr_choice,check_corr_length,rhobar_maxaverage_choice,T,bat_y,dyn_y,clean_sig,mixed_x)

if check_corr_choice==1
	temp=abs(corrcoef(dyn_y'));          disp(['corrcoef dynamic y = ',num2str(temp(1,2))]);
    temp=abs(corrcoef(bat_y'));          disp(['corrcoef batch y = ',num2str(temp(1,2))]);
    temp=abs(corrcoef(clean_sig'));  disp(['corrcoef of clean_sig = ',num2str(temp(1,2))]);
    temp=abs(corrcoef(mixed_x'));    disp(['corrcoef of mixture = ',num2str(temp(1,2))]);

	nny=size(dyn_y,2)-T;
    rho_dyn_array=zeros(2*check_corr_length+1,1);
    rho_bat_array=zeros(2*check_corr_length+1,1);
    rho_clean_array=zeros(2*check_corr_length+1,1);
    rho_mix_array=zeros(2*check_corr_length+1,1);
    for k=-check_corr_length:check_corr_length
        temp=corrcoef(dyn_y(1,T:nny)',dyn_y(2,T+k:nny+k)');           rho_dyn_array(k+1+check_corr_length)=temp(1,2);
        temp=corrcoef(bat_y(1,T:nny)',bat_y(2,T+k:nny+k)');           rho_bat_array(k+1+check_corr_length)=temp(1,2);
        temp=corrcoef(clean_sig(1,T:nny)',clean_sig(2,T+k:nny+k)');   rho_clean_array(k+1+check_corr_length)=temp(1,2);            
        temp=corrcoef(mixed_x(1,T:nny)',mixed_x(2,T+k:nny+k)');       rho_mix_array(k+1+check_corr_length)=temp(1,2);
    end
    if rhobar_maxaverage_choice==1
        rho_dynamic=sum(abs(rho_dyn_array))/(2*check_corr_length+1);
        rho_batch  =sum(abs(rho_bat_array))/(2*check_corr_length+1);
        rho_clean  =sum(abs(rho_clean_array))/(2*check_corr_length+1);
        rho_mix    =sum(abs(rho_mix_array))/(2*check_corr_length+1);
    else
        rho_dynamic=max(abs(rho_dyn_array));
        rho_batch=max(abs(rho_bat_array));
        rho_clean=max(abs(rho_clean_array));
        rho_mix=max(abs(rho_mix_array));
    end
    
    clean_dyn_array=zeros(2*check_corr_length+1,1);
    clean_bat_array=zeros(2*check_corr_length+1,1);
    clean_mix_array=zeros(2*check_corr_length+1,1);
	for i=1:2
        for j=1:2
            for k=-check_corr_length:check_corr_length
                temp=corrcoef(clean_sig(i,T:nny)',dyn_y(j,T+k:nny+k)');      clean_dyn_array(k+1+check_corr_length)=temp(1,2);
                temp=corrcoef(clean_sig(i,T:nny)',bat_y(j,T+k:nny+k)');      clean_bat_array(k+1+check_corr_length)=temp(1,2);
                temp=corrcoef(clean_sig(i,T:nny)',mixed_x(j,T+k:nny+k)');    clean_mix_array(k+1+check_corr_length)=temp(1,2);
            end
            if rhobar_maxaverage_choice==1
                clean_dynamic(i,j)=sum(abs(clean_dyn_array))/(2*check_corr_length+1);
                clean_batch(i,j)=sum(abs(clean_bat_array))/(2*check_corr_length+1);
                clean_mix(i,j)=sum(abs(clean_mix_array))/(2*check_corr_length+1);
            else
                clean_dynamic(i,j)=max(abs(clean_dyn_array));
                clean_batch(i,j)=max(abs(clean_bat_array));
                clean_mix(i,j)=max(abs(clean_mix_array));
            end
        end
	end
end

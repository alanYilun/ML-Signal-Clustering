function z = label_reconstruction(x,y,recon_vec)
   receive_signal1_t0 = x(1);
   receive_signal1_t1 = x(2);
   receive_signal2_t0 = y(1);
   receive_signal2_t1 = y(2);
   for ii = 1:16
    recon_label(:,1) = recon_vec(:,1)*receive_signal1_t0  + recon_vec(:,2)*receive_signal2_t0;
    recon_label(:,2) = recon_vec(:,3)*receive_signal1_t1  + recon_vec(:,4)*receive_signal2_t1;
   end
   z = recon_label;
end
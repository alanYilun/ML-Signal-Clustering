%% This is a QPSK siso rayleigh fading channel with 4-means clustering with label
% 1Tx 1Rx
%% parameters
function ErrorRate = siso_gmm(SNR,packet_num)
fprintf('Running EmGm for SISO channel ... \n');
sample_len = 100;
%packet_num =200;
label1 = 1+1i;
qpsk = [1+1i 1-1i -1+1i -1-1i];
recon_vec = qpsk/label1;
%% transmission
for kSNR = 1:length(SNR)
   varNoise = 2/10^(SNR(kSNR)/10);
   for  npacket = 1:packet_num
       
   % generate QPSK signal
   signal = sign(rand(sample_len-1,1)-0.5) + 1i * sign(rand(sample_len-1,1)-0.5);
   %add label
   s = [label1
       signal];
   % complex gaussian distribution for rayleigh fading channel
   h = normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1);
  
   % AWGN 
   noise = normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1);
   % transmission, received signal y
   y = h*s+noise; centroid = transpose(y(1) .* recon_vec);

   %% EmGm Algorithm
   
   % input format transfer
   ini_centroids = siso_trans_format(centroid);
   X = siso_trans_format(y);
    
   % EmGm algorithm 
   PX = gmm_try(X,ini_centroids); 
   [~,index] = max(PX'); 
   y = [y index']; 
    
    % symbol detection
    s1 = zeros(length(y),1);
    for jj = 1:length(y)
        if y(jj,2) == 1
            s1(jj,1) = 1+1j;
        elseif y(jj,2) ==2
            s1(jj,1) = 1-1j;
        elseif y(jj,2) ==3
            s1(jj,1) = -1+1j;
        else
            s1(jj,1) = -1-1j;
        end
    end
    % estimate received signal according to the 2 centroids
    Error_num = length(find((s-s1)~=0));
    Packet_ErrorRate(npacket,1) = Error_num/(sample_len);
   end
   ErrorRate(kSNR,1) = mean(Packet_ErrorRate);
end
% 
% scatter(real(cluster1),imag(cluster1),'y');
% hold on
% scatter(real(cluster2),imag(cluster2),'g');
% hold on 
% scatter(real(cluster3),imag(cluster3),'r');
% hold on 
% scatter(real(cluster4),imag(cluster4),'b');
% hold off
% figure(2)
% plot(SNR,ErrorRate);
end

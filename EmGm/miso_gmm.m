
%% QPSK; SIMO; Kmeans Clustering; Label Reconstruction with only 1 label; 
% 2 Tx; 1 Rx
function ErrorRate = miso_gmm(SNR,packet_num)
fprintf('Running gmm for MISO channel ... \n');
% set sample length/ packet number/ SNR Reigon
sample_len = 1000;
%packet_num = 10;

%% generate 16 label sets
label_gen = [1-1j 1+1j 1+1j 1+1j -1+1j 1+1j -1-1j 1+1j];
pair_vec = [1 -1 1j -1j];
[Ref,recon_vec] = pre_define(label_gen,pair_vec);

%% main function
for kSNR = 1:length(SNR)
   varNoise = 2/10^(SNR(kSNR)/10);
   for  npacket = 1:packet_num
   %generate QPSK signal
   signal = sign(rand(sample_len,1)-0.5) + 1i * sign(rand(sample_len,1)-0.5); 
   label = label_gen(1:2); % only 1 label(2 symbols) transmitted
   signal(1:2) = label;
   
   % alamouti encoder
   [tx1,tx2] = alamouti_encoder(signal);
    
   % complex gaussian distribution for rayleigh fading channel
   h1 = normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1);
   h2 = normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1);
   
   % generate noise
   noise1 = normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1);
   noise2 = normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1);  
  
   % pre-operations before combiner (used to define initial 16 centroids)
   receive_signal1 = h1 * tx1 + noise1 ;
   receive_signal2 = h2 * tx2 + noise2 ;
   
   % set initial 16 centroids
   centroid = label_reconstruction(receive_signal1,receive_signal2,recon_vec);
  
   % after combiner
   receive_signal = receive_signal1 + receive_signal2;
   r1 = receive_signal(1:2:end); r2 = receive_signal(2:2:end); 
   y = [r1 r2]; % received data pairs
   
   
   %% GMM algorithm 
   % transform the format of input
   X = miso_trans_format(y);
   ini_centroids = miso_trans_format(centroid);
   %y = gmm(y,centroid);
   PX = gmm_try(X, ini_centroids);
   [~,index] = max(PX'); 
   y = [y index']; 
    %% symbol detection
    [s1,S] = symbol_detection(y,Ref);
    
    % Error Detection
    Error_num = length(find((signal-S(:,1))~=0));
    Packet_ErrorRate(npacket,1) = Error_num/(sample_len);
   end
   ErrorRate(kSNR,1) = mean(Packet_ErrorRate);
end
r = [receive_signal S(:,2)];

%% see the clustering of received signal (there are 16 clusters here)
% for uu = 1:16
% eval(['Q',num2str(uu),'=','r(find(r(:,2)==uu))',';']);
% end 
% 
% figure(1);
% for uu = 1:16
% scatter(real(eval(['Q',num2str(uu)])),imag(eval(['Q',num2str(uu)])));
% hold on
% end
% hold off;
% figure(2);
semilogy(SNR,ErrorRate);
end

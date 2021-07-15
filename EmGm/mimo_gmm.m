%% QPSK; MIMO; Kmeans Clustering; Label Reconstruction with only 1 label; 
% 2 Tx; 2 Rx
function ErrorRate = mimo_gmm(SNR,packet_num)
fprintf('Running EmGm for MIMO channel ... \n');
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
       npacket
   %generate QPSK signal
   signal = sign(rand(sample_len,1)-0.5) + 1i * sign(rand(sample_len,1)-0.5); 
   label = label_gen(1:2); % only 1 label(2 symbols) transmitted
   signal(1:2) = label;
   
   % alamouti encoder
   [tx1,tx2] = alamouti_encoder(signal); % first transmit antenna
     
   % complex gaussian distribution for rayleigh fading channel
   for uu = 1:4
   eval(['h',num2str(uu),'=','normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1)',';']);
   end 
   
   % generate noise
   for kk = 1:4
   eval(['noise',num2str(kk),'=','normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1)',';']);
   end 
 
   % pre-operations before combiner (used to define initial 16 centroids)
   receive_signal1 = h1 * tx1 + noise1;
   receive_signal2 = h2 * tx2 + noise2;
   receive_signal3 = h3 * tx1 + noise3;
   receive_signal4 = h4 * tx2 + noise4;
   
   % set initial 16 centroids
   centroid1 = label_reconstruction(receive_signal1,receive_signal2,recon_vec);
   centroid2 = label_reconstruction(receive_signal3,receive_signal4,recon_vec);
   centroid = [centroid1 centroid2];
   % after combiner
   Rx1 = receive_signal1 + receive_signal2;
   Rx2 = receive_signal3 + receive_signal4;
   r11 = Rx1(1:2:end); r12 = Rx1(2:2:end); 
   r21 = Rx2(1:2:end); r22 = Rx2(2:2:end);
   y = [r11 r12 r21 r22]; % received data pairs
   
   %% GMM algorithm 
   % transform input format for GMM algorithm
   ini_centroids = mimo_trans_format(centroid);
   X = mimo_trans_format(y);
   % EM-GMM alogrithm
   PX = gmm_try(X,ini_centroids);
   [~,index] = max(PX'); 
   y = [y index']; 
    %% symbol detection
    [~,S] = symbol_detection(y,Ref);
    
    % Error Detection
    %Packet_ErrorRate = zeros(packet_num,1);
    Error_num = length(find((signal-S(:,1))~=0));
    Packet_ErrorRate(npacket,1) = Error_num/(sample_len);
   end
   %ErrorRate = zeros(length(SNR),1);
   ErrorRate(kSNR,1) = mean(Packet_ErrorRate);
end
%r = [y S(:,2)];

%% see the clustering of rece ived signal (there are 16 clusters here)
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
% semilogy(SNR,ErrorRate);
end

%% This is a QPSK simo rayleigh fading channel with 4-means clustering with label 
%(1Tx, 2Rx)
%% parameters
function ErrorRate = SIMO_Kmeans(SNR,packet_num)
fprintf('Running Kmeans for SIMO channel ... \n');
sample_len = 100;
%packet_num =200;

label1 = 1+1i;
qpsk = [1+1i 1-1i -1+1i -1-1i];
est_vec = qpsk/label1;
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
   h1 = normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1);
   h2 = normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1);
   % AWGN
   noise1 = normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1);
   noise2 = normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1);
   % transmission, received signal y
   Rx1 = h1 * s + noise1;
   Rx2 = h2 * s + noise2;
   R = [Rx1 Rx2];

   %% clustering
   % choose initial 4 centroids with label reconstruction
   centroid1 = transpose(Rx1(1) .* est_vec);
   centroid2 = transpose(Rx2(1) .* est_vec);
   centroid = [centroid1 centroid2];

    y = kmeans_cluster(R,centroid);
    s1 = zeros(length(y),1);
    % symbol detection
    for jj = 1:length(y)
        if y(jj,3) == 1
            s1(jj,1) = 1+1j;
        elseif y(jj,3) ==2
            s1(jj,1) = 1-1j;
        elseif y(jj,3) ==3
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


end

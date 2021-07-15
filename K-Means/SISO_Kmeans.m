%% This is a QPSK siso rayleigh fading channel with 4-means clustering with label
% 1Tx 1Rx
%% parameters
function ErrorRate = SISO_Kmeans(SNR,packet_num)
fprintf('Running Kmeans for SISO channel ... \n');
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
   h = normrnd(0,sqrt(1/2),1,1) + 1i*normrnd(0,sqrt(1/2),1,1);
  
   % AWGN 
   noise = normrnd(0,sqrt(varNoise/2),sample_len,1)+1i*normrnd(0,sqrt(varNoise/2),sample_len,1);
   % transmission, received signal y
   y = h*s+noise;

   %% clustering
   % choose initial 4 centroids with label reconstruction
   centroid = y(1) .* est_vec;
   % kmeans algorithm 
    for p = 1:20 % literation number = 20
    Distance = zeros(4,1);
    for ii = 1:length(y(:,1))
        Distance(1) = abs((y(ii,1))-(centroid(1)));
        Distance(2) = abs((y(ii,1))-(centroid(2)));
        Distance(3) = abs((y(ii,1))-(centroid(3)));
        Distance(4) = abs((y(ii,1))-(centroid(4)));
        y(ii,2) = find(Distance == min(Distance));    %cluster each symbol
    end
    cluster1=y(find(y(:,2)==1),1);
    cluster2=y(find(y(:,2)==2),1);
    cluster3=y(find(y(:,2)==3),1);
    cluster4=y(find(y(:,2)==4),1);
    centroid = [mean(cluster1) mean(cluster2) mean(cluster3) mean(cluster4)]; %update the centroids
    end
    
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

end

sample_len = 100;
%packet_num =200;
SNR = 15;
packet_num = 1;
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
   % kmeans algorithm 
%     for p = 1:20 % literation number = 20
%     Distance = zeros(4,1);
%     for ii = 1:length(y(:,1))
%         Distance(1) = abs((y(ii,:))-(centroid(1,:)));
%         Distance(2) = abs((y(ii,:))-(centroid(2)));
%         Distance(3) = abs((y(ii,:))-(centroid(3)));
%         Distance(4) = abs((y(ii,:))-(centroid(4)));
%         y(ii,2) = find(Distance == min(Distance));    %cluster each symbol
%     end
%     
%     centroid = [mean(cluster1) mean(cluster2) mean(cluster3) mean(cluster4)]; %update the centroids
%     end
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
    cluster1=y(find(y(:,3)==1),1:2);
    cluster2=y(find(y(:,3)==2),1:2);
    cluster3=y(find(y(:,3)==3),1:2);
    cluster4=y(find(y(:,3)==4),1:2);
% X = [real(y(:,1)) imag(y(:,1)) real(y(:,2)) imag(y(:,2))];
% X1 = 
% scatter3(X(:,1),X(:,2),ones(100,1));
% hold on
% scatter3(X(:,3),X(:,4),2*ones(100,1));

z = [1 2];
x = [real(cluster1(:,1)) imag(cluster1(:,1))];
y = [real(cluster1(:,2)) imag(cluster1(:,2))];
plot3(y,z,x,'Color','b','Marker','*','MarkerEdgeColor','b');
hold on;

x = [real(cluster2(:,1)) imag(cluster2(:,1))];
y = [real(cluster2(:,2)) imag(cluster2(:,2))];
plot3(y,z,x,'Color','g','Marker','*','MarkerEdgeColor','g');
hold on;

x = [real(cluster3(:,1)) imag(cluster3(:,1))];
y = [real(cluster3(:,2)) imag(cluster3(:,2))];
plot3(y,z,x,'Color','r','Marker','*','MarkerEdgeColor','r');
hold on;

x = [real(cluster4(:,1)) imag(cluster4(:,1))];
y = [real(cluster4(:,2)) imag(cluster4(:,2))];
plot3(y,z,x,'Color','c','Marker','*','MarkerEdgeColor','c');
hold on;

x = [real(centroid(:,1)) imag(centroid(:,1))];
y = [real(centroid(:,2)) imag(centroid(:,2))];
plot3(y,z,x,'Color','k','Marker','*','MarkerEdgeColor','k','LineWidth',2);
hold on;
xlabel('real');
ylabel('first &second antenna');
zlabel('imag');
grid on;
%  scatter3(real(cluster1(:,1)),ones(length(cluster1),1),imag(cluster1(:,1)),'b');
%  hold on;
% scatter3(real(cluster1(:,2)),2*ones(length(cluster1),1),imag(cluster1(:,2)),'b');
%  hold on;
% 
% x = [real(cluster2(:,1)) imag(cluster2(:,1))];
% y = [real(cluster2(:,2)) imag(cluster2(:,2))];
% plot3(x,z,y,'r');
% hold on;
% x = [real(cluster3(:,1)) imag(cluster3(:,1))];
% y = [real(cluster3(:,2)) imag(cluster3(:,2))];
% plot3(x,z,y,'g');
% hold on;
% x = [real(cluster4(:,1)) imag(cluster4(:,1))];
% y = [real(cluster4(:,2)) imag(cluster4(:,2))];
% plot3(x,z,y,'y');

% hold on;
% scatter3(real(cluster1(:,1)),ones(length(cluster1),1),imag(cluster1(:,1)),'b');
% hold on;
% scatter3(real(cluster1(:,2)),2*ones(length(cluster1),1),imag(cluster1(:,2)),'b');
% hold on;
% scatter3(real(cluster2(:,1)),ones(length(cluster2),1),imag(cluster2(:,1)),'r');
% hold on;
% scatter3(real(cluster2(:,2)),2*ones(length(cluster2),1),imag(cluster2(:,2)),'r');
% hold on;
% scatter3(real(cluster3(:,1)),ones(length(cluster3),1),imag(cluster3(:,1)),'y');
% hold on;
% scatter3(real(cluster3(:,2)),2*ones(length(cluster3),1),imag(cluster3(:,2)),'y');
% hold on;
% scatter3(real(cluster4(:,1)),ones(length(cluster4),1),imag(cluster4(:,1)),'g');
% hold on;
% scatter3(real(cluster4(:,2)),2*ones(length(cluster4),1),imag(cluster4(:,2)),'g');
% hold on;

grid on;



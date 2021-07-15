%% this is K-means algorithm 
% transmit QPSK symbols
% via Rayleigh channel
% SIMO channel is 1x2 (Receive Diversity)
% MISO channel is 2x1 using space-time coding (Alamouti scheme)
% MIMO channel is 2x2 using space-time coding (Alamouti scheme)
% both with label reconstruction: each packet contains 1 label (2 symbols)

%% 
clear;
SNR = transpose(0);
packet_num = 10000;
%% SISO channel
%ErrorRate1 = SISO_Kmeans(SNR,packet_num);
%% SIMO channel
%ErrorRate2 = SIMO_Kmeans(SNR,packet_num);
%% MISO channel 
%ErrorRate3 = MISO_Kmeans(SNR,packet_num);
%% MIMO channel
ErrorRate4 = MIMO_Kmeans(SNR,packet_num); 
%% Save the results in one document
save('Result', 'SNR', 'ErrorRate1', 'ErrorRate2','ErrorRate3','ErrorRate4');
%% Plot the Error Rate for 4 scenarios
Result = load('Result.mat');
fprintf('Plottiing...\n');
for ii = 1:4
semilogy(Result.SNR,(eval(['Result.ErrorRate',num2str(ii)])));
hold on;
end 
legend('siso channel','simo channel','miso channel','mimo channel');
grid on;
hold off;
fprintf('Completed\n');
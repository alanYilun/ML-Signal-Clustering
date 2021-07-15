function [Ref,recon_vec] = pre_define(label_gen,pair_vec)
pair = transpose(label_gen) * pair_vec;
Pair = reshape(pair,32,1);
% there are totally 16 pairs, all can be obtained by inital 1 pair (at least 1 pair of label)
%% generate reference vector for symbol detection
ref_vec = [Pair(1:2:end) Pair(2:2:end)];
[tx1,tx2] = alamouti_encoder(Pair); 
R1 = tx1 + tx2 ; 
Rx = [R1(1:2:end) R1(2:2:end)];
Ref = [Rx ref_vec];

% generate label_reconstruct vector
T = [tx1 tx2];
Tx1t0 = tx1(1:2:end); Tx1t1 = tx1(2:2:end);
Tx2t0 = tx2(1:2:end); Tx2t1 = tx2(2:2:end);
recon_vec1 = Tx1t0/Tx1t0(1); recon_vec2 = Tx2t0/Tx2t0(1);
recon_vec3 = Tx1t1/Tx1t1(1); recon_vec4 = Tx2t1/Tx2t1(1);
recon_vec = [recon_vec1 recon_vec2 recon_vec3 recon_vec4];
end
function y = kmeans_cluster(y,centroid)
 p = 1; % p is iteration times
 x = size(centroid,2);
 while p<15
   Distance = zeros(length(centroid),1);  % cluster by distance to centroid
   for ii = 1:length(y(:,1))
        for jj = 1:length(Distance)
           Distance(jj) = norm((y(ii,1:x))-(centroid(jj,:))); 
        end    
        min_dis = find(Distance == min(Distance));
        y(ii,x+1) = min_dis(randi([1 length(min_dis)],1));    %cluster each symbol
    end
    
    % clustering & update the centroids
    for kk = 1:length(centroid)
       centroid_update(kk,1:x) = mean(y(find(y(:,x+1)==kk),1:x));
    end   
    if mean(norm(centroid_update-centroid)./norm(centroid))<0.001; % to judge whether the centroid has converged
    break   % if converged, stop iteration
    end
    p = p + 1;
    centroid = centroid_update; % update centroids of 16 clusters
 end
end
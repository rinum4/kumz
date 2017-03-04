% S4min = [];
% for i=1:nsz
%     ind = find(S_merge(1,:)==S3min(i));
%     if ~isempty(ind)
%         S4min=[S4min S3min(i) S_merge(2,ind)];
%     else    
%         S4min=[S4min S3min(i)];
%     end
% end

S5min = [];
for i=1:numel(S_out)
    S5min = [S5min ones(1,19)*S_out(i)];
end
S5min = [S5min S4min];
function group_membership_weightedvote(ADJFN, GTGRPFN, thresh, OUTGRPFN)
group = textread(GTGRPFN);
adj = csvread(ADJFN);
N = size(adj,1);
%thresh = 6;

autogroup = zeros(N,1);
prob = zeros(N,1);
maxidx = zeros(N,1);
averagesim = zeros(N,1);
for i = 1:N
    sim = adj(i, :);
    averagesim(i) = mean(sim);
    sim_sort = sort( sim, 'descend' );
   
    group_selected = group(sim >= sim_sort(thresh));
    sim_selected = sim(sim >= sim_sort(thresh));
    vote_1 = sum(sim_selected(group_selected == 1));
    vote_2 = sum(sim_selected(group_selected == 2));
    prob(i) = vote_1/(vote_1+vote_2);
    
    if prob(i) >= 0.5
        autogroup(i) = 1;
    else
        autogroup(i) = 2;
    end
    
    sim(group~=autogroup(i)) = 0;
    [~,maxidx(i)] = max(sim);
    
end

% save
OUT = [autogroup, maxidx, prob, averagesim]';
fileID = fopen(OUTGRPFN,'w');
fprintf(fileID,'%d %d %2.4f %2.4f\n',OUT);
fclose(fileID);
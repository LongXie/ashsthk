function group_membership_MultiAndUnifiTemp(ADJFN, GTGRPFN, thresh, OUTGRPFN, OUTUTFN)
group = textread(GTGRPFN);
adj = csvread(ADJFN);
N = size(adj,1);
if isa(thresh, 'char') || isa(thresh, 'string')
    thresh = str2num(thresh);
end

autogroup = zeros(N,1);
prob = zeros(N,1);
maxidx = zeros(N,1);
averagesim = zeros(N,1);

autogroupUT = zeros(N,1);
UTmaxidx = zeros(N,1);

for i = 1:N
    sim = adj(i, :);
    averagesim(i) = mean(sim);
    sim_sort = sort( sim, 'descend' );
   
    group_selected = group(sim >= sim_sort(thresh));
    sim_selected = sim(sim >= sim_sort(thresh));
    vote_1 = sum(sim_selected(group_selected == 1));
    vote_2 = sum(sim_selected(group_selected == 2));
    prob(i) = vote_1/(vote_1+vote_2);
    
    if prob(i) >= 0.6
        autogroup(i) = 1;
    else
        autogroup(i) = 2;
    end
    
    sim(group~=autogroup(i)) = 0;
    [~,maxidx(i)] = max(sim);
  
    autogroupUT(i) = 1;
    [~,UTmaxidx(i)] = max(sim);
  
end

% save
OUT = [autogroup, maxidx, prob, averagesim]';
fileID = fopen(OUTGRPFN,'w');
fprintf(fileID,'%d %d %2.4f %2.4f\n',OUT);
fclose(fileID);

% save unified template info
OUT = [autogroupUT, UTmaxidx, prob, averagesim]';
fileID = fopen(OUTUTFN,'w');
fprintf(fileID,'%d %d %2.4f %2.4f\n',OUT);
fclose(fileID);

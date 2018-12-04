function accumulation = testsim_evidence_acc(N_reps,num_outcomes,num_samples,num_coherences)


cmap = cool(num_coherences);


for r = 1:N_reps
    accumulation(r,:,:) = testsim_evidence_acc_one_run(num_outcomes,num_samples,num_coherences);
end

for c = 1:num_coherences,
    plot(1:num_samples,mean(accumulation(:,:,c),1),'Color',cmap(c,:)); hold on
end
xlabel('samples');
ylabel('evidence for UP');


function accumulation = testsim_evidence_acc_one_run(num_outcomes,num_samples,num_coherences)

coherence_values = linspace(1/num_outcomes,1,num_coherences); % set of coherence values with which to parameterize sensory uncertainty 

evidence_halfcircle = sin(linspace(pi/2,-pi/2,num_outcomes/2+1));
evidence = [evidence_halfcircle fliplr(evidence_halfcircle(2:end-1))]/num_samples;

for c = 1:length(coherence_values),
    outcomes = randsample(evidence,num_samples,true,[coherence_values(c) repmat((1-coherence_values(c))/(num_outcomes - 1),1,num_outcomes - 1)]);
    accumulation(:,c) = cumsum(outcomes); 
end

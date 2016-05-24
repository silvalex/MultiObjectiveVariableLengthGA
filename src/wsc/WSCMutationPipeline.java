package wsc;

import java.util.Collections;
import java.util.LinkedList;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;

public class WSCMutationPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wscmutationpipeline");
	}

	@Override
	public int numSources() {
		return 1;
	}

	@Override
	public int produce(int min, int max, int start, int subpopulation,
			Individual[] inds, EvolutionState state, int thread) {

		int n = sources[0].produce(min, max, start, subpopulation, inds, state, thread);

        if (!(sources[0] instanceof BreedingPipeline)) {
            for(int q=start;q<n+start;q++)
                inds[q] = (Individual)(inds[q].clone());
        }

        if (!(inds[start] instanceof SequenceVectorIndividual))
            // uh oh, wrong kind of individual
            state.output.fatal("WSCMutationPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds[start]);

        WSCInitializer init = (WSCInitializer) state.initializer;

        // Perform mutation
        for(int q=start;q<n+start;q++) {
        	SequenceVectorIndividual ind = (SequenceVectorIndividual)inds[q];

        	LinkedList<Service> extras = new LinkedList<Service>(init.relevant);
        	Collections.shuffle(extras, init.random);

        	int count = 0;
        	while(count != WSCInitializer.numMutations && !extras.isEmpty()) {
        		Service next = extras.poll();
        		ind.genome.add(0, next);
        	}
        	ind.genome.addAll(extras);
            ind.evaluated=false;
        }
        return n;
	}
}

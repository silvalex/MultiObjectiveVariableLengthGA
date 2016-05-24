package wsc;

import java.util.ArrayList;
import java.util.List;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;

public class WSCCrossoverPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wsccrossoverpipeline");
	}

	@Override
	public int numSources() {
		return 2;
	}

	@Override
	public int produce(int min, int max, int start, int subpopulation,
			Individual[] inds, EvolutionState state, int thread) {

		WSCInitializer init = (WSCInitializer) state.initializer;

		Individual[] inds1 = new Individual[inds.length];
		Individual[] inds2 = new Individual[inds.length];

		int n1 = sources[0].produce(min, max, 0, subpopulation, inds1, state, thread);
		int n2 = sources[1].produce(min, max, 0, subpopulation, inds2, state, thread);

        if (!(sources[0] instanceof BreedingPipeline)) {
            for(int q=0;q<n1;q++)
                inds1[q] = (Individual)(inds1[q].clone());
        }

        if (!(sources[1] instanceof BreedingPipeline)) {
            for(int q=0;q<n2;q++)
                inds2[q] = (Individual)(inds2[q].clone());
        }

        if (!(inds1[0] instanceof SequenceVectorIndividual))
            // uh oh, wrong kind of individual
            state.output.fatal("WSCCrossoverPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds1[0]);

        if (!(inds2[0] instanceof SequenceVectorIndividual))
            // uh oh, wrong kind of individual
            state.output.fatal("WSCCrossoverPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds2[0]);

        int nMin = Math.min(n1, n2);

        // Perform crossover
        for(int q=start,x=0; q < nMin + start; q++,x++) {
    		SequenceVectorIndividual t1 = ((SequenceVectorIndividual)inds1[x]);
    		SequenceVectorIndividual t2 = ((SequenceVectorIndividual)inds2[x]);

    		// Select a random index number for each cut
    		int index1 = init.random.nextInt(t1.genome.size());
    		int index2 = init.random.nextInt(t2.genome.size());

    		// Cut the genomes
    		List<Service> genome1prefix = t1.genome.subList(0, index1);
    		List<Service> genome1suffix = t1.genome.subList(index1, t1.genome.size());
    		List<Service> genome2prefix = t2.genome.subList(0, index2);
    		List<Service> genome2suffix = t2.genome.subList(index2, t2.genome.size());

    		// Create new genomes, by appending other candidate's prefix
    		List<Service> newGenome1 = new ArrayList<Service>();
    		newGenome1.addAll(genome2prefix);
    		newGenome1.addAll(t1.genome);
    		newGenome1.addAll(genome2suffix);
    		List<Service> newGenome2 = new ArrayList<Service>();
    		newGenome2.addAll(genome1prefix);
    		newGenome2.addAll(t2.genome);
    		newGenome2.addAll(genome1suffix);

    		// Replace the old genomes with the new ones
    		t1.genome = newGenome1;
    		t2.genome = newGenome2;

	        inds[q] = t1;
	        inds[q].evaluated=false;

	        if (q+1 < inds.length) {
	        	inds[q+1] = t2;
	        	inds[q+1].evaluated=false;
	        }
        }
        return n1;
	}
}

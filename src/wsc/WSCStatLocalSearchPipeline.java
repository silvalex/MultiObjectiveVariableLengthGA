package wsc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.util.Parameter;

public class WSCStatLocalSearchPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wscstatlocalsearchpipeline");
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
            state.output.fatal("WSCStatLocalSearchPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds[start]);

        WSCInitializer init = (WSCInitializer) state.initializer;

        // Perform local search
        for(int q=start;q<n+start;q++) {
        	SequenceVectorIndividual ind = (SequenceVectorIndividual)inds[q];

        	double bestFitness = ind.fitness.fitness();
        	List<Service> bestNeighbour = ind.genome;

        	// Choose random service to swap with
        	int idx = init.random.nextInt(ind.genome.size());
        	Service chosenServ = ind.genome.get(idx);

        	// Find the index of this service in the local search list
        	int chosenIdx = 0;
        	for (int i = 0; i < init.localSearchList.size(); i++) {
        		if (init.localSearchList.get(i).name.equals(chosenServ.name)) {
        			chosenIdx = i;
        			break;
        		}
        	}

        	SequenceVectorIndividual neighbour = new SequenceVectorIndividual();
        	neighbour.genome = new ArrayList<Service>();
        	neighbour.genome.addAll(init.localSearchList);

        	for (int i = 0; i < neighbour.genome.size(); i++) {
        		// Perform swap
        		Collections.swap(neighbour.genome, chosenIdx, i);

        		// Calculate fitness, and update the best neighbour if necessary
        		neighbour.calculateSequenceFitness(init.numLayers, init.endServ, init, state, true, false);
    			if (ind.fitness.fitness() > bestFitness) {
    				bestFitness = ind.fitness.fitness();
    				bestNeighbour = new ArrayList<Service>(neighbour.genome);
    			}
    			// Swap back (on the original sequence)
    			Collections.swap(neighbour.genome, chosenIdx, i);
        	}

            // Update the tree to contain the best genome found
        	ind.genome = bestNeighbour;
            ind.evaluated = false;
        }
        return n;
	}
}

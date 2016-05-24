package wsc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ec.BreedingPipeline;
import ec.EvolutionState;
import ec.Individual;
import ec.simple.SimpleFitness;
import ec.util.Parameter;

public class WSCLinearLocalSearchPipeline extends BreedingPipeline {

	private static final long serialVersionUID = 1L;

	@Override
	public Parameter defaultBase() {
		return new Parameter("wsclinearlocalsearchpipeline");
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
            state.output.fatal("WSCLinearLocalSearchPipeline didn't get a SequenceVectorIndividual. The offending individual is in subpopulation "
            + subpopulation + " and it's:" + inds[start]);

        WSCInitializer init = (WSCInitializer) state.initializer;

        // Perform local search
        for(int q=start;q<n+start;q++) {
        	SequenceVectorIndividual ind = (SequenceVectorIndividual)inds[q];

        	double bestFitness = ind.fitness.fitness();
        	List<Service> bestNeighbour = ind.genome;

        	List<Service> extras = new ArrayList<Service>(init.relevant);
        	Collections.shuffle(extras, init.random);
        	
        	int extraStart = ind.genome.size();
        	int extraLength = extras.size();
        	int chosen = init.random.nextInt(extraStart);
        	
        	SequenceVectorIndividual neighbour = new SequenceVectorIndividual();
        	neighbour.genome.addAll(ind.genome);
        	neighbour.genome.addAll(extras);

        	for (int i = extraStart; i < extraStart + extraLength; i++) {
        		// Perform swap
        		Collections.swap(neighbour.genome, chosen, i);
        		
        		// Calculate fitness, and update the best neighbour if necessary
        		neighbour.calculateSequenceFitness(init.numLayers, init.endServ, init, state, true);
    			if (ind.fitness.fitness() > bestFitness) {
    				bestFitness = ind.fitness.fitness();
    				bestNeighbour = new ArrayList<Service>(neighbour.genome);
    			}
    			// Swap back
    			Collections.swap(neighbour.genome, chosen, i);
        	}
        	
            // Update the tree to contain the best genome found
        	ind.genome = bestNeighbour;
            ind.evaluated = false;
        }
        return n;
	}
}

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

        	List<Service> neighbour = new ArrayList<Service>();
        	
        	int extraStart = ind.genome.size();
        	int extraLength = extras.size();
        	int chosen = init.random.nextInt(extraStart);

        	for (int i = extraStart; i < extraStart + extraLength; i++) {
        		neighbour.clear();
        		
        		Collections.swap(neighbour, 0, 1); // XXX

        		Set<Service> predecessors = findPredecessors(init, s);
        		neighbour.addAll(0, predecessors);
        		Collections.shuffle(neighbour, init.random);
        		neighbour.addAll(ind.genome);
        		neighbour.addAll(extras);

        		// Calculate fitness, and update the best neighbour if necessary
        		ind.calculateSequenceFitness(init.numLayers, init.endServ, init, state, true);
    			if (ind.fitness.fitness() > bestFitness) {
    				bestFitness = ind.fitness.fitness();
    				bestNeighbour = new ArrayList<Service>(neighbour);
    			}

        	}
            // Update the tree to contain the best genome found
        	ind.genome = bestNeighbour;
            ind.evaluated = false;
        }
        return n;
	}

	public Set<Service> findPredecessors(WSCInitializer init, Service s) {
		Set<Service> predecessors = new HashSet<Service>();

		// Get only inputs that are not subsumed by the given composition inputs (i.e. the start node)
		Set<String> inputsNotSatisfied = init.getInputsNotSubsumed(s.getInputs(), init.startServ.outputs);
		Set<String> inputsToSatisfy = new HashSet<String>(inputsNotSatisfied);

		// If start node is one of the predecessors, add it to set
		if (inputsToSatisfy.size() < s.getInputs().size())
			predecessors.add(init.startServ);

		// Randomly find services to satisfy all remaining inputs
		for (String i : inputsNotSatisfied) {
			if (inputsToSatisfy.contains(i)) {
				List<Service> candidates = init.taxonomyMap.get(i).servicesWithOutput;
				Collections.shuffle(candidates, init.random);

				Service chosen = null;
				candLoop:
				for(Service cand : candidates) {
					if (init.relevant.contains(cand) && cand.layer < s.layer) {
						predecessors.add(cand);
						chosen = cand;
						break candLoop;
					}
				}

				inputsToSatisfy.remove(i);

				// Check if other outputs can also be fulfilled by the chosen candidate, and remove them also
				Set<String> subsumed = init.getInputsSubsumed(inputsToSatisfy, chosen.outputs);
				inputsToSatisfy.removeAll(subsumed);
			}
		}
		return predecessors;
	}
}

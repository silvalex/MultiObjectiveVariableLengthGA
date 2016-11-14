package wsc;

import java.util.Collections;

import ec.*;
import ec.simple.*;

public class WSCProblem extends Problem implements SimpleProblemForm {
	private static final long serialVersionUID = 1L;

	public void evaluate(final EvolutionState state, final Individual ind,
			final int subpopulation, final int threadnum) {
		if (ind.evaluated)
			return;

		if (!(ind instanceof SequenceVectorIndividual))
			state.output.fatal("Whoa!  It's not a SequenceVectorIndividual!!!",
					null);

		SequenceVectorIndividual ind2 = (SequenceVectorIndividual) ind;
		WSCInitializer init = (WSCInitializer) state.initializer;

		//if (!(ind2.fitness instanceof NSGA2MultiObjectiveFitness)) state.output.fatal("Whoa!  It's not a NSGA2MultiObjectiveFitness!!!", null);


		ind2.calculateSequenceFitness(init.numLayers, init.endServ, init, state, false, true);

	}

	@Override
	public void finishEvaluating(EvolutionState state, int threadnum) {
		WSCInitializer init = (WSCInitializer) state.initializer;

		// Update the overall local search list
		Collections.sort(init.localSearchList);


		// Get population
		Subpopulation pop = state.population.subpops[0];

		double minAvailability = 2.0;
		double maxAvailability = -1.0;
		double minReliability = 2.0;
		double maxReliability = -1.0;
		double minTime = Double.MAX_VALUE;
		double maxTime = -1.0;
		double minCost = Double.MAX_VALUE;
		double maxCost = -1.0;

		// Keep track of means
		double meanAvailability = 0.0;
		double meanReliability = 0.0;
		double meanTime = 0.0;
		double meanCost = 0.0;

		// Find the normalisation bounds
		for (Individual ind : pop.individuals) {
			SequenceVectorIndividual wscInd = (SequenceVectorIndividual) ind;
			double a = wscInd.getAvailability();
			double r = wscInd.getReliability();
			double t = wscInd.getTime();
			double c = wscInd.getCost();

			meanAvailability += a;
			meanReliability += r;
			meanTime += t;
			meanCost += c;

			if (WSCInitializer.dynamicNormalisation) {
				if (a < minAvailability)
					minAvailability = a;
				if (a > maxAvailability)
					maxAvailability = a;
				if (r < minReliability)
					minReliability = r;
				if (r > maxReliability)
					maxReliability = r;
				if (t < minTime)
					minTime = t;
				if (t > maxTime)
					maxTime = t;
				if (c < minCost)
					minCost = c;
				if (c > maxCost)
					maxCost = c;
			}
		}

		WSCInitializer.meanAvailPerGen[WSCInitializer.availIdx++] = meanAvailability / pop.individuals.length;
		WSCInitializer.meanReliaPerGen[WSCInitializer.reliaIdx++] = meanReliability / pop.individuals.length;
		WSCInitializer.meanTimePerGen[WSCInitializer.timeIdx++] = meanTime / pop.individuals.length;
		WSCInitializer.meanCostPerGen[WSCInitializer.costIdx++] = meanCost / pop.individuals.length;

		if (WSCInitializer.dynamicNormalisation) {
			// Update the normalisation bounds with the newly found values
			init.minAvailability = minAvailability;
			init.maxAvailability = maxAvailability;
			init.minReliability = minReliability;
			init.maxReliability = maxReliability;
			init.minCost = minCost;
			init.maxCost = maxCost;
			init.minTime = minTime;
			init.maxTime = maxTime;

			// Finish calculating the fitness of each candidate
			for (Individual ind : pop.individuals) {
				((SequenceVectorIndividual)ind).finishCalculatingSequenceFitness(init, state);
			}
		}

	}
}
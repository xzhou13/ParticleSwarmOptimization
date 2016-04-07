#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <sstream>
using namespace std;

class PSO {
	public:
		// constructor that intializes the initial positions, velocities, personal bests,
		// global best and other protected parameters
		PSO (char lbm, char di, int pn, int ne, int np, double in, double cognition, double social, int w, int h, double mv, double me, double localpar);
		// void function that prints out the initial positions, velocities and personal bests
		void printInitial();
		// function that returns the fitness according to problem 1's Q() equation
		double fitness1(double xp, double yp);
		// function that returns the fitness value based on problem 2's Q() equation
		double fitness2(double xp, double yp);
		// calculates the fitness values of all inidividuals and outputs the index of
		// the individual with the highest fitness
		int bestfit();
		// runs the particle swarm optimization (PSO)
		void run();
		// find the index of neighbors for a particle
		vector <int> findLocal(int index);
	protected:
		char localBest;				// lobal best mode
		char decInertia;			// decreasing inertia mode
		int problemNum;				// number of problem with different Q() eqns
		int numEpochs;				// number of maximum iterations
		int numParticles;			// number of particles
		double inertia;				// inertia (range 0-1, typically close to 1)
		double cognitionPar;		// cognition parameter (range 0-4, typically ~2)
		double socialPar;			// social parameter (range 0-4, typically ~2)
		int worldWidth;				// range in x direction [-worldWidth, worldWidth]
		int worldHeight;			// range in y direction
		double maxVelocity;			// maximum velocity
		double maxError;			// max error to run the pso 
		double localPar;			// local best parameter

		FILE * fout;				// file to write errors to
		FILE * fout1;				// file to write position vectors to
		FILE * fout2;				// file to write personal best vectors to
		FILE * fout3;				// file to write number of particles that converged to w/in max error
		stringstream ss;				
		double mdist;							// sqrt(worldWidth^2 + worldHeight^2)/2
		vector < vector <double> > position;	// position vectors
		vector < vector <double> > velocity;	// velocity vectors
		vector < vector <double> > personalBest;// personal best position for each particle
		vector <double> globalBest;				// global best position
};

// constructor that intializes the initial positions, velocities, personal bests,
// global best and other protected parameters
PSO::PSO(char lbm, char di, int pn, int ne, int np, double in, double cognition, double social, int w, int h, double mv, double me, double localpar){
	int bestFitIndex, nc;
	double error_x, error_y, error;
	string str, str1, str2, str3;
	
	// files to output to
	ss << pn;
	str = ss.str()+"errors.csv";
	str1 = ss.str()+"positions.csv";
	str2 = ss.str()+"pbest.csv";
	str3 = ss.str()+"num_convergence.csv";
	fout = fopen(str.c_str(), "w+");
	fout1 = fopen(str1.c_str(), "w+");
	fout2 = fopen(str2.c_str(), "w+");
	fout3 = fopen(str3.c_str(), "w+");

	// initialize parameters
	localBest = lbm;
	decInertia = di;
	problemNum = pn;
	numEpochs = ne;
	numParticles = np;
	inertia = in;
	cognitionPar = cognition;
	socialPar = social;
	worldWidth = w;
	worldHeight = h;
	maxVelocity = mv;
	maxError = me;
	localPar = localpar;
	
	// initialize size of data (position, velocity, personal best) vectors
	position.resize(numParticles);
	velocity.resize(numParticles);
	personalBest.resize(numParticles);
	globalBest.resize(2);
	for (int i=0; i < numParticles; i++){
		position[i].resize(2);
		velocity[i].resize(2);
		personalBest[i].resize(2);
	}
	
	// initialize positions to any point in the world (world width x world height)
	// write to file positions for each individual (iteration 0)
	srand(time(0));
	fprintf(fout1, "%s%c%s\n", "Local best mode,", localBest, ",");
	fprintf(fout1, "%s%c%s\n", "Decreasing mode,", decInertia, ",");
	fprintf(fout1, "%s%d%s\n", "Number of problem,", problemNum, ",");
	fprintf(fout1, "%s%d%s\n", "Number of epochs,", numEpochs, ",");
	fprintf(fout1, "%s%d%s\n", "Number of particles,", numParticles, ",");
	fprintf(fout1, "%s%f%s\n", "Inertia,", inertia, ",");
	fprintf(fout1, "%s%f%s\n", "Cognition parameter,", cognitionPar, ",");
	fprintf(fout1, "%s%f%s\n", "Social parameter,", socialPar, ",");
	fprintf(fout1, "%s%d%s\n", "World width,", worldWidth, ",");
	fprintf(fout1, "%s%d%s\n", "World height,", worldHeight, ",");
	fprintf(fout1, "%s%f%s\n", "Max velocity,", maxVelocity, ",");
	fprintf(fout1, "%s%f%s\n", "Max error,", maxError, ",");
	fprintf(fout1, "%s%f%s\n", "Local parameter,", localPar, ",");
	fprintf(fout1, "%s\n", " ");
	fprintf(fout1, "%s\n", "Iteration 0:,");
	fprintf(fout1, "%s%s%s\n", "Individuals,", "x,", "y,");
	for (int i=0; i < numParticles; i++){
		position[i][0] = rand() % (2*worldWidth) - worldWidth;
		position[i][1] = rand() % (2*worldHeight) - worldHeight;
		fprintf(fout1, "%d%s%f%s%f%s\n", i+1, ",", position[i][0], ",", position[i][1], ",");
		velocity[i][0] = 0.00;
		velocity[i][1] = 0.00;
	}
	fprintf(fout1, "%s\n", " "); 

	fprintf(fout2, "%s%c%s\n", "Local best mode,", localBest, ",");
	fprintf(fout2, "%s%c%s\n", "Decreasing mode,", decInertia, ",");
	fprintf(fout2, "%s%d%s\n", "Number of problem,", problemNum, ",");
	fprintf(fout2, "%s%d%s\n", "Number of epochs,", numEpochs, ",");
	fprintf(fout2, "%s%d%s\n", "Number of particles,", numParticles, ",");
	fprintf(fout2, "%s%f%s\n", "Inertia,", inertia, ",");
	fprintf(fout2, "%s%f%s\n", "Cognition,", cognitionPar, ",");
	fprintf(fout2, "%s%f%s\n", "Social,", socialPar, ",");
	fprintf(fout2, "%s%d%s\n", "World width,", worldWidth, ",");
	fprintf(fout2, "%s%d%s\n", "World height,", worldHeight, ",");
	fprintf(fout2, "%s%f%s\n", "Max velocity,", maxVelocity, ",");
	fprintf(fout2, "%s%f%s\n", "Max error,", maxError, ",");
	fprintf(fout2, "%s%f%s\n", "Local parameter,", localPar, ",");
	fprintf(fout2, "%s\n", " ");
	fprintf(fout2, "%s\n", "Iteration 0:,");
	fprintf(fout2, "%s%s%s\n", "Individuals,", "pbest_x,", "pbest_y,");
	// initialize personal best to its initial position
	for (int i=0; i < numParticles; i++){
		personalBest[i][0] = position[i][0];
		personalBest[i][1] = position[i][1];
		fprintf(fout2, "%d%s%f%s%f%s\n", i+1, ",", personalBest[i][0], ",", personalBest[i][1], ",");
	}
	fprintf(fout2, "%s\n", " "); 
	
	mdist = sqrt(worldWidth*worldWidth + worldHeight*worldHeight)/2;
	
	// find best fit position and set as global best, output it to std output
	bestFitIndex = bestfit();
	globalBest[0] = position[bestFitIndex][0];
	globalBest[1] = position[bestFitIndex][1];

	// print out position, velocity, personal best after initialization
	// to standard output
	cout << "Iteration 0" << endl;
	printInitial();
	
	// calculate/output error (x and y components)
	error_x = error_y = 0;
	for (int i=0; i < numParticles; i++){
		error_x += (position[i][0]-globalBest[0])*(position[i][0]-globalBest[0]);
		error_y += (position[i][1]-globalBest[1])*(position[i][1]-globalBest[1]);
	}
	error_x = sqrt((1.00/(2.00*numParticles))*error_x);
	error_y = sqrt((1.00/(2.00*numParticles))*error_y); 
	cout << "error_x is " << error_x << endl;
	cout << "error_y is " << error_y << endl;
	
	// write to file errors after initialization
	fprintf(fout, "%s%c%s\n", "Local best mode,", localBest, ",");
	fprintf(fout, "%s%c%s\n", "Decreasing mode,", decInertia, ",");
	fprintf(fout, "%s%d%s\n", "Number of problem,", problemNum, ",");
	fprintf(fout, "%s%d%s\n", "Number of epochs,", numEpochs, ",");
	fprintf(fout, "%s%d%s\n", "Number of particles,", numParticles, ",");
	fprintf(fout, "%s%f%s\n", "Inertia,", inertia, ",");
	fprintf(fout, "%s%f%s\n", "Cognition,", cognitionPar, ",");
	fprintf(fout, "%s%f%s\n", "Social,", socialPar, ",");
	fprintf(fout, "%s%d%s\n", "World width,", worldWidth, ",");
	fprintf(fout, "%s%d%s\n", "World height,", worldHeight, ",");
	fprintf(fout, "%s%f%s\n", "Max velocity,", maxVelocity, ",");
	fprintf(fout, "%s%f%s\n", "Max error,", maxError, ",");
	fprintf(fout, "%s%f%s\n", "Local parameter,", localPar, ",");
	fprintf(fout, "%s\n", " ");
	fprintf(fout, "%s%s%s\n", "Iteration,", "Error_x,", "Error_y,");
	fprintf(fout, "%d%s%f%s%f%s\n", 0, ",", error_x, ",", error_y, ",");

	// calculate the scalar value of the error and output it to standard error
	error = sqrt(error_x*error_x + error_y*error_y);
	cout << "Iteration 0: error = " << error << endl << endl;
	cout << endl;

	// calculate and write number of convergence to file
	fprintf(fout3, "%s%c%s\n", "Local best mode,", localBest, ",");
	fprintf(fout3, "%s%c%s\n", "Decreasing mode,", decInertia, ",");
	fprintf(fout3, "%s%d%s\n", "Number of problem,", problemNum, ",");
	fprintf(fout3, "%s%d%s\n", "Number of epochs,", numEpochs, ",");
	fprintf(fout3, "%s%d%s\n", "Number of particles,", numParticles, ",");
	fprintf(fout3, "%s%f%s\n", "Inertia,", inertia, ",");
	fprintf(fout3, "%s%f%s\n", "Cognition,", cognitionPar, ",");
	fprintf(fout3, "%s%f%s\n", "Social,", socialPar, ",");
	fprintf(fout3, "%s%d%s\n", "World width,", worldWidth, ",");
	fprintf(fout3, "%s%d%s\n", "World height,", worldHeight, ",");
	fprintf(fout3, "%s%f%s\n", "Max velocity,", maxVelocity, ",");
	fprintf(fout3, "%s%f%s\n", "Max error,", maxError, ",");
	fprintf(fout3, "%s%f%s\n", "Local parameter,", localPar, ",");
	fprintf(fout3, "%s\n", " ");
	fprintf(fout3, "%s%s\n", "Iterations,", "Number of Convergence,");
	nc = 0;
	for (int i=0; i < numParticles; i++){
		if ((abs(position[i][0]-20.00) <= maxError) && (abs(position[i][1]-7.00) <= maxError)){
			nc++;
		}
	}
	fprintf(fout3, "%d%s%d%s\n", 0, "," , nc, ",");
}

// void function that prints out the initial positions, velocities, and personal
// best positions to standard output
void PSO::printInitial(){
	cout << "Initial positions, fitness, velocities, personal best positions, and global best position:" << endl;
	printf("%-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %10s\n", "Individual", "x_position", "y_position", "fitness", "velocity_x", "velocity_y", "pbest_x", "pbest_y", "gbest_x", "gbest_y");
	for (int i=0; i < numParticles; i++){
		if (problemNum == 1){
		printf("%10d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n", i+1, position[i][0], position[i][1], fitness1(position[i][0], position[i][1]), velocity[i][0], velocity[i][1], personalBest[i][0], personalBest[i][1], globalBest[0], globalBest[1]);
		}
		else{
		printf("%10d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n", i+1, position[i][0], position[i][1], fitness2(position[i][0], position[i][1]), velocity[i][0], velocity[i][1], personalBest[i][0], personalBest[i][1], globalBest[0], globalBest[1]);
		}
	}

	cout << "Iteration 0: global best positions are " << "(" << globalBest[0] << ", " << globalBest[1] << ")" << endl;
	cout << "mdist = " << mdist << endl;

}

// function that output the fitness according to problem 1's fitness function
double PSO::fitness1(double xp, double yp){
	return 100.00*(1.00 - (sqrt((xp-20.00)*(xp-20.00) + (yp-7.00)*(yp-7.00))/(sqrt(worldWidth*worldWidth + worldHeight*worldHeight)/2)));
}

// function that output the fitness according to problem 2's fitness function
double PSO::fitness2(double xp, double yp){
	double pdist, ndist;
	pdist = sqrt((xp - 20)*(xp - 20) + (yp - 7)*(yp - 7));
	ndist = sqrt((xp + 20)*(xp + 20) + (yp + 7)*(yp + 7));
	return 9.00*max(0.00,10.00-pdist*pdist)+10.00*(1.00-pdist/mdist)+70.00*(1.00-ndist/mdist);
}

// function that calculates the fitness of all the positions and output the index
// of the individual with the highest fitness
int PSO::bestfit(){
	double bestFitness;
	int bestIndex;
	double currentFitness;
	
	// problem 1
//	printf("%-10s  %-10s\n", "Individual", "Fitness");
	if (problemNum == 1){	
		bestFitness = fitness1(position[0][0], position[0][1]);
		bestIndex = 0;
//		printf("%10d  %10f\n", 1, bestFitness);
		for (int i=1; i < numParticles; i++){
//			printf("%10d  %10f\n", i+1, fitness1(position[i][0], position[i][1]));
			if (fitness1(position[i][0], position[i][1]) > bestFitness){
				bestFitness = fitness1(position[i][0], position[i][1]);
				bestIndex = i;
			}
		}
	}
	
	// problem 2
	else{
		bestFitness = fitness2(position[0][0], position[0][1]);
		bestIndex = 0;
//		printf("%10d  %10f\n", 1, bestFitness);
		for (int i=1; i < numParticles; i++){
			currentFitness = fitness2(position[i][0], position[i][1]);
//			printf("%10d  %10f\n", i+1, currentFitness);
			if (currentFitness > bestFitness){
				bestFitness = fitness2(position[i][0], position[i][1]);
				bestIndex = i;
			}
		}
	}
	cout << "best fitness value is " << bestFitness << endl;
	cout << "most fit inidividual is " << bestIndex+1 << endl;

	return bestIndex;
}

// run the pso
void PSO::run(){
	double fullInertia, pdec, error, r1, r2, r3, sumVelocity, error_x, error_y, currentFitness, highQ;
	int bestFitIndex, numConverged, ln, hn, high;
	int ni = 0;
	vector <int> locals;
	
	fullInertia = inertia;

	pdec = 1.00;
	// do the following while number of epoches is < numEpoches and error > maxError
	do{
		cout << "Iteration " << ni+1 << endl;
		//cout << "Inertia is " <<  inertia << endl;
	
		// update velocity (includes local best and non-local-best modes)
		for (int i=0; i < numParticles; i++){	
			
			// local best mode to update velocity
			if (localBest == 'y'){
				locals = findLocal(i);
				ln = locals[0];
				hn = locals[1];
				//cout << "low neighbor of " << i+1 << " is " << ln+1 << endl;
				//cout << "high neighbor of " << i+1 << " is " << hn+1 << endl;
				// find neighbor with highest fitness
				if (problemNum == 1){
					highQ = fitness1(position[i][0], position[i][1]);
					high = i;
					if (fitness1(position[ln][0], position[ln][1]) > highQ){
						//cout << "individual " << ln+1 << " has a higher Q than individual " << i+1 << endl;
						high = ln;
						highQ = fitness1(position[ln][0], position[ln][1]);
					}
					if (fitness1(position[hn][0], position[hn][1]) > highQ){
						//cout << "individual " << hn+1 << " has a higher Q than individual " << i+1 << " and/or individual " << ln+1 <<  endl;
						high = hn;
						highQ = fitness1(position[hn][0], position[hn][1]);
					}
					//cout << "particle " << i+1 << " at (" << position[i][0] << ", " << position[i][1] << ") has the most fit neighbor as individual " << high+1 << " with fitness " << highQ << " at (" << position[high][0] << ", " << position[high][1] << ")" << endl;
				}else{
					highQ = fitness2(position[i][0], position[i][1]);
					high = i;
					if (fitness2(position[ln][0], position[ln][1]) > highQ){
						high = ln;
						highQ = fitness2(position[ln][0], position[ln][1]);
					}
					if (fitness2(position[hn][0], position[hn][1]) > highQ){
						high = hn;
						highQ = fitness2(position[hn][0], position[hn][1]);
					}
					//cout << "particle " << i+1 << " at (" << position[i][0] << ", " << position[i][1] << ") has the most fit neighbor at " << high+1 << "th position with fitness value " << highQ << " at (" << position[high][0] << ", " << position[high][1] << ")" << endl;
				}
				
				// update velocity (local best mode)
				r1 = drand48();
				r2 = drand48();	
				r3 = drand48();
				velocity[i][0] = inertia*velocity[i][0] + cognitionPar*r1*(personalBest[i][0]-position[i][0]) + socialPar*r2*(globalBest[0]-position[i][0]) +localPar*r3*(position[high][0]-position[i][0]);
				velocity[i][1] = inertia*velocity[i][1] + cognitionPar*r1*(personalBest[i][1]-position[i][1]) + socialPar*r2*(globalBest[1]-position[i][1]) + localPar*r3*(position[high][1]-position[i][1]);
			}

			// non local best mode to update velocity
			else{
				r1 = drand48();
				r2 = drand48();	
				velocity[i][0] = inertia*velocity[i][0] + cognitionPar*r1*(personalBest[i][0]-position[i][0]) + socialPar*r2*(globalBest[0]-position[i][0]);
				velocity[i][1] = inertia*velocity[i][1] + cognitionPar*r1*(personalBest[i][1]-position[i][1]) + socialPar*r2*(globalBest[1]-position[i][1]);
			}

			// rescale velocity if necessary
			sumVelocity = sqrt(velocity[i][0]*velocity[i][0] + velocity[i][1]*velocity[i][1]);
			if (abs(velocity[i][0]) > (maxVelocity*maxVelocity)){
				velocity[i][0] = (maxVelocity/sumVelocity)*velocity[i][0];
			}
			if (abs(velocity[i][1]) > (maxVelocity*maxVelocity)){
				velocity[i][1] = (maxVelocity/sumVelocity)*velocity[i][1];
			}

		}
			
		// decrease inertia if in that mode
		if (decInertia == 'y'){
			pdec = fullInertia/(numEpochs-1);
			inertia -= pdec;
			cout << "new inertia is " << inertia << endl;
		}
		
		fprintf(fout1, "%s%d\n", "Iteration ", ni+1);
		fprintf(fout1, "%s%s%s\n", "Individuals,", "x,", "y,");
		for (int i=0; i < numParticles; i++){
			// update position
			position[i][0] += velocity[i][0];
			position[i][1] += velocity[i][1];
			fprintf(fout1, "%d%s%f%s%f%s\n", i+1, ",", position[i][0], ",", position[i][1], ",");
		}

		// update personal best and global best
		for (int i=0; i < numParticles; i++){
			if (problemNum == 1){
				currentFitness = fitness1(position[i][0], position[i][1]);
				// update personal best
				if (currentFitness > fitness1(personalBest[i][0], personalBest[i][1])){
					personalBest[i] = position[i];
				}
				// update global best
				if (currentFitness > fitness1(globalBest[0], globalBest[1])){
					globalBest[0] = position[i][0];
					globalBest[1] = position[i][1];
				}
			}
			else{
				currentFitness = fitness2(position[i][0], position[i][1]);
				// update personal best
				if (currentFitness > fitness2(personalBest[i][0], personalBest[i][1])){
					personalBest[i] = position[i];
				}
				// update global best
				if (currentFitness > fitness2(globalBest[0], globalBest[1])){
					globalBest[0] = position[i][0];
					globalBest[1] = position[i][1];
				}
			}
		}
			
		// 1) calculate number of particles converged to within the maximum error
		// 2) print updated position, velocity, personal best, and global best to
		// 3) standard output; write personal best and number of convergence to file
		fprintf(fout2, "%s%d\n", "Iteration ", ni+1);
		fprintf(fout2, "%s%s%s\n", "Individuals,", "gbest_x,", "gbest_y,");
		printf("%-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n", "Individual", "x_position", "y_position", "fitness", "velocity_x", "velocity_y", "pbest_x", "pbest_y", "gbest_x", "gbest_y");
		numConverged = 0;
		for (int i=0; i < numParticles; i++){
			if ((abs(position[i][0]-20.00) <= maxError) && (abs(position[i][1]-7.00) <= maxError)){
				numConverged++;
			}
			printf("%10d  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f  %10f\n", i+1, position[i][0], position[i][1], fitness1(position[i][0], position[i][1]), velocity[i][0], velocity[i][1], personalBest[i][0], personalBest[i][1], globalBest[0], globalBest[1]);
			// write personal best to file
			fprintf(fout2, "%d%s%f%s%f%s\n", i+1, ",", personalBest[i][0], ",", personalBest[i][1], ",");
		}
		fprintf(fout3, "%d%s%d%s\n", ni+1, ",", numConverged, ",");

		// calculate errors and print to standard output
		error_x = error_y = 0;
		for (int i=0; i < numParticles; i++){
			error_x += (position[i][0]-globalBest[0])*(position[i][0]-globalBest[0]);
			error_y += (position[i][1]-globalBest[1])*(position[i][1]-globalBest[1]);
		}
		error_x = sqrt((1.00/(2.00*numParticles))*error_x);
		error_y = sqrt((1.00/(2.00*numParticles))*error_y); 
		cout << "error_x is " << error_x << endl;
		cout << "error_y is " << error_y << endl;

		// write iteration number and errors to file
		fprintf(fout, "%d%s%f%s%f%s\n", ni+1, ",", error_x, ",", error_y, ",");
		
		// calculate scalar error and print it on standard output
		error = sqrt(error_x*error_x + error_y*error_y);
		cout << "Iteration " << ni+1 << ": scalar error = " << error << endl;
		
		// call bestfit() to print out global best to standard output
		bestfit();
		cout << "iteration " << ni+1 << ": global best positions are " << "(" << globalBest[0] << ", " << globalBest[1] << ")" << endl << endl;
		
		// write global best to file with position vectors
		fprintf(fout1, "%s%f%s%f%s\n", "global best,", globalBest[0], ",", globalBest[1], ",");
		fprintf(fout1, "%s\n", " "); 
		
		// increment number of iterations
		ni++;

	}while((error > maxError) &&  (ni < numEpochs));
	
	// close files
	fclose(fout);
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
}

vector <int> PSO::findLocal(int index){
	vector <int> local;
	int ln, hn, numNeighbors;
	
	numNeighbors = 2;
	local.resize(2);

	if ((index - (numNeighbors-1)) < 0){
		ln = numParticles - 1;	
	}else{
		ln = index - 1;
	}

	if ((index + (numNeighbors-1)) > (numParticles-1)){
		hn = 0;
	}else{
		hn = index + 1;
	}

	local[0] = ln;
	local[1] = hn;
	return local;
}

int main(){
	char lm, di;
	int nn, pn, ne, np, w, h;
	double in, cognition, social, mv, me, localpar;
	
	// read in if it's local best mode, and number of neighbors
	cin >> lm;

	// read in if it's decreasing inertia mode
	cin >> di;
	
	// read in problem number
	cin >> pn;

	// read in other parameters and output to standard output
	cin >> ne >> np >> in >> cognition >> social >> w >> h >> mv >> me >> localpar;
	if (lm == 'y'){
		cout << "Local Best Mode" << endl;
	}

	if (di == 'y'){
		cout << "Decrease Inertia Mode" << endl;
	}
	cout << "Problem number: " << pn << endl;
	cout << "Number of epochs (iterations): " << ne << endl;
	cout << "Number of particles: " << np << endl;
	cout << "Inertia: " << in << endl;
	cout << "Cognition parameter: " << cognition << endl;
	cout << "Social parameter: " << social << endl;
	cout << "Local parameter: " << localpar << endl;
	cout << "World width: " << w << endl;
	cout << "World height: " << h << endl;
	cout << "Maximum velocity: " << mv << endl;
	cout << "Minimum error: " << me << endl << endl;
	
	// initialize instance of PSO and run
	PSO pso(lm, di, pn, ne, np, in, cognition, social, w, h, mv, me, localpar);
	pso.run();
	
	return 0;
}


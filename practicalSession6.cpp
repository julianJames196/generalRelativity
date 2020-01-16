//A programme to calculate trategory in spacetime according to General Relativity
//Based of the paper http://old.hanno-rein.de/download/schwarzschild.pdf
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

//A global time variable
double time = 0.0;

struct planet {
	double mass;
	double position[2];
	double velocity[2];
};

double centerMass = 0.1;
double schwarzRadius = centerMass * 2;

planet point;

double dist(double x, double y) {
	return sqrt(pow(x, 2.0) + pow(y, 2.0));
}

double angularVelocity(double x, double y, double vx, double vy) {
	double theta = atan2(y, x);
	double distance = dist(x, y);
	return ( vy*cos(theta) - vx*sin(theta) ) / distance;
}

double schwarzPotential(double mass, double x, double y, planet p) {
	double thetaV = angularVelocity(x,y,p.velocity[0],p.velocity[1]);
	if (x == 0 || y == 0) {
		cout << "A collision with the central mass has occured";
		return 0.0;
	}
	else {
		double distance = dist(x, y);
		//The l^2/2m^2 given in the paper i.e the rotational energy factor
		double factor = 0.5 * pow(distance,4) * pow(thetaV, 2); 
		double potential = (-1 * centerMass / distance) 
			+ (factor / pow(distance,2)) 
			- (factor *schwarzRadius / pow(distance, 2));
		return potential;
	}
}

//Calculate derviatives of the potential to get the force
double forceX(planet p) {
	double diffStep = 0.0001;
	return (schwarzPotential(p.mass, p.position[0]+diffStep, p.position[1], p) 
		- schwarzPotential(p.mass, p.position[0], p.position[1], p) ) / diffStep;
}

double forceY(planet p) {
	double diffStep = 0.0001;
	return (schwarzPotential(p.mass, p.position[0], p.position[1] + diffStep, p)
		- schwarzPotential(p.mass, p.position[0], p.position[1], p) ) / diffStep;
}

planet update(planet p, double timeStep) {
	//Takes a planet p and updates it by 1 timestep using the leapfrog method
	double increament = 0.5 * timeStep;
	//Go through each planet and calculate the corresponding gravitional force from it
	p.position[0] += (increament * p.velocity[0]);
	p.position[1] += (increament * p.velocity[1]);
	time += increament;
	//Using the force to calculate the new velocities. 
	p.velocity[0] += (timeStep * -1* forceX(p) / p.mass); //The minus is because the Force = - grad(F)
	p.velocity[1] += (timeStep * -1 *forceY(p) / p.mass);

	p.position[0] += (increament * p.velocity[0]);
	p.position[1] += (increament * p.velocity[1]);
	time += increament;

	return p;
}

int main(){
	cout << "This programme calculates spacetime trategory using the the Schwarzchld Metric\n";
	cout << "Natural Units are used so G=1 and c=1\n";
	cout << "The center mass is: " + to_string(centerMass) + " hence Schwarzchild radius = " + to_string(schwarzRadius) + "\n";
	double x, y, vx, vy , mass;
	cout << "Input x-position: ";
	cin >> x;
	cout << "Input y-position: ";
	cin >> y;
	cout << "Input x-velocity: ";
	cin >> vx;
	cout << "Input y-velocity: ";
	cin >> vy;
	cout << "Input mass: ";
	cin >> mass;
	//Declaring
	point = { mass, {x, y}, {vx, vy} };
	double timeStep = 0.0001;
	double closestApproach = dist(x, y);
	//Running the simulation initally
	while (time < 10) {
		point = update(point, timeStep);

		if (dist(point.position[0], point.position[1]) < closestApproach) {
			closestApproach = dist(point.position[0], point.position[1]); 
		}
	}
	//Checks for stability with initally conditions
	double checkerVar[2] = { x, y };

	while (time <= 100) {
		point = update(point, timeStep);
		//Calculate the closest approach
		if (dist(point.position[0], point.position[1]) < closestApproach) {
			closestApproach = dist(point.position[0], point.position[1]);
		}
		//Compare with checker var
		//X-Position within 4%
		if (abs((checkerVar[0] - point.position[0]) / checkerVar[0]) < 0.04) {
			//Y-Position within 4%
			if (abs((checkerVar[1] - point.position[1]) < 0.04) / checkerVar[1]) {
				cout << "\nThe planet has returned to the same position at time:" + to_string(time) + "\n";
			}
		}

		//Has the object flang out to space?
		if (abs(point.position[0]) > 100000 && abs(point.velocity[0]) > 1000) {
			cout << "Will likely continue on to infinity due to high velocity and large position"; break;
		}

		if (abs(point.position[1]) > 1000000 && abs(point.velocity[1]) > 1000) {
			cout << "Will likely continue on to infinity due to high velocity and large position"; break;
		}
	}
	//Final output
	cout << "\nAfter time = " + to_string(time) + " which is " + to_string(time / timeStep) + " iterations";
	cout << "\nX = " + to_string(point.position[0]);
	cout << "\nY = " + to_string(point.position[1]);
	cout << "\nVx = " + to_string(point.velocity[0]);
	cout << "\nVy = " + to_string(point.velocity[1]);
	cout << "\nClosest Approach was " + to_string(closestApproach);
}
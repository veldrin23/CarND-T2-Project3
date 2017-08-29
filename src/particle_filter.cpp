/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1500;

	std::default_random_engine gen; //black magic random number generator

	// gaussian noise for x, y and yaw
	std::normal_distribution<double> x_gaus(x, std[0]);
	std::normal_distribution<double> y_gaus(y, std[1]);
	std::normal_distribution<double> theta_gaus(theta, std[2]);


	// create particles
	for (int i = 0; i < num_particles; ++i) 
	{
		Particle p;
		p.id = i;
		p.weight = 1;
		p.x = x_gaus(gen);
		p.y = y_gaus(gen);
		p.theta = theta_gaus(gen);

		particles.push_back(p);
	}

	weights.resize(num_particles);
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	for (int i = 0; i < num_particles; ++i) {

		std::default_random_engine gen;
		std::normal_distribution<double> x_gaus(particles[i].x, std_pos[0]);
		std::normal_distribution<double> y_gaus(particles[i].y, std_pos[1]);
		std::normal_distribution<double> theta_gaus(particles[i].theta, std_pos[2]);

		particles[i].x = x_gaus(gen);
		particles[i].y = y_gaus(gen);


		if (fabs(yaw_rate) > 1e-6) {
			// yaw rate != 0 
			particles[i].x = particles[i].x + velocity/yaw_rate *
						(sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y = particles[i].y + velocity/yaw_rate *
						(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
			particles[i].theta = particles[i].theta + yaw_rate*delta_t;
		} else {  
			// (yaw rate == 0 or veleeelllly velly close to)
			particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			particles[i].theta = theta_gaus(gen);
		}
	}
} 


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(int i = 0; i < observations.size(); ++i) {

		double c_dist = 0.0f; // current measurement
		double m_dist = 0.0f; // min measurement
		int m_dist_ind = 0; // index of min measurement 

		// from heler_functions.h: dist(double x1, double y1, double x2, double y2) 
		for(int j = 0; i < predicted.size(); ++j) {

			c_dist = dist(observations[i].x, observations[i].y,
				predicted[j].x, predicted[j].y);

			if(j == 0) {
				m_dist = c_dist;
				}
			else {
				if(c_dist < m_dist) m_dist = c_dist;
			}

}
	observations[i].id = m_dist_ind;
	}



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	std::vector<LandmarkObs> trans_observations(observations.size());
	std::vector<LandmarkObs> landmark_observations(map_landmarks.landmark_list.size());

	for (int i = 0; i < landmark_observations.size(); ++i) 
	{
		landmark_observations[i].id = 0;
		landmark_observations[i].x = map_landmarks.landmark_list[i].x_f;
		landmark_observations[i].y = map_landmarks.landmark_list[i].y_f;
	}


	double p_x, p_y, l_x, l_y, norm_x, norm_y, new_w, denom;

	// translate from vehicle's coordinates to map coordinates
	for (int i = 0; i < particles.size(); ++i) 
	{
		p_x = particles[i].x;
		p_y = particles[i].y;
		for(int j = 0; j < observations.size(); ++j)
		{
			trans_observations[j].x = observations[j].x*cos(particles[i].theta) - observations[i].y*sin(particles[i].theta) + p_x;
			trans_observations[j].y = observations[j].x*sin(particles[i].theta) + observations[i].y*cos(particles[i].theta) + p_y;
		}
			double new_weight = 1.0f;
	for(auto& trans: trans_observations)
	{
		l_x = landmark_observations[trans.id].x;
		l_y = landmark_observations[trans.id].y;

		denom = 1.0/(2.0*M_PI*std_landmark[0]*std_landmark[1]);
		norm_x = (pow(trans.x-l_x, 2)) / (2*pow(std_landmark[0],2)); 
		norm_y = (pow(trans.y-l_y, 2)) / (2*pow(std_landmark[1],2));
		new_w = denom * exp(-1*(norm_x + norm_y));
		new_weight *= new_w;
	}
	particles[i].weight = new_weight;
	weights[i] = new_weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::vector<Particle> new_particles;
	std::default_random_engine gen;
	std::discrete_distribution<> dist(weights.begin(), weights.end());

	int new_index;

	for(int i = 0; i < num_particles; ++i) 
	{
		new_index = dist(gen);
		new_particles.push_back(particles[new_index]);

	}
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

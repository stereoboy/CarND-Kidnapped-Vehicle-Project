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
	const double std_x = std[0];
	const double std_y = std[1];
	const double std_theta = std[2];
	num_particles = 100; // Controllable Factor

	std::default_random_engine gen;
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);


	for (int i = 0; i < num_particles; ++i) {
		// TODO: Sample  and from these normal distrubtions like this:
		//	 sample_x = dist_x(gen);
		//	 where "gen" is the random engine initialized earlier.
		struct Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = fmod(dist_theta(gen) + (2*M_PI), (2*M_PI));
		particle.weight = 1;
		particles.push_back(particle);
		//cout << particle.x << "," << particle.y << endl;
	}

  weights.resize(num_particles);

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	const double std_x = std_pos[0];
	const double std_y = std_pos[1];
	const double std_theta = std_pos[2];

	std::default_random_engine gen;
	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

  double x, y, theta;
	// Lesson 12 Motion Model
	if (yaw_rate == 0)
	//if (abs(yaw_rate) <= e-6 )
  {
    // x = x + v(dt)*cos(theta)
    // y = y + v(dt)*sin(theta)
    // theta = theta
    for (int i = 0; i < particles.size(); i++)
    {
      x = particles[i].x;
      y = particles[i].y;
      theta = particles[i].theta;

      particles[i].x = x + velocity*delta_t*cos(theta);
      particles[i].x += dist_x(gen);
      particles[i].y = y + velocity*delta_t*sin(theta);
      particles[i].y += dist_y(gen);
      particles[i].theta = theta;
      particles[i].theta += dist_theta(gen);
      particles[i].theta = fmod(particles[i].theta + (2*M_PI), (2*M_PI));
    }
  }
  else
  {
    // x = x + v/yaw_rate*(sin(theta + yaw_rate(dt)) - sin(theta))
    // y = y + v/yaw_rate*(cos(theta)  - cos(theta + yaw_rate(dt)))
    // theta = theta  + yaw_rate(dt)
    for (int i = 0; i < particles.size(); i++)
    {
      x = particles[i].x;
      y = particles[i].y;
      theta = particles[i].theta;

      particles[i].x = x + velocity/yaw_rate*(sin(theta + yaw_rate*delta_t) - sin(theta));
      particles[i].x += dist_x(gen);

      particles[i].y = y + velocity/yaw_rate*(-cos(theta + yaw_rate*delta_t) + cos(theta));
      particles[i].y += dist_y(gen);

      particles[i].theta = theta + yaw_rate*delta_t;
      particles[i].theta += dist_theta(gen);
      particles[i].theta = fmod(particles[i].theta + (2*M_PI), (2*M_PI));
    }
  }

  weights.resize(num_particles);
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  for (int i = 0; i < observations.size(); i++)
  {
    double min_dist = 0xFFFFFFF;
    for (int j = 0; j < predicted.size(); j++)
    {
      double len = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if (len < min_dist)
      {
        min_dist = len;
        //observations[i].id = predicted[j].id;
        observations[i].id = j;
      }
    }
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

  const double std_x = std_landmark[0];
  const double std_y = std_landmark[1];

  for (int i = 0; i < particles.size(); i++)
  {
    const double x = particles[i].x;
    const double y = particles[i].y;
    const double theta = particles[i].theta;

    vector<LandmarkObs> predicted;
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      const int m_id = map_landmarks.landmark_list[j].id_i;
      const double m_x  = map_landmarks.landmark_list[j].x_f;
      const double m_y  = map_landmarks.landmark_list[j].y_f;
      const double len = dist(x, y, m_x, m_y);

      if ( len <= sensor_range)
      {
        LandmarkObs landmark;
        landmark.id = m_id;
        landmark.x = m_x;
        landmark.y = m_y;
        predicted.push_back(landmark);
      }

    }
    if (predicted.size() == 0)
    {
      weights[i] = 0;
    }
    else
    {
      vector<LandmarkObs> observations_map;
      observations_map.resize(observations.size());
      for (int j = 0; j < observations_map.size(); j++)
      {
        observations_map[j].x = x + observations[j].x*cos(theta) - observations[j].y*sin(theta);
        observations_map[j].y = y + observations[j].x*sin(theta) + observations[j].y*cos(theta);
      }

      dataAssociation(predicted, observations_map);

      double w = 0;
      for (int j = 0; j < observations_map.size(); j++)
      {
        int id_o = observations_map[j].id;
        const double x_o = observations_map[j].x;
        const double y_o = observations_map[j].y;
        const double x_p = predicted[id_o].x;
        const double y_p = predicted[id_o].y;
        double temp = (x_o - x_p)*(x_o - x_p)/(std_x*std_x) + (y_o - y_p)*(y_o - y_p)/(std_y*std_y);
        w += temp;
      }

      double weight  = exp(-1.0/2.0*w);
      weights[i] = weight;
      particles[i].weight = weight;
      //cout << weights[i] <<endl;
    }
  }

  double total_w = 0;
  for (int i = 0; i < weights.size(); i++)
  {
    total_w += weights[i];
  }

  for (int i = 0; i < weights.size(); i++)
  {
    weights[i] /= total_w;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::default_random_engine gen;
  std::discrete_distribution<> d(weights.begin(), weights.end());

  vector<Particle> particles_new;
  particles_new.resize(particles.size());

  for (int i = 0; i < particles.size(); i++)
  {
    particles_new[i] = particles[d(gen)];
  }

  particles = particles_new;
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

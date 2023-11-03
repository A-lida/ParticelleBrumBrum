#include <array>
#include <cmath>
#include <iostream>

#include "Headers.hpp"
#include "TApplication.h"
#include "TDirectory.h"
#include "TRandom.h"
#include "TRint.h"
#include "particles/Particle.hpp"
#include "particles/ParticleType.hpp"
#include "particles/ResonanceType.hpp"

int main(int argc, char** argv) {
// object ownership setup
#ifdef PROGRAM_USE_LOCAL_OWNERSHIP
  TDirectory::AddDirectory(kFALSE);
#endif
  // create root application
  APP_TYPE app(APP_NAME, &argc, argv);

  Particle::AddParticleType("Pion+", 0.13957, 1, 0.);
  Particle::AddParticleType("Pion-", 0.13957, -1, 0.);
  Particle::AddParticleType("Kaon+", 0.49367, 1, 0.);
  Particle::AddParticleType("Kaon-", 0.49367, -1, 0.);
  Particle::AddParticleType("Proton+", 0.93827, +1, 0.);
  Particle::AddParticleType("Proton-", 0.93827, -1, 0.);
  Particle::AddParticleType("Kaon*", 0.89166, 0, 0.050);
  // Particle::DeallocateMemory();

  gRandom->SetSeed(56);
  int Events = 100000;
  int N = 120;
  Particle* EventParticles[N];
  double phi;
  double theta;
  double p;

  void ConversionPoltoCart(double r, double t, double p, Particle& P) {
    double Px = r * sin(t) * cos(p);
    double Py = r * sin(t) * sin(p);
    double Pz = r * cos(t);
    P.set_P(Px, Py, Pz);
  };

  double probabilityPionPlus = 0.4;
  double probabilityPionMinus = 0.4;
  double probabilityKaonPlus = 0.05;
  double probabilityKaonMinus = 0.05;
  double probabilityProtonPlus = 0.045;
  double probabilityProtonMinus = 0.045;
  double probabilityKStar = 0.01;

  for (int i = 0; i < Events; i++) {
    double randomValue = randomGenerator.Rndm();

    double charge = 0.0;
    double mass = 0.0;

    for (int j = 0; j < 100; ++j) {
      phi = TRandom::Uniform(1 / (2 * M_PI));
      theta = TRandom::Uniform(1 / M_PI);
      p = TRandom::Exp(1.);
      ConversionPoltoCart(p, theta, phi, EventParticles[j]);

      if (randomValue < probabilityPionPlus) {
        EventParticles[j]->set_fIndex("Pion+");

      } else if (randomValue < probabilityPionPlus + probabilityPionMinus) {
        EventParticles[j]->set_fIndex("Pion-");

      } else if (randomValue < probabilityPionPlus + probabilityPionMinus +
                                   probabilityKaonPlus) {
        EventParticles[j]->set_fIndex("Kaon+");

      } else if (randomValue < probabilityPionPlus + probabilityPionMinus +
                                   probabilityKaonPlus + probabilityKaonMinus) {
        EventParticles[j]->set_fIndex("Kaon-");

      } else if (randomValue < probabilityPionPlus + probabilityPionMinus +
                                   probabilityKaonPlus + probabilityKaonMinus +
                                   probabilityProtonPlus) {
        EventParticles[j]->set_fIndex("Proton+");

      } else if (randomValue < probabilityPionPlus + probabilityPionMinus +
                                   probabilityKaonPlus + probabilityKaonMinus +
                                   probabilityProtonPlus +
                                   probabilityProtonMinus) {
        EventParticles[j]->set_fIndex("Proton-");

      } else {
        EventParticles[j]->set_fIndex("Kaon*");
        double x->TRandom::Uniform(1);
        if (x < 0.5) {
          EventParticles[j].Particle::Decay2body(
              Particle("Kaon+", 0.49367, 1, 0),
              Particle("Pion-", 0.13957, -1, 0));
        } else {
          EventParticles[j].Particle::Decay2body(
              Particle("Kaon-", 0.49367, -1, 0),
              Particle("Pion+", 0.13957, 1, 0));
        }
      }
    }
  }
  // do stuff

  // run application
  app.Run();

  return 0;
}

#ifndef PARTICLE_HPP
#define PARTICLE_HPP
//#include <array>
#include <cmath>
#include <iostream>

#include "ParticleType.hpp"
#include "ResonanceType.hpp"

class Particle
{
private:
  const char *f_Name;
  static const int fMaxNumParticleType = 10;
  static ParticleType *fParticleType[fMaxNumParticleType];
  static int fNParticleType;

  int fIndex;

  double fPx;
  double fPy;
  double fPz;

  static int FindParticle(const char *Name);

  void Boost(double bx, double by, double bz);

public:
  Particle();
  Particle(const char *, double{0}, double{0}, double{0});

  int get_fIndex();

  double const get_Xmomentum();
  double const get_Ymomentum();
  double const get_Zmomentum();
  double get_Mass() const;

  void set_fIndex(int index);
  void set_fIndex(const char *pname);
  void set_P(double px, double py, double pz);

  double const ParticleEnergy();
  double InvMass(Particle &p);

  static void ArrayStatus(); // può essere statico? si
  void ParticleStatus();

  static void AddParticleType(const char *, const double, const int,
                              const double);
  static void DeallocateMemory();
  int Decay2body(Particle &dau1, Particle &dau2) const;
};
#endif

#include "particles/Particle.hpp"

//#include <array> non serve, noi usiamo gli array nativi di c non la classe in c++
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "particles/ParticleType.hpp"
#include "particles/ResonanceType.hpp"

int const Particle::fMaxNumParticleType;
int Particle::fNParticleType = 0;
ParticleType *Particle::fParticleType[fMaxNumParticleType];

int Particle::FindParticle(const char *Name) {
  for (int i = 0; i < fMaxNumParticleType; ++i) {
    if (fParticleType[i]->get_ParticleName() == Name) {
      return i;
    }
  }
  std::cout << "No correspondence, sorry";
  return -1;
}
Particle::Particle(){};
Particle::Particle(const char *parName, double Px, double Py, double Pz)
    : f_Name(parName), fPx(Px = 0), fPy(Py = 0), fPz(Pz = 0) {
  fIndex = FindParticle(parName);
};

int Particle::get_fIndex() { return fIndex; }

double const Particle::get_Xmomentum() { return fPx; }
double const Particle::get_Ymomentum() { return fPy; }
double const Particle::get_Zmomentum() { return fPz; }

double Particle::get_Mass() const {
  return fParticleType[fIndex]->get_ParticleMass();
}

double const Particle::ParticleEnergy() {
  double const m = fParticleType[fIndex]->get_ParticleMass();
  double P = fPx * fPx + fPy * fPy + fPz * fPz;
  return sqrt(m * m + P);
}

double Particle::InvMass(Particle &p) {
  const double E1 = ParticleEnergy();
  const double E2 = p.ParticleEnergy();
  const double PX = fPx + p.get_Xmomentum();
  const double PY = fPy + p.get_Ymomentum();
  const double PZ = fPz + p.get_Zmomentum();

  return sqrt((E1 + E2) * (E1 + E2) - (PX + PY + PZ) * (PX + PY + PZ));
}

void Particle::set_P(double px, double py, double pz) {
  fPx = px;
  fPy = py;
  fPz = pz;
}

void Particle::ArrayStatus()  // può essere statico? non lo so
{
  for (int i = 0; i < fNParticleType; ++i) {
    ParticleType p = *fParticleType[i];
    p.print();
    std::cout << std::endl;
  }
}

void Particle::ParticleStatus() {
  std::cout << "Particle code: " << get_fIndex() << std::endl << std::endl;
  std::cout << "Name: " << f_Name << std::endl;
  std::cout << "Momentum: "
            << "(" << get_Xmomentum() << ", " << get_Ymomentum() << " ,"
            << get_Zmomentum() << ")" << std::endl;
}

void Particle::AddParticleType(const char *name, const double Mass,
                               const int Charge, const double Width) {
  if (FindParticle(name) == -1 && fNParticleType < fMaxNumParticleType) {
    for (int i = fNParticleType; i < fMaxNumParticleType; ++i) {
      fParticleType[i] = new ResonanceType(name, Mass, Charge, Width);
      ++fNParticleType;
    }
  }

  else {
    std::cout << "The array is already full, sorry bro" << std::endl;
  }
}

void Particle::set_fIndex(int index) {
  if (get_fIndex() == index) {
    std::cout << "The particle is already here." << std::endl;
  } else {
    fIndex = index;
  }
}

void Particle::set_fIndex(const char *pname) {
  if (FindParticle(pname) != -1) {
    std::cout << "The particle is already here." << std::endl;
  } else {
    fIndex = FindParticle(pname);
  }
}

void Particle::DeallocateMemory() {
  for (int i = 0; i < fNParticleType; ++i) {
    delete fParticleType[i];
  }
}

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (get_Mass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = get_Mass();
  double massDau1 = dau1.get_Mass();
  double massDau2 = dau2.get_Mass();

  if (fIndex > -1) {  // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += fParticleType[fIndex]->get_ParticleWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf(
        "Decayment cannot be preformed because mass is too low in this "
        "channel\n");
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.set_P(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
             pout * cos(theta));
  dau2.set_P(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
             -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz) {
  double energy = ParticleEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}

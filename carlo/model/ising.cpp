#include "ising.hpp"
#include <array>
#include <iostream>

namespace carlo::model
{

ising::ising(std::size_t dimension) :
  _dimension{dimension}
{
  data.resize(dimension * dimension);
}

ising::spin ising::operator[](std::size_t x, std::size_t y) const
{
  return static_cast<spin>( ((int)data[y * _dimension + x] * 2) - 1 );
}

void ising::set(std::size_t x, std::size_t y, spin val)
{
  data[y * _dimension + x] = static_cast<bool>( ((int)val + 1) / 2 );
}

void ising::flip(std::size_t x, std::size_t y)
{
  data[y * _dimension + x] = !data[y * _dimension + x];
}

void ising::scatter() 
{
  for (std::size_t i = 0; i < data.size(); i++)
    data[i] = rd();
}

ising::random_generator::random_generator() :
  dist(0, 1)
{ }

ising::spin ising::random_generator::operator() () 
{
  return static_cast<spin>( dist(generator) );
}

}

namespace carlo::method
{

metropolis::metropolis(double J, double T, std::size_t dimension) :
  _J{J},
  beta{1.0 / (/*1.380649e-23 */ T)},
  base(dimension),
  rd(0, 3),
  rfd(0, 1.0)
{ 
  _model.scatter();

  // Calculate the energy
  double E = 0.0;
  for (std::size_t y = 0; y < _model.size(); y++)
    for (std::size_t x = 0; x < _model.size(); x++)
    {
      const std::array<std::pair<int, int>, 4> offsets = {
        { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } }
      };

      const auto this_spin = _model[x, y];

      for (auto& off : offsets)
      {
        const std::pair<std::size_t, std::size_t> pos = 
        {
          ( x + off.first  ) % _model.size(),
          ( y + off.second ) % _model.size()
        };

        E += this_spin * _model[ pos.first, pos.second ];
      }
    }

  _E = -0.5 * _J * E;
}

void metropolis::method_step() 
{
  // Metropolis step logic
  // Go through each site and randomly flip a neighbor
  const auto inner = [this](std::size_t x, std::size_t y)
  {
    const std::array<std::pair<int, int>, 4> offsets = {
      { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } }
    };

    // We're trying to flip at (x, y) so we need to sum up all the surrounding values
    double s  = _model[x, y];
    double nb = 0.0;
    for (const auto& off : offsets)
    {
      const std::pair<std::size_t, std::size_t> pos = 
      {
        ( (int)x + off.first  ) % _model.size(),
        ( (int)y + off.second ) % _model.size()
      };
      nb += (int)_model[pos.first, pos.second];
    } 

    const auto dE = 2.0 * s * nb;
    if (dE <= 0.0)
      _model.flip(x, y);
    else
    {
      // If rand pass then flip
      if (rfd(gen) < exp( -1.0 * dE * beta ))
        _model.flip(x, y);
    }
  };

  // Written like this so that we can then parallelize it
  for (std::size_t y = 0; y < _model.size(); y++)
    for (std::size_t x = 0; x < _model.size(); x++)
      inner(x, y);

  /*
  for (std::size_t y = 1; y < _model.size(); y += 2)
    for (std::size_t x = 1; x < _model.size(); x += 2)
      inner(x, y);*/
}

}

std::string operator*(carlo::model::ising::spin s)
{
  switch (s)
  {
  case carlo::model::ising::up:    return "u";
  case carlo::model::ising::down:  return "d";
  }
}
#include <iostream>

#include "carlo/model/ising.hpp"
#include "carlo/image.hpp"

int main()
{
  using namespace carlo;
  
  gif_export<method::metropolis> mtr("test.gif", 2.0, 1.8, 256);
  mtr->add_attachment<model::total_mag>();

  for (int i = 0; i < 100; i++)
  {
    //if (i % 100 == 0)
      //draw_state("test" + std::to_string(i) + ".bmp", mtr.model());
    mtr.step();
  }
  mtr.write();

  mtr->get_attachment<model::total_mag>()->show();

  /*
  gif_export<method::metropolis> exp("test.gif", 5.0, 0.05, 256);
  for (int i = 0; i < 1000; i++) exp.step();
  exp.write();*/
}
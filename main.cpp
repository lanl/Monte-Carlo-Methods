#include <iostream>

#include "carlo/model/ising.hpp"
#include "carlo/image.hpp"

int main()
{
  using namespace carlo;
  
  // Construct the metropolis method for use with the ising model
  // wrapped in a gif_export that will write the state out to a .gif file
  // at the end
  auto mtr = gif_export<model::ising>::make<method::metropolis>("test.gif", 2.0, 1.8, 256);

  // Add a total_magentization attachment which writes the data 
  // out to total_mag.png
  mtr->add_attachment<model::total_mag>();

  // Go through the method
  for (int i = 0; i < 100; i++) mtr.step();

  // Write the .gif file
  mtr.write();

  // Export the attachment data
  mtr->get_attachment<model::total_mag>()->write("total_mag.png");
}
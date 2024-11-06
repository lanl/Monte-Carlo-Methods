#include <iostream>

#include "carlo/model/ising.hpp"
#include "carlo/image.hpp"

#include <sstream>
#include <iostream>

static void draw_progress_bar(int length, int current, int max_iterations)
{
  std::stringstream ss;
  const auto completion = (float)current / (float)max_iterations;

  ss << "[";
  for (int i = 0; i < length; i++)
    ss << ((float)i / (float)length <= completion ? '=' : ' ');
  ss << "] " << current << "/" << max_iterations << " " << (int)(completion * 100.f) << "%";

  std::cout << ss.str() << "\r";
  std::cout.flush();
}

int main()
{
  using namespace carlo;
  
  // Construct the metropolis method for use with the ising model
  // wrapped in a gif_export that will write the state out to a .gif file
  // at the end
  auto mtr = gif_export<model::ising>::make<method::metropolis>("test.gif", 2.0, 1.8, 256);

  // Add a total_magentization attachment which writes the data 
  // out to total_mag.png
  mtr->add_attachment<model::characteristic_length>("l.gif");

  // Go through the method
  for (int i = 0; i < 100; i++) 
  {
    draw_progress_bar(40, i, 100);
    mtr.step();
  }
  std::cout << "\n";

  // Write the .gif file
  mtr.write();

  // Export the attachment data
  //mtr->get_attachment<model::characteristic_length>()->write("total_mag.png");
}
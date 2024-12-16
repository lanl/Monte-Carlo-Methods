#include <iostream>

#include "carlo/model/ising.hpp"
#include "carlo/image.hpp"

#include <sstream>
#include <chrono>
#include <iomanip>

static void draw_progress_bar(int length, int current, int max_iterations, std::string extra = "")
{
  std::stringstream ss;
  const auto completion = (float)current / (float)max_iterations;

  ss << "[";
  for (int i = 0; i < length; i++)
    ss << ((float)i / (float)length <= completion ? '=' : ' ');
  ss << "] " << current << "/" << max_iterations << " " << (int)(completion * 100.f) << "%";
  if (extra.length()) ss << " (" << extra << ")                    ";

  std::cout << ss.str() << "\r";
  std::cout.flush();
}

int main()
{
  using namespace carlo;
  using namespace std::chrono;
  
  // Construct the metropolis method for use with the ising model
  // wrapped in a gif_export that will write the state out to a .gif file
  // at the end
  std::vector<std::shared_ptr<method::metropolis>> methods;

  for (int i = 0; i < 10; i++)
    methods.emplace_back(std::make_shared<method::metropolis>(2.0, 1.8, 256));

  // Add a total_magentization attachment which writes the data 
  // out to total_mag.png
  for (auto& method : methods)
  {
    method->add_attachment<model::total_mag>();
    method->add_attachment<model::characteristic_length>();
  }
  // Go through the method
  constexpr auto iterations = 200U;

  std::vector<float> iteration_times;

  const auto now = high_resolution_clock::now();
  auto last = high_resolution_clock::now();
  for (int i = 0; i < iterations; i++) 
  {
    std::stringstream extra;
    if (iteration_times.size())
    {
      const auto elapsed = duration<float, seconds::period>(high_resolution_clock::now() - now).count();
      // remaining equals average dt * iterations - elapsed
      float avg = 0.f;
      for (const auto& it : iteration_times)
        avg += it;
      avg /= iteration_times.size();

      extra << std::fixed << std::setprecision(2) << avg * iterations - elapsed << " seconds remaining";
    }

    draw_progress_bar(40, i + 1, iterations, extra.str());
    for (auto& method : methods)
      method->step();
      
    iteration_times.push_back(duration<float, seconds::period>(high_resolution_clock::now() - last).count());
    last = high_resolution_clock::now();
  }

  for (int i = 0; i < methods.size(); i++)
  {
    draw_progress_bar(40, i + 1, methods.size(), "Writing outputs...");
    methods[i]->get_attachment<model::total_mag>()->write("total_mag-" + std::to_string(i) + ".png");
    methods[i]->get_attachment<model::characteristic_length>()->write("c_length-" + std::to_string(i) + ".png");
  }
  std::cout << "\n";

  reduce<model::total_mag>(methods).write("total_mag_temporal_avg.png");
  reduce<model::characteristic_length>(methods).write("c_length_temporal_avg.png");
}
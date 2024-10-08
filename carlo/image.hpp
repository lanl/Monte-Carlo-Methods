#pragma once

#include <memory>
#include <string>
#include <vector>

#include "model/ising.hpp"

namespace carlo
{

template<typename T>
void draw_state(const std::string& file, const T& state);

struct digit
{
  digit();

  void blit_number(uint32_t number, std::size_t x, std::size_t y, std::vector<std::uint8_t>& frame, uint32_t frame_width);

private:
  std::size_t width, height, digit_width;
  std::vector<uint8_t> image_data;
};

template<typename Model>
struct gif_export
{
  method::base<Model>* operator->() const
  {
    return method.get();
  }

  template<typename Method, typename... Args>
  inline static gif_export make(const std::string& filename, Args&&... args)
  {
    static_assert(std::is_base_of_v<method::base<Model>, Method>);
    gif_export g;
    g.method = std::make_shared<Method>(std::forward<Args>(args)...);
    g._construct_writer(filename);
    return g;
  } 

  void write();
  void step();

private:
  gif_export() = default;
  /*
  template<typename... Args>
  gif_export(const std::string& filename, Args&&... args) :
    method( std::make_unique<M>(std::forward<Args>(args)...) )
  { _construct_writer(filename); }*/

  void _construct_writer(const std::string& filename);

  std::size_t frame_number = 0;
  digit _digit;
  std::vector<std::uint8_t> frame;
  std::shared_ptr<method::base<Model>> method;
  std::shared_ptr<void> gif_writer;
};

}
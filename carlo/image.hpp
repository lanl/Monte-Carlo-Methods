#pragma once

#include <memory>
#include <string>
#include <vector>

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

template<typename M>
struct gif_export
{
  template<typename... Args>
  gif_export(const std::string& filename, Args&&... args) :
    method( std::make_unique<M>(std::forward<Args>(args)...) )
  { _construct_writer(filename); }

  M* operator->() const
  {
    return method.get();
  }

  void write();
  void step();

private:
  void _construct_writer(const std::string& filename);

  std::size_t frame_number = 0;
  digit _digit;
  std::vector<std::uint8_t> frame;
  std::unique_ptr<M> method;
  std::shared_ptr<void> gif_writer;
};

}
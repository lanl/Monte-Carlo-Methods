#include "image.hpp"
#include "model/ising.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image_write.h>
#include <stb_image.h>
#include <cassert>
#include <cstring>

#include <gif.h>

namespace carlo
{

template<>
void draw_state(const std::string& file, const model::ising& state)
{
  std::vector<uint8_t> data(state.size() * state.size());
  for (std::size_t y = 0; y < state.size(); y++)
    for (std::size_t x = 0; x < state.size(); x++)
      data[y * state.size() + x] = 255 * static_cast<int>(state[x, y]);

  assert(stbi_write_bmp(file.c_str(), state.size(), state.size(), 1, data.data()));
}

digit::digit()
{
  // Load in the image
  int x, y, n;
  uint8_t* data = stbi_load("../res/images/digits.bmp", &x, &y, &n, 1);
  assert(data);

  width = x;
  digit_width = width / 10;
  height = y;
  image_data.resize(x * y);
  std::memcpy(image_data.data(), data, x * y);
  STBI_FREE(data);
}

void digit::blit_number(uint32_t number, std::size_t x, std::size_t y, std::vector<std::uint8_t>& frame, uint32_t frame_width)
{
  const auto rendering = std::to_string(number);

  for (const auto& c : rendering)
  {
    std::size_t c_offset = c - '0';
    for (int d_x = 0; d_x < digit_width; d_x++)
      for (int d_y = 0; d_y < height; d_y++)
      {
        const auto pixel_data = image_data[d_y * width + (c_offset * digit_width + d_x)];

        if (!pixel_data)
        {
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 0] = 255;
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 1] = 0;
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 2] = 0;
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 3] = 255;
        }
        else
        {
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 0] = 0;
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 1] = 255;
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 2] = 255;
          frame[((y + d_y) * frame_width + (x + d_x)) * 4 + 3] = 255;
        }
      }

    x += digit_width + 1;
  }
}

// Do the gif writing, then execute the step method
template<>
void gif_export<method::metropolis>::step()
{
  if (!frame.size())
    frame.resize(method->model().size() * method->model().size() * 4);
  
  for (std::size_t y = 0; y < method->model().size(); y++)
    for (std::size_t x = 0; x < method->model().size(); x++)
    {
      const auto val = 255 * static_cast<int>( method->model()[x, y] );
      frame[(y * method->model().size() + x) * 4 + 0] = val;
      frame[(y * method->model().size() + x) * 4 + 1] = val;
      frame[(y * method->model().size() + x) * 4 + 2] = val;
      frame[(y * method->model().size() + x) * 4 + 3] = 255;
    }

  // Draw frame number on top
  _digit.blit_number(frame_number, 2, 2, frame, method->model().size());

  GifWriteFrame(
    reinterpret_cast<GifWriter*>(gif_writer.get()), 
    frame.data(), 
    method->model().size(), 
    method->model().size(), 
    15
  );
  
  method->step();

  frame_number++;
}

template<>
void gif_export<method::metropolis>::write()
{
  GifEnd(reinterpret_cast<GifWriter*>(gif_writer.get()));
}

// Construct the gif writer
template<>
void gif_export<method::metropolis>::_construct_writer(const std::string& filename)
{
  auto* writer = new GifWriter();
  gif_writer = std::shared_ptr<void>(writer, [](void* ptr) { delete reinterpret_cast<GifWriter*>(ptr); });
  GifBegin(writer, filename.c_str(), method->model().size(), method->model().size(), 15);
}

}
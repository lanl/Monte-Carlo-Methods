#pragma once

#include <gif.h>

#include <matplotlibcpp.h>

#include <fftw3.h>
#include <random>
#include <vector>
#include <optional>
#include <cassert>
#include <iostream>
#include <memory>
#include <typeindex>
#include <unordered_map>

namespace carlo
{
namespace data
{
  template<typename Model>
  struct attachment
  {
    attachment(const Model&) {}
    virtual void step(const Model&) = 0;
  };
}

namespace model
{
  struct ising
  {
    enum spin : int { down = -1, up = 1 };

    ising(std::size_t dimension);

    spin operator()(std::size_t x, std::size_t y) const;
    void set(std::size_t x, std::size_t y, spin val);
    void flip(std::size_t x, std::size_t y);
    void scatter();

    constexpr auto size() const { return _dimension; }

  private:
    const std::size_t _dimension;
    std::vector<bool> data;

    struct random_generator
    { 
      random_generator();
      spin operator() ();

    private:
      std::mt19937 generator;
      std::uniform_int_distribution<int> dist;
    };

    inline static random_generator rd;
  };

  struct characteristic_length : data::attachment<ising>
  {
    characteristic_length(const ising& model, const std::string& gif_name) :
      attachment(model),
      gif(new GifWriter())
    { 
      GifBegin(gif.get(), gif_name.c_str(), model.size(), model.size(), 15);
    }

    ~characteristic_length()
    {
      GifEnd(gif.get());
    }

    void step(const ising& model) override
    {
      // Construct a representation of the model
      std::vector<double> model_rep(model.size() * model.size() * 2);
      for (int y = 0; y < model.size(); y++)
        for (int x = 0; x < model.size(); x++)
        {
          model_rep[(y * model.size() + x * 2) + 0] = static_cast<int>(model(x, y));
          model_rep[(y * model.size() + x * 2) + 1] = 0;
        }
        
      // Take the fourier transform
      std::vector<double> out(model_rep.size() * 2);
      {
        const auto plan = fftw_plan_dft_2d(
          model.size(), 
          model.size(), 
          reinterpret_cast<fftw_complex*>(model_rep.data()), 
          reinterpret_cast<fftw_complex*>(out.data()), 
          FFTW_FORWARD, 
          0
        );

        fftw_execute(plan);

        // Take the absolute value squared of all of the entries 
        for (int y = 0; y < model.size(); y++)
          for (int x = 0; x < model.size(); x++)
          {
            const auto& real = out[(y * model.size() + x) * 2 + 0];
            const auto& imag = out[(y * model.size() + x) * 2 + 1];
            model_rep[(y * model.size() + x) * 2 + 0] = real*real + imag*imag;
            model_rep[(y * model.size() + x) * 2 + 1] = 0;
          }
      }

      std::vector<double> values(model.size() * model.size());
      {
        // Compute the inverse fourier transform
        const auto plan = fftw_plan_dft_2d(
          model.size(), 
          model.size(), 
          reinterpret_cast<fftw_complex*>(model_rep.data()), 
          reinterpret_cast<fftw_complex*>(out.data()), 
          FFTW_BACKWARD, 
          0
        );
        
        fftw_execute(plan);
        
        // Keep only the real values
        for (int y = 0; y < model.size(); y++)
          for (int x = 0; x < model.size(); x++)
          {
            const auto& real = out[(y * model.size() + x) * 2 + 0];
            const auto& imag = out[(y * model.size() + x) * 2 + 1];
            values[y * model.size() + x] = real;// + imag*imag;
          }
      }

      const auto write_gif_slide = [this, &model](const auto& data)
      {
        const auto max_val = *std::max_element(data.begin(), data.end());
        const auto min_val = *std::min_element(data.begin(), data.end());

        // for debug, convert to gif slide
        std::vector<uint8_t> gif_slide(data.size() * 4);
        for (int y = 0; y < model.size(); y++)
          for (int x = 0; x < model.size(); x++)
          {
            const auto pixel_val = (data[y * model.size() + x] - min_val) / (max_val - min_val) * 255;
            gif_slide[(y * model.size() + x) * 4 + 0] = pixel_val;
            gif_slide[(y * model.size() + x) * 4 + 1] = pixel_val;
            gif_slide[(y * model.size() + x) * 4 + 2] = pixel_val;
            gif_slide[(y * model.size() + x) * 4 + 3] = 255;
          }

        GifWriteFrame(gif.get(), gif_slide.data(), model.size(), model.size(), 15);
      };

      write_gif_slide(values);


      // Calculate characteristic length
      // Take the inverse fourier transform of this
      // Take the real part of this
      // Display this data

      // Cr = real.(ifft(abs.(fft(sys)) .^ 2))
      // integrate_theta(normalize(circshift(Cr, (L / 2, L / 2))), Int(L / 2), Int(L / 2), L / 2, 0.25)

      current_step++;
    }

    void write(const std::string& filename)
    {
      // Do matplotlib with steps and lengths
    }

  private:
    std::size_t current_step = 0;

    std::unique_ptr<GifWriter> gif;

    std::vector<std::size_t> steps;
    std::vector<float> lengths;
  };
  
  struct total_mag : data::attachment<ising>
  {
    void step(const ising& model) override
    {
      int total_mag = 0;
      for (int y = 0; y < model.size(); y++)
        for (int x = 0; x < model.size(); x++)
          total_mag += static_cast<int>(model(x, y));
        
      steps.push_back(current_step);
      total_mag_steps.push_back(total_mag);

      current_step++;
    }

    void write(const std::string& filename)
    {
      namespace plt = matplotlibcpp;

      // Generate matplotlib image
      plt::plot(steps, total_mag_steps);
      plt::title("Total Magnetization");
      plt::xlabel("Iteration");
      plt::ylabel("Magnetization");
      plt::save(filename);
    }

  private:
    std::size_t current_step = 0;

    std::vector<std::size_t> steps;
    std::vector<int> total_mag_steps;
  };
}

// The methods own the above models
namespace method
{
  template<typename Model>
  struct base
  {
    template<typename... Args>
    base(Args&&... args) : 
      _model(std::forward<Args>(args)...)
    { }

    virtual ~base() = default;

    template<typename Att, typename... Args>
    void add_attachment(Args&&... args)
    {
      static_assert(std::is_base_of_v<data::attachment<Model>, Att>);
      assert(!attachments.count(std::type_index(typeid(Att))));
      attachments[std::type_index(typeid(Att))] = std::make_shared<Att>(_model, std::forward<Args>(args)...);
    }

    template<typename Att>
    std::shared_ptr<Att> get_attachment()
    {
      assert(attachments.count(std::type_index(typeid(Att))));
      return std::reinterpret_pointer_cast<Att>(attachments[std::type_index(typeid(Att))]);
    }

    const Model& model() const { return _model; }

    void step()
    {
      method_step();
      for (const auto& [ _, a ] : attachments)
        a->step(model());
    }

  protected:
    virtual void method_step() = 0;

    Model _model;

  private:
    std::unordered_map<std::type_index, std::shared_ptr<data::attachment<Model>>> attachments;
  };

  struct metropolis : base<model::ising>
  {
    metropolis(double J, double T, std::size_t dimension);

  protected:
    void method_step() override;

  private:
    double _E;
    const double _J, beta;

    std::mt19937 gen;
    std::uniform_real_distribution<double> rfd;
    std::uniform_int_distribution<int> rd;
  };

  struct block_gibbs : base<model::ising>
  {
    block_gibbs(double beta, std::size_t dimension);
  };
}

}

std::string operator*(carlo::model::ising::spin s);
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
#include <thread>

namespace carlo
{
namespace data
{
  template<typename Model>
  struct attachment
  {
    attachment() = default;
    attachment(const Model&) {}
    virtual void step(const Model&) = 0;
  };
}

namespace model
{
  // We can make this IIsing, which is used all below and is simply an interface with
  // pure virtual operator(), set, flip, scatter, and size
  // Then we can create two different implementations: one that is CPU (like it is now)
  // and another that is for GPU as the memory layout will have to be fundamentally different
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
    using attachment::attachment;

    void step(const ising& model) override
    {
      // calculate the exponent!
      // sigma(\vec{r}) -FFT-> sigma_hat(\vec{k})
      // S(\vec{k}) = average of |sigma_hat(\vec{k})|^2
      // then do ifft(S(\vec{k})) to get g(\vec{r}) then find 
      // then plot g(\vec{r}, t)

      // Construct a representation of the model
      std::vector<fftw_complex> model_rep(model.size() * model.size());
      for (int y = 0; y < model.size(); y++)
        for (int x = 0; x < model.size(); x++)
        {
          model_rep[y * model.size() + x][0] = static_cast<int>(model(x, y));
          model_rep[y * model.size() + x][1] = 0;
        }
        
      // Take the fourier transform
      std::vector<fftw_complex> out(model_rep.size());
      {
        mutex.lock();
        const auto plan = fftw_plan_dft_2d(
          model.size(), 
          model.size(), 
          model_rep.data(), 
          out.data(), 
          FFTW_FORWARD, 
          0
        );
        mutex.unlock();

        fftw_execute(plan);

        // Take the absolute value squared of all of the entries 
        for (int y = 0; y < model.size(); y++)
          for (int x = 0; x < model.size(); x++)
          {
            const auto& real = out[y * model.size() + x][0];
            const auto& imag = out[y * model.size() + x][1];
            model_rep[y * model.size() + x][0] = real * real + imag * imag;
            model_rep[y * model.size() + x][1] = 0;
          }
      }

      std::vector<double> values(model.size() * model.size());
      {
        // Compute the inverse fourier transform
        mutex.lock();
        const auto plan = fftw_plan_dft_2d(
          model.size(), 
          model.size(), 
          model_rep.data(), 
          out.data(), 
          FFTW_BACKWARD, 
          0
        );
        mutex.unlock();
        
        fftw_execute(plan);
        
        // Keep only the real values
        for (int y = 0; y < model.size(); y++)
          for (int x = 0; x < model.size(); x++)
          {
            const auto& real = out[y * model.size() + x][0];
            const auto& imag = out[y * model.size() + x][1];
            values[y * model.size() + x] = real / (model.size() * model.size());
          }
      }

      std::vector<float> rs, vs;

      int dr = 1;
      int final_r = model.size() / 2;

      for (int r = 0; r < final_r; r += dr)
      {
        rs.push_back(r);
        std::vector<float> radial;

        const int bounds = model.size() / 2;
        for (int y = -bounds; y < bounds; y++)
          for (int x = -bounds; x < bounds; x++)
          {
            const auto distance = x*x + y*y;
            if (distance >= r*r && distance <= (r + dr)*(r + dr))
              radial.push_back(values[(y % model.size()) * model.size() + (x % model.size())]);
          }
        
        vs.push_back(std::accumulate(radial.begin(), radial.end(), 0.f) / radial.size());
      }
      
      // Find the first value that drops below zero
      float length = -1.f;
      for (int i = 1; i < vs.size(); i++)
        if (vs[i] < 0.f && vs[i - 1] >= 0.f)
        {
          length = (rs[i - 1] + rs[i]) / 2.f;
          break;
        }

      if (length >= 0.f)
      {
        steps.push_back(current_step);
        lengths.push_back(length);
      }

      current_step++;
    }

    void write(const std::string& filename)
    {
      namespace plt = matplotlibcpp;
      plt::plot(steps, lengths);
      plt::title("Characteristic Length");
      plt::xlabel("Iteration");
      plt::ylabel("l");
      plt::save(filename);
      plt::close();
    }

    std::vector<float> steps, lengths;

  private:
    std::size_t current_step = 0;
    
    inline static std::mutex mutex;
  };
  
  struct total_mag : data::attachment<ising>
  {
    using attachment::attachment;

    void step(const ising& model) override
    {
      int total_mag = 0;
      for (int y = 0; y < model.size(); y++)
        for (int x = 0; x < model.size(); x++)
          total_mag += static_cast<int>(model(x, y));
        
      steps.push_back(current_step);
      total_mag_steps.push_back((float)total_mag / (float)(model.size() * model.size()));

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
      plt::close();
    }

    std::vector<std::size_t> steps;
    std::vector<float> total_mag_steps;
  
  private:
    std::size_t current_step = 0;
  };
}

// The methods own the above models
namespace method
{
  template<typename Model>
  struct base
  {
    using BaseModel = Model;

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
      assert(has_attachment<Att>());
      return std::reinterpret_pointer_cast<Att>(attachments[std::type_index(typeid(Att))]);
    }

    template<typename Att>
    bool has_attachment() const
    {
      return attachments.count(std::type_index(typeid(Att)));
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

  template<typename M>
  struct multirun
  {
    template<typename... Args>
    multirun(std::size_t runs, Args&&... args)
    {
      methods.reserve(runs);
      for (std::size_t i = 0; i < runs; i++)
        methods.emplace_back(std::forward<Args>(args)...);
    }

    template<typename Att, typename... Args>
    void add_attachment_to_each(Args&&... args)
    {
      for (auto& m : methods) m.add_attachment<Att>(std::forward<Args>(args)...);
    }

    void step()
    {
      for (auto& m : methods) m.step();
    }

  private:
    std::vector<M> methods;
  };

  template<typename Att>
  static Att
  reduce_attachment(const std::vector<std::shared_ptr<Att>>& attachments);
  
  template<>
  model::total_mag
  reduce_attachment(const std::vector<std::shared_ptr<model::total_mag>>& attachments)
  {
    model::total_mag att;

    // All of the instances' attachments must have the same dimension and step values
    std::optional<std::size_t> iterations;
    const auto check = [&attachments, &iterations, &att]() -> bool
    {
      for (const auto& base_att : attachments)
      {
        if (!iterations)
        {
          iterations = base_att->steps.size();

          att.steps.reserve(*iterations);
          for (const auto& step : base_att->steps)
            att.steps.push_back(step);
        }
        else if (base_att->steps.size() != iterations)
          return false;
        
        for (std::size_t i = 0; i < *iterations; i++)
          if (att.steps[i] != base_att->steps.at(i))
            return false;
      }
      return true;
    }();
    assert(check);

    att.total_mag_steps.resize(*iterations, 0.f);
    for (std::size_t i = 0; i < *iterations; i++)
    {
      auto& step = att.total_mag_steps[i];
      for (const auto& a : attachments)
        step += a->total_mag_steps.at(i);
      step /= attachments.size();
    }

    return att;
  }

  
  template<>
  model::characteristic_length
  reduce_attachment(const std::vector<std::shared_ptr<model::characteristic_length>>& attachments)
  {
    model::characteristic_length att;

    struct TimePoint
    {
      float value;
    };

    std::vector<TimePoint> time_points;
    std::map<std::size_t, std::unordered_map<std::size_t, float>> data_points;
    const auto get_time_index = [&time_points](float time_value)
    {
      for (std::size_t i = 0; i < time_points.size(); i++)
        if (std::abs(time_points[i].value - time_value) <= 0.0001f)
          return i;

      TimePoint tp{ .value = time_value };
      const auto it = time_points.insert( 
          std::upper_bound( time_points.begin(), time_points.end(), tp, 
            [](const auto& a, const auto& b){ return a.value < b.value; } ),
          tp 
      );
      return static_cast<std::size_t>( std::distance(time_points.begin(), it) );
    };

    for (std::size_t att_i = 0; att_i < attachments.size(); att_i++)
      for (std::size_t i = 0; i < attachments[att_i]->steps.size(); i++)
        data_points[get_time_index(attachments[att_i]->steps[i])][att_i] = attachments[att_i]->lengths[i];

    std::vector<std::size_t> to_remove;
    for (const auto& [ time_index, map ] : data_points)
      if (map.size() != attachments.size())
        to_remove.push_back(time_index);
    
    for (const auto& r : to_remove)
      data_points.erase(r);
    
    att.steps.reserve(data_points.size());
    att.lengths.reserve(data_points.size());
    for (const auto& [ time_index, map ] : data_points)
    {
      att.steps.push_back(time_points[time_index].value);
      auto& length = att.lengths.emplace_back(0);
      for (const auto& [ index, value ] : map)
        length += value;
      length /= map.size();
    }

    // we want to find the first and last instances of an unordered_map that contains all of the values
    // then delete all the instances before the first and after the last

    // data_points is ragged, to make it not, we need to go through the unordered_map
    // and determine which indices are missing, then we take the previous value for this index
    // and the next value for this index 

    return att;
  }

  template<typename Att, typename Method>
  static Att 
  reduce(const std::vector<std::shared_ptr<Method>>& instances)
  {
    static_assert(std::is_base_of_v<base<typename Method::BaseModel>, Method>);
    
    using Model    = typename Method::BaseModel;
    using BaseType = base<Model>;

    // Convert to vector of pointers to base type
    std::vector<std::shared_ptr<Att>> casted;
    casted.reserve(instances.size());
    for (const auto& method : instances)
      if (method->template has_attachment<Att>())
        casted.emplace_back(method->template get_attachment<Att>());

    return reduce_attachment<Att>(casted);
  }
}

}

std::string operator*(carlo::model::ising::spin s);
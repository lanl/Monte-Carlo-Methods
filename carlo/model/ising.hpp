#pragma once

#include <matplotlibcpp.h>

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
    virtual void step(const Model&) = 0;
  };
}

namespace model
{
  struct ising
  {
    enum spin : int { down = -1, up = 1 };

    ising(std::size_t dimension);

    spin operator[](std::size_t x, std::size_t y) const;
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
  
  struct total_mag : data::attachment<ising>
  {
    void step(const ising& model)
    {
      int total_mag = 0;
      for (int y = 0; y < model.size(); y++)
        for (int x = 0; x < model.size(); x++)
          total_mag += static_cast<int>(model[x, y]);
        
      steps.push_back(current_step);
      total_mag_steps.push_back(total_mag);

      current_step++;
    }

    void show()
    {
      namespace plt = matplotlibcpp;

      // Generate matplotlib image
      plt::plot(steps, total_mag_steps);
       plt::save("./total_mag.png");
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
  template<typename Model, typename T>
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
      attachments[std::type_index(typeid(Att))] = std::make_shared<Att>(std::forward<Args>(args)...);
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

  struct metropolis : base<model::ising, metropolis>
  {
    metropolis(double J, double T, std::size_t dimension);

    const model::ising& model() const { return _model; }

  protected:
    void method_step() override;

  private:
    double _E;
    const double _J, beta;

    std::mt19937 gen;
    std::uniform_real_distribution<double> rfd;
    std::uniform_int_distribution<int> rd;
  };
}

}

std::string operator*(carlo::model::ising::spin s);
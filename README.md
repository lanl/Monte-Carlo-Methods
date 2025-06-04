# Overview
The purpose of this framework is to separate *models* from *methods* when solving basic iterative systems. This way, we can write statistical measures on the *model*, while replacing the *method* and keeping most of the functionality invariant.

## Models and Methods
The model represents the underlying data. In the case of `model::ising`, the underlying data is a (-1, 1) value arranged in a 2x2 grid. 

The method is a layer above this. It operates on the underlying model to advance its state. Because attachments only require models, this means we can use the same attachments to extract important quantities and then exchange the method out to compare the difference.

## Attachments
The measures come in the form of `data::attachment`, which takes a model as a template parameter. It can directly access the current state of the model (in the case of the Ising model, the 2x2 matrix of [-1, 1] spins) from which it can derive its statistical measure. In the case of `data::total_mag`, it simply sums the value at all the sites. For `data::characteristic_length`, it does a more complicated procedure involving Fourier transforms to determine the current characteristic length.

# Copyright
© 2025. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

(Copyright request O#: O4859).

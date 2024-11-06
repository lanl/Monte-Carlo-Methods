# Overview
The purpose of this framework is to separate *models* from *methods* when solving basic iterative systems. This way, we can write statistical measures on the *model*, while replacing the *method* and keeping most of the functionality invariant.

## Models and Methods
The model represents the underlying data. In the case of `model::ising`, the underlying data is a (-1, 1) value arranged in a 2x2 grid. 

The method is a layer above this. It operates on the underlying model to advance its state. Because attachments only require models, this means we can use the same attachments to extract important quantities and then exchange the method out to compare the difference.

## Attachments
The measures come in the form of `data::attachment`, which takes a model as a template parameter. It can directly access the current state of the model (in the case of the Ising model, the 2x2 matrix of [-1, 1] spins) from which it can derive its statistical measure. In the case of `data::total_mag`, it simply sums the value at all the sites. For `data::characteristic_length`, it does a more complicated procedure involving Fourier transforms to determine the current characteristic length.
Based on pybind11, I created a Python module in C++. Possible minor changes and enhancements can be expected.

This module contains:

- A function returning the n-th element of the Fibonacci sequence.

- Given the following ordinary differential equation (ODE) of order k:
 - x^{k}(t) = f(t, x(t),..., x^{k-1}(t))
   Subject to the initial conditions:
   - x(t_{0}) = x_{0}
   - x'(t_{0}) = x'_{0}
   - ...
   - x^{k-1}(t_{0}) = x^{k-1}_{0}

- **Vector class:** This class is a 1-dimensional  C++ container for template type data. Based on Vector class there is created a class Vector_d which is Python embedding of Vector<double>. Below is a list of functionalities available in the C++ class:
  - **Template:**
    - Default constructor
    - Constructor depending on the size of the Vector
    - Constructor based on the size of the Vector and a given coordinate
    - Constructor based on the `std::vector` data container
    - Constructor based on the object Vector
    - Destructor
  - **Operators:**
    - Assignation `=`
    - Comparison `==`
    - Addition `+`
    - Coordinate-wise multiplication `*`
    - Multiplication by scalar `*`
    - Inner product operator `^`
    - Input stream `>>`
    - Output stream `<<`
  - Insertion and appendation
  - Removal of coordinate by index or by value
  - Merge of two Vectors
  - Type conversion
  - Action of an external function in a given rolling window on Vector's coordinate
  - Sorting
  - Action of an external function on rows of several Vectors
  - Getting length and coordinate value by index
  - Iteration functionality (for Python usage only)
  - Modification of chosen coordinate
  - Returning of indices with provided values
  - Verification if a Vector has a given value stored in coords

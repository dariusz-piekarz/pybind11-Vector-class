Based on pybind11, I created a Python module in C++.

This module contains:

- A function returning the n-th element of the Fibonacci sequence.
- **Vector_d class:** This class is a 1-dimensional container for double data. Below is a list of functionalities available in the C++ class:
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

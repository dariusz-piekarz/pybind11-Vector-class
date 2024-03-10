Based on pybind11, I created a Python module in C++. 


This module contains:

1) function returning n-th element of the Fibonacci sequence.

2) class Vector_d. This class is in fact 1 dimensional double-data container.
  Below list of functionalities available in the cpp class:
    a) Template: 
		- default constructor
   		- constructor depending on the size of the Vector
        - constructor based on the size of the Vector and coordinate given
        - constructor based on std::vector<T> data container
        - constructor based on the object Vector<T>
        -destructor
    b) operators:
        - assignation =
        - comparison ==
        - addition +
        - coordinatewise multplication *
        - multiplication by scalar *
        - inner product operator ^
        - istream >>
        - ostream <<
    c) insertation and appendation
    d) removal of coordinate by index or by value
    d) merge of two Vectors
    e) type conversion
    f) action of an external function in a given rolling windows on Vector's coordinate
    g) sorting
    h) action of an external function on rows of several Vectors 
    i) getting length and coordinate value by index
    j) interation functionality (for python usage only)
    k) modification of choosen coordinate
    l) returning of indices with provided values
    m) verfication if a Vector has a given value stored in coords


    

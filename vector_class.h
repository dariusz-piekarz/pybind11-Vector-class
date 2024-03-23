#ifndef VECTOR_CLASS_H
#define VECTOR_CLASS_H

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <stdexcept>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

template<typename T>
class Vector {
protected:
    size_t n;
    T* coord;
    size_t current_iter;

public:
    Vector();
    ~Vector();
    explicit Vector(size_t );
    Vector(size_t , T* );
    explicit Vector(const std::vector<T>&);
    Vector(const Vector<T>& vect);

    Vector<T>& operator=(const Vector<T>& );
    Vector<T> operator+(const Vector<T>& ) const;
    Vector<T> operator*(const Vector<T>& ) const;
    template <typename U>
    Vector<T> operator*(U ) const;
    T operator^(const Vector<T>& ) const;
    bool operator==(const Vector<T>& ) const;

    template <typename U>
    friend std::istream& operator>>(std::istream&, Vector<U>&);

    template <typename U>
    friend std::ostream& operator<<(std::ostream&, Vector<U>&);

    template <typename U>
    friend size_t len(const Vector<U>& );

    void modify_coord(size_t, T);
    void insert(T, int);
    void append(T);
    bool is_in(T);
    void remove_item(T);
    void del_item_by_index(size_t j);
    size_t* get_pos(T);
    std::vector<size_t> get_position(T);
    void delete_all(T);
    std::string vstr(Vector<T>& );
    Vector<T> merge(Vector<T>& );
    Vector<T> sort(bool);
    Vector<double> vrolling(size_t, T(*)(T*, size_t));
    Vector<double> rolling(size_t, T(*)(std::vector<T>));
    Vector<double> pyrolling(size_t, py::function);

    Vector<double> roll(size_t, std::string);

    template<typename U>
    Vector<U> convert();

    template<typename U, typename V>
    friend Vector<U> act(std::vector<Vector<V>> , U(*)(T*, size_t));

    template<typename U, typename V>
    friend Vector<U> actv(std::vector<Vector<V>>, U(*)(std::vector<V>));

    template<typename U, typename V>
    friend Vector<U> actpy(py::list , py::function);

    template<typename U, typename V>
    friend Vector<U> actpy2(std::vector<Vector<V>>, py::function );
    
    template <typename U>
    friend double (*vchoice(std::string))(std::vector<U>);
    
    template <typename U>
    friend U(*choice(std::string))(U*, size_t);
    
    py::iterator iter();
    T next();
    T get_coordinate(size_t );
};

std::tuple<std::vector<double>, std::vector<double>> general_ode(double (*)(std::vector<double>), double, double, std::vector<double>, int, std::string);
std::tuple<std::vector<double>, std::vector<double>> general_ode_py(py::function, double, double, std::vector<double>, int, std::string);
#endif // VECTOR_CLASS_H
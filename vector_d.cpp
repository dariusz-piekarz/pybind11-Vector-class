#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include <unordered_map>
#include <functional>
#include <tuple>
#include <thread>
#include <mutex>


#include <exprtk.hpp> //https://www.partow.net/programming/exprtk/index.html

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "vector_class.h"

using namespace std;
namespace py = pybind11;


/**
 * @brief Calculates the nth Fibonacci number
 *
 * @param n the index of the Fibonacci number to calculate
 * @return long long the nth Fibonacci number
 */
long long fib(int n)
{
    if (n == 0 || n == 1)
        return 1;
    else if (n > 1)
        return fib(n - 1) + fib(n - 2);
    else
        throw std::runtime_error("n cannot be negative!");
}


template<typename T>
py::list vector_to_list(const std::vector<T>& vec) {
    py::list lst;
    for (const auto& elem : vec) {
        lst.append(elem);
    }
    return lst;
}


template <typename T>
size_t Partition(T* arr, size_t start, size_t end) {

    size_t pivot = end;
    size_t j = start;
    for (size_t i = start; i < end; ++i) {
        if (arr[i] < arr[pivot]) {
            swap(arr[i], arr[j]);
            ++j;
        }
    }
    swap(arr[j], arr[pivot]);
    return j;

}

template <typename T>
void Quicksort(T arr[], size_t start, size_t end) {

    if (start < end) {
        size_t p = Partition(arr, start, end);
        Quicksort(arr, start, p - 1);
        Quicksort(arr, p + 1, end);
    }

}


template <typename T>
Vector<T>::Vector() : n(-1), current_iter(0), coord(nullptr) {}

template <typename T>
Vector<T>::~Vector()
{
    delete[] coord;
}

template <typename T>
Vector<T>::Vector(size_t m) : n(m), current_iter(0) {
    coord = new T[n];
}

template <typename T>
Vector<T>::Vector(size_t m, T* other) : n(m), current_iter(0) {
    coord = new T[n];
    for (size_t i = 0; i < n; i++)
        coord[i] = other[i];
}

template <typename T>
Vector<T>::Vector(const std::vector<T>& vect) : current_iter(0)
{
    n = vect.size();
    coord = new T[n];
    for (size_t i = 0; i < n; i++)
        coord[i] = vect[i];
}

template <typename T>
Vector<T>::Vector(const Vector<T>& vect) : n(vect.n) {
    current_iter = 0;
    coord = new T[n];
    for (size_t i = 0; i < n; i++)
        coord[i] = vect.coord[i];
}

template <typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& vect) {
    delete[] coord;
    n = vect.n;
    coord = new T[n];
    for (size_t i = 0; i < n; i++)
        coord[i] = vect.coord[i];
    return *this;
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& vect) const {
    if (n == vect.n) {
        Vector<T> result(n);
        for (size_t i = 0; i < n; i++)
            result.coord[i] = coord[i] + vect.coord[i];
        return result;
    }
    else
        throw std::runtime_error("Cannot add two vectors of different size!");

}

template <typename T>
Vector<T> Vector<T>::operator*(const Vector<T>& other) const
{
    if (n == other.n) {
        Vector<T> result(n);
        for (size_t i = 0; i < n; i++)
            result.coord[i] = coord[i] * other.coord[i];
        return result;
    }
    else
        throw std::runtime_error("Cannot add two vectors of different size!");
}

template <typename T>
template <typename U>
Vector<T> Vector<T>::operator*(U scalar) const {
    Vector<T> result(n);
    for (size_t i = 0; i < n; i++)
        result.coord[i] = coord[i] * scalar;
    return result;
}

template <typename T>
T Vector<T>::operator^(const Vector<T>& vect) const {
    if (n == vect.n) {
        T inner = static_cast<T>(0);
        for (size_t i = 0; i < n; i++)
            inner += coord[i] * vect.coord[i];
        return inner;
    }
    else
        throw std::runtime_error("Two vectors of different size!");
}

template<typename T>
bool Vector<T>::operator==(const Vector<T>& v) const
{
    if (n != v.n) {
        return false;
    }
    for (size_t i = 0; i < n; i++) {
        if (coord[i] != v.coord[i]) {
            return false;
        }
    }
    return true;
}

template <typename U>
std::istream& operator>>(std::istream& istr, Vector<U>& vect) {
    if (vect.n == -1) {
        std::cout << "Enter the vector's dimension:" << std::endl;
        istr >> vect.n;
    }
    vect.coord = new U[vect.n];
    for (size_t i = 0; i < vect.n; i++)
        istr >> vect.coord[i];
    return istr;
}

template <typename U>
std::ostream& operator<<(std::ostream& ostr, Vector<U>& vect) {
    for (size_t i = 0; i < vect.n; i++)
        if (i == 0)
            ostr << "[" << vect.coord[i] << ", ";
        else if (i == vect.n - 1)
            ostr << vect.coord[i] << "]";
        else
            ostr << vect.coord[i] << ", ";
    return ostr;
}

template <typename T>
Vector<T> inputv(std::string str)
{
    Vector <T> output;
    if (!str.empty())
        cout << str << endl;
    cin >> output;
    return output;
}

template<typename T>
py::iterator Vector<T>::iter()
{
    return py::make_iterator(coord, coord + n);
}

template<typename T>
T Vector<T>::next()
{
    if (current_iter >= n) 
        throw py::stop_iteration();
    return coord[current_iter++];
}

template <typename T>
template <typename U>
Vector<U> Vector<T>::convert()
{
    U* temp;
    temp = new U[n];
    for (size_t i = 0; i < n; i++)
        temp[i] = static_cast<U>(coord[i]);
    Vector <U> output_vector(n, temp);
    delete []temp;
    return output_vector;
}

template <typename T>
T Vector<T>::get_coordinate(size_t i) {
    if (i < n)
        return coord[i];
    else
        throw std::runtime_error("Index out of bounds!");
}

template <typename T>
size_t len(const Vector<T>& vect)
{
    return vect.n;
}

template<typename T>
void Vector<T>::modify_coord(size_t i, T value_to_set)
{
    if (i < 0 || i>n - 1)
        throw std::runtime_error("i parameter must be in vector's dimension range!");
    coord[i] = value_to_set;
}

template<typename T>
void Vector<T>::insert(T element, int j)
{
    T* temp = new T[n + 1];
    for (size_t i = 0; i < n + 1; i++)
        if (i == j)
            temp[i] = element;
        else if (i < j)
            temp[i] = coord[i];
        else
            temp[i] = coord[i - 1];
    delete[]coord;
    coord = temp;
    n += 1;
}

template<typename T>
void Vector<T>::append(T element)
{
    T* temp = new T[n + 1];
    for (size_t i = 0; i < n; i++)
        temp[i] = coord[i];
    temp[n] = element;
    delete[]coord;
    coord = temp;
    n += 1;
}

template<typename T>
bool Vector<T>::is_in(T element)
{
    for (size_t i = 0; i < n; i++)
        if (coord[i] == element)
            return true;
    return false;
}

template<typename T>
void Vector<T>::remove_item(T element)
{
    size_t element_coordinate = NULL;
    for (size_t i = 0; i < n; i++)
        if (coord[i] == element)
        {
            element_coordinate = i;
            break;
        }
    if (element_coordinate != NULL)
        (*this).del_item_by_index(element_coordinate);
}

template<typename T>
void Vector<T>::del_item_by_index(size_t j)
{
    T* temp = new T[n - 1];
    for (size_t i = 0; i < n ; i++)
        if (i < j)
            temp[i] = coord[i];
        else
            temp[i] = coord[i + 1];
    delete[]coord;
    coord = temp;
    n -= 1;
}

template<typename T>
size_t* Vector<T>::get_pos(T element)
{
    size_t couter = 0;
    for (size_t i = 0; i < n; i++)
        if (coord[i] == element)
            couter += 1;
    size_t* indices;
    if (couter > 0)
        indices = new size_t[couter];
    else
        return NULL;
    couter = 0;
    for (size_t i = 0; i < n; i++)
        if (coord[i] == element)
        {
            indices[couter] = i;
            couter += 1;
        }
    return indices;
}

template<typename T>
std::vector<size_t> Vector<T>::get_position(T element)
{
    vector<size_t> indices;
    for (size_t i = 0; i < n; i++)
        if (coord[i] == element)
            indices.push_back(i);
    return indices;
}

template<typename T>
void Vector<T>::delete_all(T element)
{
    vector<size_t> ind_to_del = get_position(element);
    size_t size = ind_to_del.size();
 
    if (size > 0)
    {
        size_t counter = 0;
        T* temp = new T[n - size];
        for (size_t i = 0; i < n; i++)
        {
            while (coord[i + counter] == element)
                counter += 1;
            temp[i] = coord[i + counter];
        }
        delete[]coord;
        coord = temp;
        n -= size;
    }
}

template <typename T>
std::string vstr(Vector<T>& v) {
    std::string output = "[";
    for (size_t i = 0; i < len(v) - 1; i++) {
        output += std::to_string(v.get_coordinate(i));
        output += ", ";
    }
    if (len(v) > 0) {
        output += std::to_string(v.get_coordinate(len(v) - 1));
    }
    output += "]";
    return output;
}

template <typename T>
Vector<T> Vector<T>::merge(Vector <T>& vect)
{
    
    size_t m = n + vect.n;
    T* temp = new T[m];
    for (size_t i = 0; i < m; i++)
        if (i < n)
            temp[i] = coord[i];
        else
            temp[i] = vect.coord[i - n];
    Vector<double> output(m, temp);
    return output;
}

template <typename T>
Vector<T> Vector<T>::sort(bool ascending)
{
    Quicksort(coord, 0, n - 1);
    if (!ascending)
        for (size_t i = 0; i < n / 2; ++i)
            std::swap(coord[i], coord[n - i - 1]);
    return *this;
}

template <typename T>
Vector<double> Vector<T>::vrolling(size_t window, T(*fun)(T*, size_t))
{
    double* out;
    out = new double[n];
    if (window > n)
        throw std::length_error("window should be shorter than length of the vector!");
    for (size_t i = 0; i < n; i++)
    {

        if (i < window - 1)
            out[i] = std::numeric_limits<double>::quiet_NaN();
        else
        {
            T* temp;
            temp = new T[window];
            for (size_t j = 0; j < window; j++)
                temp[j] = coord[i - window + 1 + j];
            out[i] = fun(temp, window).cast<double>();
            delete[] temp;
        }
    }
    Vector<double> output(n, out);
    delete[]out;
    return output;
}

template <typename T>
Vector<double> Vector<T>::rolling(size_t window, T(*fun)(std::vector<T>))
{
    double* out;
    out = new double[n];
    if (window > n)
        throw std::length_error("window should be shorter than length of the vector!");
    for (size_t i = 0; i < n; i++)
    {

        if (i < window - 1)
            out[i] = std::numeric_limits<double>::quiet_NaN();
        else
        {
            std::vector<T> temp;
            temp.reserve(window);
            for (size_t j = i - window + 1; j < i + 1; j++)
                temp.push_back(coord[j]);
            out[i] = (double) fun(temp);
        }
    }
    Vector<double> output(n, out);
    delete[]out;
    return output;
}

template <typename T>
Vector<double> Vector<T>::pyrolling(size_t window, py::function fun)
{
    double* out;
    out = new double[n];
    if (window > n)
        throw std::length_error("window should be shorter than length of the vector!");
    for (size_t i = 0; i < n; i++)
    {

        if (i < window - 1)
            out[i] = std::numeric_limits<double>::quiet_NaN();
        else
        {
            std::vector<T> temp;
            temp.reserve(window);
            for (size_t j = i - window + 1; j < i + 1; j++)
                temp.push_back(coord[j]);
            out[i] = fun(vector_to_list(temp)).cast<double>();
        }
    }
    Vector<double> output(n, out);
    delete[]out;
    return output;
        
}

template<typename T>
T add(std::vector<T> v)
{
    T sum = (T)0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += v[i];
    }
    return sum;
}

template<typename T>
T mult(std::vector<T> v)
{
    T pr = (T)1;
    for (size_t i = 0; i < v.size(); ++i) {
        pr *= v[i];
    }
    return pr;
}

template<typename T>
T minimal(std::vector<T> v)
{
    T pr = v[0];
    for (size_t i = 1; i < v.size(); ++i)
        if (v[i] < pr)
            pr = v[i];
    return pr;
}

template<typename T>
T maximal(std::vector<T> v)
{
    T pr = v[0];
    for (size_t i = 1; i < v.size(); ++i)
        if (v[i] > pr)
            pr = v[i];
    return pr;
}

template<typename T>
double mean(std::vector<T> v)
{
    T pr = add<T>(v);
    return pr / v.size();
}

template<typename T>
double stdev(std::vector<T> v)
{
    T pr=0;
    T avg = mean<T>(v);
    for (size_t i = 0; i < v.size(); ++i)
        pr += pow(v[i] - avg, 2);
    return sqrt(pr);
}

template<typename T>
T add1(T* v, size_t size)
{
    T sum = (T)0;
    for (size_t i = 0; i < size; ++i)
        sum += v[i];
    return sum;
}

template<typename T>
T prod1(T* v, size_t size)
{
    T pr = (T)1;
    for (size_t i = 0; i < size; ++i)
        pr *= v[i];
    return pr;
}

template<typename T>
double p_norm(std::vector<T> v, int p)
{
    T sum = (T)0;
    for (size_t i = 0; i < v.size(); ++i) {
        sum += pow(abs(v[i]), p);
    }
    return pow(sum, (double)1/p);
}

template<typename T>
double (*vchoice(std::string name))(std::vector<T>)
{
    static std::unordered_map<std::string, double (*)(std::vector<T>)> functions;
    functions["add"] = reinterpret_cast<double (*)(std::vector<T>)>(add<T>);
    functions["mul"] = reinterpret_cast<double (*)(std::vector<T>)>(mult<T>);
    functions["mean"] = mean<T>;
    functions["std"] = stdev<T>;
    functions["min"] = reinterpret_cast<double (*)(std::vector<T>)>(minimal<T>);
    functions["max"] = reinterpret_cast<double (*)(std::vector<T>)>(maximal<T>);
    //functions["p_norm"] = p_norm<T>;
    return functions[name];
}

template<typename T>
T (*choice(std::string name))(T*, size_t)
{
    static std::unordered_map<std::string, T (*)(T*, size_t)> functions;
    functions["add"] = add1<T>;
    functions["mul"] = prod1<T>;
    return functions[name];
}

template <typename T>
Vector<double> Vector<T>::roll(size_t window, std::string name)
{
    double* out;
    out = new double[n];
    if (window > n)
        throw std::length_error("window should be shorter than length of the vector!");
    for (size_t i = 0; i < n; i++)
    {
        if (i < window - 1)
            out[i] = std::numeric_limits<double>::quiet_NaN();
        else
        {
            std::vector<T> temp;
            temp.reserve(window );
            for (size_t j = i - window + 1; j < i + 1; j++)
                temp.push_back(coord[j]);
            out[i] = (*vchoice<T>(name))(temp);
        }
    }
    Vector<double> output(n, out);
    delete[]out;
    return output;

}

template<typename U, typename T>
Vector<U> actv(std::vector<Vector<T>> collection, U(*fun)(std::vector<T>))
{
    size_t dim = len(collection[0]);
    for (size_t i = 0; i < collection.size(); i++)
        if (len(collection[i]) != dim)
            throw std::runtime_error("Inconsisten dimensions!");

    std::vector<U> out;
    out.reserve(dim);
    static vector<T> temp;
    temp.reserve(collection.size());

    for (size_t i = 0; i < dim; i++)
    {

        for (size_t j = 0; j < collection.size(); j++)
            temp.push_back(collection[j].get_coordinate(i));
        out.push_back(fun(temp));
        temp.clear();
    }
    Vector<U> output(out);
    return output;
}

template<typename U, typename T>
Vector<U> act(std::vector<Vector<T>> collection, U(*fun)(T*, size_t))
{
    size_t dim = len(collection[0]);
    for (size_t i = 0; i < collection.size(); i++)
        if (len(collection[i]) != dim)
            throw std::runtime_error("Inconsisten dimensions!");

    U* out;
    out = new U[dim];
    T* temp;
    temp = new T[dim];

    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j < collection.size(); j++)
            temp[j] = collection[j].get_coordinate(i);
        out[i] = fun(temp, collection.size());
        delete[]temp;
    }
    Vector<U> output(dim, out);
    delete[] out;
    return output;
}

template<typename U, typename T>
Vector<U> actpy(py::list col, py::function fun)
{
    std::vector<Vector<T>> collection;
    collection.reserve(len(col));
    for (size_t i = 0; i < len(col); i++)
        collection.push_back(col[i].cast<Vector<T>>());
    size_t dim = len(collection[0]);
    for (size_t i = 0; i < collection.size(); i++)
        if (len(collection[i]) != dim)
            throw std::runtime_error("Inconsisten dimensions!");

    std::vector<U> out;
    out.reserve(dim);
    static vector<T> temp;
    temp.reserve(collection.size());

    for (size_t i = 0; i < dim; i++)
    {

        for (size_t j = 0; j < collection.size(); j++)
            temp.push_back(collection[j].get_coordinate(i));
        out.push_back(fun(vector_to_list(temp)).cast<U>());
        temp.clear();
    }
    Vector<U> output(out);
    return output;
}

template<typename U, typename T>
Vector<U> actpy2(std::vector<Vector<T>> collection, py::function fun)
{
    size_t dim = len(collection[0]);
    for (size_t i = 0; i < collection.size(); i++)
        if (len(collection[i]) != dim)
            throw std::runtime_error("Inconsisten dimensions!");

    std::vector<U> out;
    out.reserve(dim);
    static vector<T> temp;
    temp.reserve(collection.size());

    for (size_t i = 0; i < dim; i++)
    {
        for (size_t j = 0; j < collection.size(); j++)
            temp.push_back(collection[j].get_coordinate(i));
        out.push_back(fun(vector_to_list(temp)).cast<U>());
        temp.clear();
    }
    Vector<U> output(out);
    return output;
}


tuple<vector<double>, vector<double>> general_ode(
    double (*fun)(vector<double>),
    double t0,
    double tmax,
    vector<double> x0,
    int no_nodes,
    string plot_expansion) {
    if (plot_expansion == "right" || plot_expansion == "left")
    {
        vector<double> t(no_nodes);
        for (size_t i = 0; i < no_nodes; i++)
            t[i] = t0 + i * (tmax - t0) / no_nodes;

        double dt = t[1] - t[0];
        size_t k = x0.size();
        vector<vector<double>> x(k + 1, vector<double>(no_nodes));

        for (size_t i = 0; i < k; i++)
            x[i][0] = x0[i];

        for (size_t i = k; i > 0; i--)
            for (size_t j = 0; j < count(x[i].begin(), x[i].end(), double()); j++)
                x[i - 1][j + 1] = dt * x[i][j] + x[i - 1][j];

        for (size_t j = 0; j < no_nodes - k; j++)
            for (size_t i = k; i > -1; i--)
                if (i == k)
                {
                    vector<double> arg(k + 1);
                    for (size_t l = 0; l < k + 1; l++)
                        arg[l] = (l == 0) ? t[j] : x[l - 1][j];
                    x[i][j] = fun(arg);
                }
                else
                    x[i][k - i + j] = x[i + 1][k - i + j - 1] * dt + x[i][k - i + j - 1];

        if (plot_expansion == "left")
        {
            reverse(x[0].begin(), x[0].end());
            reverse(t.begin(), t.end());
            return make_tuple(t, x[0]);
        }
        else
            return make_tuple(t, x[0]);
    }
    else if (plot_expansion == "both")
    {
        vector<double> t_left, t_right;
        vector<double> x_left, x_right;
        if (tmax > t0)
        {
            tie(t_left, x_left) = general_ode(fun, t0, 2 * t0 - 1 * tmax, x0, no_nodes, "left");
            tie(t_right, x_right) = general_ode(fun, t0, tmax, x0, no_nodes, "right");
        }
        else if (tmax < t0)
        {
            tie(t_left, x_left) = general_ode(fun, t0, tmax, x0, no_nodes, "left");
            tie(t_right, x_right) = general_ode(fun, t0, 2 * t0 - 1 * tmax, x0, no_nodes, "right");
        }
        t_left.insert(t_left.end(), t_right.begin() + 1, t_right.end());
        x_left.insert(x_left.end(), x_right.begin() + 1, x_right.end());

        return make_tuple(t_left, x_left);
    }
    else
        return make_tuple(vector<double>(), vector<double>());
}

tuple<vector<double>, vector<double>> general_ode_py(
    py::function fun,
    double t0,
    double tmax,
    vector<double> x0,
    int no_nodes,
    string plot_expansion) {
    if (plot_expansion == "right" || plot_expansion == "left")
    {
        vector<double> t(no_nodes);
        for (size_t i = 0; i < no_nodes; i++)
            t[i] = t0 + i * (tmax - t0) / no_nodes;

        double dt = t[1] - t[0];
        size_t k = x0.size();
        vector<vector<double>> x(k + 1, vector<double>(no_nodes));

        for (size_t i = 0; i < k; i++)
            x[i][0] = x0[i];

        for (size_t i = k; i > 0; i--)
            for (size_t j = 0; j < count(x[i].begin(), x[i].end(), double()); j++)
                x[i - 1][j + 1] = dt * x[i][j] + x[i - 1][j];

        for (size_t j = 0; j < no_nodes - k; j++)
            for (size_t i = k; i > -1; i--)
                if (i == k)
                {
                    py::list arg;
                    for (size_t l = 0; l < k + 1; l++)
                        arg.append((l == 0) ? t[j] : x[l - 1][j]);
                    x[i][j] = fun(*py::tuple(arg)).cast<double>();
                }
                else
                    x[i][k - i + j] = x[i + 1][k - i + j - 1] * dt + x[i][k - i + j - 1];

        if (plot_expansion == "left")
        {
            reverse(x[0].begin(), x[0].end());
            reverse(t.begin(), t.end());
            return make_tuple(t, x[0]);
        }
        else
            return make_tuple(t, x[0]);
    }
    else if (plot_expansion == "both")
    {
        vector<double> t_left, t_right;
        vector<double> x_left, x_right;
        if (tmax > t0)
        {
            tie(t_left, x_left) = general_ode_py(fun, t0, 2 * t0 - 1 * tmax, x0, no_nodes, "left");
            tie(t_right, x_right) = general_ode_py(fun, t0, tmax, x0, no_nodes, "right");
        }
        else if (tmax < t0)
        {
            tie(t_left, x_left) = general_ode_py(fun, t0, tmax, x0, no_nodes, "left");
            tie(t_right, x_right) = general_ode_py(fun, t0, 2 * t0 - 1 * tmax, x0, no_nodes, "right");
        }
        t_left.insert(t_left.end(), t_right.begin() + 1, t_right.end());
        x_left.insert(x_left.end(), x_right.begin() + 1, x_right.end());

        return make_tuple(t_left, x_left);
    }
    else
        return make_tuple(vector<double>(), vector<double>());
}

inline double evaluateExpression(const string& expression, double x, double y, double z)
{
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    typedef exprtk::symbol_table<double> symbol_table_t;

    expression_t expr;
    symbol_table_t symbol_table;
    parser_t parser;

    symbol_table.add_variable("x", x);
    symbol_table.add_variable("y", y);
    symbol_table.add_variable("z", z);

    expr.register_symbol_table(symbol_table);

    if (!parser.compile(expression, expr))
        throw runtime_error("Cannot compile...");

    return expr.value();
}


tuple<vector<double>,
      vector<double>, 
      vector<double>> implicit3D(
        string equation,
        vector<double> x,
        vector<double> y,
        vector<double>z,
        double dt,
        double precision,
        int interval_div)
{
   size_t no_nodes_x = ceil((size_t)((x[1] - x[0]) / dt));
   size_t no_nodes_y = ceil((size_t)((y[1] - y[0]) / dt));
   size_t no_nodes_z = ceil((size_t)((z[1] - z[0]) / dt));

   vector<double> x_interval(no_nodes_x);
   vector<double> y_interval(no_nodes_y);
   vector<double> z_interval(no_nodes_z);

   for (size_t i = 0; i < no_nodes_x; i++)
       x_interval[i] = x[0] + i * dt;
   for (size_t i = 0; i < no_nodes_y; i++)
       y_interval[i] = y[0] + i * dt;
   for (size_t i = 0; i < no_nodes_z; i++)
       z_interval[i] = z[0] + i * dt;

    vector<double> x_solution;
    vector<double> y_solution;
    vector<double> z_solution;

    size_t no_break_x = (size_t) floor(x_interval.size() / interval_div);
    size_t no_break_y = (size_t) floor(y_interval.size() / interval_div);
    size_t no_break_z = (size_t) floor(z_interval.size() / interval_div);
    vector<thread> threads(pow(interval_div, 3));
    mutex mtx;

    for (int i=0; i < interval_div; i++)
    {
        vector<size_t> data_ind_x(2);
        if (i < interval_div - 1)
            data_ind_x = { i * no_break_x, (i + 1) * no_break_x };
        else
            data_ind_x = { i * no_break_x, x_interval.size() };
        for (int j = 0; j < interval_div; j++)
        {
            vector<size_t> data_ind_y(2);
            if (j < interval_div - 1)
                data_ind_y = { j * no_break_y, (j + 1) * no_break_y };
            else
                data_ind_y = { j * no_break_y, y_interval.size() };

            for (int k = 0; k < interval_div; k++)
            {
                vector<size_t> data_ind_z(2);
                if (k < interval_div - 1)
                    data_ind_z = {k * no_break_z, (k + 1) * no_break_z};
                else
                    data_ind_z = {k * no_break_z, z_interval.size()};
                
                int thread_index = i * interval_div * interval_div + j * interval_div + k;
                threads[thread_index] = thread(loop_atomic, equation, ref(x_solution), ref(y_solution), ref(z_solution),
                    ref(x_interval), ref(y_interval), ref(z_interval), data_ind_x, data_ind_y, data_ind_z, precision, ref(mtx));
            }
        }
    }

    for (thread & thr : threads)
        thr.join();

    return make_tuple(x_solution, y_solution, z_solution);
}

void loop_atomic(string equation, 
    vector<double>& x_solution,
    vector<double>& y_solution,
    vector<double>& z_solution,
    vector<double>& x_interval,
    vector<double>& y_interval, 
    vector<double>& z_interval,
    vector<size_t> data_ind_x,
    vector<size_t> data_ind_y,
    vector<size_t> data_ind_z,
    double precision,
    mutex& mtx)
{
    vector<double> partial_solution_x;
    vector<double> partial_solution_y;
    vector<double> partial_solution_z;

    for (size_t i = data_ind_x[0]; i < data_ind_x[1]; i++)
        for (size_t j = data_ind_y[0]; j < data_ind_y[1]; j++)
            for (size_t k = data_ind_z[0]; k < data_ind_z[1]; k++)
                if (abs(evaluateExpression(equation, x_interval[i], y_interval[j], z_interval[k])) <= precision)
                {
                    partial_solution_x.push_back(x_interval[i]);
                    partial_solution_y.push_back(y_interval[j]);
                    partial_solution_z.push_back(z_interval[k]);
                }
            
    if (!partial_solution_x.empty()) 
    {
        mtx.lock();
            x_solution.insert(x_solution.end(), partial_solution_x.begin(), partial_solution_x.end());
            y_solution.insert(y_solution.end(), partial_solution_y.begin(), partial_solution_y.end());
            z_solution.insert(z_solution.end(), partial_solution_z.begin(), partial_solution_z.end());
        mtx.unlock();
    }
}

PYBIND11_MODULE(math_module, m)
{
    m.def("implicit3D", &implicit3D);
    m.def("ode", &general_ode_py);
    m.def("add", &add<double>, py::arg("v"));
    m.def("mul", &mult<double>, py::arg("v"));
    m.def("mean", &mean<double>, py::arg("v"));
    m.def("min", &minimal<double>, py::arg("v"));
    m.def("max", &maximal<double>, py::arg("v"));
    m.def("stdev", &stdev<double>, py::arg("v"));
    m.def("p_norm", &p_norm<double>, py::arg("v"), py::arg("p"));
    m.def("fib", &fib);
    m.def("mult", &prod1<double>, py::arg("v"), py::arg("size"));
    m.def("inp", &inputv<double>, py::arg("str"));
    m.def("act", &actpy<double, double>, py::arg("collection"), py::arg("function"));
    m.def("act2", &actpy2<double, double>, py::arg("collection"), py::arg("function"));

    py::class_<Vector<double>>(m, "Vector_d")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<const std::vector<double>&>(), py::return_value_policy::take_ownership)
        .def(py::init<const Vector<double>&>())
        .def("__iassign__", [](Vector<double>& self, const Vector<double>& other)
            {
                self = other; 
                return self;
            })
        .def("__add__", [](Vector<double>& self, Vector<double>& other)
            {
                Vector <double> output;
                output = self + other;
                return output;
            })
        .def("__sub__", [](Vector<double>& self, Vector<double>& other)
            {
                Vector <double> output;
                output = self + other*(-1);
                return output;
            })
        .def("__mul__", [](Vector<double>& self, double num)
            {
                Vector <double> output;
                output = self * num;
                return output;
            })
        .def("__mul__", [](Vector<double>& self, Vector<double>& other)
            {
                return self*other;
            })
        .def("__setitem__", &Vector<double>::modify_coord, py::arg("i"), py::arg("value_to_set"))
        .def("__getitem__", &Vector<double>::get_coordinate, py::arg("i"))
        .def("__delitem__", &Vector<double>::del_item_by_index, py::arg("j"))
        .def("__str__", [](Vector<double>& vect)
            {
                std::string out;
                out = vstr(vect);
                return out;
            })
        .def("__len__", [](Vector <double>& vect)
            {
                return len(vect);
            })
        .def("__eq__", [](const Vector<double>& v1, const Vector<double>& v2) 
            {
                return v1 == v2;
            })
        .def("__ne__", [](const Vector<double>& v1, const Vector<double>& v2)
            {
                return !(v1 == v2);
            })
        .def("__iter__", &Vector<double>::iter)
        .def("__next__", &Vector<double>::next)
        .def("__contains__", &Vector<double>::is_in)
        .def("append", &Vector<double>::append, py::arg("element"))
        .def("insert", &Vector<double>::insert, py::arg("element"), py::arg("j"))
        .def("del_all",&Vector<double>::delete_all, py::arg("element"))
        .def("merge", &Vector<double>::merge, py::arg("vect"))
        .def("remove", &Vector<double>::remove_item, py::arg("element"))
        .def("sort", &Vector<double>::sort, py::arg("ascending"))
        .def("rolling", &Vector<double>::pyrolling, py::arg("window"), py::arg("function"))
        .def("roll", &Vector<double>::roll, py::arg("window"), py::arg("fun_name"))
        .def("get_pos", &Vector<double>::get_position, py::arg("element"));    
}


int main()
{
    return 0;
}

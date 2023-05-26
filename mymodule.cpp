#include <comp.hpp>
#include <python_comp.hpp>

// using namespace ngcomp;

#include "hdivsym.cpp"
#include "hdivsymfespace.cpp"
   
extern "C" void mymodule(py::object & res) {
  cout << "called mymodule" << endl;
  // import ngsolve such that python base classes are defined    
  auto ngs = py::module::import("ngsolve");    

  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    
    
  ExportHDivSymFESpace(m);
  res = m;    
}    


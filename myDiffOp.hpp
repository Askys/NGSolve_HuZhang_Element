#ifndef FILE_MYDIFFOP_HPP
#define FILE_MYDIFFOP_HPP

#include <fem.hpp>

#include "hdivsym.hpp"

namespace ngfem
{
  /* 
     DiffOps provide the link between function evaluation, and finite elements.
     Different DiffOps evaluate either shape functions, or derivatives.
     Typically, a FiniteElementSpace creates a DiffOp to for function evaluation, 
     and for the canonical derivative, as well as DiffOps for evaluation at the
     boundary. These DiffOps are used when evaluating GridFunctions, and setting 
     up element-matrices from trial- and test-functions.
     DiffOps use static polymorphism, aka Curiously Recurring Template Pattern (CRTP).
   */

  // 
   
  class DiffOpIdHDivSym : public DiffOp<DiffOpIdHDivSym>
  {
  public:
    // some constants for the diffop:
    
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 2*2 };
    enum { DIFFORDER = 0 };
    enum { DIM_STRESS = 2*2 };
    
//     static string Name() { return "id"; }

    static Array<int> GetDimensions() { return Array<int> ({2,2}); }

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      dynamic_cast<const HDivSymFE&> (fel).CalcMappedShape(mip, Trans(mat));
    }

    // can overload more functionality for performance optimization,
    // like evaluation in the whole integration rule
  };
    
  class DiffOpDivHDivSym : public DiffOp<DiffOpDivHDivSym>
  {
  public:
    // some constants for the diffop:    
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 2 };
    enum { DIFFORDER = 1 };
    enum { DIM_STRESS = 2 };
    
    static string Name() { return "div"; }

    //static Array<int> GetDimensions() { return Array<int> ({2,1}); }

    // fill evaluation matrix of dimension 1 times fel.ndof 
    // the input mip is a mapped integration point, which has also
    // access to the integration point on the reference element
    template<typename MIP, typename MAT>
    static void GenerateMatrix(const FiniteElement & fel, const MIP & mip,
                               MAT & mat, LocalHeap & lh)
    {
      dynamic_cast<const HDivSymFE&> (fel).CalcMappedDShape(mip, Trans(mat));
    }

    // can overload more functionality for performance optimization,
    // like evaluation in the whole integration rule
  };
}
#endif // FILE_MYDIFFOP_HPP

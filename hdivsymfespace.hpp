#ifndef FILE_HDIVSYMFESPACE_HPP
#define FILE_HDIVSYMFESPACE_HPP


namespace ngcomp
{
  class HDivSymFESpace : public FESpace
  {
    int order;
    
    Array<int> first_edge_dof;
    Array<int> first_cell_dof;
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    HDivSymFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    string GetClassName () const override
    {
      return "HDivSymFESpace";
    }

    void Update() override;

    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
  };
}    

void ExportHDivSymFESpace(py::module m);
#endif

#ifndef FILE_HDIVSYM_HPP
#define FILE_HDIVSYM_HPP


namespace ngfem
{

  class HDivSymFE : public FiniteElement, public VertexOrientedFE<ET_TRIG>
  {
  public:
    enum { DIM = 2 };
    //enum { DIM_STRESS = (DIM*(DIM+1))/2 };
    enum { DIM_STRESS = (DIM*DIM) };
    
  public:
    using VertexOrientedFE<ET_TRIG>::SetVertexNumbers;
    
    HDivSymFE (int order);
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const
    {
        throw Exception("HDivSymFE: can't evaluate from ip");
    }
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             BareSliceVector<> dshape) const
    {
        throw Exception("HDivSymFE: can't evaluate from ip");
    }
    
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip, 
                          BareSliceMatrix<> shape) const;
                          
    virtual void CalcMappedDShape (const MappedIntegrationPoint<2,2> & mip, 
                   BareSliceMatrix<> dshape) const;
                   
    template <class T>
    void T_CalcShape (const T & x, const T & y, const Mat<2,2,T> & F, BareSliceMatrix<T> shape) const;
  };


}

#endif


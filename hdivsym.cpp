#include <fem.hpp>
#include "hdivsym.hpp"


namespace ngfem
{

   // dyadic product of two vectors
  template <int H, int W, typename T>
  Mat<H,W,T> DyadProd(Vec<H,T> a, Vec<W,T> b)
  {
    Mat<H,W,T> m;
    for (int i = 0; i < H; i++)
      for (int j = 0; j < W; j++)
        m(i,j) = a(i)*b(j);
    return m;
  }
    
  // symmetric dyadic product of two vectors
  template <int S, typename T>
  Mat<S,S,T> SymDyadProd(Vec<S,T> a, Vec<S,T> b)
  {
    Mat<S,S,T> m;
    for (int i = 0; i < S; i++)
      for (int j = 0; j < S; j++)
        m(i,j) = a(i)*b(j)+a(j)*b(i);
    return m;
  }
        
  HDivSymFE :: HDivSymFE (int order)
      // (order+1)*(order+2)/2-> #dof/trig
    : FiniteElement (3*(order+1)*(order+2)/2, order)
  {
       if (order < 3)
           throw Exception("In HDivSymFE: order must be 3 or higher");
  }
    
  template <class T>
  void HDivSymFE :: T_CalcShape (const T & x, const T & y, const Mat<2,2,T> & F, BareSliceMatrix<T> shape) const {

    // Barycentric base functions
    T lam[3] = { x, y, 1-x-y };
    Vec<2> pnts[3] = { { 1, 0 }, { 0, 1 } , { 0, 0 } };
      
    //sort vertices of triangle
    //INT<4> f = ET_trait<ET_TRIG>::GetFaceSort(0, vnums); 
    // no sorting
    INT<4> f = {0,1,2,3}; 
     
    // symmetric cartesian basis
    Mat<2,2> e1e1 =  {{1,0},{0,0}};
    Mat<2,2> sym_e1e2 =  {{0,0.5},{0.5,0}};
    Mat<2,2> e2e2 =  {{0,0},{0,1}};  
      
    // vertex shapes functions
    for (int i = 0; i < 3; i++) {
        shape.Row(3*i).Range(0,DIM_STRESS) = lam[f[i]]*e1e1;
        shape.Row(3*i+1).Range(0,DIM_STRESS) = lam[f[i]]*sym_e1e2;
        shape.Row(3*i+2).Range(0,DIM_STRESS) = lam[f[i]]*e2e2;
    }
      
    // there are 9 vertex base functions in total (3 per vertex)
    int ii = 9;
    
    ArrayMem<T, 20> polx(order+1), poly(order+1);
      
    // edge-based shapes
    for (int i = 0; i < 3; i++)
      if (order >= 2)
	{ 
          auto edge = GetVertexOrientedEdge(i);
          auto ls = lam[edge[0]];
          auto le = lam[edge[1]];
          
          // (unnormalized) tangential vector on reference triangle
          Vec<2> tau = pnts[edge[0]] - pnts[edge[1]];
          // (unnormalized) tangential vector on physical triangle
          Vec<2,T> t = F * tau;
          // (unnormalized) normal vector on physical triangle
          Vec<2,T> n = Vec<2,T>(t[1],-t[0]);
          
          // Dyadic products on the physical edge
          Mat<2,2,T> sym_tn = SymDyadProd(t,n);
          Mat<2,2,T> nn = DyadProd(n,n);

          // Li ((le-ls)/(le+ls)) * (le+ls)
          // * coincides with Li(le-ls) on edge
          // * vanishes on other edges with ls=0 and le=0
          // * is a polynomial of order i
          ScaledIntegratedLegendrePolynomial (order, le-ls, le+ls, polx);
          for (int j = 2; j <= order; j++) {
            // on each edge there are two base functions with connectivity 
            shape.Row(ii++).Range(0,DIM_STRESS) = polx[j] * sym_tn;
            shape.Row(ii++).Range(0,DIM_STRESS) = polx[j] * nn;
          }
	}
      
        // edge-cell shapes (no connectivity across interfaces!)
    for (int i = 0; i < 3; i++)
      if (order >= 2)
	{ 
          auto edge = GetVertexOrientedEdge(i);
          auto ls = lam[edge[0]];
          auto le = lam[edge[1]];
          
          // (unnormalized) tangential vector on reference triangle
          Vec<2> tau = pnts[edge[0]] - pnts[edge[1]];
          // (unnormalized) tangential vector on physical triangle
          Vec<2,T> t = F * tau;
          
          // Dyadic products on the physical edge
          Mat<2,2,T> tt = DyadProd(t,t);

          // Li ((le-ls)/(le+ls)) * (le+ls)
          // * coincides with Li(le-ls) on edge
          // * vanishes on other edges with ls=0 and le=0
          // * is a polynomial of order i
          ScaledIntegratedLegendrePolynomial (order, le-ls, le+ls, polx);
          for (int j = 2; j <= order; j++) {
            // on each edge there is one cell function 
            shape.Row(ii++).Range(0,DIM_STRESS) = polx[j] * tt;
          }
	}

    // inner shapes
    if (order >= 3)
      {
        //INT<4> f = GetVertexOrientedFace (0);
        T bub = lam[0]*lam[1]*lam[2];
        ScaledLegendrePolynomial (order-2, lam[1]-lam[0], lam[1]+lam[0], polx);
        LegendrePolynomial (order-1, 2*lam[2]-1, poly);

        for (int i = 0; i <= order-3; i++) {
            for (int j = 0; j <= order-3-i; j++) {
            // for each scalar base functions define 3 tensor base functions
            shape.Row(ii++).Range(0,DIM_STRESS) = bub*polx[j]*poly[i]*e1e1;
            shape.Row(ii++).Range(0,DIM_STRESS) = bub*polx[j]*poly[i]*sym_e1e2;
            shape.Row(ii++).Range(0,DIM_STRESS) = bub*polx[j]*poly[i]*e2e2;
            }
        }
      }
  }
  
  void HDivSymFE :: CalcMappedShape (const BaseMappedIntegrationPoint & bmip, 
                          BareSliceMatrix<> shape) const
  {
    auto mip = static_cast< const MappedIntegrationPoint<2,2> &>(bmip);
    IntegrationPoint ip = mip.IP();
    double x = ip(0);
    double y = ip(1);
      
    // Jacobian from mapping evaluated at integration point
    Mat<2,2> F = mip.GetJacobian();
      
    T_CalcShape (x, y, F, shape);
  }


  void HDivSymFE :: CalcMappedDShape (const MappedIntegrationPoint<2,2> & mip, 
                   BareSliceMatrix<> dshape) const
  {
    if (false)
    {
    // Jacobian from mapping evaluated at integration point
    Mat<2,2> F = mip.GetJacobian();
    // Hessian from mapping evaluated at integration point
    Mat<2,2> H[2]; 
    mip.CalcHesse (H[0], H[1]);
    
    // AutoDiff Jacobian
    Mat<2,2, AutoDiff<2>> F2;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            F2(i,j).Value()= F(i,j);
            for (int k = 0; k < 2; k++) {
                F2(i,j).DValue(k)= H[k](i,j);
            }    
        }
    }
    
    // Integration points
    AutoDiff<2> x(mip(0), 0);
    AutoDiff<2> y(mip(1), 1);
    Mat<2,2> Finv = mip.GetJacobianInverse();

    for (int j=0; j < 2; j++)
    {
      x.DValue(j) = Finv(0,j);
      y.DValue(j) = Finv(1,j);
    }
      
    Matrix <AutoDiff<2>> shapearray(ndof,DIM_STRESS);
    T_CalcShape<AutoDiff <2>> (x, y, F2, shapearray);
    for (int i = 0; i < ndof; i++) {
         dshape(i, 0) = shapearray(i,0).DValue(0) + shapearray(i,1).DValue(1);
         dshape(i, 1) = shapearray(i,2).DValue(0) + shapearray(i,3).DValue(1);
     }
    }
    else
    {
        LocalHeapMem<30000> lh("lh-hsymdiv");
        Matrix<double> gradshape(ndof,2*2*2);
        CalcDShapeFE<HDivSymFE,2,2,2*2>(*this, mip, gradshape, lh, 1e-4);
        /*
        for (int j = 0; j < DIM_STRESS; j++)
         for (int k = 0; k < nd_u; k++)
          for (int l = 0; l < DIMSPACE; l++)
            gradshape(k, l*DIM_STRESS+j) = dshape_u(k,l);
            */
        for (int i = 0; i < ndof; i++) {
         dshape(i, 0) = gradshape(i,0*4+0)+gradshape(i,1*4+1);
         dshape(i, 1) = gradshape(i,0*4+2)+gradshape(i,1*4+3);
     }
    }
          
  }
}

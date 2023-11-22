#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "ngbem.hpp"


namespace ngbem
{


  
  SingleLayerPotential :: SingleLayerPotential(shared_ptr<FESpace> aspace, int _intorder)
    : space(aspace), intorder(_intorder)
  {
    auto mesh = space->GetMeshAccess();
    // setup global-2-boundary mappings;
    BitArray bnddofs(space->GetNDof());
    bnddofs.Clear();
    for (int i = 0; i < mesh->GetNE(BND); i++)
      {
        Array<DofId> dnums;
        space->GetDofNrs(ElementId(BND, i), dnums);
        for (auto d : dnums) 
          bnddofs.SetBit(d);
      }
    
    mapglob2bnd.SetSize(space->GetNDof());
    mapglob2bnd = -1;
    for (int i = 0; i < space->GetNDof(); i++)
      if (bnddofs.Test(i))
        {
          mapglob2bnd[i] = mapbnd2glob.Size();
          mapbnd2glob.Append(i);
        }
    // cout << "dim: " << space->GetSpatialDimension() << endl;
    // cout << "bnddofs: " << bnddofs << endl;
  }

  void SingleLayerPotential ::
  CalcElementMatrix(FlatMatrix<double> matrix,  // matrix dim = ndof_bnd x ndof_bnd
                    LocalHeap &lh) const
  {
    matrix = 0;

    const double Pi = M_PI;

    int numGaussPoints = 5;
    Array<double> xi, wi;
    ComputeGaussRule(numGaussPoints, xi, wi);  

    
    IntegrationRule irsegm(ET_SEGM, intorder);
    IntegrationRule irhex (ET_HEX, intorder);    
    IntegrationRule irtrig(ET_TRIG, intorder); // order=4

    
    Array<Vec<4>> identic_Duffies;
    Array<double> identic_weights;

    static Timer tall("SingleLayer - all");
    static Timer t_identic("SingleLayer - identic panel");
    static Timer t_common_vertex("SingleLayer - common vertex");        
    static Timer t_common_edge("SingleLayer - common edge");        
    static Timer t_disjoint("SingleLayer - disjoint");
    static Timer t_disjoint2("SingleLayer - disjoint2");        
    RegionTimer reg(tall);

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          identic_Duffies.Append (xi*Vec<4>(1, 1-e1+e1*e2, 1-e1*e2*e3, 1-e1));
          identic_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, 1-e1, 1, 1-e1+e1*e2));
          identic_Duffies.Append (xi*Vec<4>(1, e1*(1-e2+e2*e3), 1-e1*e2, e1*(1-e2) ));
          identic_Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*(1-e2+e2*e3) ));
          identic_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*(1-e2) ));
          identic_Duffies.Append (xi*Vec<4>(1, e1*(1-e2), 1-e1*e2*e3, e1*(1-e2*e3) ));
          for (int j = 0; j < 6; j++)
            identic_weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());
        }


    Array<Vec<4>> common_edge_Duffies;
    Array<double> common_edge_weights;
    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          common_edge_Duffies.Append (xi*Vec<4>(1, e1*e3, 1-e1*e2, e1*(1-e2)));
          common_edge_Duffies.Append (xi*Vec<4>(1, e1, 1-e1*e2*e3, e1*e2*(1-e3)));
          common_edge_Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*e2*e3));
          common_edge_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*e2*(1-e3), 1, e1));
          common_edge_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*e2));
          
          common_edge_weights.Append (xi*xi*xi*e1*e1    * ipeta.Weight()*ipxi.Weight());
          for (int j = 0; j < 4; j++)
            common_edge_weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());          
        }


    Array<Vec<4>> common_vertex_Duffies;
    Array<double> common_vertex_weights;
    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          common_vertex_Duffies.Append (xi*Vec<4>(1, e1, e2, e2*e3 ));
          common_vertex_Duffies.Append (xi*Vec<4>(e2, e2*e3, 1, e1 ));
          for (int j = 0; j < 2; j++)
            common_vertex_weights.Append (xi*xi*xi*e2  * ipeta.Weight()*ipxi.Weight());
        }



    auto mesh = space->GetMeshAccess();
    auto evaluator = space->GetEvaluator(BND);

    
    for (int i = 0; i < mesh->GetNSE(); i++)
      for (int j = 0; j < mesh->GetNSE(); j++)
        {
          HeapReset hr(lh);
          ElementId ei(BND, i);
          ElementId ej(BND, j);


          auto verti = mesh->GetElement(ei).Vertices();
          auto vertj = mesh->GetElement(ej).Vertices();          
          
          BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ei, lh));
          BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ej, lh));

          ElementTransformation &trafoi = mesh->GetTrafo(ei, lh);
          ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
          
          Array<DofId> dnumsi, dnumsj;
          space->GetDofNrs(ei, dnumsi); // mapping to global dof
          space->GetDofNrs(ej, dnumsj);

          FlatVector<> shapei(feli.GetNDof(), lh);
          FlatVector<> shapej(felj.GetNDof(), lh);

          FlatMatrix elmat(feli.GetNDof(), felj.GetNDof(), lh); // e.g. 3 x 3
          elmat = 0;

          
          int n_common_vertices = 0;
          for (auto vi : verti)
            if (vertj.Contains(vi))
              n_common_vertices++;

          
          // Sauter-Schwab, page 240 German edition
          switch (n_common_vertices)
            {
            case 3: //identical panel
              {
                RegionTimer reg(t_identic);    
                
                elmat = 0.0;
                for (int k = 0; k < identic_weights.Size(); k++)
                  {
                    Vec<4> xy = identic_Duffies[k];
                    IntegrationPoint xhat(xy(0)-xy(1), xy(1), 0, 0);
                    IntegrationPoint yhat(xy(2)-xy(3), xy(3), 0, 0);
                      
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    // feli.CalcShape (xhat, shapei);
                    // felj.CalcShape (yhat, shapej);
                    evaluator->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
                    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_weights[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }

                // elmat = 0.0;
                // cout << "single panel elmat = " << endl << elmat << endl;
                break;
              }
            case 2: //common edge
              {
                RegionTimer reg(t_common_edge);    
                
                const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
                int cex, cey;
                for (int cx = 0; cx < 3; cx++)
                  for (int cy = 0; cy < 3; cy++)
                    {
                      INT<2> ex (verti[edges[cx][0]], verti[edges[cx][1]]);
                      INT<2> ey (vertj[edges[cy][0]], vertj[edges[cy][1]]);
                      if (ex.Sort() == ey.Sort())
                        {
                          cex = cx;
                          cey = cy;
                          break;
                      }
                    }
                int vpermx[3] = { edges[cex][0], edges[cex][1], -1 };
                vpermx[2] = 3-vpermx[0]-vpermx[1];
                int vpermy[3] = { edges[cey][1], edges[cey][0], -1 };
                vpermy[2] = 3-vpermy[0]-vpermy[1];
                
                int ivpermx[3], ivpermy[3];
                for (int i = 0; i < 3; i++)
                  {
                    ivpermx[vpermx[i]] = i;
                    ivpermy[vpermy[i]] = i;
                  }
                
                
                elmat = 0.0;
                for (int k = 0; k < common_edge_weights.Size(); k++)
                  {
                    Vec<4> xy = common_edge_Duffies[k];
                    
                    Vec<3> lamx (xy(0)-xy(1), xy(1), 1-xy(0));   // other ref-triangle
                    Vec<3> lamy (xy(2)-xy(3), xy(3), 1-xy(2));
                    // lamx0, lamx1 ... common edge
                    
                    IntegrationPoint xhat(lamx(ivpermx[0]), lamx(ivpermx[1]), 0, 0);
                    IntegrationPoint yhat(lamy(ivpermy[0]), lamy(ivpermy[1]), 0, 0);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    // feli.CalcShape (xhat, shapei);
                    // felj.CalcShape (yhat, shapej);
                    evaluator->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
                    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weights[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }

                // cout.precision(12);                
                // cout << "common edge: " << endl << elmat << endl;
                // elmat = 0.0;
                break;
              }

            case 1: //common vertex
              {
                RegionTimer reg(t_common_vertex);    
                
                int cvx=-1, cvy=-1;
                for (int cx = 0; cx < 3; cx++)
                  for (int cy = 0; cy < 3; cy++)
                    {
                      if (verti[cx] == vertj[cy])
                        {
                          cvx = cx;
                          cvy = cy;
                          break;
                        }
                    }
                
                int vpermx[3] = { cvx, (cvx+1)%3, (cvx+2)%3 };
                int vpermy[3] = { cvy, (cvy+1)%3, (cvy+2)%3 };
                
                int ivpermx[3], ivpermy[3];
                for (int i = 0; i < 3; i++)
                  {
                    ivpermx[vpermx[i]] = i;
                    ivpermy[vpermy[i]] = i;
                  }
                
                elmat = 0.0;
                for (int k = 0; k < common_vertex_weights.Size(); k++)
                  {
                    Vec<4> xy = common_vertex_Duffies[k];
                    
                    Vec<3> lamx (xy(0)-xy(1), xy(1), 1-xy(0));   // other ref-triangle
                    Vec<3> lamy (xy(2)-xy(3), xy(3), 1-xy(2));
                    // lam2 .. singular vertex
                    
                    IntegrationPoint xhat(lamx(ivpermx[2]), lamx(ivpermx[0]), 0, 0);
                    IntegrationPoint yhat(lamy(ivpermy[2]), lamy(ivpermy[0]), 0, 0);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    // feli.CalcShape (xhat, shapei);
                    // felj.CalcShape (yhat, shapej);
                    evaluator->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
                    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weights[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }

                // cout.precision(12);
                // cout << "common vertex: " << endl << elmat << endl;
                // elmat = 0.0;
                break;
              }

            case 0: //disjoint panels
              {
                RegionTimer r(t_disjoint);    
                
                elmat = 0.0;

              /*
                // naive version
              for (auto ipx : irtrig)              
                for (auto ipy : irtrig)              
                  {
                    MappedIntegrationPoint<2,3> mipx(ipx, trafoi);
                    MappedIntegrationPoint<2,3> mipy(ipy, trafoj);

                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();

                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    feli.CalcShape (ipx, shapei);
                    felj.CalcShape (ipy, shapej);
                    double fac = mipx.GetMeasure()*mipy.GetMeasure() * ipx.Weight()*ipy.Weight();
                    elmat += fac*kernel* shapei * Trans(shapej);
                  }
              */

              /*
              // shapes+geom out of loop
              auto & mirx = trafoi(irtrig, lh);
              auto & miry = trafoj(irtrig, lh);
              FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
              FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
              feli.CalcShape (irtrig, shapesi);
              felj.CalcShape (irtrig, shapesj);

              for (int ix = 0; ix < irtrig.Size(); ix++)
                for (int iy = 0; iy < irtrig.Size(); iy++)
                  {
                    Vec<3> x = mirx[ix].GetPoint();
                    Vec<3> y = miry[iy].GetPoint();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                    elmat += fac*kernel* shapesi.Col(ix) * Trans(shapesj.Col(iy));
                  }
              */

              // shapes+geom out of loop, matrix multiplication
                auto & mirx = trafoi(irtrig, lh);
                auto & miry = trafoj(irtrig, lh);
                FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
                FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
                FlatMatrix<> kernel_shapesj(felj.GetNDof(), irtrig.Size(), lh);
                
                
                // feli.CalcShape (irtrig, shapesi);
                // felj.CalcShape (irtrig, shapesj);
                evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
                evaluator -> CalcMatrix(felj, miry, Trans(shapesj), lh);
                
                RegionTimer r2(t_disjoint2);
                kernel_shapesj = 0;
                for (int ix = 0; ix < irtrig.Size(); ix++)
                  for (int iy = 0; iy < irtrig.Size(); iy++)
                    {
                      Vec<3> x = mirx[ix].GetPoint();
                      Vec<3> y = miry[iy].GetPoint();
                    
                      double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                      double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                      kernel_shapesj.Col(ix) += fac*kernel*shapesj.Col(iy);
                    }
                
                // elmat += shapesi * Trans(kernel_shapesj);
                elmat += kernel_shapesj * Trans(shapesi);
                // cout << "disjoint panel " << endl << elmat << endl;
                // elmat = 0.0;
                break;
              }
            default:
              throw Exception ("not possible");
            }
          
          for (int ii = 0; ii < dnumsi.Size(); ii++)
            for (int jj = 0; jj < dnumsj.Size(); jj++)
              matrix(mapglob2bnd[dnumsj[jj]], mapglob2bnd[dnumsi[ii]]) += elmat(jj, ii);
        }
  }

  void SingleLayerPotential :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
    //cout << "GetDofNrsSL: " << dnums << endl;
  }



  DoubleLayerPotential :: DoubleLayerPotential (shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace,
                                                int _intorder)
    : space(aspace), space2(bspace), intorder(_intorder)
  {

    auto mesh = space->GetMeshAccess();
    auto mesh2 = space2->GetMeshAccess();
    
    // setup global-2-boundary mappings;
    BitArray bnddofs(space->GetNDof());
    bnddofs.Clear();
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
        Array<DofId> dnums;
        space->GetDofNrs(ElementId(BND, i), dnums);
        for (auto d : dnums)
          bnddofs.SetBit(d);
      }
    mapglob2bnd.SetSize(space->GetNDof());
    mapglob2bnd = -1;
    for (int i = 0; i < space->GetNDof(); i++)
      if (bnddofs.Test(i))
        {
          mapglob2bnd[i] = mapbnd2glob.Size();
          mapbnd2glob.Append(i);
        }
    
    // cout << "dim1: " << space->GetSpatialDimension() << endl;
    // cout << "bnddofs1: " << bnddofs << endl;
    // cout << "mapglob2bnd: " << mapglob2bnd << endl;
    // cout << "mapbnd2glob: " << mapbnd2glob << endl;

    BitArray bnddofs2(space2->GetNDof());
    bnddofs2.Clear();
    for (int i = 0; i < mesh2->GetNSE(); i++)
      {
        Array<DofId> dnums;
        space2->GetDofNrs(ElementId(BND, i), dnums);
        for (auto d : dnums)
          bnddofs2.SetBit(d);
      }
    
    mapglob2bnd2.SetSize(space2->GetNDof());
    mapglob2bnd2 = -1;
    for (int i = 0; i < space2->GetNDof(); i++)
      if (bnddofs2.Test(i))
        {
          mapglob2bnd2[i] = mapbnd2glob2.Size();
          mapbnd2glob2.Append(i);
        }
    
    // cout << "dim2: " << space2->GetSpatialDimension() << endl;
    // cout << "bnddofs2: " << bnddofs2 << endl;
    // cout << "mapglob2bnd2: " << mapglob2bnd2 << endl;
    // cout << "mapbnd2glob2: " << mapbnd2glob2 << endl;
  }

  void DoubleLayerPotential ::
  CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd2
                    LocalHeap &lh) const
  {
    auto mesh = space->GetMeshAccess();    // trialspace = H1
    auto mesh2 = space2->GetMeshAccess();  // testspace = L2
    
    //cout << "CalcElementMatrix: " << endl;
    matrix = 0; 


    IntegrationRule irsegm(ET_SEGM, intorder);
    IntegrationRule irhex (ET_HEX, intorder);    
    IntegrationRule irtrig(ET_TRIG, intorder); // order=4

    
    Array<Vec<4>> identic_Duffies;
    Array<double> identic_weights;

    static Timer tall("SingleLayer - all");
    static Timer t_identic("SingleLayer - identic panel");
    static Timer t_common_vertex("SingleLayer - common vertex");        
    static Timer t_common_edge("SingleLayer - common edge");        
    static Timer t_disjoint("SingleLayer - disjoint");
    static Timer t_disjoint2("SingleLayer - disjoint2");        
    RegionTimer reg(tall);

    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          identic_Duffies.Append (xi*Vec<4>(1, 1-e1+e1*e2, 1-e1*e2*e3, 1-e1));
          identic_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, 1-e1, 1, 1-e1+e1*e2));
          identic_Duffies.Append (xi*Vec<4>(1, e1*(1-e2+e2*e3), 1-e1*e2, e1*(1-e2) ));
          identic_Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*(1-e2+e2*e3) ));
          identic_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*(1-e2) ));
          identic_Duffies.Append (xi*Vec<4>(1, e1*(1-e2), 1-e1*e2*e3, e1*(1-e2*e3) ));
          for (int j = 0; j < 6; j++)
            identic_weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());
        }


    Array<Vec<4>> common_edge_Duffies;
    Array<double> common_edge_weights;
    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          common_edge_Duffies.Append (xi*Vec<4>(1, e1*e3, 1-e1*e2, e1*(1-e2)));
          common_edge_Duffies.Append (xi*Vec<4>(1, e1, 1-e1*e2*e3, e1*e2*(1-e3)));
          common_edge_Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*e2*e3));
          common_edge_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*e2*(1-e3), 1, e1));
          common_edge_Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*e2));
          
          common_edge_weights.Append (xi*xi*xi*e1*e1    * ipeta.Weight()*ipxi.Weight());
          for (int j = 0; j < 4; j++)
            common_edge_weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());          
        }


    Array<Vec<4>> common_vertex_Duffies;
    Array<double> common_vertex_weights;
    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          common_vertex_Duffies.Append (xi*Vec<4>(1, e1, e2, e2*e3 ));
          common_vertex_Duffies.Append (xi*Vec<4>(e2, e2*e3, 1, e1 ));
          for (int j = 0; j < 2; j++)
            common_vertex_weights.Append (xi*xi*xi*e2  * ipeta.Weight()*ipxi.Weight());
        }


    auto evaluator = space->GetEvaluator(BND);
    auto evaluator2 = space2->GetEvaluator(BND);
    // cout << "type(eval2) = " << typeid(*evaluator2).name() << endl
    
    for (int i = 0; i < mesh->GetNSE(); i++)
      for (int j = 0; j < mesh2->GetNSE(); j++)
        {
          HeapReset hr(lh);
          ElementId ei(BND, i);
          ElementId ej(BND, j);
          //cout << "bem elements: " << ei << ", " << ej << endl;

          auto verti = mesh->GetElement(ei).Vertices();
          auto vertj = mesh->GetElement(ej).Vertices();          

          
          BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ei, lh));
          BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(space2->GetFE(ej, lh));
          
          ElementTransformation &trafoi = mesh->GetTrafo(ei, lh);
          ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
          
          Array<DofId> dnumsi, dnumsj;
          space->GetDofNrs(ei, dnumsi); // mapping to global dof
          space2->GetDofNrs(ej, dnumsj);

          FlatVector<> shapei(feli.GetNDof(), lh);
          FlatVector<> shapej(felj.GetNDof(), lh);

          FlatMatrix elmat(felj.GetNDof(), feli.GetNDof(), lh); 
          elmat = 0;


          
          int n_common_vertices = 0;
          for (auto vi : verti)
            if (vertj.Contains(vi))
              n_common_vertices++;


          
          // Sauter-Schwab, page 240 German edition
          switch (n_common_vertices)
            {
            case 3: //identical panel
              {
                RegionTimer reg(t_identic);    
                
                elmat = 0.0;
                for (int k = 0; k < identic_weights.Size(); k++)
                  {
                    Vec<4> xy = identic_Duffies[k];
                    IntegrationPoint xhat(xy(0)-xy(1), xy(1), 0, 0);
                    IntegrationPoint yhat(xy(2)-xy(3), xy(3), 0, 0);
                      
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();

                    Vec<3> nx = mipx.GetNV();
                    double nxy = InnerProduct(nx, (x-y));
                    double normxy = L2Norm(x-y);
                    double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                    
                    // feli.CalcShape (xhat, shapei);
                    // felj.CalcShape (yhat, shapej);
                    evaluator->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
                    evaluator2->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_weights[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }
                
                // cout << "new elmat = " << endl << elmat << endl;
                break;
              }
            case 2: //common edge
              {
                RegionTimer reg(t_common_edge);    
                
                const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
                int cex, cey;
                for (int cx = 0; cx < 3; cx++)
                  for (int cy = 0; cy < 3; cy++)
                    {
                      INT<2> ex (verti[edges[cx][0]], verti[edges[cx][1]]);
                      INT<2> ey (vertj[edges[cy][0]], vertj[edges[cy][1]]);
                      if (ex.Sort() == ey.Sort())
                        {
                          cex = cx;
                          cey = cy;
                          break;
                      }
                    }
                int vpermx[3] = { edges[cex][0], edges[cex][1], -1 };
                vpermx[2] = 3-vpermx[0]-vpermx[1];
                int vpermy[3] = { edges[cey][0], edges[cey][1], -1 };
                vpermy[2] = 3-vpermy[0]-vpermy[1];
                
                int ivpermx[3], ivpermy[3];
                for (int i = 0; i < 3; i++)
                  {
                    ivpermx[vpermx[i]] = i;
                    ivpermy[vpermy[i]] = i;
                  }
                
                
                elmat = 0.0;
                for (int k = 0; k < common_edge_weights.Size(); k++)
                  {
                    Vec<4> xy = common_edge_Duffies[k];
                    
                    Vec<3> lamx (xy(0)-xy(1), xy(1), 1-xy(0));   // other ref-triangle
                    Vec<3> lamy (xy(2)-xy(3), xy(3), 1-xy(2));
                    // lamx0, lamx1 ... common edge
                    
                    IntegrationPoint xhat(lamx(ivpermx[0]), lamx(ivpermx[1]), 0, 0);
                    IntegrationPoint yhat(lamy(ivpermy[0]), lamy(ivpermy[1]), 0, 0);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    // double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    Vec<3> nx = mipx.GetNV();
                    double nxy = InnerProduct(nx, (x-y));
                    double normxy = L2Norm(x-y);
                    double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                    
                    // feli.CalcShape (xhat, shapei);
                    // felj.CalcShape (yhat, shapej);
                    evaluator->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
                    evaluator2->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weights[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }
                break;
              }

            case 1: //common vertex
              {
                RegionTimer reg(t_common_vertex);    
                
                int cvx=-1, cvy=-1;
                for (int cx = 0; cx < 3; cx++)
                  for (int cy = 0; cy < 3; cy++)
                    {
                      if (verti[cx] == vertj[cy])
                        {
                          cvx = cx;
                          cvy = cy;
                          break;
                        }
                    }
                // cout << "cvx = " << cvx << ", cvy = " << cvy << endl;
                int vpermx[3] = { cvx, (cvx+1)%3, (cvx+2)%3 };
                vpermx[2] = 3-vpermx[0]-vpermx[1];
                int vpermy[3] = { cvy, (cvy+1)%3, (cvy+2)%3 };
                vpermy[2] = 3-vpermy[0]-vpermy[1];
                
                int ivpermx[3], ivpermy[3];
                for (int i = 0; i < 3; i++)
                  {
                    ivpermx[vpermx[i]] = i;
                    ivpermy[vpermy[i]] = i;
                  }
                
                elmat = 0.0;
                for (int k = 0; k < common_vertex_weights.Size(); k++)
                  {
                    Vec<4> xy = common_vertex_Duffies[k];
                    
                    Vec<3> lamx (xy(0)-xy(1), xy(1), 1-xy(0));   // other ref-triangle
                    Vec<3> lamy (xy(2)-xy(3), xy(3), 1-xy(2));
                    // lamx0, lamx1 ... common edge
                    
                    IntegrationPoint xhat(lamx(ivpermx[0]), lamx(ivpermx[1]), 0, 0);
                    IntegrationPoint yhat(lamy(ivpermy[0]), lamy(ivpermy[1]), 0, 0);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    // double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    Vec<3> nx = mipx.GetNV();
                    double nxy = InnerProduct(nx, (x-y));
                    double normxy = L2Norm(x-y);
                    double kernel = nxy / (4*M_PI*normxy*normxy*normxy);

                    // feli.CalcShape (xhat, shapei);
                    // felj.CalcShape (yhat, shapej);
                    evaluator->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
                    evaluator2->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weights[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }
                // cout << "new: common vertex elmat = " << elmat << endl;
                break;
              }

            case 0: //disjoint panels
              {
                RegionTimer r(t_disjoint);    
                
                elmat = 0.0;

                /*
                // naive version
                for (auto ipx : irtrig)              
                  for (auto ipy : irtrig)              
                    {
                      MappedIntegrationPoint<2,3> mipx(ipx, trafoi);
                      MappedIntegrationPoint<2,3> mipy(ipy, trafoj);
                      
                      Vec<3> x = mipx.Point();
                      Vec<3> y = mipy.Point();

                      // double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                      Vec<3> nx = mipx.GetNV();
                      double nxy = InnerProduct(nx, (x-y));
                      double normxy = L2Norm(x-y);
                      double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                      
                      feli.CalcShape (ipx, shapei);
                      felj.CalcShape (ipy, shapej);
                      double fac = mipx.GetMeasure()*mipy.GetMeasure() * ipx.Weight()*ipy.Weight();
                      elmat += fac*kernel* shapej * Trans(shapei);
                    }
                */

              /*
              // shapes+geom out of loop
              auto & mirx = trafoi(irtrig, lh);
              auto & miry = trafoj(irtrig, lh);
              FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
              FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
              feli.CalcShape (irtrig, shapesi);
              felj.CalcShape (irtrig, shapesj);

              for (int ix = 0; ix < irtrig.Size(); ix++)
                for (int iy = 0; iy < irtrig.Size(); iy++)
                  {
                    Vec<3> x = mirx[ix].GetPoint();
                    Vec<3> y = miry[iy].GetPoint();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                    elmat += fac*kernel* shapesi.Col(ix) * Trans(shapesj.Col(iy));
                  }
              */
                
              // shapes+geom out of loop, matrix multiplication
                // auto & mirx = trafoi(irtrig, lh);
                // auto & miry = trafoj(irtrig, lh);
                MappedIntegrationRule<2,3> mirx(irtrig, trafoi, lh);
                MappedIntegrationRule<2,3> miry(irtrig, trafoj, lh);
                
                FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
                FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
                FlatMatrix<> kernel_shapesj(felj.GetNDof(), irtrig.Size(), lh);
                
                
                evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
                evaluator2 -> CalcMatrix(felj, miry, Trans(shapesj), lh);
                
                RegionTimer r2(t_disjoint2);
                kernel_shapesj = 0;
                for (int ix = 0; ix < irtrig.Size(); ix++)
                  for (int iy = 0; iy < irtrig.Size(); iy++)
                    {
                      Vec<3> x = mirx[ix].GetPoint();
                      Vec<3> y = miry[iy].GetPoint();
                    
                      // double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                      Vec<3> nx = mirx[ix].GetNV();
                      double nxy = InnerProduct(nx, (x-y));
                      double normxy = L2Norm(x-y);
                      double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                      
                      double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                      kernel_shapesj.Col(ix) += fac*kernel*shapesj.Col(iy);
                    }
                
                // elmat += shapesi * Trans(kernel_shapesj);
                elmat += kernel_shapesj * Trans(shapesi);

                // cout << "new: disjoint elmat = " << elmat << endl;
                // cout << "dnumsj = " << dnumsj[0] << ", rowdof = " << mapglob2bnd2[dnumsj[0]] << endl;
                // cout << "dnumsi = " << dnumsi << endl;
                // elmat = 0.0;
                break;
              }
            default:
              throw Exception ("not possible");
            }

          for (int ii = 0; ii < dnumsi.Size(); ii++)
            for (int jj = 0; jj < dnumsj.Size(); jj++)
              matrix(mapglob2bnd2[dnumsj[jj]], mapglob2bnd[dnumsi[ii]]) += elmat(jj,ii);

        }
  }

  void DoubleLayerPotential :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
  }

  void DoubleLayerPotential ::  GetDofNrs2(Array<int> &dnums) const   
  {
    dnums = mapbnd2glob2;
  }

}

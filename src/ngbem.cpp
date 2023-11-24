#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "ngbem.hpp"


namespace ngbem
{

  // x, y in triangle [(0,0), (1,0), (0,1)]
  tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>> IdenticPanelIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);    

    Array<Vec<4>> Duffies;
    Array<double> weights;

    // Sauter-Schwab integration points in triangle [(0,0), (1,0), (1,1)]
    // page 240 German edition    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          Duffies.Append (xi*Vec<4>(1, 1-e1+e1*e2, 1-e1*e2*e3, 1-e1));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, 1-e1, 1, 1-e1+e1*e2));
          Duffies.Append (xi*Vec<4>(1, e1*(1-e2+e2*e3), 1-e1*e2, e1*(1-e2) ));
          Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*(1-e2+e2*e3) ));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*(1-e2) ));
          Duffies.Append (xi*Vec<4>(1, e1*(1-e2), 1-e1*e2*e3, e1*(1-e2*e3) ));
          for (int j = 0; j < 6; j++)
            weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());
        }

	
    // trafo to [(0,0), (1,0), (0,1)]
    Array<Vec<2>> ipx, ipy;
    for (auto ip : Duffies)
      {
        ipx += Vec<2>(ip(0)-ip(1), ip(1));
        ipy += Vec<2>(ip(2)-ip(3), ip(3));
      }

    return tuple { ipx, ipy, weights };
  }


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x=(0,0) and y=(0,0) are common vertices
  tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>> CommonVertexIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);    

    Array<Vec<4>> Duffies;
    Array<double> weights;

    // Sauter-Schwab integration points: [(0,0), (1,0), (1,1)]
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          Duffies.Append (xi*Vec<4>(1, e1, e2, e2*e3 ));
          Duffies.Append (xi*Vec<4>(e2, e2*e3, 1, e1 ));
          for (int j = 0; j < 2; j++)
            weights.Append (xi*xi*xi*e2  * ipeta.Weight()*ipxi.Weight());
        }

    // trafo to [(0,0), (1,0), (0,1)]
    Array<Vec<2>> ipx, ipy;
    for (auto ip : Duffies)
      {
        ipx += Vec<2>(ip(0)-ip(1), ip(1));
        ipy += Vec<2>(ip(2)-ip(3), ip(3));
      }

    return tuple { ipx, ipy, weights };
  }


  // x, y in triangle [(0,0), (1,0), (0,1)]
  // x in [(0,0),(1,0)] and y in [(0,0),(1,0)] are common edges
  tuple<Array<Vec<2>>, Array<Vec<2>>, Array<double>> CommonEdgeIntegrationRule (int order)
  {
    IntegrationRule irsegm(ET_SEGM, order);
    IntegrationRule irhex (ET_HEX, order);    

    Array<Vec<4>> Duffies;
    Array<double> weights;

    // Sauter-Schwab integration points: [(0,0), (1,0), (1,1)]
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);
          
          
          Duffies.Append (xi*Vec<4>(1, e1*e3, 1-e1*e2, e1*(1-e2)));
          Duffies.Append (xi*Vec<4>(1, e1, 1-e1*e2*e3, e1*e2*(1-e3)));
          Duffies.Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*e2*e3));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*e2*(1-e3), 1, e1));
          Duffies.Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*e2));
          
          weights.Append (xi*xi*xi*e1*e1    * ipeta.Weight()*ipxi.Weight());
          for (int j = 0; j < 4; j++)
            weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());          
        }

    Array<Vec<2>> ipx, ipy;
    for (auto ip : Duffies)
      {
        ipx += Vec<2>(ip(0)-ip(1), ip(1));
        ipy += Vec<2>(ip(2)-ip(3), ip(3));
      }

    return tuple { ipx, ipy, weights };
  }


  
  
  SingleLayerPotentialOperator :: SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, int _intorder)
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

  void SingleLayerPotentialOperator ::
  CalcElementMatrix(FlatMatrix<double> matrix,  // matrix dim = ndof_bnd x ndof_bnd
                    LocalHeap &lh) const
  {
    static Timer tall("SingleLayer - all");
    static Timer t_identic("SingleLayer - identic panel");
    static Timer t_common_vertex("SingleLayer - common vertex");        
    static Timer t_common_edge("SingleLayer - common edge");        
    static Timer t_disjoint("SingleLayer - disjoint");
    static Timer t_disjoint2("SingleLayer - disjoint2");        
    RegionTimer reg(tall);


    IntegrationRule irtrig(ET_TRIG, intorder); // order=4
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(intorder);

    auto mesh = space->GetMeshAccess();
    auto evaluator = space->GetEvaluator(BND);

    matrix = 0;
    
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

          FlatMatrix<double,ColMajor> mshapei(1, feli.GetNDof(), lh);
          FlatMatrix<double,ColMajor> mshapej(1, felj.GetNDof(), lh);

          FlatMatrix elmat(feli.GetNDof(), felj.GetNDof(), lh); // e.g. 3 x 3
          elmat = 0;

          
          // get number of common nodes for integration rule
          int n_common_vertices = 0;
          for (auto vi : verti)
            if (vertj.Contains(vi))
              n_common_vertices++;

          

          switch (n_common_vertices)
            {
            case 3: //identical panel
              {
                RegionTimer reg(t_identic);    
                
                elmat = 0.0;
                for (int k = 0; k < identic_panel_weight.Size(); k++)
                  {
                    IntegrationPoint xhat (identic_panel_x[k]);
                    IntegrationPoint yhat (identic_panel_y[k]);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));

                    evaluator->CalcMatrix(feli, mipx, mshapei, lh);
                    evaluator->CalcMatrix(felj, mipy, mshapej, lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_panel_weight[k];
                    elmat += fac*kernel* Trans(mshapej) * mshapei;
                  }

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

                // permute vertices st common edge comes first
                int vpermx[3] = { edges[cex][0], edges[cex][1], -1 };
                vpermx[2] = 3-vpermx[0]-vpermx[1];
                int vpermy[3] = { edges[cey][1], edges[cey][0], -1 };
                vpermy[2] = 3-vpermy[0]-vpermy[1];

                /*
                cout << "px = " << verti[vpermx[0]] << ", "
                     << verti[vpermx[1]] << ", "
                     << verti[vpermx[2]] << endl;
                cout << "py = " << vertj[vpermy[0]] << ", "
                     << vertj[vpermy[1]] << ", "
                     << vertj[vpermy[2]] << endl;
                */
                  
                
                elmat = 0.0;
                for (int k = 0; k < common_edge_weight.Size(); k++)
                  {
                    Vec<2> xk = common_edge_x[k];
                    Vec<2> yk = common_edge_y[k];


                    Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                    Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );

			
                    // consider permutation 
                    Vec<3> plamx, plamy;
                    for (int i = 0; i < 3; i++)
                      {
                        plamx(vpermx[i]) = lamx(i);
                        plamy(vpermy[i]) = lamy(i);
                      }
                    
                    IntegrationPoint xhat(plamx(0), plamx(1), 0, 0);
                    IntegrationPoint yhat(plamy(0), plamy(1), 0, 0);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                    
                    evaluator->CalcMatrix(feli, mipx, mshapei, lh);
                    evaluator->CalcMatrix(felj, mipy, mshapej, lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weight[k];
                    elmat += fac*kernel* Trans(mshapej) * mshapei;
                  }

                // cout.precision(12);                
                // cout << "common edge: " << endl << elmat << endl;
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
                
                elmat = 0.0;
                for (int k = 0; k < common_vertex_weight.Size(); k++)
                  {
                    Vec<2> xk = common_vertex_x[k];
                    Vec<2> yk = common_vertex_y[k];

                    Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                    Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );

                    Vec<3> plamx, plamy;
                    for (int i = 0; i < 3; i++)
                      {
                        plamx(vpermx[i]) = lamx(i);
                        plamy(vpermy[i]) = lamy(i);
                      }
                    
                    IntegrationPoint xhat(plamx(0), plamx(1), 0, 0);
                    IntegrationPoint yhat(plamy(0), plamy(1), 0, 0);
                    
                    MappedIntegrationPoint<2,3> mipx(xhat, trafoi);
                    MappedIntegrationPoint<2,3> mipy(yhat, trafoj);
                    
                    Vec<3> x = mipx.Point();
                    Vec<3> y = mipy.Point();
                    
                    double kernel = 1.0 / (4*M_PI*L2Norm(x-y));

                    evaluator->CalcMatrix(feli, mipx, mshapei, lh);
                    evaluator->CalcMatrix(felj, mipy, mshapej, lh);
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weight[k];
                    elmat += fac*kernel* Trans(mshapej) * mshapei;
                  }

                // cout.precision(12);
                // cout << "common vertex: " << endl << elmat << endl;
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
                
                elmat += kernel_shapesj * Trans(shapesi);
                */


                auto & mirx = trafoi(irtrig, lh);
                auto & miry = trafoj(irtrig, lh);
                FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
                FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
                
                evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
                evaluator -> CalcMatrix(felj, miry, Trans(shapesj), lh);
                
                RegionTimer r2(t_disjoint2);

                FlatMatrix<> kernel_ixiy(irtrig.Size(), irtrig.Size(), lh);
                for (int ix = 0; ix < irtrig.Size(); ix++)
                  {
                    for (int iy = 0; iy < irtrig.Size(); iy++)
                      {
                        Vec<3> x = mirx[ix].GetPoint();
                        Vec<3> y = miry[iy].GetPoint();
                        
                        double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
                        double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                        kernel_ixiy(ix, iy) = fac*kernel;
                      }
                  }

                FlatMatrix<double> kernel_shapesj(felj.GetNDof(), irtrig.Size(), lh);
                kernel_shapesj = shapesj * Trans(kernel_ixiy);
                elmat += kernel_shapesj * Trans(shapesi);

                
                // cout << "disjoint panel " << endl << elmat << endl;
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

  
  void SingleLayerPotentialOperator :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
  }



  DoubleLayerPotentialOperator :: DoubleLayerPotentialOperator (shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace,
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

  void DoubleLayerPotentialOperator ::
  CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd_H1 x ndof_bnd_L2
                    LocalHeap &lh) const
  {
    auto mesh = space->GetMeshAccess();    // trialspace = H1
    auto mesh2 = space2->GetMeshAccess();  // testspace = L2
    
    static Timer tall("DoubleLayer - all");
    static Timer t_identic("DoubleLayer - identic panel");
    static Timer t_common_vertex("DoubleLayer - common vertex");        
    static Timer t_common_edge("DoubleLayer - common edge");        
    static Timer t_disjoint("DoubleLayer - disjoint");
    static Timer t_disjoint2("DoubleLayer - disjoint2");        
    RegionTimer reg(tall);


    IntegrationRule irtrig(ET_TRIG, intorder); // order=4
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(intorder);

    //cout << "CalcElementMatrix: " << endl;
    matrix = 0; 

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


          switch (n_common_vertices)
            {
            case 3: //identical panel
              {
                RegionTimer reg(t_identic);    
                
                elmat = 0.0;
                for (int k = 0; k < identic_panel_weight.Size(); k++)
                  {
                    IntegrationPoint xhat (identic_panel_x[k]);
                    IntegrationPoint yhat (identic_panel_y[k]);
			
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
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_panel_weight[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }
                
                // cout << "new elmat = " << endl << elmat << endl;
                break;
              }
            case 2: //common edge
              {
                RegionTimer reg(t_common_edge);    
                
                const EDGE * edges = ElementTopology::GetEdges (ET_TRIG); // 0 1 | 1 2 | 2 0 
                int cex, cey;
                for (int cx = 0; cx < 3; cx++)
                  for (int cy = 0; cy < 3; cy++)
                    {
                      INT<2> ex (verti[edges[cx][0]], verti[edges[cx][1]]); 
                      INT<2> ey (vertj[edges[cy][0]], vertj[edges[cy][1]]); 
                      if (ex.Sort() == ey.Sort()) 
                        {
                          cex = cx;  // -> "common" edge number triangle i
                          cey = cy;  // -> "common" edge number triangle j
                          break;
                      }
                    }
                int vpermx[3] = { edges[cex][0], edges[cex][1], -1 }; // common edge gets first
                vpermx[2] = 3-vpermx[0]-vpermx[1]; 
                int vpermy[3] = { edges[cey][1], edges[cey][0], -1 }; // common edge gets first
                vpermy[2] = 3-vpermy[0]-vpermy[1];
                
                
                elmat = 0.0;
                for (int k = 0; k < common_edge_weight.Size(); k++)
                  {

                    Vec<2> xk = common_edge_x[k]; // int point in [(0,0), (1,0), (0,1)] 
                    Vec<2> yk = common_edge_y[k]; // int point in [(0,0), (1,0), (0,1)] 

                    Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) ); 
                    Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );

			
                    // consider permuation 
                    Vec<3> plamx, plamy;
                    for (int i = 0; i < 3; i++)
                      {
                        plamx(vpermx[i]) = lamx(i);
                        plamy(vpermy[i]) = lamy(i);
                      }

                    IntegrationPoint xhat(plamx(0), plamx(1), 0, 0); // int point in [(0,1),(1,0),(0,0)]
                    IntegrationPoint yhat(plamy(0), plamy(1), 0, 0);
                    
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
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weight[k];
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
                
                elmat = 0.0;
                for (int k = 0; k < common_vertex_weight.Size(); k++)
                  {
                    
                    Vec<2> xk = common_vertex_x[k];
                    Vec<2> yk = common_vertex_y[k];

                    Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                    Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );

                    Vec<3> plamx, plamy;
                    for (int i = 0; i < 3; i++)
                      {
                        plamx(vpermx[i]) = lamx(i);
                        plamy(vpermy[i]) = lamy(i);
                      }
                    
                    IntegrationPoint xhat(plamx(0), plamx(1), 0, 0);
                    IntegrationPoint yhat(plamy(0), plamy(1), 0, 0);
                    
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
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weight[k];
                    elmat += fac*kernel* shapej * Trans(shapei);
                  }
                // cout << "new: common vertex elmat = " << elmat << endl;
                break;
              }

            case 0: //disjoint panels
              {
                RegionTimer r(t_disjoint);    
                
                elmat = 0.0;
                
              // shapes+geom out of loop, matrix multiplication
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
                    
                      Vec<3> nx = mirx[ix].GetNV();
                      double nxy = InnerProduct(nx, (x-y));
                      double normxy = L2Norm(x-y);
                      double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                      
                      double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
                      kernel_shapesj.Col(ix) += fac*kernel*shapesj.Col(iy);
                    }
                
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

  void DoubleLayerPotentialOperator :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
  }

  void DoubleLayerPotentialOperator ::  GetDofNrs2(Array<int> &dnums) const   
  {
    dnums = mapbnd2glob2;
  }

}

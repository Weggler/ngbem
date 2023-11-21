#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "ngbem.hpp"


namespace ngbem
{

  bool commonEdge(Vec<3,int> verti, Vec<3,int> vertj) {
    bool test = false;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        test = test || (verti[i] == vertj[j] && verti[(i+1)%3] == vertj[(j+1)%3]);
        test = test || (verti[i] == vertj[j] && verti[(i+2)%3] == vertj[(j+2)%3]);
        test = test || (verti[i] == vertj[j] && verti[(i+1)%3] == vertj[(j+2)%3]);
        test = test || (verti[i] == vertj[j] && verti[(i+2)%3] == vertj[(j+1)%3]);
      }
    }
    return test;
  }

  bool commonVertex(Vec<3> verti, Vec<3> vertj) {
    bool test = false;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        test = test || (verti[i] == vertj[j]);
      }
    }
    return test;
  }

  void commonEdgeVec(Vec<3> &u, Vec<3> &v, Vec<3> &w, Vec<3> verti, Vec<3> vertj, const std::shared_ptr<MeshAccess>& mesh) {
    Vec<3> cords1, cords2;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        if (verti[i] == vertj[j] && verti[(i+1)%3] == vertj[(j+1)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+1)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+2)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+2)%3], cords1);
          w = cords1-cords2;
          return;
        }
        if (verti[i] == vertj[j] && verti[(i+2)%3] == vertj[(j+2)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+2)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+1)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+1)%3], cords1);
          w = cords1-cords2;
          return;
        }
        if (verti[i] == vertj[j] && verti[(i+1)%3] == vertj[(j+2)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+1)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+2)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+1)%3], cords1);
          w = cords1-cords2;
          return;
        }
        if (verti[i] == vertj[j] && verti[(i+2)%3] == vertj[(j+1)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+2)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+1)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+2)%3], cords1);
          w = cords1-cords2;
          return;
        }
      }
    }
  }

  void commonEdgeVec2(Vec<3> &u, Vec<3> &v, Vec<3> &w, Vec<3> &nt, int* idmap, Vec<3> verti, Vec<3> vertj, const std::shared_ptr<MeshAccess>& mesh) {
    Vec<3> cords1, cords2;
    mesh->GetPoint(vertj[0], cords1);
    mesh->GetPoint(vertj[1], cords2);
    nt = cords2 - cords1;
    mesh->GetPoint(vertj[2], cords2);
    nt = Cross(nt, cords2 - cords1);
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        if (verti[i] == vertj[j] && verti[(i+1)%3] == vertj[(j+1)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+1)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+2)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+2)%3], cords1);
          w = cords1-cords2;
          idmap[0] = j;
          idmap[1] = (j+1)%3;
          idmap[2] = (j+2)%3;
          return;
        }
        if (verti[i] == vertj[j] && verti[(i+2)%3] == vertj[(j+2)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+2)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+1)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+1)%3], cords1);
          w = cords1-cords2;
          idmap[0] = j;
          idmap[1] = (j+2)%3;
          idmap[2] = (j+1)%3;
          return;
        }
        if (verti[i] == vertj[j] && verti[(i+1)%3] == vertj[(j+2)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+1)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+2)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+1)%3], cords1);
          w = cords1-cords2;
          idmap[0] = j;
          idmap[1] = (j+2)%3;
          idmap[2] = (j+1)%3;
          return;
        }
        if (verti[i] == vertj[j] && verti[(i+2)%3] == vertj[(j+1)%3]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+2)%3], cords2);
          u = cords2-cords1;
          mesh->GetPoint(verti[(i+1)%3], cords1);
          v = cords1-cords2;
          mesh->GetPoint(vertj[(j+2)%3], cords1);
          w = cords1-cords2;
          idmap[0] = j;
          idmap[1] = (j+1)%3;
          idmap[2] = (j+2)%3;
          return;
        }
      }
    }
  }

  void commonVertexVec(Vec<3> &r1, Vec<3> &r2, Vec<3> &rt1, Vec<3> &rt2, Vec<3> verti, Vec<3> vertj, const std::shared_ptr<MeshAccess>& mesh) {
    Vec<3> cords1, cords2;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        if (verti[i] == vertj[j]) {
          mesh->GetPoint(verti[i], cords1);
          mesh->GetPoint(verti[(i+1)%3], cords2);
          r1 = cords2 - cords1;
          mesh->GetPoint(verti[(i+2)%3], cords1);
          r2 = cords1 - cords2;
          mesh->GetPoint(vertj[j], cords1);
          mesh->GetPoint(vertj[(j+1)%3], cords2);
          rt1 = cords2 - cords1;
          mesh->GetPoint(vertj[(j+2)%3], cords1);
          rt2 = cords1 - cords2;
          return;
        }
      }
    }
  }

  void commonVertexVec2(Vec<3> &r1, Vec<3> &r2, Vec<3> &rt1, Vec<3> &rt2, Vec<3> &nt, int* idmap, Vec<3> verti, Vec<3> vertj, const std::shared_ptr<MeshAccess>& mesh) {
    Vec<3> cords1, cords2;
    mesh->GetPoint(vertj[0], cords1);
    mesh->GetPoint(vertj[1], cords2);
    nt = cords2 - cords1;
    mesh->GetPoint(vertj[2], cords2);
    nt = Cross(nt, cords2 - cords1);
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        if (verti[i] == vertj[j]) {
          idmap[0] = i;
          if (verti[(i+1)%3] < verti[(i+2)%3]) {
            idmap[1] = (i+1)%3;
            idmap[2] = (i+2)%3;
          }
          else{
            idmap[1] = (i+2)%3;
            idmap[2] = (i+1)%3;
          }
          mesh->GetPoint(verti[idmap[0]], cords1);
          mesh->GetPoint(verti[idmap[1]], cords2);
          r1 = cords2 - cords1;
          mesh->GetPoint(verti[idmap[2]], cords1);
          r2 = cords1 - cords2;

          idmap[0] = j;
          if (vertj[(j+1)%3] < vertj[(j+2)%3]) {
            idmap[1] = (j+1)%3;
            idmap[2] = (j+2)%3;
          }
          else{
            idmap[1] = (j+2)%3;
            idmap[2] = (j+1)%3;
          }
          mesh->GetPoint(vertj[idmap[0]], cords1);
          mesh->GetPoint(vertj[idmap[1]], cords2);
          rt1 = cords2 - cords1;
          mesh->GetPoint(vertj[idmap[2]], cords1);
          rt2 = cords1 - cords2;
          return;
        }
      }
    }
  }




  //static bool commonEdge(Vec<3> verti, Vec<3> vertj);

  SingleLayerPotential :: SingleLayerPotential(shared_ptr<FESpace> aspace)
    : space(aspace) 
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
    
    IntegrationRule irsegm(ET_SEGM, 10);
    IntegrationRule irhex (ET_HEX, 10);    
    IntegrationRule irtrig(ET_TRIG, 10); // order=4

    
    Array<Vec<4>> identic_Duffies[6];
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

          identic_Duffies[0].Append (xi*Vec<4>(1, 1-e1+e1*e2, 1-e1*e2*e3, 1-e1));
          identic_Duffies[1].Append (xi*Vec<4>(1-e1*e2*e3, 1-e1, 1, 1-e1+e1*e2));
          identic_Duffies[2].Append (xi*Vec<4>(1, e1*(1-e2+e2*e3), 1-e1*e2, e1*(1-e2) ));
          identic_Duffies[3].Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*(1-e2+e2*e3) ));
          identic_Duffies[4].Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*(1-e2) ));
          identic_Duffies[5].Append (xi*Vec<4>(1, e1*(1-e2), 1-e1*e2*e3, e1*(1-e2*e3) ));
          identic_weights.Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());
        }


    Array<Vec<4>> common_edge_Duffies[5];
    Array<double> common_edge_weights[5];
    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          common_edge_Duffies[0].Append (xi*Vec<4>(1, e1*e3, 1-e1*e2, e1*(1-e2)));
          common_edge_Duffies[1].Append (xi*Vec<4>(1, e1, 1-e1*e2*e3, e1*e2*(1-e3)));
          common_edge_Duffies[2].Append (xi*Vec<4>(1-e1*e2, e1*(1-e2), 1, e1*e2*e3));
          common_edge_Duffies[3].Append (xi*Vec<4>(1-e1*e2*e3, e1*e2*(1-e3), 1, e1));
          common_edge_Duffies[4].Append (xi*Vec<4>(1-e1*e2*e3, e1*(1-e2*e3), 1, e1*e2));
          
          common_edge_weights[0].Append (xi*xi*xi*e1*e1    * ipeta.Weight()*ipxi.Weight());
          common_edge_weights[1].Append (xi*xi*xi*e1*e1*e2 * ipeta.Weight()*ipxi.Weight());          
        }
    common_edge_weights[2] = common_edge_weights[1];
    common_edge_weights[3] = common_edge_weights[1];
    common_edge_weights[4] = common_edge_weights[1];




    Array<Vec<4>> common_vertex_Duffies[2];
    Array<double> common_vertex_weights;
    
    for (auto ipeta : irhex)
      for (auto ipxi : irsegm)
        {
          double e1 = ipeta(0);
          double e2 = ipeta(1);
          double e3 = ipeta(2);
          double xi = ipxi(0);

          common_vertex_Duffies[0].Append (xi*Vec<4>(1, e1, e2, e2*e3 ));
          common_vertex_Duffies[1].Append (xi*Vec<4>(e2, e2*e3, 1, e1 ));
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
          //cout << "bem elements: " << ei << ", " << ej << endl;

          Ngs_Element ngsEli = mesh->GetElement(ei);
          auto tmp_verti = ngsEli.Vertices();
          
          Vec<3,int> verti;
          verti[0] = tmp_verti[0];
          verti[1] = tmp_verti[1];
          verti[2] = tmp_verti[2];
          
          //cout << verti << endl;
          Ngs_Element ngsElj = mesh->GetElement(ej);
          auto tmp_vertj = ngsElj.Vertices();
          Vec<3,int> vertj;
          vertj[0] = tmp_vertj[0];
          vertj[1] = tmp_vertj[1];
          vertj[2] = tmp_vertj[2];
          //cout << vertj << endl;


          BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ei, lh));
          BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ej, lh));

          ElementTransformation &trafoi = mesh->GetTrafo(ei, lh);
          ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
          
          Array<DofId> dnumsi, dnumsj;
          space->GetDofNrs(ei, dnumsi); // mapping to global dof
          space->GetDofNrs(ej, dnumsj);

          FlatVector<> shapei(feli.GetNDof(), lh);
          FlatVector<> shapej(felj.GetNDof(), lh);
          // cout << dnumsi.length() << endl;
          //  return;

          FlatMatrix elmat(feli.GetNDof(), felj.GetNDof(), lh); // e.g. 3 x 3
          elmat = 0;

          
          int intCase = -1;
          // 2d cases
          if (space->GetSpatialDimension() == 2) {
            if (i == j)
              intCase = 0;
            else if (verti[0] == vertj[0] || verti[0] == vertj[1] ||
                     verti[1] == vertj[0] || verti[1] == vertj[1])
              intCase = 1;
            else
              intCase = 2;
          }

          // 3d cases
          if (space->GetSpatialDimension() == 3) {
            if (i == j)
              intCase = 0;
            else if (commonEdge(verti, vertj))
              intCase = 1;
            else if (commonVertex(verti, vertj)) 
              intCase = 2;
            else
              intCase = 3;
          }

          //cout << "Case: " << intCase << endl;

          // feli.GetNDof() == dnumsi.Size()

          double I = 0;
          Matrix tmpmat(feli.GetNDof(), felj.GetNDof());

          if (space->GetSpatialDimension() == 2) {
            // 2d cases

            if (intCase == 0) {
              Vec<3> int_start, int_end, start, end;
              int_start = 0;
              int_end = 0;
              int_end[0] = 1;

              IntegrationPoint ip_start(int_start, wi[0]);
              IntegrationPoint ip_end(int_end, wi[0]);

              trafoi.CalcPoint(ip_start, start);
              trafoi.CalcPoint(ip_end, end);

              double seg_length = L2Norm(end-start);
              I = seg_length * seg_length * ( log(seg_length) - 1.5);
              elmat = I;
            }
            else {
              //case 1 + 2, to do Guass-Log
              for (int ki = 0; ki < xi.Size(); ki++)
                {
                  for (int kj = 0; kj < xi.Size(); kj++)
                    {
                      Vec<3> xi3, xj3, pxi, pxj;
                      xi3 = 0;
                      xj3 = 0;
                      xi3[0] = xi[ki];
                      xj3[0] = xi[kj];

                      IntegrationPoint ipi(xi3, wi[ki]);
                      IntegrationPoint ipj(xj3, wi[kj]);

                      trafoi.CalcPoint(ipi, pxi);
                      trafoj.CalcPoint(ipj, pxj);

                      I += wi[ki] * wi[kj] * log(L2Norm(pxi-pxj));
                    }
                }

              Vec<3> int_start, int_end, start, end;
              int_start = 0;
              int_end = 0;
              int_end[0] = 1;

              IntegrationPoint p_start(int_start, wi[0]);
              IntegrationPoint p_end(int_end, wi[0]);

              trafoi.CalcPoint(p_start, start);
              trafoi.CalcPoint(p_end, end);

              double segi_length = L2Norm(end-start);

              trafoj.CalcPoint(p_start, start);
              trafoj.CalcPoint(p_end, end);

              double segj_length = L2Norm(end-start);

              I = I * segi_length * segj_length;
              elmat = I;
            }
          }
          else
            {
              // 3d cases
              
              if (intCase == 0) {
                //identical panel
                RegionTimer reg(t_identic);    
    
                

                // Sauter-Schwab, page 240 German edition
                elmat = 0.0;
                for (int k = 0; k < identic_weights.Size(); k++)
                  for (int l = 0; l < 6; l++)
                    {
                      Vec<4> xy = identic_Duffies[l][k];
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
                      elmat += fac*kernel* shapei * Trans(shapej);
                    }
                
                // cout << "new elmat = " << endl << elmat << endl;
            }
            else if (intCase == 1) {
              //common edge
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

              
              // Sauter-Schwab, page 240 German edition
              elmat = 0.0;
              for (int k = 0; k < common_edge_weights[0].Size(); k++)
                for (int l = 0; l < 5; l++)
                  {
                    Vec<4> xy = common_edge_Duffies[l][k];

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
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weights[l][k];
                    elmat += fac*kernel* shapei * Trans(shapej);
                  }
            }
            else if (intCase == 2) {
              //common vertex
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
              vpermx[2] = 3-vpermx[0]-vpermx[1];
              int vpermy[3] = { cvy, (cvy+1)%3, (cvy+2)%3 };
              vpermy[2] = 3-vpermy[0]-vpermy[1];

              int ivpermx[3], ivpermy[3];
              for (int i = 0; i < 3; i++)
                {
                  ivpermx[vpermx[i]] = i;
                  ivpermy[vpermy[i]] = i;
                }

              // Sauter-Schwab, page 240 German edition
              elmat = 0.0;
              for (int k = 0; k < common_vertex_weights.Size(); k++)
                for (int l = 0; l < 2; l++)
                  {
                    Vec<4> xy = common_vertex_Duffies[l][k];

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
                    
                    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weights[k];
                    elmat += fac*kernel* shapei * Trans(shapej);
                  }
            }
            else {
              //disjoint panels
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
                    // double fac = irtrig[ix].Weight()*irtrig[iy].Weight();
                    kernel_shapesj.Col(ix) += fac*kernel*shapesj.Col(iy);
                  }

              elmat += shapesi * Trans(kernel_shapesj);
            }
          }

          for (int ii = 0; ii < dnumsi.Size(); ii++)
            for (int jj = 0; jj < dnumsj.Size(); jj++)
              matrix(mapglob2bnd[dnumsi[ii]], mapglob2bnd[dnumsj[jj]]) += elmat(ii, jj);
        }
  }

  void SingleLayerPotential :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
    //cout << "GetDofNrsSL: " << dnums << endl;
  }



  DoubleLayerPotential :: DoubleLayerPotential (shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace)
    : space(aspace), space2(bspace)
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
    cout << "dim1: " << space->GetSpatialDimension() << endl;

    cout << "bnddofs1: " << bnddofs << endl;
    //cout << "mapglob2bnd: " << mapglob2bnd << endl;
    //cout << "mapbnd2glob: " << mapbnd2glob << endl;

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
    cout << "dim2: " << space2->GetSpatialDimension() << endl;

    cout << "bnddofs2: " << bnddofs2 << endl;
    //cout << "mapglob2bnd2: " << mapglob2bnd2 << endl;
    //cout << "mapbnd2glob2: " << mapbnd2glob2 << endl;
  }

  void DoubleLayerPotential :: CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd2
                                         LocalHeap &lh) const
  {
    auto mesh = space->GetMeshAccess();
    auto mesh2 = space2->GetMeshAccess();

    
    //cout << "CalcElementMatrix: " << endl;
    matrix = 0; 

    double Pi = 3.14159265358979323;

    int intOrder = 5;
    int intCase;

    Array<DofId> dnumsi, dnumsj;
    Array<double> xi, wi, xi2, wi2;

    ComputeGaussRule(intOrder, xi, wi);
    ComputeGaussRule(2, xi2, wi2);

    Array<double> quadTriPoints1, quadTriPoints2, quadTriWeights;
    quadTriPoints1.SetSize(7);
    quadTriPoints2.SetSize(7);
    quadTriWeights.SetSize(7);
    quadTriPoints1[0] = 0.333333333333333;
    quadTriPoints2[0] = 0.333333333333333;

    quadTriPoints1[1] = 0.470142064105115;
    quadTriPoints2[1] = 0.470142064105115;

    quadTriPoints1[2] = 0.470142064105115;
    quadTriPoints2[2] = 0.05971587178977;

    quadTriPoints1[3] = 0.05971587178977;
    quadTriPoints2[3] = 0.470142064105115;

    quadTriPoints1[4] = 0.101286507323456;
    quadTriPoints2[4] = 0.101286507323456;

    quadTriPoints1[5] = 0.101286507323456;
    quadTriPoints2[5] = 0.797426985353087;

    quadTriPoints1[6] = 0.797426985353087;
    quadTriPoints2[6] = 0.101286507323456;

    quadTriWeights[0] = 0.1125;
    quadTriWeights[1] = 0.066197076394253;
    quadTriWeights[2] = 0.066197076394253;
    quadTriWeights[3] = 0.066197076394253;
    quadTriWeights[4] = 0.0629695902724135;
    quadTriWeights[5] = 0.0629695902724135;
    quadTriWeights[6] = 0.0629695902724135;

    for (int i = 0; i < mesh->GetNSE(); i++)
      for (int j = 0; j < mesh2->GetNSE(); j++)
        {
          HeapReset hr(lh);
          ElementId ei(BND, i);
          ElementId ej(BND, j);
          //cout << "bem elements: " << ei << ", " << ej << endl;

          Ngs_Element ngsEli = mesh->GetElement(ei);
          auto tmp_verti = ngsEli.Vertices();
          Vec<3> verti;
          verti[0] = tmp_verti[0];
          verti[1] = tmp_verti[1];
          verti[2] = tmp_verti[2];

          //cout << verti << endl;
          Ngs_Element ngsElj = mesh2->GetElement(ej);
          auto tmp_vertj = ngsElj.Vertices();
          Vec<3> vertj;
          vertj[0] = tmp_vertj[0];
          vertj[1] = tmp_vertj[1];
          vertj[2] = tmp_vertj[2];
          //cout << vertj << endl;

          BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ei, lh));
          //cout << feli << endl;
          BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(space2->GetFE(ej, lh));
          //cout << felj << endl;

          ElementTransformation &trafoi = mesh->GetTrafo(ei, lh);
          ElementTransformation &trafoj = mesh2->GetTrafo(ej, lh);
          Array<DofId> dnumsi, dnumsj;
          space->GetDofNrs(ei, dnumsi); // mapping to global dof
          space2->GetDofNrs(ej, dnumsj);

          //cout << dnumsi.Size() << endl;
          //cout << dnumsj.Size() << endl;
          //  return;

          Matrix elmat(feli.GetNDof(), felj.GetNDof()); // e.g. 3 x 1
          elmat = 0;
          /*
            for (auto d : dnumsi)
            cout << d;

            cout << " x ";

            for (auto d : dnumsj)
            cout << d;

            cout << endl;
          */
          intCase = -1;
          // 2d cases
          if (space->GetSpatialDimension() == 2) {
            if (i == j)
              intCase = 0;
            else if (verti[0] == vertj[0] || verti[0] == vertj[1] ||
                     verti[1] == vertj[0] || verti[1] == vertj[1])
              intCase = 1;
            else
              intCase = 2;
          }

          // 3d cases
          if (space->GetSpatialDimension() == 3) {
            if (i == j)
              intCase = 0;
            else if (commonEdge(verti, vertj))
              intCase = 1;
            else if (commonVertex(verti, vertj)) 
              intCase = 2;
            else
              intCase = 3;
          }

          //cout << "Case: " << intCase << endl;

          // 3d cases
          if (intCase == 0) {

          }
          else if (intCase == 1) {
            //common edge
          
            double elmat1, elmat2, elmat3;
            int idmap[3];
            Vec<3> u, v, w, Ji;
            commonEdgeVec2(u, v, w, Ji, idmap, verti, vertj, mesh);
            /*cout << u << endl << v << endl << w << endl;
              cout << "u" << u << endl;
              cout << "v" << v << endl;
              cout << "w" << w << endl;
            */
            double J = L2Norm(Cross(u, v));
            //Vec<3> Ji = Cross(u, w);
            double JT = L2Norm(Ji);

            double n1 = Ji[0]/JT, n2 = Ji[1]/JT, n3 = Ji[2]/JT;
            //cout << "nt " << n1 << " " << n2 << " " << " " << n3 << endl;
            J *= JT;

            Vec<3> xy;
            elmat1 = 0;
            elmat2 = 0;
            elmat3 = 0;

            double eta1 = 0.5;
            for (int k=0; k < xi.Size(); k++)
              {
                double eta2 = xi[k];
                double weta2 = wi[k];
                double eta12 = eta1*eta2;
                for (int l=0; l < xi.Size(); l++)
                  {
                    double eta3 = xi[l];
                    double weta3 = wi[l];
                    double eta123 = eta12*eta3;
                    xy = eta2*(eta3*(u+w)-w)+v;
                    double normxy = L2Norm(xy);
                    double xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                    double H = xyn/normxy/normxy/normxy;
                    double I1ij = wi2[0]*xi2[0]*(1.0-xi2[0]*(1.0-eta123));
                    I1ij+= wi2[1]*xi2[1]*(1.0-xi2[1]*(1.0-eta123));
                    I1ij *= H;
                    double I2ij = wi2[0]*xi2[0] * xi2[0]*(1.0-eta12);
                    I2ij += wi2[1]*xi2[1] * xi2[1]*(1.0-eta12);
                    I2ij *= H;
                    double I3ij = wi2[0]*xi2[0] * xi2[0]*eta12*(1.0-eta3);
                    I3ij += wi2[1]*xi2[1] * xi2[1]*eta12*(1.0-eta3);
                    I3ij *= H;
              
                    xy = -eta2*(u+v+eta3*w)+v;
                    normxy = L2Norm(xy);
                    xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                    H = xyn/normxy/normxy/normxy;
                    I1ij += wi2[0]*xi2[0] * (1.0-xi2[0]) * H;
                    I1ij += wi2[1]*xi2[1] * (1.0-xi2[1]) * H;
                    I2ij += wi2[0]*xi2[0] * xi2[0]*(1.0-eta123) * H;
                    I2ij += wi2[1]*xi2[1] * xi2[1]*(1.0-eta123) * H;
                    I3ij += wi2[0]*xi2[0] * xi2[0]*eta123 * H;
                    I3ij += wi2[1]*xi2[1] * xi2[1]*eta123 * H;
              
                    xy = eta2*(v-eta3*(u+v))-w;
                    normxy = L2Norm(xy);
                    xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                    H = xyn/normxy/normxy/normxy;
                    I1ij += wi2[0]*xi2[0] * (1.0-xi2[0]) * H;
                    I1ij += wi2[1]*xi2[1] * (1.0-xi2[1]) * H;
                    I2ij += wi2[0]*xi2[0] * xi2[0]*(1.0-eta1) * H;
                    I2ij += wi2[1]*xi2[1] * xi2[1]*(1.0-eta1) * H;
                    I3ij += wi2[0]*xi2[0] * xi2[0]*eta1 * H;
                    I3ij += wi2[1]*xi2[1] * xi2[1]*eta1 * H;
              
                    xy = -eta2*(w+eta3*(u+v))+v;
                    normxy = L2Norm(xy);
                    xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                    H = xyn/normxy/normxy/normxy;
                    I1ij += wi2[0]*xi2[0] * (1.0-xi2[0]) * H;
                    I1ij += wi2[1]*xi2[1] * (1.0-xi2[1]) * H;
                    I2ij += wi2[0]*xi2[0] * xi2[0]*(1.0-eta12) * H;
                    I2ij += wi2[1]*xi2[1] * xi2[1]*(1.0-eta12) * H;
                    I3ij += wi2[0]*xi2[0] * xi2[0]*eta12 * H;
                    I3ij += wi2[1]*xi2[1] * xi2[1]*eta12 * H;
              
                    I1ij *= eta2;
                    I2ij *= eta2;
                    I3ij *= eta2;
              
                    xy = eta2*(u+w)+eta3*v-w;
                    normxy = L2Norm(xy);
                    xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                    H = xyn/normxy/normxy/normxy;
                    I1ij += wi2[0]*xi2[0]*(1.0-xi2[0]*(1.0-eta12)) * H;
                    I1ij += wi2[1]*xi2[1]*(1.0-xi2[1]*(1.0-eta12)) * H;
                    I2ij += wi2[0]*xi2[0] * xi2[0]*(1.0-eta1) * H;
                    I2ij += wi2[1]*xi2[1] * xi2[1]*(1.0-eta1) * H;
                    I3ij += wi2[0]*xi2[0] * xi2[0]*eta1*(1.0-eta2) * H;
                    I3ij += wi2[1]*xi2[1] * xi2[1]*eta1*(1.0-eta2) * H;

                    elmat1 += weta2*weta3*I1ij;
                    elmat2 += weta2*weta3*I2ij;
                    elmat3 += weta2*weta3*I3ij;
                  }
              }

            elmat(idmap[0],0) = elmat1 * J / 4 / Pi;
            elmat(idmap[1],0) = elmat2 * J / 4 / Pi;
            elmat(idmap[2],0) = elmat3 * J / 4 / Pi;
          
            //cout << I << endl;
            //cout << 3 << "-d , P: " << p1 << p2 << p3 << endl;
          }
          else if (intCase == 2) {
            //common vertex
          
            Vec<3> u, v, w, z, J1;
            double elmat1, elmat2, elmat3;
            int idmap[3];
            commonVertexVec2(u, v, w, z, J1, idmap, verti, vertj, mesh);
            /*cout << "u" << u << endl;
              cout << "v" << v << endl;
              cout << "w" << w << endl;
              cout << "z" << z << endl;
              cout << r1 << endl << r2 << endl << rt1 << endl << rt2 << endl;*/
            double J = L2Norm(Cross(u, v));
            double JT = L2Norm(J1);
            double n1 = J1[0]/JT, n2 = J1[1]/JT, n3 = J1[2]/JT;
            //cout << "nt" << n1 << " " << n2 << " " << " " << n3 << endl;
            J *= JT;

            Vec<3> xy;
            elmat1 = 0;
            elmat2 = 0;
            elmat3 = 0;

            for (int k = 0; k < xi.Size(); k++)
              {
                double eta1 = xi[k];
                double weta1 = wi[k];
                for (int l = 0; l < xi.Size(); l++)
                  {
                    double eta2 = xi[l];
                    double weta2 = wi[l];
                    for (int m = 0; m < xi.Size(); m++)
                      {
                        double eta3 = xi[m];
                        double weta3 = wi[m];
                        xy = u+eta1*v-eta2*(w+eta3*z);
                        double normxy = L2Norm(xy);
                        double xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                        double H = xyn/normxy/normxy/normxy;
                        double I1ijk = wi2[0]*xi2[0]*(1.0-xi2[0]*eta2);
                        I1ijk += wi2[1]*xi2[1]*(1.0-xi2[1]*eta2);
                        I1ijk *= H;
                        double I2ijk = wi2[0]*xi2[0]*xi2[0]*eta2*(1.0-eta3);
                        I2ijk += wi2[1]*xi2[1]*xi2[1]*eta2*(1.0-eta3);
                        I2ijk *= H;
                        double I3ijk = wi2[0]*xi2[0]*xi2[0]*eta2*eta3;
                        I3ijk += wi2[1]*xi2[1]*xi2[1]*eta2*eta3;
                        I3ijk *= H;
                
                        xy = eta2*(u+eta3*v)-(w+eta1*z);
                        normxy = L2Norm(xy);
                        xyn = xy[0]*n1+xy[1]*n2+xy[2]*n3;
                        H = xyn/normxy/normxy/normxy;
                        I1ijk += wi2[0]*xi2[0]*(1.0-xi2[0]) * H;
                        I1ijk += wi2[1]*xi2[1]*(1.0-xi2[1]) * H;
                        I2ijk += wi2[0]*xi2[0]*xi2[0]*(1.0-eta1) * H;
                        I2ijk += wi2[1]*xi2[1]*xi2[1]*(1.0-eta1) * H;
                        I3ijk += wi2[0]*xi2[0]*xi2[0]*eta1 * H;
                        I3ijk += wi2[1]*xi2[1]*xi2[1]*eta1 * H;
                
                        elmat1 += weta1*weta2*weta3 * eta2 * I1ijk;
                        elmat2 += weta1*weta2*weta3 * eta2 * I2ijk;
                        elmat3 += weta1*weta2*weta3 * eta2 * I3ijk;
                      }
                  }
              }

            elmat(idmap[0],0) = J * elmat1 / 4.0 / Pi;
            elmat(idmap[1],0) = J * elmat2 / 4.0 / Pi;
            elmat(idmap[2],0) = J * elmat3 / 4.0 / Pi;
            //cout << "idmap:" << idmap[0] << idmap[1] << idmap[2] << endl;
          
          
            //cout << I << endl;
            //cout << 3 << "-d , P: " << p1 << p2 << p3 << endl;
          }
          else {
            //disjoint panels
          
            Vec<3> A,u,v,B,w,z,nt;
            mesh->GetPoint(verti[0], A);
            mesh->GetPoint(verti[1], u);
            mesh->GetPoint(verti[2], v);
            mesh->GetPoint(vertj[0], B);
            mesh->GetPoint(vertj[1], w);
            mesh->GetPoint(vertj[2], z);
            u = u - A;
            v = v - A;
            w = w - B;
            z = z - B;
          
            nt = Cross(w,z);
            double J = L2Norm(nt);
            nt = 1/J * nt;
            J *= L2Norm(Cross(u, v));

            for (int k=0; k<7; k++)
              {
                double ksi1 = quadTriPoints1[k];
                double ksi2 = quadTriPoints2[k];
                double wksi = quadTriWeights[k];
                double x1 = A[0]+ksi1*u[0]+ksi2*v[0];
                double x2 = A[1]+ksi1*u[1]+ksi2*v[1];
                double x3 = A[2]+ksi1*u[2]+ksi2*v[2];
                for (int l=0; l<7; l++)
                  {
                    double eta1 = quadTriPoints1[l];
                    double eta2 = quadTriPoints2[l];
                    double weta = quadTriWeights[l];
                    double xy1 = x1-B[0]-eta1*w[0]-eta2*z[0];
                    double xy2 = x2-B[1]-eta1*w[1]-eta2*z[1];
                    double xy3 = x3-B[2]-eta1*w[2]-eta2*z[2];
                    double normxy = sqrt(xy1*xy1+xy2*xy2+xy3*xy3);
                    double xyn = xy1*nt[0]+xy2*nt[1]+xy3*nt[2];
                    double H = xyn/normxy/normxy/normxy;
                    elmat(0,0) += wksi*weta * (1.0f-eta1-eta2) * H;
                    elmat(1,0) += wksi*weta * eta1 * H;
                    elmat(2,0) += wksi*weta * eta2 * H;
                  }
              }

            elmat(0,0) = J * elmat(0,0) / 4.0 / Pi;
            elmat(1,0) = J * elmat(1,0) / 4.0 / Pi;
            elmat(2,0) = J * elmat(2,0) / 4.0 / Pi;
          
          }

          // feli.GetNDof() == dnumsi.Size()
          /*
            cout << elmat << endl;
            cout << i << endl << "x" << endl;
            cout << vertj[0] << endl;
            cout << vertj[1] << endl;
            cout << vertj[2] << endl;
          */

          for(int jj = 0; jj < 3; jj++)
            matrix(i, vertj[jj]) += elmat(jj, 0);

          /*for (int ii = 0; ii < dnumsi.Size(); ii++)    // 3
            for (int jj = 0; jj < dnumsj.Size(); jj++)  // 1
            matrix(mapglob2bnd2[dnumsj[jj]], mapglob2bnd[dnumsi[ii]]) += elmat(ii, jj);
          */
        }                                           
  }

  void DoubleLayerPotential :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
    //cout << "GetDOFNrs: " << dnums << endl;
  }

  void DoubleLayerPotential ::  GetDofNrs2(Array<int> &dnums) const   
  {
    dnums = mapbnd2glob2;
    //cout << "GetDOFNrs2: " << dnums << endl;
  }

}

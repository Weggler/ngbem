#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "ngbem.hpp"
#include "hmat.hpp"


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
  
  SingleLayerPotentialOperator :: SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, struct BEMParameters _param)
    : space(aspace), param(_param)
  {
    auto mesh = space->GetMeshAccess();
    cluster_tree = make_shared<ClusterTree>(space, param.leafsize);
    
    // setup global-2-boundary mappings:
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

    //cout << "dim: " << space->GetSpatialDimension() << endl;
    //cout << "dofs: " << mapbnd2glob << endl;

    // run through bnd elements, for each element get its dofs and add the elem
    /*
    elems4dof.SetSize(space->GetNDof());
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
	ElementId ei(BND, i);
      
	Array<DofId> dnumsi;
	space->GetDofNrs(ei, dnumsi); 
	for (int ii = 0; ii < dnumsi.Size(); ii++)
	  {
	    elems4dof[dnumsi[ii]].Append(i);
	  }
      }
    */
    // Table-creator creates table with one big block of memory,
    // avoids memory fragmentation 
    TableCreator<int> creator;
    Array<DofId> dnumsi;
    for ( ; !creator.Done(); creator++)
      for (int i = 0; i < mesh->GetNSE(); i++)
        {
          space->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator.Add (d, i);
        }
    elems4dof = creator.MoveTable();
    //cout << "elems4dof: " << elems4dof << endl;

    /*START: TEST hmatrix: compare approximation with dense matrix. */   

    // create hmatrix
    hmatrix = make_shared<HMatrix>(cluster_tree, cluster_tree,
                                   param.eta, space->GetNDof(), space->GetNDof());
    // compute all its blocks
    LocalHeap lh(10000000);
    CalcHMatrix(*hmatrix, lh, param);
    cout << "HMatrix done " << endl;


    if (param.testhmatrix)
      {
        HeapReset hr(lh);    
        
        Matrix<double> dense(space->GetNDof(), space->GetNDof());
        CalcElementMatrix(dense, lh);
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // Test with vector
        Vector<double> x(space->GetNDof()), y(space->GetNDof());
        x = 1.;
        y = 0.;
        y = dense * x;
        
        S_BaseVectorPtr<> x_base(space->GetNDof(), 1, x.Data());
        S_BaseVectorPtr<> y_base(space->GetNDof(), 1, y.Data());    
        hmatrix->MultAdd(-1., x_base, y_base);
        
        double err = 0.;
        for (int i = 0; i < y.Size(); i++)
          err += y(i) * y(i);
        cout << "error " << sqrt(err) << endl;
      }
    /*END: TEST hmatrix: compare approximation with dense matrix. */   
  }

  void SingleLayerPotentialOperator ::
  CalcElementMatrix(FlatMatrix<double> matrix,  // matrix dim = ndof_bnd x ndof_bnd
		    LocalHeap &lh) const
  {
    Array<int> range;
    for (int i = 0; i < space->GetNDof(); i++)
      {
	range.Append(i);
      }
    CalcBlockMatrix(matrix, range, range, lh);
  }

  /* compute single layer matrix for given dofs - dofs are global dof numbers*/
  void SingleLayerPotentialOperator ::
  CalcBlockMatrix(FlatMatrix<double> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
		  LocalHeap &lh) const
  {
    static Timer t("ngbem - SLP::CalcBlockMatrix"); RegionTimer regt(t);
    
    static Timer tall("SingleLayer - block");
    static Timer t_identic("SingleLayer - identic panel");
    static Timer t_common_vertex("SingleLayer - common vertex");        
    static Timer t_common_edge("SingleLayer - common edge");        
    static Timer t_disjoint("SingleLayer - disjoint");
    static Timer t_disjoint2("SingleLayer - disjoint2");        

    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder); // order=4
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(param.intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(param.intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(param.intorder);

    auto mesh = space->GetMeshAccess();
    auto evaluator = space->GetEvaluator(BND);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv;
    Array<int> testdofsinv;
    trialdofsinv.SetSize(space->GetNDof()); 
    testdofsinv.SetSize(space->GetNDof());

    trialdofsinv = -1;
    testdofsinv = -1;
  
    // determine support (patch) for each trial function 
    for (int j = 0; j < trialdofs.Size(); j++)
      {
	trialdofsinv[trialdofs[j]] = j;
	tmp.Append(elems4dof[ trialdofs[j] ]);
      }
    QuickSort( tmp );
    for (int j = 0; j < tmp.Size(); j++)
      {
	patchj.Append(tmp[j]);
	int tmpj = tmp[j];
	while (j < tmp.Size() && tmp[j] == tmpj)
	  j++;
	j--;
      }
    int nj = patchj.Size(); 
  
    // determine support (patch) for each test function 
    for (int i = 0; i < testdofs.Size(); i++)
      {
	testdofsinv[testdofs[i]] = i; 
	tmp2.Append(elems4dof[ testdofs[i] ]);
      }
    QuickSort( tmp2 );
    for (int i = 0; i < tmp2.Size(); i++)
      {
	patchi.Append(tmp2[i]);
	int tmpi = tmp2[i];
	while (i < tmp2.Size() && tmp2[i] == tmpi)
	  i++;
	i--;
      }
    int ni = patchi.Size(); 

    matrix = 0; 

    // compute element matrices for all patch elements - this is more than we actually need
    for (int i = 0; i < ni; i++) 
      for (int j = 0; j < nj; j++)
	{
	  HeapReset hr(lh);
	  ElementId ei(BND, patchi[i]);
	  ElementId ej(BND, patchj[j]);
            
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
	        // RegionTimer reg(t_identic);    
                    
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
		    elmat += fac*kernel* Trans(mshapei) * mshapej;
		  }
              
	        // cout << "single panel elmat = " << endl << elmat << endl;
	        break;
	      }
                  
	    case 2: //common edge
	      {
		// RegionTimer reg(t_common_edge);    
                        
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
		    elmat += fac*kernel* Trans(mshapei) * mshapej;
		  }
                
	        // cout.precision(12);                
	        // cout << "common edge: " << endl << elmat << endl;
	        break;
	      }
              
	    case 1: //common vertex
	      {
		// RegionTimer reg(t_common_vertex);    
                    
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
		    elmat += fac*kernel* Trans(mshapei) * mshapej;
		  }
                
		// cout.precision(12);
		// cout << "common vertex: " << endl << elmat << endl;
		break;
	      }
            
	    case 0: //disjoint panels
	      {
		// RegionTimer r(t_disjoint);    
                  
		elmat = 0.0;
          
		auto & mirx = trafoi(irtrig, lh);
		auto & miry = trafoj(irtrig, lh);
		FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
		FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
                  
		evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
		evaluator -> CalcMatrix(felj, miry, Trans(shapesj), lh);
                  
		// RegionTimer r2(t_disjoint2);
          
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
          
		FlatMatrix<double> kernel_shapesj(irtrig.Size(), felj.GetNDof(), lh);
		kernel_shapesj = kernel_ixiy * Trans(shapesj);
		elmat += shapesi * kernel_shapesj;
          
		// cout << "disjoint panel " << endl << elmat << endl;
		break;
	      }
	    default:
	      throw Exception ("not possible");
	    }
          
	  // add element matrix entries to global matrix entries  
	  for (int ii = 0; ii < dnumsi.Size(); ii++)
	    for (int jj = 0; jj < dnumsj.Size(); jj++)
	      if(testdofsinv[dnumsi[ii]] != -1 && trialdofsinv[dnumsj[jj]] != -1 ) // consider only trial and test function
		matrix(testdofsinv[dnumsi[ii]], trialdofsinv[dnumsj[jj]]) += elmat(ii, jj);
	}
  }

  void StochasticTSVD1 (MatrixView<> A, MatrixView<> U, MatrixView<> V, VectorView<> S)
  {
    for (int i = 0; i < V.Height(); i++)
      for (int j = 0; j < V.Width(); j++)
        V(i,j) = double (rand()) / RAND_MAX;

    Matrix AVt = A * Trans(V);
    Matrix tmp(V.Height(), V.Height());
    LapackSVD (AVt, U, tmp, S, false);    
    Matrix UtA = Trans(U) * A;
    LapackSVD (UtA, tmp, V, S, false);
    U = Matrix<> (U * tmp);
  }

  size_t StochasticTSVD (MatrixView<> A, MatrixView<> U, MatrixView<> V, VectorView<> S, double eps)
  {
    int rank = 5;
    int p = min(A.Height(), A.Width());
    while (rank < p)
      {
        StochasticTSVD1 (A, U.Cols(rank), V.Rows(rank), S.Range(rank));
        if (S[rank-1] < eps)
          {
            for (int j = 1; j < p; j++)
              if (S[j] < eps)
                return j-1;
          }

        rank = int (rank*1.5);
      }
    
    LapackSVD(A, U, V, S, false);
    for (int j = 1; j < p; j++)
      if (S[j] < eps)
        return j-1;
    return p;
  }
  

  unique_ptr<LowRankMatrix> SingleLayerPotentialOperator ::
  CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, LocalHeap &lh) const
  {
    static Timer t("ngbem - SLP::CalcFarFieldBlock"); RegionTimer reg(t);
    static Timer tblock("ngbem - SLP::CalcFarFieldBlock, block");
    static Timer tsvd("ngbem - SLP::CalcFarFieldBlock, SVD");    
    
    int m = testdofs.Size();
    int n = trialdofs.Size();
    int p = min(n, m);
  
    Matrix<double> A(m, n);
    {
      RegionTimer rblock(tblock);
      CalcBlockMatrix(A, trialdofs, testdofs, lh);
    }
    // Calculate SVD for A^\top = V S U^\top
    Matrix<double> V(p, n), Ut(m, p);
    Vector<> S(p);
    // Array<double> work(n * m + 100);
    // integer info;
    // char jobu = 'S', jobv = 'S';
    // integer lda = Trans(A).Dist(), ldu = Ut.Dist(), ldv = V.Dist();
    // integer lwork = work.Size();
    int k;
    {
      RegionTimer rsvd(tsvd);
      // dgesvd_(&jobv, &jobu, &n, &m, Trans(A).Data(), &lda, S.Data(), V.Data(), &ldv,
      // Ut.Data(), &ldu, work.Data(), &lwork, &info);
      // needs ngsolve update, prefer to call NGSolve wrapper function:
      // LapackSVD (A, Trans(Ut), Trans(V), S, false);
      k = StochasticTSVD (A, Ut, V, S, param.eps);
    }
    //Truncate according to eps. k is the rank
    
    // if (n > 200 || m > 200)
    // *testout << "svd, n,m = " << n << ", " << m << ", k = " << k << endl;
    
    // Low-rank approximation from truncated svd
    Matrix<double> U_trc(m, k), Vt_trc(k, n);
    for (size_t j = 0; j < k; j++)
      U_trc.Col(j) = sqrt(S(j)) * Ut.Col(j);

    for (size_t i = 0; i < k; i++)    
      Vt_trc.Row(i) = sqrt(S(i)) * V.Row(i);

    return make_unique<LowRankMatrix> (std::move(U_trc), std::move(Vt_trc));
  }

  void SingleLayerPotentialOperator :: CalcHMatrix(HMatrix & hmatrix, LocalHeap & clh, struct BEMParameters &param) const
  {
    static Timer t("ngbem - SLP::CalcHMatrix"); RegionTimer reg(t);    
    auto & matList = hmatrix.GetMatList();

    // run through all blocks of the #HMatrix
    // for (int k = 0; k < matList.Size(); k++)

    ParallelForRange (matList.Size(), [&](IntRange r)
    {
      LocalHeap lh = clh.Split();
      for (int k : r)
        {
          HeapReset hr(lh);
          BEMBlock & block = matList[k];
          auto trialdofs = block.GetTrialDofs();
          auto testdofs = block.GetTestDofs();
          if(block.IsNearField())
            {
              // Compute dense block
              Matrix<> near(testdofs.Size(), trialdofs.Size());
              CalcBlockMatrix(near, trialdofs, testdofs, lh);	    
              block.SetMat(make_unique<BaseMatrixFromMatrix>(std::move(near)));
            }
          else
            {
              // Compute low-rank block
              block.SetMat(CalcFarFieldBlock(trialdofs, testdofs, lh));	    
            }
        }
    }, TasksPerThread(4));
  }

  void SingleLayerPotentialOperator :: Apply(FlatVector<double> elx, FlatVector<double> ely, 
					     LocalHeap & lh) const
  {
    static Timer t("ngbem - SLP::Apply"); RegionTimer reg(t);
    /*
    for (int i = 0; i < ely.Size(); i++)
      ely(i) = 0.;
    S_BaseVectorPtr<> xp_base(elx.Size(), 1, elx.Data());
    S_BaseVectorPtr<> yp_base(ely.Size(), 1, ely.Data());
    hmatrix->MultAdd(1., xp_base, yp_base);
    */
    ely = 0.0;
    VFlatVector<> xp_base(elx);
    VFlatVector<> yp_base(ely);
    hmatrix->MultAdd(1., xp_base, yp_base);
  }
  
  void SingleLayerPotentialOperator :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
  }

  // aspace, space == H1 trialspace, bspace, space2 ==  L2 testspace
  DoubleLayerPotentialOperator ::
  DoubleLayerPotentialOperator (shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace,
                                struct BEMParameters _param)
    : space(aspace), space2(bspace), param(_param)
  {

    auto mesh = space->GetMeshAccess(); // trialspace
    auto mesh2 = space2->GetMeshAccess(); // testspace


    cluster_tree = make_shared<ClusterTree>(aspace, param.leafsize);
    cluster_tree2 = make_shared<ClusterTree>(bspace, param.leafsize);
    
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

    // run through surface elements and add elem to related trial dofs
    /*
    elems4dof.SetSize(space->GetNDof());
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
        ElementId ei(BND, i);
          
        Array<DofId> dnumsi;
        space->GetDofNrs(ei, dnumsi); 
        for (int ii = 0; ii < dnumsi.Size(); ii++)
          elems4dof[dnumsi[ii]].Append(i);
      }
    */
    // Table-creator creates table with one big block of memory,
    // avoids memory fragmentation 
    TableCreator<int> creator;
    Array<DofId> dnumsi;
    for ( ; !creator.Done(); creator++)
      for (int i = 0; i < mesh->GetNSE(); i++)
        {
          space->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator.Add (d, i);
        }
    elems4dof = creator.MoveTable();
    //cout << "elems4dof: " << elems4dof << endl;

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

    // run through surface elements and add elem to related test dofs
/*
    elems4dof2.SetSize(space2->GetNDof());
    for (int i = 0; i < mesh2->GetNSE(); i++)
      {
        ElementId ei(BND, i);
            
        Array<DofId> dnumsi;
        space2->GetDofNrs(ei, dnumsi); 
        for (int ii = 0; ii < dnumsi.Size(); ii++)
          {
            elems4dof2[dnumsi[ii]].Append(i);
          }
      }
*/
    TableCreator<int> creator2;
    for ( ; !creator2.Done(); creator2++)
      for (int i = 0; i < mesh2->GetNSE(); i++)
        {
          space2->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator2.Add (d, i);
        }
    elems4dof2 = creator2.MoveTable();
    //cout << "elems4dof: " << elems4dof2 << endl;

    // create hmatrix
    hmatrix = make_shared<HMatrix>(cluster_tree, // trial space H1
                                   cluster_tree2, // test space L2
                                   param.eta, space->GetNDof(), space2->GetNDof());

    LocalHeap lh(100000000);
    CalcHMatrix(*hmatrix, lh, param);

    /*START: TEST hmatrix: compare approximation with dense matrix. */
    if (param.testhmatrix)
      {
        Matrix<double> dense(space2->GetNDof(), space->GetNDof()); // ndof(L2) x ndof(H1)
        CalcElementMatrix(dense, lh);
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // compute all its blocks
        HeapReset hr(lh);    
        
        // Test with vector
        Vector<double> x(space->GetNDof()), y(space2->GetNDof());
        x = 1.;
        y = 0.;
        y = dense * x;
        
        S_BaseVectorPtr<> x_base(space->GetNDof(), 1, x.Data());
        S_BaseVectorPtr<> y_base(space2->GetNDof(), 1, y.Data());    
        hmatrix->MultAdd(-1., x_base, y_base);
  
        cout << "error " << L2Norm (y) << endl;
    }
    /*END: TEST hmatrix: compare approximation with dense matrix. */   
  }

  /* compute dense double layer matrix,  dim = ndof_bnd_L2 x ndof_bnd_H1 */
  void DoubleLayerPotentialOperator ::
  CalcElementMatrix(FlatMatrix<double> matrix, LocalHeap &lh) const
  {
    Array<int> range, range2;
    for (int i = 0; i < space->GetNDof(); i++) // trial H1
      range.Append(i);
    for (int j = 0; j < space2->GetNDof(); j++) // test L2
      range2.Append(j);
    CalcBlockMatrix(matrix, range, range2, lh);
  }

  /* compute double layer matrix for given dof sets  - global boundry dof numbers*/
  void DoubleLayerPotentialOperator ::
  CalcBlockMatrix(FlatMatrix<double> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
		  LocalHeap &lh) const 
  {
    auto mesh = space->GetMeshAccess();    // trialspace = H1
    auto mesh2 = space2->GetMeshAccess();  // testspace = L2
    
    static Timer tall("ngbem DLP - all");
    static Timer tloops("ngbem DLP - loops");    
    static Timer t_identic("ngbem DLP - identic panel");
    static Timer t_common_vertex("ngbem DLP - common vertex");        
    static Timer t_common_edge("ngbem DLP - common edge");        
     static Timer t_disjoint("ngbem DLP - disjoint");
    // static Timer t_disjoint2("ngbem DLP - disjoint2");        

    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder); // order=4
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(param.intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(param.intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(param.intorder);

    //cout << "CalcElementMatrix: " << endl;
    matrix = 0; 

    auto evaluator = space->GetEvaluator(BND);
    auto evaluator2 = space2->GetEvaluator(BND);
    // cout << "type(eval2) = " << typeid(*evaluator2).name() << endl

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> testdofsinv;
    Array<int> trialdofsinv;

    trialdofsinv.SetSize(space->GetNDof()); 
    testdofsinv.SetSize(space2->GetNDof());

    trialdofsinv = -1;
    testdofsinv = -1;
  
    for (int i = 0; i < testdofs.Size(); i++)
      {
	testdofsinv[testdofs[i]] = i;
	tmp2.Append(elems4dof2[ testdofs[i] ]);
      }
    QuickSort( tmp2 );
    for (int i = 0; i < tmp2.Size(); i++)
      {
	patchi.Append(tmp2[i]);
	int tmpi = tmp2[i];
	while (i < tmp2.Size() && tmp2[i] == tmpi)
	  i++;
	i--;
      }
    int ni = patchi.Size(); 
    
    for (int j = 0; j < trialdofs.Size(); j++)
      {
	trialdofsinv[trialdofs[j]] = j;
	tmp.Append(elems4dof[ trialdofs[j] ]);
      }
    QuickSort( tmp );
    for (int j = 0; j < tmp.Size(); j++)
      {
	patchj.Append(tmp[j]);
	int tmpj = tmp[j];
	while (j < tmp.Size() && tmp[j] == tmpj)
	  j++;
	j--;
      }
    int nj = patchj.Size(); 


    RegionTimer regloops(tloops);    
    for (int i = 0; i < ni; i++) // test
      for (int j = 0; j < nj; j++) // trial
	{
	  HeapReset hr(lh);
	  ElementId ei(BND, patchi[i]);
	  ElementId ej(BND, patchj[j]);
	  //cout << "bem elements: " << ei << ", " << ej << endl;
          
	  auto verti = mesh2->GetElement(ei).Vertices();
	  auto vertj = mesh->GetElement(ej).Vertices();          
            
	  BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(space2->GetFE(ei, lh));
	  BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(space->GetFE(ej, lh));
              
	  ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
	  ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
              
	  Array<DofId> dnumsi, dnumsj;
	  space2->GetDofNrs(ei, dnumsi); // mapping to global dof
	  space->GetDofNrs(ej, dnumsj);
        
	  FlatVector<> shapei(feli.GetNDof(), lh);
	  FlatVector<> shapej(felj.GetNDof(), lh);
          
	  FlatMatrix elmat(feli.GetNDof(), felj.GetNDof(), lh); 
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
                
		    Vec<3> ny = mipy.GetNV();
		    double nxy = InnerProduct(ny, (x-y));
		    double normxy = L2Norm(x-y);
		    double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                        
		    evaluator2->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
		    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                        
		    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_panel_weight[k];
		    elmat += fac*kernel* shapei * Trans(shapej);
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

		    Vec<3> ny = mipy.GetNV();
		    double nxy = InnerProduct(ny, (x-y));
		    double normxy = L2Norm(x-y);
		    double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                        
		    evaluator2->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
		    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                        
		    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weight[k];
		    elmat += fac*kernel* shapei * Trans(shapej);
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
                  
		    Vec<3> ny = mipy.GetNV();
		    double nxy = InnerProduct(ny, (x-y));
		    double normxy = L2Norm(x-y);
		    double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
          
		    evaluator2->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
		    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                  
		    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weight[k];
		    elmat += fac*kernel* shapei * Trans(shapej);
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
		// FlatMatrix<> kernel_shapesj(felj.GetNDof(), irtrig.Size(), lh);
                
		evaluator2 -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
		evaluator-> CalcMatrix(felj, miry, Trans(shapesj), lh);

                /*
		// RegionTimer r2(t_disjoint2);
		kernel_shapesj = 0;
		for (int ix = 0; ix < irtrig.Size(); ix++)
		  for (int iy = 0; iy < irtrig.Size(); iy++)
		    {
		      Vec<3> x = mirx[ix].GetPoint();
		      Vec<3> y = miry[iy].GetPoint();
                  
		      Vec<3> ny = miry[iy].GetNV();
		      double nxy = InnerProduct(ny, (x-y));
		      double normxy = L2Norm(x-y);
		      double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                    
		      double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
		      kernel_shapesj.Col(ix) += fac*kernel*shapesj.Col(iy);
		    }
                  
		elmat += shapesi * Trans(kernel_shapesj);
                */

		FlatMatrix<> kernel_ixiy(irtrig.Size(), irtrig.Size(), lh);
		for (int ix = 0; ix < irtrig.Size(); ix++)
		  {
		    for (int iy = 0; iy < irtrig.Size(); iy++)
		      {
			Vec<3> x = mirx[ix].GetPoint();
			Vec<3> y = miry[iy].GetPoint();

                        Vec<3> ny = miry[iy].GetNV();
                        double nxy = InnerProduct(ny, (x-y));
                        double normxy = L2Norm(x-y);
                        double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                        
			// double kernel = 1.0 / (4*M_PI*L2Norm(x-y));
			double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
			kernel_ixiy(ix, iy) = fac*kernel;
		      }
		  }
          
		FlatMatrix<double> kernel_shapesj(irtrig.Size(), felj.GetNDof(), lh);
		kernel_shapesj = kernel_ixiy * Trans(shapesj);
		elmat += shapesi * kernel_shapesj;


                
		// cout << "new: disjoint elmat = " << elmat << endl;
		// cout << "dnumsj = " << dnumsj[0] << ", rowdof = " << mapglob2bnd2[dnumsj[0]] << endl;
		// cout << "dnumsi = " << dnumsi << endl;
		break;
	      }
	    default:
	      throw Exception ("not possible");
	    }
	  
	  for (int ii = 0; ii < dnumsi.Size(); ii++) // test
	    for (int jj = 0; jj < dnumsj.Size(); jj++) // trial
	      if(trialdofsinv[dnumsj[jj]] != -1 && testdofsinv[dnumsi[ii]] != -1)
		matrix(testdofsinv[dnumsi[ii]], trialdofsinv[dnumsj[jj]]) += elmat(ii, jj);

	}
  }

  unique_ptr<LowRankMatrix> DoubleLayerPotentialOperator ::
  CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, LocalHeap &lh) const
  {
    static Timer t("ngbem - DLP::CalcFarFieldBlock"); RegionTimer reg(t);
    static Timer tsvd("ngbem - DLP::CalcFarFieldBlock svd"); 
    int m = testdofs.Size();
    int n = trialdofs.Size();
    int p = min(n, m);
    
    Matrix<double> A(m, n);
    CalcBlockMatrix(A, trialdofs, testdofs, lh);
    
    // Calculate SVD for A^\top = V S U^\top
    Matrix<double, ColMajor> V(n, p), Ut(p, m);
    Vector<> S(p);
    Array<double> work(n * m + 100);
    integer info;
    char jobu = 'S', jobv = 'S';
    integer lda = Trans(A).Dist(), ldu = Ut.Dist(), ldv = V.Dist();
    integer lwork = work.Size();
    int k;
    {
      RegionTimer reg(tsvd);          
      // dgesvd_(&jobv, &jobu, &n, &m, Trans(A).Data(), &lda, S.Data(), V.Data(), &ldv,
      // Ut.Data(), &ldu, work.Data(), &lwork, &info);
      
      k = StochasticTSVD (A, Trans(Ut), Trans(V), S, param.eps);      
    }
    //Truncate according to eps. k is the rank
    
    // Low-rank approximation from truncated svd
    // if (n > 200 || m > 200)
    // *testout << "svd, n,m = " << n << ", " << m << ", k = " << k << endl;
    
    Matrix<double> U_trc(m, k), Vt_trc(k, n);
    for (size_t j = 0; j < k; j++)
      U_trc.Col(j) = sqrt(S(j)) * Ut.Row(j);

    for (size_t i = 0; i < k; i++)    
      Vt_trc.Row(i) = sqrt(S(i)) * V.Col(i);

    return make_unique<LowRankMatrix> (Matrix<>(U_trc), Matrix<>(Vt_trc));
  }

  void DoubleLayerPotentialOperator :: CalcHMatrix(HMatrix & hmatrix, LocalHeap &clh, struct BEMParameters &param) const
  {
    static Timer t("ngbem - DLP::CalcHMatrix"); RegionTimer reg(t);    
    auto & matList = hmatrix.GetMatList();

    // run through all blocks of the #HMatrix    
    // for (int k = 0; k < matList.Size(); k++)

    ParallelForRange (matList.Size(), [&](IntRange r)
    {
      LocalHeap lh = clh.Split();
      for (int k : r)
        {
	HeapReset hr(lh);
          
	BEMBlock & block = matList[k];
	auto trialdofs = block.GetTrialDofs();
	auto testdofs = block.GetTestDofs();
	if(block.IsNearField())
	  {
	    //// Compute dense block
	    Matrix<> near(testdofs.Size(), trialdofs.Size());
	    CalcBlockMatrix(near, trialdofs, testdofs, lh);
	    block.SetMat(make_unique<BaseMatrixFromMatrix>(near));
	  }
	else
	  {
	    // Compute low-rank block
	    block.SetMat(CalcFarFieldBlock(trialdofs, testdofs, lh));	    
	  }
      }
    }, TasksPerThread(4));
  }

  void DoubleLayerPotentialOperator :: Apply(FlatVector<double> elx, FlatVector<double> ely, 
					     LocalHeap & lh) const
  {
    static Timer t("ngbem - SLP::Apply"); RegionTimer reg(t);
    /*
    for (int i = 0; i < ely.Size(); i++)
      ely(i) = 0.;
    S_BaseVectorPtr<> xp_base(elx.Size(), 1, elx.Data());
    S_BaseVectorPtr<> yp_base(ely.Size(), 1, ely.Data());
    hmatrix->MultAdd(1., xp_base, yp_base);
    */
    VVector<> vx(hmatrix->Width());
    VVector<> vy(hmatrix->Height());
    *testout << "Apply, hmat.h = " << hmatrix->Height() << ", w = " << hmatrix->Width() << endl;
    vy = 0.0;
    for (int i = 0; i < mapbnd2glob.Size(); i++)
      vx(mapbnd2glob[i]) = elx[i];
    hmatrix -> MultAdd (1, vx, vy);
    for (int i = 0; i < mapbnd2glob2.Size(); i++)
      ely[i] = vy(mapbnd2glob2[i]);
    
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

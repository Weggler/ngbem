#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "intrules.hpp"
#include "ngbem.hpp"
#include "hmat.hpp"


namespace ngbem
{

  template <typename T>  
  IntegralOperator<T> ::
  IntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                   BEMParameters _param)
    : trial_space(_trial_space), test_space(_test_space), param(_param)
  {
    if (!test_space)
      test_space = trial_space;
    
    trial_ct = make_shared<ClusterTree>(trial_space, param.leafsize);
    if (trial_space == test_space)
      test_ct = trial_ct; // the same
    else
      test_ct = make_shared<ClusterTree>(test_space, param.leafsize);


    auto mesh = trial_space->GetMeshAccess(); // trialspace
    auto mesh2 = test_space->GetMeshAccess(); // testspace

    
    // setup global-2-boundary mappings;
    BitArray bnddofs(trial_space->GetNDof());
    bnddofs.Clear();
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
	Array<DofId> dnums;
	trial_space->GetDofNrs(ElementId(BND, i), dnums);
	for (auto d : dnums)
	  bnddofs.SetBit(d);
      }
    mapglob2bnd.SetSize(trial_space->GetNDof());
    mapglob2bnd = -1;
    for (int i = 0; i < trial_space->GetNDof(); i++)
      if (bnddofs.Test(i))
	{
	  mapglob2bnd[i] = mapbnd2glob.Size();
	  mapbnd2glob.Append(i);
	}

    // run through surface elements and add elem to related trial dofs
    // Table-creator creates table with one big block of memory,
    // avoids memory fragmentation 
    TableCreator<int> creator;
    Array<DofId> dnumsi;
    for ( ; !creator.Done(); creator++)
      for (int i = 0; i < mesh->GetNSE(); i++)
        {
          trial_space->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator.Add (d, i);
        }
    elems4dof = creator.MoveTable();

    BitArray bnddofs2(test_space->GetNDof());
    bnddofs2.Clear();
    for (int i = 0; i < mesh2->GetNSE(); i++)
      {
	Array<DofId> dnums;
	test_space->GetDofNrs(ElementId(BND, i), dnums);
	for (auto d : dnums)
	  bnddofs2.SetBit(d);
      }
    
    mapglob2bnd2.SetSize(test_space->GetNDof());
    mapglob2bnd2 = -1;
    for (int i = 0; i < test_space->GetNDof(); i++)
      if (bnddofs2.Test(i))
	{
	  mapglob2bnd2[i] = mapbnd2glob2.Size();
	  mapbnd2glob2.Append(i);
	}

    // run through surface elements and add elem to related test dofs
    TableCreator<int> creator2;
    for ( ; !creator2.Done(); creator2++)
      for (int i = 0; i < mesh2->GetNSE(); i++)
        {
          test_space->GetDofNrs( ElementId(BND, i), dnumsi); 
          for (auto d : dnumsi)
            creator2.Add (d, i);
        }
    elems4dof2 = creator2.MoveTable();
  }

  template <typename T>    
  void IntegralOperator<T> :: GetDofNrs(Array<int> &dnums) const
  {
    dnums = mapbnd2glob;
  }

  template <typename T>    
  void IntegralOperator<T> ::  GetDofNrs2(Array<int> &dnums) const   
  {
    dnums = mapbnd2glob2;
  }


  /* compute dense double layer matrix,  dim = ndof_bnd_L2 x ndof_bnd_H1 */
  template <typename T>    
  void IntegralOperator<T> ::
  CalcElementMatrix(FlatMatrix<double> matrix, LocalHeap &lh) const
  {
    if constexpr (is_same<T,double>())
      CalcBlockMatrix(matrix, mapbnd2glob, mapbnd2glob2, lh);    
  }

  template <typename T>    
  void IntegralOperator<T> :: Apply(FlatVector<double> elx, FlatVector<double> ely, 
                                    LocalHeap & lh) const
  {
    static Timer t("ngbem - IntegralOperator::Apply"); RegionTimer reg(t);
    VVector<> vx(hmatrix->Width());
    VVector<> vy(hmatrix->Height());

    vy = 0.0;
    for (int i = 0; i < mapbnd2glob.Size(); i++)
      vx(mapbnd2glob[i]) = elx[i];
    hmatrix -> MultAdd (1, vx, vy);
    for (int i = 0; i < mapbnd2glob2.Size(); i++)
      ely[i] = vy(mapbnd2glob2[i]);
  }



  
  template <typename T>    
  void IntegralOperator<T> :: CalcHMatrix(HMatrix<T> & hmatrix, LocalHeap & clh, struct BEMParameters &param) const
  {
    static Timer t("ngbem - BaseClass::CalcHMatrix"); RegionTimer reg(t);    
    auto & matList = hmatrix.GetMatList();

    ParallelForRange (matList.Size(), [&](IntRange r)
    {
      LocalHeap lh = clh.Split();
      for (int k : r)
        {
          HeapReset hr(lh);
          BEMBlock<T> & block = matList[k];
          auto trialdofs = block.GetTrialDofs();
          auto testdofs = block.GetTestDofs();
          if(block.IsNearField())
            {
              // Compute dense block
              Matrix<T> near(testdofs.Size(), trialdofs.Size());
              CalcBlockMatrix(near, trialdofs, testdofs, lh);	    
              block.SetMat(make_unique<BaseMatrixFromMatrix<T>>(std::move(near)));
            }
          else
            {
              // Compute low-rank block
              try
                {
                  block.SetMat(CalcFarFieldBlock(trialdofs, testdofs, lh));
                }
              catch (netgen::NgException & e)
                {
                  // cout << "not seperated, size = " << testdofs.Size() << " x " << trialdofs.Size() << endl;
                  Matrix<T> near(testdofs.Size(), trialdofs.Size());
                  CalcBlockMatrix(near, trialdofs, testdofs, lh);	    
                  block.SetMat(make_unique<BaseMatrixFromMatrix<T>>(std::move(near)));
                }
            }
        }
    }, TasksPerThread(4));
  }



  template <typename T>
  void StochasticTSVD1 (MatrixView<T> A, MatrixView<T> U, MatrixView<T> V, VectorView<> S)
  {
    for (int i = 0; i < V.Height(); i++)
      for (int j = 0; j < V.Width(); j++)
        V(i,j) = double (rand()) / RAND_MAX;

    Matrix<T> AVt = A * Trans(V);
    Matrix<T> tmp(V.Height(), V.Height());
    LapackSVD (AVt, U, tmp, S, false);    
    Matrix<T> UtA = Trans(U) * A;
    LapackSVD (UtA, tmp, V, S, false);
    U = Matrix<T> (U * tmp);
  }

  void StochasticTSVD1 (MatrixView<Complex> A, MatrixView<Complex> U, MatrixView<Complex> V, VectorView<> S)
  {
    for (int i = 0; i < V.Height(); i++)
      for (int j = 0; j < V.Width(); j++)
        V(i,j) = double (rand()) / RAND_MAX;

    // Matrix<Complex> AVt = A * Conj(Trans(V));
    Matrix<Complex> AVt = A * Trans(V) | Lapack;
    
    Matrix<Complex> tmp(V.Height(), V.Height());
    LapackSVD (AVt, U, tmp, S, false);    
    // Matrix<Complex> UtA = Conj(Trans(U)) * A;
    Matrix<Complex> UtA = Matrix<Complex>(Conj(Trans(U))) * A | Lapack;
    LapackSVD (UtA, tmp, V, S, false);
    U = Matrix<Complex> (U * tmp | Lapack);
  }

  
  template <typename T>
  size_t StochasticTSVD (MatrixView<T> A, MatrixView<T> U, MatrixView<T> V, VectorView<> S, double eps)
  {
    static Timer tsvd("ngbem - StochasticTSVD"); 
    RegionTimer reg(tsvd);
    
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

  template <typename T>    
  unique_ptr<LowRankMatrix<T>> IntegralOperator<T> ::
  CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, LocalHeap &lh) const
  {
    static Timer t("ngbem - IntegralOperator::CalcFarFieldBlock"); RegionTimer reg(t);
    int m = testdofs.Size();
    int n = trialdofs.Size();
    int p = min(n, m);
    
    Matrix<T> A(m, n);
    CalcBlockMatrix(A, trialdofs, testdofs, lh);

    // Calculate SVD for A^\top = V S U^\top
    Matrix<T, ColMajor> V(n, p), Ut(p, m);
    Vector<> S(p);

    int k = StochasticTSVD<T> (A, Trans(Ut), Trans(V), S, param.eps);

    // Low-rank approximation from truncated svd
    
    Matrix<T> U_trc(m, k), Vt_trc(k, n);
    for (size_t j = 0; j < k; j++)
      U_trc.Col(j) = sqrt(S(j)) * Ut.Row(j);

    for (size_t i = 0; i < k; i++)    
      Vt_trc.Row(i) = sqrt(S(i)) * V.Col(i);

    return make_unique<LowRankMatrix<T>> (Matrix<T>(U_trc), Matrix<T>(Vt_trc));
  }

  

  
  
  SingleLayerPotentialOperator ::
  SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, struct BEMParameters _param)
    : IntegralOperator(aspace, aspace, _param)
  {
    /*START: TEST hmatrix: compare approximation with dense matrix. */   

    // create hmatrix
    hmatrix = make_shared<HMatrix<double>>(trial_ct, test_ct,
                                           param.eta, trial_space->GetNDof(), trial_space->GetNDof());
    // compute all its blocks
    LocalHeap lh(10000000);
    CalcHMatrix(*hmatrix, lh, param);
    cout << "HMatrix done " << endl;


    if (param.testhmatrix)
      {
        HeapReset hr(lh);    
        
        Matrix<double> dense(trial_space->GetNDof(), trial_space->GetNDof());
        CalcElementMatrix(dense, lh);
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // Test with vector
        Vector<double> x(trial_space->GetNDof()), y(trial_space->GetNDof());
        x = 1.;
        y = 0.;
        y = dense * x;
        
        S_BaseVectorPtr<> x_base(trial_space->GetNDof(), 1, x.Data());
        S_BaseVectorPtr<> y_base(trial_space->GetNDof(), 1, y.Data());    
        hmatrix->MultAdd(-1., x_base, y_base);
        
        double err = 0.;
        for (int i = 0; i < y.Size(); i++)
          err += y(i) * y(i);
        cout << "error " << sqrt(err) << endl;
      }
    /*END: TEST hmatrix: compare approximation with dense matrix. */   
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

    auto mesh = trial_space->GetMeshAccess();
    auto evaluator = trial_space->GetEvaluator(BND);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv;
    Array<int> testdofsinv;
    trialdofsinv.SetSize(trial_space->GetNDof()); 
    testdofsinv.SetSize(trial_space->GetNDof());

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
              
	  BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(trial_space->GetFE(ei, lh));
	  BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(trial_space->GetFE(ej, lh));
            
	  ElementTransformation &trafoi = mesh->GetTrafo(ei, lh);
	  ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
              
	  Array<DofId> dnumsi, dnumsj;
	  trial_space->GetDofNrs(ei, dnumsi); // mapping to global dof 
	  trial_space->GetDofNrs(ej, dnumsj);
            
	  FlatMatrix<double,ColMajor> mshapei(1, feli.GetNDof(), lh);
	  FlatMatrix<double,ColMajor> mshapej(1, felj.GetNDof(), lh);
            
	  FlatMatrix<> elmat(feli.GetNDof(), felj.GetNDof(), lh); // e.g. 3 x 3
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


  // aspace, space == H1 trialspace, bspace, space2 ==  L2 testspace
  DoubleLayerPotentialOperator ::
  DoubleLayerPotentialOperator (shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace,
                                BEMParameters _param)
    : IntegralOperator(aspace, bspace, _param)
  {
    // create hmatrix
    hmatrix = make_shared<HMatrix<double>>(trial_ct, test_ct, 
                                           param.eta, trial_space->GetNDof(), test_space->GetNDof());

    LocalHeap lh(100000000);
    CalcHMatrix(*hmatrix, lh, param);

    /*START: TEST hmatrix: compare approximation with dense matrix. */
    if (param.testhmatrix)
      {
        Matrix<double> dense(test_space->GetNDof(), trial_space->GetNDof()); // ndof(L2) x ndof(H1)
        CalcElementMatrix(dense, lh);
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // compute all its blocks
        HeapReset hr(lh);    
        
        // Test with vector
        Vector<double> x(trial_space->GetNDof()), y(test_space->GetNDof());
        x = 1.;
        y = 0.;
        y = dense * x;
        
        S_BaseVectorPtr<> x_base(trial_space->GetNDof(), 1, x.Data());
        S_BaseVectorPtr<> y_base(test_space->GetNDof(), 1, y.Data());    
        hmatrix->MultAdd(-1., x_base, y_base);
  
        cout << "error " << L2Norm (y) << endl;
      }
    /*END: TEST hmatrix: compare approximation with dense matrix. */   
  }


  /* compute double layer matrix for given dof sets  - global boundry dof numbers*/
  void DoubleLayerPotentialOperator ::
  CalcBlockMatrix(FlatMatrix<double> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
		  LocalHeap &lh) const 
  {
    auto mesh = trial_space->GetMeshAccess();    // trialspace = H1
    auto mesh2 = test_space->GetMeshAccess();  // testspace = L2
    
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

    auto evaluator = trial_space->GetEvaluator(BND);
    auto evaluator2 = test_space->GetEvaluator(BND);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv(trial_space->GetNDof()); 
    Array<int> testdofsinv(test_space->GetNDof());

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

    RegionTimer regloops(tloops);    
    for (int i = 0; i < patchi.Size(); i++) // test
      for (int j = 0; j < patchj.Size(); j++) // trial
	{
	  HeapReset hr(lh);
	  ElementId ei(BND, patchi[i]);
	  ElementId ej(BND, patchj[j]);
          
	  auto verti = mesh2->GetElement(ei).Vertices();
	  auto vertj = mesh->GetElement(ej).Vertices();          
            
	  BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(test_space->GetFE(ei, lh));
	  BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(trial_space->GetFE(ej, lh));
              
	  ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
	  ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
              
	  Array<DofId> dnumsi, dnumsj;
	  test_space->GetDofNrs(ei, dnumsi); // mapping to global dof
	  trial_space->GetDofNrs(ej, dnumsj);
        
	  FlatVector<> shapei(feli.GetNDof(), lh);
	  FlatVector<> shapej(felj.GetNDof(), lh);
          
	  FlatMatrix<> elmat(feli.GetNDof(), felj.GetNDof(), lh); 
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
                
		evaluator2 -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
		evaluator-> CalcMatrix(felj, miry, Trans(shapesj), lh);


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
                        
			double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
			kernel_ixiy(ix, iy) = fac*kernel;
		      }
		  }
          
		FlatMatrix<double> kernel_shapesj(irtrig.Size(), felj.GetNDof(), lh);
		kernel_shapesj = kernel_ixiy * Trans(shapesj);
		elmat += shapesi * kernel_shapesj;
                
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




  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          KERNEL _kernel,
                          BEMParameters _param)
    : IntegralOperator<value_type>(_trial_space, _test_space, _param), kernel(_kernel)
  {
    hmatrix =
      make_shared<HMatrix<value_type>>(trial_ct, test_ct, 
                                       param.eta, trial_space->GetNDof(), test_space->GetNDof());
    
    LocalHeap lh(100000000);
    this->CalcHMatrix(*hmatrix, lh, param);
  }
                                  
  template <typename KERNEL>
  void GenericIntegralOperator<KERNEL> ::
  CalcBlockMatrix(FlatMatrix<value_type> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
                  LocalHeap &lh) const
  {
    auto mesh = this->trial_space->GetMeshAccess();  
    auto mesh2 = this->test_space->GetMeshAccess();  
    
    static Timer tall("ngbem DLP - all");
    static Timer tloops("ngbem DLP - loops");    
    static Timer t_identic("ngbem DLP - identic panel");
    static Timer t_common_vertex("ngbem DLP - common vertex");        
    static Timer t_common_edge("ngbem DLP - common edge");        
    static Timer t_disjoint("ngbem DLP - disjoint");

    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder);
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(param.intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(param.intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(param.intorder);

    matrix = 0; 

    auto evaluator = trial_space->GetEvaluator(BND);
    auto evaluator2 = test_space->GetEvaluator(BND);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv(trial_space->GetNDof()); 
    Array<int> testdofsinv(test_space->GetNDof());

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

    RegionTimer regloops(tloops);    
    for (int i = 0; i < patchi.Size(); i++) // test
      for (int j = 0; j < patchj.Size(); j++) // trial
	{
	  HeapReset hr(lh);
	  ElementId ei(BND, patchi[i]);
	  ElementId ej(BND, patchj[j]);
          
	  auto verti = mesh2->GetElement(ei).Vertices();
	  auto vertj = mesh->GetElement(ej).Vertices();          
            
	  BaseScalarFiniteElement &feli = dynamic_cast<BaseScalarFiniteElement &>(test_space->GetFE(ei, lh));
	  BaseScalarFiniteElement &felj = dynamic_cast<BaseScalarFiniteElement &>(trial_space->GetFE(ej, lh));
              
	  ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
	  ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
              
	  Array<DofId> dnumsi, dnumsj;
	  test_space->GetDofNrs(ei, dnumsi); // mapping to global dof
	  trial_space->GetDofNrs(ej, dnumsj);
        
	  FlatVector<> shapei(feli.GetNDof(), lh);
	  FlatVector<> shapej(felj.GetNDof(), lh);
          
	  FlatMatrix<value_type> elmat(feli.GetNDof(), felj.GetNDof(), lh); 
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
		    Vec<3> ny = mipy.GetNV();                    
		    double nxy = InnerProduct(ny, (x-y));
		    double normxy = L2Norm(x-y);
		    // double kernel = nxy / (4*M_PI*normxy*normxy*normxy);
                    value_type kernel_ = kernel.Evaluate(x, y, nx, ny)(0,0);
                        
		    evaluator2->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
		    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                        
		    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_panel_weight[k];
		    elmat += fac*kernel_* shapei * Trans(shapej);
		  }

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
		    Vec<3> ny = mipy.GetNV();
		    // double nxy = InnerProduct(ny, (x-y));
		    // double normxy = L2Norm(x-y);
		    value_type kernel_ = kernel.Evaluate(x,y,nx,ny)(0,0);
                        
		    evaluator2->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
		    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                        
		    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weight[k];
		    elmat += fac*kernel_* shapei * Trans(shapej);
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
		    Vec<3> nx = mipy.GetNV();
		    Vec<3> ny = mipy.GetNV();
                    
		    value_type kernel_ = kernel.Evaluate(x, y, nx, ny)(0,0);
          
		    evaluator2->CalcMatrix(feli, mipx, Trans(shapei.AsMatrix(feli.GetNDof(),1)), lh);
		    evaluator->CalcMatrix(felj, mipy, Trans(shapej.AsMatrix(felj.GetNDof(),1)), lh);
                  
		    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weight[k];
		    elmat += fac*kernel_* shapei * Trans(shapej);
		  }

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
                
		evaluator2 -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
		evaluator-> CalcMatrix(felj, miry, Trans(shapesj), lh);


		FlatMatrix<value_type> kernel_ixiy(irtrig.Size(), irtrig.Size(), lh);
		for (int ix = 0; ix < irtrig.Size(); ix++)
		  {
		    for (int iy = 0; iy < irtrig.Size(); iy++)
		      {
			Vec<3> x = mirx[ix].GetPoint();
			Vec<3> y = miry[iy].GetPoint();

                        Vec<3> nx = miry[iy].GetNV();
                        Vec<3> ny = miry[iy].GetNV();
                        value_type kernel_ = kernel.Evaluate(x, y, nx, ny)(0,0);
                        
			double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
			kernel_ixiy(ix, iy) = fac*kernel_;
		      }
		  }
          
		FlatMatrix<value_type> kernel_shapesj(irtrig.Size(), felj.GetNDof(), lh);
		kernel_shapesj = kernel_ixiy * Trans(shapesj);
		elmat += shapesi * kernel_shapesj;
                
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




  template <typename KERNEL>
  unique_ptr<LowRankMatrix<typename KERNEL::value_type>> GenericIntegralOperator<KERNEL> ::    
  CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                    LocalHeap &lh) const 
  
  {
    auto mesh = this->trial_space->GetMeshAccess();  
    auto mesh2 = this->test_space->GetMeshAccess();  
    
    static Timer tall("ngbem generic FarFieldBlock");
    static Timer tkernel("ngbem generic FarFieldBlock - kernel");    
    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder);
    
    auto evaluator = trial_space->GetEvaluator(BND);
    auto evaluator2 = test_space->GetEvaluator(BND);

    Array<int> tmp, tmp2;
    Array<int> patchi, patchj;
    Array<int> trialdofsinv(trial_space->GetNDof()); 
    Array<int> testdofsinv(test_space->GetNDof());

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

    Matrix<value_type> matrix(testdofs.Size(), trialdofs.Size());
    matrix = value_type(0.0);


    // new code
    Array<Vec<3>> xi, yj, nxi, nyj;  // i..test, j... trial
    Array<double> wxi, wyj;
    
    BitArray test_vertices(mesh->GetNV());
    test_vertices.Clear();

    // test patches
    for (int i = 0; i < patchi.Size(); i++)
      {
        HeapReset hr(lh);
        ElementId ei(BND, patchi[i]);
          
        for (auto v : mesh2->GetElement(ei).Vertices())
          test_vertices.SetBit(v);
        
        ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
        MappedIntegrationRule<2,3> mirx(irtrig, trafoi, lh);
        
        for (int ix = 0; ix < irtrig.Size(); ix++)
          {
            xi.Append (mirx[ix].GetPoint());
            nxi.Append (mirx[ix].GetNV());
            wxi.Append (mirx[ix].GetWeight());
          }
      }

    // trial patches
    for (int j = 0; j < patchj.Size(); j++)
      {
        HeapReset hr(lh);
        ElementId ej(BND, patchj[j]);

        if (mesh == mesh2)
          for (auto v : mesh->GetElement(ej).Vertices())
            if (test_vertices.Test(v))
              throw Exception("far field block must not have common vertices");              
        
        FiniteElement &felj = trial_space->GetFE(ej, lh);
        ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);

        MappedIntegrationRule<2,3> miry(irtrig, trafoj, lh);
            
        for (int iy = 0; iy < irtrig.Size(); iy++)
          {
            yj.Append (miry[iy].GetPoint());
            nyj.Append (miry[iy].GetNV());
            wyj.Append (miry[iy].GetWeight());
          }
      }

    tkernel.Start();
    /*
    Matrix<value_type> kernel_matrix(xi.Size(), yj.Size());
    for (int i = 0; i < xi.Size(); i++)
      for (int j = 0; j < yj.Size(); j++)
        kernel_matrix(i,j) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j])(0,0);

    size_t p = 200;
    Matrix<value_type> Umax(kernel_matrix.Height(), p), Vmax(p, kernel_matrix.Width());
    Vector<> S(p);

    int k = StochasticTSVD<value_type> (kernel_matrix, Umax, Vmax, S, param.eps);


    for (size_t j = 0; j < k; j++)
      {
        Umax.Col(j) *= sqrt(S(j));
        Vmax.Row(j) *= sqrt(S(j));
      }
    
    */

    size_t p = 300;
    int rank = p;
    auto GetRow = [&](int i, SliceVector<value_type> row)
    {
      for (int j = 0; j < yj.Size(); j++)
        row(j) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j])(0,0);
    };
    auto GetCol = [&](int j, SliceVector<value_type> col)
    {
      for (int i = 0; i < xi.Size(); i++)
        col(i) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j])(0,0);
    };

    Matrix<value_type> Umax(xi.Size(), p), Vmax(p, yj.Size());

    // for quasi random sequence of pivot indices
    size_t primes[] = { 71, 73, 79, 83, 89, 97 };
    size_t prime;
    for (auto tp : primes)
      {
        if (xi.Size()%tp != 0)
          {
            prime = tp;
            break;
          }
      }
    // ACA compression 
    for (size_t k = 0; k < p; k++)
      {
        // int ik = k;  // what else ?
        size_t ik = (k*prime)%xi.Size(); 
        
        GetRow(ik, Vmax.Row(k));
        Vmax.Row(k) -= Trans(Vmax.Rows(0,k)) * Umax.Row(ik).Range(0,k);
         
        double err = L2Norm(Vmax.Row(k));
        // cout << "Norm vk = " << err << endl;
        if (err < param.eps)
          {
            rank = k;
            break;
          }
        
        int jmax = 0;
        for (int j = 0; j < Vmax.Width(); j++)
          if (fabs (Vmax(k,j)) > fabs(Vmax(k,jmax)))
            jmax = j;
        Vmax.Row(k) *= 1.0 / Vmax(k,jmax);

        GetCol(jmax, Umax.Col(k));
        Umax.Col(k) -= Umax.Cols(0,k) * Vmax.Col(jmax).Range(0,k);
      }
    // *testout << "rank = " << rank << endl;
    int k = rank;
    tkernel.Stop();
    
    auto U = Umax.Cols(0,k);
    auto V = Vmax.Rows(0,k);
    // cout << "k = " << k << ", err = " << L2Norm(help-U*V) << endl;
    
    for (int i = 0; i < U.Height(); i++)
      U.Row(i) *= wxi[i];
    for (int j = 0; j < V.Width(); j++)
      V.Col(j) *= wyj[j];

    Matrix<value_type> U2(testdofs.Size(), k);
    Matrix<value_type> V2(k, trialdofs.Size());
    U2 = value_type(0.0);
    V2 = value_type(0.0);
    
    int cnt = 0;
    for (int i = 0; i < patchi.Size(); i++) // test
      {
        HeapReset hr(lh);
        ElementId ei(BND, patchi[i]);
            
        FiniteElement &feli = test_space->GetFE(ei, lh);
        ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
              
        Array<DofId> dnumsi;
        test_space->GetDofNrs(ei, dnumsi);
          
        MappedIntegrationRule<2,3> mirx(irtrig, trafoi, lh);
        FlatMatrix<> shapesi(feli.GetNDof(), irtrig.Size(), lh);
            
        evaluator2 -> CalcMatrix(feli, mirx, Trans(shapesi), lh);

        Matrix<value_type> tmp = shapesi * U.Rows(cnt, cnt+irtrig.Size());
        cnt += irtrig.Size();
        for (int ii = 0; ii < dnumsi.Size(); ii++) // test
          if (testdofsinv[dnumsi[ii]] != -1)
            U2.Row(testdofsinv[dnumsi[ii]]) += tmp.Row(ii);
      }

    cnt = 0;
    for (int j = 0; j < patchj.Size(); j++) // test
      {
        HeapReset hr(lh);
        ElementId ej(BND, patchj[j]);
            
        FiniteElement &felj = trial_space->GetFE(ej, lh);
        ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
        
        Array<DofId> dnumsj;
        trial_space->GetDofNrs(ej, dnumsj);
        
        MappedIntegrationRule<2,3> miry(irtrig, trafoj, lh);
        FlatMatrix<> shapesj(felj.GetNDof(), irtrig.Size(), lh);
            
        evaluator -> CalcMatrix(felj, miry, Trans(shapesj), lh);

        Matrix<value_type> tmp = V.Cols(cnt, cnt+irtrig.Size()) * Trans(shapesj);
        cnt += irtrig.Size();
                     
        for (int jj = 0; jj < dnumsj.Size(); jj++) // trial
          if(trialdofsinv[dnumsj[jj]] != -1)
            V2.Col(trialdofsinv[dnumsj[jj]]) += tmp.Col(jj);
      }

    return make_unique<LowRankMatrix<value_type>> (std::move(U2), std::move(V2));
  }

  
  template class GenericIntegralOperator<LaplaceSLKernel<3>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3>>;
  
  template class GenericIntegralOperator<HelmholtzSLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3>>;    
}

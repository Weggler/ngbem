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

  
  


  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          KERNEL _kernel,
                          BEMParameters _param)
    : IntegralOperator<value_type>(_trial_space, _test_space, _param), kernel(_kernel)
  {
    trial_evaluator = trial_space -> GetEvaluator(BND);
    test_evaluator = test_space -> GetEvaluator(BND);
    hmatrix =
      make_shared<HMatrix<value_type>>(trial_ct, test_ct, 
                                       param.eta, trial_space->GetNDof(), test_space->GetNDof());
    
    LocalHeap lh(100000000);
    this->CalcHMatrix(*hmatrix, lh, param);
	
    if (param.testhmatrix)
      {
        Matrix<value_type> dense(test_space->GetNDof(), trial_space->GetNDof());
        CalcBlockMatrix(dense, mapbnd2glob, mapbnd2glob2, lh);            
        cout << "dense: " << dense.Height() << " x " << dense.Width() << endl;
        
        // compute all its blocks
        HeapReset hr(lh);    
        
        // Test with vector
        Vector<value_type> x(trial_space->GetNDof()), y(test_space->GetNDof());
        x = 1.;
        y = dense * x;
        
        VFlatVector<value_type> x_base(x);
        VFlatVector<value_type> y_base(y);
        y_base -= (*hmatrix) * x_base;

	cout << "error " << L2Norm (y) << endl;
      }
  }


  template <typename KERNEL>
  GenericIntegralOperator<KERNEL> ::
  GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                          shared_ptr<DifferentialOperator> _trial_evaluator, 
                          shared_ptr<DifferentialOperator> _test_evaluator, 
                          KERNEL _kernel,
                          BEMParameters _param)
    : IntegralOperator<value_type>(_trial_space, _test_space, _param), kernel(_kernel),
      trial_evaluator(_trial_evaluator), test_evaluator(_test_evaluator)
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
    
    static Timer tall("ngbem - all " + KERNEL::Name());
    static Timer tloops("ngbem - loops " + KERNEL::Name());
    static Timer t_identic("ngbem identic panel " + KERNEL::Name());
    static Timer t_common_vertex("ngbem common vertex " + KERNEL::Name());

    static Timer t_common_edge("ngbem common edge " + KERNEL::Name());
    static Timer t_disjoint("ngbem disjoint " + KERNEL::Name());

    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder);
    
    auto [ identic_panel_x, identic_panel_y, identic_panel_weight ] =
      IdenticPanelIntegrationRule(param.intorder);

    auto [ common_vertex_x, common_vertex_y, common_vertex_weight ] =
      CommonVertexIntegrationRule(param.intorder);
    
    auto [ common_edge_x, common_edge_y, common_edge_weight ] =
      CommonEdgeIntegrationRule(param.intorder);

    matrix = 0; 

    // auto evaluator = trial_space->GetEvaluator(BND);
    // auto evaluator2 = test_space->GetEvaluator(BND);

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
            
	  FiniteElement &feli = test_space->GetFE(ei, lh);
	  FiniteElement &felj = trial_space->GetFE(ej, lh);
              
	  ElementTransformation &trafoi = mesh2->GetTrafo(ei, lh);
	  ElementTransformation &trafoj = mesh->GetTrafo(ej, lh);
              
	  Array<DofId> dnumsi, dnumsj;
	  test_space->GetDofNrs(ei, dnumsi); // mapping to global dof
	  trial_space->GetDofNrs(ej, dnumsj);
        
	  FlatMatrix<> shapei(feli.GetNDof(), test_evaluator->Dim(), lh);
	  FlatMatrix<> shapej(felj.GetNDof(), trial_evaluator->Dim(), lh);
          
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
                /*
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

                    value_type kernel_ = kernel.Evaluate(x, y, nx, ny)(0,0);
                        
		    test_evaluator->CalcMatrix(feli, mipx, Trans(shapei), lh);
		    trial_evaluator->CalcMatrix(felj, mipy, Trans(shapej), lh);
                        
		    double fac = mipx.GetMeasure()*mipy.GetMeasure()*identic_panel_weight[k];
		    elmat += fac*kernel_* shapei * Trans(shapej);
		  }
                */

                // vectorized version:
                constexpr int BS = 128;
                for (int k = 0; k < identic_panel_weight.Size(); k+=BS)
                  {
                    int num = std::min(size_t(BS), identic_panel_weight.Size()-k);
                    
                    HeapReset hr(lh);
                    
                    IntegrationRule irx(num, lh);
                    IntegrationRule iry(num, lh);

                    for (int k2 = 0; k2 < num; k2++)
                      {
                        Vec<2> xk = identic_panel_x[k+k2];
                        Vec<2> yk = identic_panel_y[k+k2];
                        
                        irx[k2] = IntegrationPoint(xk(0), xk(1), 0,
                                                   identic_panel_weight[k+k2]);
                        iry[k2] = IntegrationPoint(yk(0), yk(1), 0, 0);
                      }

                    SIMD_IntegrationRule simd_irx(irx);
                    SIMD_IntegrationRule simd_iry(iry);
                    
                    SIMD_MappedIntegrationRule<2,3> mirx(simd_irx, trafoi, lh);
                    SIMD_MappedIntegrationRule<2,3> miry(simd_iry, trafoj, lh);


                    FlatMatrix<SIMD<double>> mshapesi(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);
                    FlatMatrix<SIMD<value_type>> mshapesi_kern(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);                    
                    FlatMatrix<SIMD<double>> mshapesj(felj.GetNDof()*trial_evaluator->Dim(), miry.Size(), lh);

                    test_evaluator->CalcMatrix(feli, mirx, mshapesi);
                    trial_evaluator->CalcMatrix(felj, miry, mshapesj);

                    for (int k2 = 0; k2 < mirx.Size(); k2++)
                      {
                        Vec<3,SIMD<double>> x = mirx[k2].Point();
                        Vec<3,SIMD<double>> y = miry[k2].Point();
                        Vec<3,SIMD<double>> nx = mirx[k2].GetNV();
                        Vec<3,SIMD<double>> ny = miry[k2].GetNV();
                    
                        SIMD<value_type> kernel_ = kernel.Evaluate(x, y, nx, ny);                        
                        auto fac = mirx[k2].GetMeasure()*miry[k2].GetMeasure()*simd_irx[k2].Weight(); 
                        mshapesi_kern.Col(k2) = fac*kernel_ * mshapesi.Col(k2);
                      }

                    AddABt (mshapesi_kern.Reshape(feli.GetNDof(), test_evaluator->Dim()*mirx.Size()),
                            mshapesj.Reshape(felj.GetNDof(), trial_evaluator->Dim()*miry.Size()),
                            elmat);
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
                /*
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
                    
		    value_type kernel_ = kernel.Evaluate(x,y,nx,ny)(0,0);
                        
		    test_evaluator->CalcMatrix(feli, mipx, Trans(shapei), lh);
		    trial_evaluator->CalcMatrix(felj, mipy, Trans(shapej), lh);
                        
		    double fac = mipx.GetMeasure()*mipy.GetMeasure() * common_edge_weight[k];
		    elmat += fac*kernel_* shapei * Trans(shapej);
		  }
                */

                // vectorized version:
                constexpr int BS = 128;
                for (int k = 0; k < common_edge_weight.Size(); k+=BS)
                  {
                    int num = std::min(size_t(BS), common_edge_weight.Size()-k);
                    
                    HeapReset hr(lh);
                    
                    IntegrationRule irx(num, lh);
                    IntegrationRule iry(num, lh);

                    for (int k2 = 0; k2 < num; k2++)
                      {
                        Vec<2> xk = common_edge_x[k+k2];
                        Vec<2> yk = common_edge_y[k+k2];

                        Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                        Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );

                        Vec<3> plamx, plamy;
                        for (int i = 0; i < 3; i++)
                          {
                            plamx(vpermx[i]) = lamx(i);
                            plamy(vpermy[i]) = lamy(i);
                          }
                        
                        irx[k2] = IntegrationPoint(plamx(0), plamx(1), 0, common_edge_weight[k+k2]);
                        iry[k2] = IntegrationPoint(plamy(0), plamy(1), 0, 0);
                      }

                    SIMD_IntegrationRule simd_irx(irx);
                    SIMD_IntegrationRule simd_iry(iry);
                    
                    SIMD_MappedIntegrationRule<2,3> mirx(simd_irx, trafoi, lh);
                    SIMD_MappedIntegrationRule<2,3> miry(simd_iry, trafoj, lh);


                    FlatMatrix<SIMD<double>> mshapesi(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);
                    FlatMatrix<SIMD<value_type>> mshapesi_kern(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);                    
                    FlatMatrix<SIMD<double>> mshapesj(felj.GetNDof()*trial_evaluator->Dim(), miry.Size(), lh);

                    test_evaluator->CalcMatrix(feli, mirx, mshapesi);
                    trial_evaluator->CalcMatrix(felj, miry, mshapesj);

                    for (int k2 = 0; k2 < mirx.Size(); k2++)
                      {
                        Vec<3,SIMD<double>> x = mirx[k2].Point();
                        Vec<3,SIMD<double>> y = miry[k2].Point();
                        Vec<3,SIMD<double>> nx = mirx[k2].GetNV();
                        Vec<3,SIMD<double>> ny = miry[k2].GetNV();
                    
                        SIMD<value_type> kernel_ = kernel.Evaluate(x, y, nx, ny);                        
                        auto fac = mirx[k2].GetMeasure()*miry[k2].GetMeasure()*simd_irx[k2].Weight(); 
                        mshapesi_kern.Col(k2) = fac*kernel_ * mshapesi.Col(k2);
                      }

                    AddABt (mshapesi_kern.Reshape(feli.GetNDof(), test_evaluator->Dim()*mirx.Size()),
                            mshapesj.Reshape(felj.GetNDof(), trial_evaluator->Dim()*miry.Size()),
                            elmat);
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
                /*
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
		    Vec<3> ny = mipy.GetNV();
                    
		    value_type kernel_ = kernel.Evaluate(x, y, nx, ny)(0,0);
          
		    test_evaluator->CalcMatrix(feli, mipx, Trans(shapei), lh);
		    trial_evaluator->CalcMatrix(felj, mipy, Trans(shapej), lh);
                  
		    double fac = mipx.GetMeasure()*mipy.GetMeasure()*common_vertex_weight[k];
		    elmat += fac*kernel_* shapei * Trans(shapej);
		  }
                */

                // vectorized version:
                constexpr int BS = 128;
                for (int k = 0; k < common_vertex_weight.Size(); k+=BS)
                  {
                    int num = std::min(size_t(BS), common_vertex_weight.Size()-k);
                    
                    HeapReset hr(lh);
                    
                    IntegrationRule irx(num, lh);
                    IntegrationRule iry(num, lh);

                    for (int k2 = 0; k2 < num; k2++)
                      {
                        Vec<2> xk = common_vertex_x[k+k2];
                        Vec<2> yk = common_vertex_y[k+k2];

                        Vec<3> lamx (1-xk(0)-xk(1), xk(0), xk(1) );
                        Vec<3> lamy (1-yk(0)-yk(1), yk(0), yk(1) );

                        Vec<3> plamx, plamy;
                        for (int i = 0; i < 3; i++)
                          {
                            plamx(vpermx[i]) = lamx(i);
                            plamy(vpermy[i]) = lamy(i);
                          }
                        
                        irx[k2] = IntegrationPoint(plamx(0), plamx(1), 0, common_vertex_weight[k+k2]);
                        iry[k2] = IntegrationPoint(plamy(0), plamy(1), 0, 0);
                      }

                    SIMD_IntegrationRule simd_irx(irx);
                    SIMD_IntegrationRule simd_iry(iry);

                    SIMD_MappedIntegrationRule<2,3> mirx(simd_irx, trafoi, lh);
                    SIMD_MappedIntegrationRule<2,3> miry(simd_iry, trafoj, lh);

                    FlatMatrix<SIMD<double>> mshapesi(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);
                    FlatMatrix<SIMD<value_type>> mshapesi_kern(feli.GetNDof()*test_evaluator->Dim(), mirx.Size(), lh);                    
                    FlatMatrix<SIMD<double>> mshapesj(felj.GetNDof()*trial_evaluator->Dim(), miry.Size(), lh);

                    test_evaluator->CalcMatrix(feli, mirx, mshapesi);
                    trial_evaluator->CalcMatrix(felj, miry, mshapesj);

                    for (int k2 = 0; k2 < mirx.Size(); k2++)
                      {
                        Vec<3,SIMD<double>> x = mirx[k2].Point();
                        Vec<3,SIMD<double>> y = miry[k2].Point();
                        Vec<3,SIMD<double>> nx = mirx[k2].GetNV();
                        Vec<3,SIMD<double>> ny = miry[k2].GetNV();
                    
                        SIMD<value_type> kernel_ = kernel.Evaluate(x, y, nx, ny);
                        auto fac = mirx[k2].GetMeasure()*miry[k2].GetMeasure()*simd_irx[k2].Weight(); 
                        mshapesi_kern.Col(k2) = fac*kernel_ * mshapesi.Col(k2);
                      }

                    AddABt (mshapesi_kern.Reshape(feli.GetNDof(), test_evaluator->Dim()*mirx.Size()),
                            mshapesj.Reshape(felj.GetNDof(), trial_evaluator->Dim()*miry.Size()),
                            elmat);
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

		FlatMatrix<> shapesi(feli.GetNDof(), test_evaluator->Dim()*irtrig.Size(), lh);
		FlatMatrix<> shapesj(felj.GetNDof(), trial_evaluator->Dim()*irtrig.Size(), lh);
                
		test_evaluator -> CalcMatrix(feli, mirx, Trans(shapesi), lh);
		trial_evaluator-> CalcMatrix(felj, miry, Trans(shapesj), lh);


		FlatMatrix<value_type> kernel_ixiy(irtrig.Size(), irtrig.Size(), lh);
		for (int ix = 0; ix < irtrig.Size(); ix++)
		  {
		    for (int iy = 0; iy < irtrig.Size(); iy++)
		      {
			Vec<3> x = mirx[ix].GetPoint();
			Vec<3> y = miry[iy].GetPoint();

                        Vec<3> nx = miry[iy].GetNV();
                        Vec<3> ny = miry[iy].GetNV();
                        value_type kernel_ = kernel.Evaluate(x, y, nx, ny);
                        
			double fac = mirx[ix].GetWeight()*miry[iy].GetWeight();
			kernel_ixiy(ix, iy) = fac*kernel_;
		      }
		  }


		FlatMatrix<value_type> kernel_shapesj(irtrig.Size(), felj.GetNDof(), lh);
                if (test_evaluator->Dim()==1)
                  {
                    kernel_shapesj = kernel_ixiy * Trans(shapesj);
                    elmat += shapesi * kernel_shapesj;
                  }
                else
                  {
                    FlatMatrix<> shapesi1(feli.GetNDof(), irtrig.Size(), lh);
                    FlatMatrix<> shapesj1(felj.GetNDof(), irtrig.Size(), lh);

                    for (int k = 0; k < test_evaluator->Dim(); k++)
                      {
                        for (int j = 0; j < irtrig.Size(); j++)
                          {
                            shapesi1.Col(j) = shapesi.Col(test_evaluator->Dim()*j+k);
                            shapesj1.Col(j) = shapesj.Col(test_evaluator->Dim()*j+k);
                          }
                        kernel_shapesj = kernel_ixiy * Trans(shapesj1);
                        elmat += shapesi1 * kernel_shapesj;
                      }
                  }
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
    // if (trial_evaluator->Dim() > 1)
    // throw Exception("ACA not supported for vectorial evaluators");              
    
    auto mesh = this->trial_space->GetMeshAccess();  
    auto mesh2 = this->test_space->GetMeshAccess();  
    
    static Timer tall("ngbem FarFieldBlock " + KERNEL::Name());
    static Timer tACA("ngbem FarFieldBlock - ACA " + KERNEL::Name());
    static Timer tkernel("ngbem FarFieldBlock - kernel " + KERNEL::Name());      
    RegionTimer reg(tall);

    IntegrationRule irtrig(ET_TRIG, param.intorder);
    SIMD_IntegrationRule simd_irtrig(irtrig);

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

    tACA.Start();
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

    size_t p = min(xi.Size(), yj.Size());
    //int rank = p;
    auto GetRow = [&](int i, SliceVector<value_type> row)
    {
      RegionTimer reg(tkernel);
      tkernel.AddFlops (yj.Size());
      for (int j = 0; j < yj.Size(); j++)
        row(j) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j]);
    };
    auto GetCol = [&](int j, SliceVector<value_type> col)
    {
      RegionTimer reg(tkernel);
      tkernel.AddFlops (xi.Size());
      for (int i = 0; i < xi.Size(); i++)
        col(i) = kernel.Evaluate(xi[i], yj[j], nxi[i], nyj[j]);
    };

    Matrix<value_type> Umax(xi.Size(), p);
    Matrix<value_type> Vmax(p, yj.Size());

    size_t ik = 0, jk = 0, ikm1 = 0, jkm1 = yj.Size() + 1, rank = p;
    
    // // for quasi random sequence of pivot indices
    // size_t primes[] = { 71, 73, 79, 83, 89, 97 };
    // size_t prime;
    // for (auto tp : primes)
    //   {
    //     if (xi.Size()%tp != 0)
    //       {
    //         prime = tp;
    //         break;
    //       }
    //   }
    // // ACA compression 
    // for (size_t k = 0; k < p; k++)
    //   {
    //     // int ik = k;  // what else ?
    //     size_t ik = (k*prime)%xi.Size(); 
        
    //     GetRow(ik, Vmax.Row(k));
    //     Vmax.Row(k) -= Trans(Vmax.Rows(0,k)) * Umax.Row(ik).Range(0,k);
         
    //     double err = L2Norm(Vmax.Row(k));
    //     // cout << "Norm vk = " << err << endl;
    //     if (err < param.eps)
    //       {
    //         rank = k;
    //         break;
    //       }
        
    //     int jmax = 0;
    //     for (int j = 0; j < Vmax.Width(); j++)
    //       if (fabs (Vmax(k,j)) > fabs(Vmax(k,jmax)))
    //         jmax = j;
    //     Vmax.Row(k) *= 1.0 / Vmax(k,jmax);

    //     GetCol(jmax, Umax.Col(k));
    //     Umax.Col(k) -= Umax.Cols(0,k) * Vmax.Col(jmax).Range(0,k);
    //   }
    
    // Scale eps appropriately (see Bebendorf, Hierarchical Matrices  p. 126 & 135
    double eps = 2. / 3. * param.eps / sqrt(xi.Size() * yj.Size());
    // The Frobenius norm squared of the approximant U * V^H
    double norm2 = 0.;

    for (size_t k = 0; k < p; k++) {
	// Get the ik-th row
	GetRow(ik, Vmax.Row(k));
        Vmax.Row(k) -= Trans(Vmax.Rows(0,k)) * Umax.Row(ik).Range(0,k);

	// Find the new column pivot position jk in the new row
	double vkj = 0.;
        for (int j = 0; j < Vmax.Width(); j++)
          if (fabs (Vmax(k, j)) > vkj && j != jkm1) {
	      vkj = fabs (Vmax(k, j));
	      jk = j;
	    }

	// If the pivot element is close to zero, exit
	if (vkj == 0.) {
	  rank = k;
	  break;
	}
	
	// Scale with inverse of the pivot entry at (ik, jk)
        Vmax.Row(k) *= 1.0 / Vmax(k, jk);

	// Get the jk-th column
	GetCol(jk, Umax.Col(k));
        Umax.Col(k) -= Umax.Cols(0,k) * Vmax.Col(jk).Range(0,k);
	
	// Find the new row pivot position ik in the new column
	double uik = 0.;
	for (int i = 0; i < Umax.Height(); i++)
          if (fabs (Umax(i, k)) > uik && i != ikm1) {
	      uik = fabs (Umax(i, k));
	      ik = i;
	    }

	// Update the Frobenius norm
	double norm_k = L2Norm(Vmax.Row(k)) * L2Norm(Umax.Col(k));
	norm2 += norm_k * norm_k;
	for (int l = 0; l < k; l++)
	  norm2 += 2. * std::real(InnerProduct(Vmax.Row(k), Vmax.Row(l)) *
				  InnerProduct(Umax.Col(k), Umax.Col(l)));

	// New pivots become old pivots
	ikm1 = ik;
	jkm1 = jk;

	// Stop if the new update is relatively small, see Bebendorf pp. 141-142
	if (norm_k < eps * sqrt(norm2)) {
	    rank = k + 1;
	    break;
	  }
      }

    // *testout << "rank = " << rank << endl;
    size_t k = rank;
    tACA.Stop();
    
    auto U = Umax.Cols(0,k);
    auto V = Vmax.Rows(0,k);
    // cout << "k = " << k << ", err = " << L2Norm(help-U*V) << endl;

    for (int i = 0; i < U.Height(); i++)
      U.Row(i) *= wxi[i];
    for (int j = 0; j < V.Width(); j++)
      V.Col(j) *= wyj[j];

    Matrix<value_type> U2(testdofs.Size(), k*test_evaluator->Dim());
    Matrix<value_type> V2(k*test_evaluator->Dim(), trialdofs.Size());
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

        int dim = test_evaluator->Dim();
        SIMD_MappedIntegrationRule<2,3> mirx(simd_irtrig, trafoi, lh);
        FlatMatrix<SIMD<double>> shapesi(dim*feli.GetNDof(), simd_irtrig.Size(), lh);
        SliceMatrix<double> dshapesi(dim*feli.GetNDof(), simd_irtrig.GetNIP(), simd_irtrig.Size()*SIMD<double>::Size(),
                                     (double*)shapesi.Data());
        
        test_evaluator -> CalcMatrix(feli, mirx, shapesi);
        Matrix<value_type> tmp1 = dshapesi * U.Rows(cnt, cnt+irtrig.Size());        
        auto tmp = tmp1.Reshape(feli.GetNDof(), dim*U.Width());
                                
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


        int dim = trial_evaluator->Dim();        
        SIMD_MappedIntegrationRule<2,3> miry(simd_irtrig, trafoj, lh);
        FlatMatrix<SIMD<double>> shapesj(dim*felj.GetNDof(), simd_irtrig.Size(), lh);
        SliceMatrix<double> dshapesj(dim*felj.GetNDof(), simd_irtrig.GetNIP(), simd_irtrig.Size()*SIMD<double>::Size(),
                                     (double*)shapesj.Data());
        
        trial_evaluator -> CalcMatrix(felj, miry, shapesj);
        Matrix<value_type> tmp1 = dshapesj * Trans(V).Rows(cnt, cnt+irtrig.Size());        
        auto tmp = tmp1.Reshape(felj.GetNDof(), dim*V.Height());
        
        cnt += irtrig.Size();
        for (int jj = 0; jj < dnumsj.Size(); jj++) // trial
          if(trialdofsinv[dnumsj[jj]] != -1)
            V2.Col(trialdofsinv[dnumsj[jj]]) += tmp.Row(jj);
      }

    return make_unique<LowRankMatrix<value_type>> (std::move(U2), std::move(V2));
  }

  
  template class GenericIntegralOperator<LaplaceSLKernel<3>>;
  template class GenericIntegralOperator<LaplaceDLKernel<3>>;
  
  template class GenericIntegralOperator<HelmholtzSLKernel<3>>;
  template class GenericIntegralOperator<HelmholtzDLKernel<3>>;    
  template class GenericIntegralOperator<CombinedFieldKernel<3>>;    
}

#include "ModifiedMeanCurvatureFlow.h"
#include "MeshAttribute.h"
#include "Laplacian.h"
#include "Mass.h"

#include <iostream>
using namespace std;

namespace DDG
{
   ModifiedMeanCurvatureFlow::ModifiedMeanCurvatureFlow( Mesh& mesh_ )
   : SurfaceFlow( mesh_ )
   {
      L0 = Laplacian<Real>::build( mesh );
      A0 = mesh.area();
   }

   void ModifiedMeanCurvatureFlow::integrate( double t )
   {
      int t0 = clock();
      int nV = mesh.vertices.size();
      DenseMatrix<Real> u( nV, 1 );

      updateFlowOperator( t );

      for( int k = 0; k < 3; k++ )
      {
         getAttribute( mesh, VertexPosition( k ), u );
         u = M*u;

         solvePositiveDefinite( H, u, u );

         setAttribute( mesh, VertexPosition( k ), u );
      }

      mesh.center();
      mesh.normalizeArea();

      int t1 = clock();
      double dt = seconds( t0, t1 );
      wallTime += dt;
      nIter++;
      cerr << "[ModifiedMeanCurvatureFlow] time: " << dt << "s" << endl;
   }

   void ModifiedMeanCurvatureFlow::updateFlowOperator( double t )
   {
      M = Mass<Real>::build( mesh );
      H = M + Real(t)*L0;
   }
}


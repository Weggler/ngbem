#ifndef TEST_COMPRESSION
#define TEST_COMPRESSION

namespace ngbem
{
  
  tuple<int, double, double> TestCompressionSVD (int nx, int ny, double eta, double eps);
  tuple<int, double, double> TestCompressionTSVD (int nx, int ny, double eta, double eps);
  tuple<int, double, double> TestCompressionACA (int nx, int ny, double eta, double eps);  
  
}


#endif

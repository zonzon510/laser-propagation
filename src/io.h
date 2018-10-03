#ifndef IO_H_
#define IO_H_

#include <vector>
#include <complex>


namespace IO {
  void read(const std::string& filename, std::vector<double>& x, std::vector<double>& y);
  
  void write(const std::string& filename, const std::vector<double>& data);
  void write(const std::string& filename, const std::vector<std::complex<double>>& data,
	     int Nrows, int Ncols);
  void write_binary(const std::string& filename,
		    const std::vector<std::complex<double>>& data);

  void write_binary(const std::string& filename,
		    const std::vector<double>& data);

  std::string enumerate_filename(const std::string& filename, int i);


  void clear_contents(const std::string& filename);   // clear contents of file
  void write_append(const std::string& filename, double value);
  void write_append(const std::string& filename, const std::string& text);
}

#endif // IO_H_

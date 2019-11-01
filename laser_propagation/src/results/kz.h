#ifndef KZ_H_
#define KZ_H_

#include "result.h"

namespace Results {

  class Kz : public Result {
  public:
    Kz(const std::string& filename)
      :filename(filename) {}
    void notify(int, double, const SimulationData& data) override {
      IO::write("kz.dat", data.propagator.kz.vec(), data.field.Nkperp, data.field.Nomega);
    }
    void finalize() override {}
    
  private:
    std::string filename;
  };
  
}

#endif // KZ_H_

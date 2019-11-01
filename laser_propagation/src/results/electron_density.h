#ifndef ELECTRON_DENSITY_H_
#define ELECTRON_DENSITY_H_

#include "result.h"

namespace Results {

  class ElectronDensity : public Result {
  public:
    ElectronDensity(const std::string& filename);
    void notify(int current_step, double current_distance,
		const SimulationData& data) override;
    void finalize() override;
    
  private:
    int number;
    std::string filename;
  };
  
}

#endif // ELECTRON_DENSITY_H_

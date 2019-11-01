#ifndef MAX_DENSITY_H_
#define MAX_DENSITY_H_

#include "result.h"

namespace Results {

class MaxDensity : public Result {
public:
  MaxDensity(const std::string& filename);
  void notify(int current_step, double current_distance,
	      const SimulationData& data) override;
  void finalize() override;

private:
  std::string filename;
};

}

#endif // MAX_DENSITY_H_

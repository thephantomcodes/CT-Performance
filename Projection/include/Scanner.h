#ifndef SCANNER_H
#define SCANNER_H

#include <vector>
#include "ProjectionParameters.h"

namespace Scanner
{
	class Scanner
	{
		public:
			Scanner(Projection::ProjectionParameters &params);
			Projection::PointValue source;
			double detector_half_width;
			std::vector<Projection::PointValue> detectors;
	};
	
	void RotateScanner(Scanner *scanner, double theta);
}

#endif

#include <vector>
#include <iostream>
#include "Scanner.h"

namespace Scanner
{
	Scanner::Scanner(Projection::ProjectionParameters &params)
		: source(params.scanning_radius, 0, 0)
		, detector_half_width(0.5*params.detector_length/(double)params.num_detectors)
	{
		detectors.reserve(params.num_detectors+1);
		for(int i=0; i<detectors.capacity(); i++)
			detectors.push_back(Projection::PointValue(-params.scanning_radius, (i - (double)params.num_detectors/2)*2*detector_half_width, 0));
	}
	
	void RotateScanner(Scanner *scanner, double theta)
	{
		Projection::RotatePointValue(&scanner->source, theta);
		for(int i=0; i < scanner->detectors.size(); i++)
			Projection::RotatePointValue(&scanner->detectors[i], theta);
	}
}

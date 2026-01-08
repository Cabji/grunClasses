#include <vector>
#include "GrunObject.h"

#ifndef CLASS_NAME
#define CLASS_NAME "GrunStage"
#endif

class GrunStage
{
	private:

	std::vector<GrunObject>		m_objects;								// a vector of GrunObjects this GrunStage pwns
	std::string					m_name		= "";						// the name of this GrunStage

	public:

	bool			addGrunObject(GrunObject& gObject);
};
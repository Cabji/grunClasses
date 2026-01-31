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

	size_t		addGrunObject(	const	GrunObject& 	gObject);
	size_t		rmGrunObject(	const	int&			objIndex);
	void		setName(		const	std::string&	name);

	public:
	// ctr
	GrunStage(const std::string& stageName);
};
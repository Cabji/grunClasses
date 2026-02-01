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
	// ctr
	GrunStage(const std::string& stageName);

	size_t				addGrunObject(		const	GrunObject& 	gObject);
	size_t				createGrunObject(	const 	std::string&	shapeType,
											const 	std::string&	name,
													double			x,
													double			y,
													double			z,
											const	std::string&	areaType	= "horizontal",
											const	std::string&	stage		= "");
	const
		std::string		getName();
	size_t				rmGrunObject(		const	int&			objIndex);
	void				setName(			const	std::string&	name);
};